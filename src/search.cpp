#include "search.h"
#include "costFunctions.h"
#include "costFunctionsGPU.cuh"
#include "stdio.h"
#include "stdlib.h"
#include "memory.h"
#include "time.h"

int globalCount=0;

int coarseAngleCount = (int)((FLIRT_ANGLE_MAX - FLIRT_ANGLE_MIN) / FLIRT_COARSE_ANGLE_STEP) + 1;
int finerAngleCount = (int)((FLIRT_ANGLE_MAX - FLIRT_ANGLE_MIN) / FLIRT_FINER_ANGLE_STEP) + 1;

Transform12 registerVolume(Volume* baseVol1, Volume* regVol1)
{
	Volume* baseVol2 = halfSample(baseVol1);
	Volume* baseVol4 = halfSample(baseVol2);
	Volume* baseVol8 = halfSample(baseVol4);

	Volume* regVol2 = halfSample(regVol1);
	Volume* regVol4 = halfSample(regVol2);
	Volume* regVol8 = halfSample(regVol4);

	int t1=clock();

	Transform12 result;

	float medianScale=0;

	//8mm stage start
	initVolumeData(baseVol8,regVol8);

	float4* coarseGrid=new float4[coarseAngleCount*coarseAngleCount*coarseAngleCount];
	coarseRegister(baseVol8,regVol8,coarseGrid,medianScale);

	float4* finerGrid=new float4[finerAngleCount*finerAngleCount*finerAngleCount];
	float* finerResults=new float[finerAngleCount*finerAngleCount*finerAngleCount];

	buildFinerGrid(coarseGrid,finerGrid);

	finerRegister( baseVol8,  regVol8, finerGrid, finerResults, medianScale);

	freeVolumeData();

	int* minima=new int[finerAngleCount*finerAngleCount*finerAngleCount / 8];

	int minimaCount=0;

	findMinima(finerResults,minima,minimaCount);

	Transform7 globalMinima[FLIRT_MINIMA_COUNT];

	sortMinima(finerGrid,finerResults,minima,minimaCount,globalMinima);

	printf("MedianScale is:%f,minimaCount is%d\n",medianScale,minimaCount);

	//4mm stage start
	initVolumeData(baseVol4,regVol4);

	Transform7 higherResStart[FLIRT_HIGHER_RESOLUTION_START];
	Transform12 t12[FLIRT_HIGHER_RESOLUTION_START];
	float finalResult[FLIRT_HIGHER_RESOLUTION_START];

	perturbationStage(baseVol4,regVol4,globalMinima,higherResStart);

	freeVolumeData();

	//2mm stage start
	initVolumeData(baseVol2,regVol2);

	for(int i=0;i<FLIRT_HIGHER_RESOLUTION_START;i++)
		t12[i]=optimize2mm(baseVol2,regVol2,higherResStart[i]);

	freeVolumeData();

	//1mm stage start
	initVolumeData(baseVol1,regVol1);

	int t2=clock();

	printf("\n\ntime used:%dms,global count is:%d\n",t2-t1,globalCount);

	float minValue=0;

	for(int i=0;i<FLIRT_HIGHER_RESOLUTION_START;i++)
	{
		Transform12 T=t12[i];
		finalResult[i]=optimize1mm(baseVol1,regVol1,T);
		t12[i]=T;

		if(minValue>finalResult[i])result=t12[i];

		printf("12 DOF 1mm:rotation(%f,%f,%f),trans(%f,%f,%f),scale(%f,%f,%f),skew(%f,%f,%f),value is:%f\n",T.rotation.x,T.rotation.y,T.rotation.z,T.translation.x,T.translation.y,T.translation.z,T.scale.x,T.scale.y,T.scale.z,T.skew.x,T.skew.y,T.skew.z,finalResult[i]);
	}

	freeVolumeData();

	t2=clock();

	printf("\n\ntime used:%dms,global count is:%d\n",t2-t1,globalCount);

	releaseVolume(baseVol2);
	releaseVolume(baseVol4);
	releaseVolume(baseVol8);

	releaseVolume(regVol2);
	releaseVolume(regVol4);
	releaseVolume(regVol8);

	return result;
}

float optimize1mm(Volume* baseVol,Volume* regVol,Transform12& t12)
{
	float value;

	int pass=0;

	while(pass < FLIRT_SEARCH_PASS)

	{
		pass++;
		
		value=goldenSearch12dof(baseVol,regVol,t12,"rz",FLIRT_FINER_ANGLE_STEP / 8.0f);
		value=goldenSearch12dof(baseVol,regVol,t12,"rx",FLIRT_FINER_ANGLE_STEP / 8.0f);
		value=goldenSearch12dof(baseVol,regVol,t12,"ry",FLIRT_FINER_ANGLE_STEP / 8.0f);

		value=goldenSearch12dof(baseVol,regVol,t12,"tz",FLIRT_TRANS_Z_STEP / 8.0f);
		value=goldenSearch12dof(baseVol,regVol,t12,"ty",FLIRT_TRANS_XY_STEP / 8.0f);
		value=goldenSearch12dof(baseVol,regVol,t12,"tx",FLIRT_TRANS_XY_STEP / 8.0f);

		value=goldenSearch12dof(baseVol,regVol,t12,"sx",FLIRT_SCALE_STEP / 8.0f);
		value=goldenSearch12dof(baseVol,regVol,t12,"sy",FLIRT_SCALE_STEP / 8.0f);
		value=goldenSearch12dof(baseVol,regVol,t12,"sz",FLIRT_SCALE_STEP / 8.0f);

		value=goldenSearch12dof(baseVol,regVol,t12,"gx",FLIRT_SKEW_STEP/2.0f );
		value=goldenSearch12dof(baseVol,regVol,t12,"gy",FLIRT_SKEW_STEP/2.0f );
		value=goldenSearch12dof(baseVol,regVol,t12,"gz",FLIRT_SKEW_STEP/2.0f );
	}

	return value;
}

Transform12 optimize2mm(Volume* baseVol,Volume* regVol,Transform7 higherResStart)
{
	Transform12 t12;

	perturbationZ(baseVol,regVol,higherResStart);

	t12.rotation=higherResStart.rotation;
	t12.translation.x=higherResStart.translateAndScale.x;
	t12.translation.y=higherResStart.translateAndScale.y;
	t12.translation.z=higherResStart.translateAndScale.z;
	t12.scale.x=t12.scale.y=t12.scale.z=higherResStart.translateAndScale.w;
	t12.skew.x=t12.skew.y=t12.skew.z=0;

	float value = optimize2mm9dof(baseVol,regVol,t12);
	

	printf("9 DOF 2mm:rotation(%f,%f,%f),trans(%f,%f,%f),scale(%f,%f,%f),value is:%f\n",t12.rotation.x,t12.rotation.y,t12.rotation.z,t12.translation.x,t12.translation.y,t12.translation.z,t12.scale.x,t12.scale.y,t12.scale.z,value);
	
	optimize2mm12dof(baseVol,regVol,t12);

	printf("12 DOF 2mm:rotation(%f,%f,%f),trans(%f,%f,%f),scale(%f,%f,%f),skew(%f,%f,%f),value is:%f\n",t12.rotation.x,t12.rotation.y,t12.rotation.z,t12.translation.x,t12.translation.y,t12.translation.z,t12.scale.x,t12.scale.y,t12.scale.z,t12.skew.x,t12.skew.y,t12.skew.z,value);
	
	return t12;
}

float optimize2mm9dof(Volume* baseVol,Volume* regVol,Transform12& t12)
{
	float value;
	
	
	value=goldenSearch9dof(baseVol,regVol,t12,"rz",FLIRT_FINER_ANGLE_STEP / 4.0);
	value=goldenSearch9dof(baseVol,regVol,t12,"rx",FLIRT_FINER_ANGLE_STEP / 4.0);
	value=goldenSearch9dof(baseVol,regVol,t12,"ry",FLIRT_FINER_ANGLE_STEP / 4.0);

	value=goldenSearch9dof(baseVol,regVol,t12,"tx",FLIRT_TRANS_XY_STEP / 4.0);
	value=goldenSearch9dof(baseVol,regVol,t12,"ty",FLIRT_TRANS_XY_STEP / 4.0);
	value=goldenSearch9dof(baseVol,regVol,t12,"tz",FLIRT_TRANS_Z_STEP / 4.0);

	value=goldenSearch9dof(baseVol,regVol,t12,"sx",FLIRT_SCALE_STEP / 4.0);
	value=goldenSearch9dof(baseVol,regVol,t12,"sy",FLIRT_SCALE_STEP / 4.0);
	value=goldenSearch9dof(baseVol,regVol,t12,"sz",FLIRT_SCALE_STEP / 4.0);

	return value;
}

float optimize2mm12dof(Volume* baseVol,Volume* regVol,Transform12& t12)
{
	float value;
	
	value=goldenSearch12dof(baseVol,regVol,t12,"rz",FLIRT_FINER_ANGLE_STEP / 4.0);
	value=goldenSearch12dof(baseVol,regVol,t12,"rx",FLIRT_FINER_ANGLE_STEP / 4.0);
	value=goldenSearch12dof(baseVol,regVol,t12,"ry",FLIRT_FINER_ANGLE_STEP / 4.0);

	value=goldenSearch12dof(baseVol,regVol,t12,"tz",FLIRT_TRANS_Z_STEP / 4.0);
	value=goldenSearch12dof(baseVol,regVol,t12,"ty",FLIRT_TRANS_XY_STEP / 4.0);
	value=goldenSearch12dof(baseVol,regVol,t12,"tx",FLIRT_TRANS_XY_STEP / 4.0);
	
	value=goldenSearch12dof(baseVol,regVol,t12,"sx",FLIRT_SCALE_STEP / 4.0);
	value=goldenSearch12dof(baseVol,regVol,t12,"sy",FLIRT_SCALE_STEP / 4.0);
	value=goldenSearch12dof(baseVol,regVol,t12,"sz",FLIRT_SCALE_STEP / 4.0);

	value=goldenSearch12dof(baseVol,regVol,t12,"gx",FLIRT_SKEW_STEP );
	value=goldenSearch12dof(baseVol,regVol,t12,"gy",FLIRT_SKEW_STEP );
	value=goldenSearch12dof(baseVol,regVol,t12,"gz",FLIRT_SKEW_STEP );

	return value;
}


void perturbationZ(Volume* baseVol,Volume* regVol,Transform7& higherResStart)
{
	float result2mm[7];

	Transform7 grid2mm[7];

	float minValue=1000;
	int minIndex=-1;

	for(int i=0;i<7;i++)
	{
		grid2mm[i]=higherResStart;
		grid2mm[i].translateAndScale.z += (i-3);

		result2mm[i]=singlePerturbationZ(baseVol,regVol,grid2mm[i]);

		if(result2mm[i] < minValue)
		{
			minValue=result2mm[i];
			minIndex = i;
		}

		printf("%d:rotation(%f,%f,%f),trans(%f,%f,%f),scale(%f),value is:%f\n",i,grid2mm[i].rotation.x,grid2mm[i].rotation.y,grid2mm[i].rotation.z,grid2mm[i].translateAndScale.x,grid2mm[i].translateAndScale.y,grid2mm[i].translateAndScale.z,grid2mm[i].translateAndScale.w,result2mm[i]);
	}

	higherResStart=grid2mm[minIndex];
}

float singlePerturbationZ(Volume* baseVol,Volume* regVol,Transform7& rotTransScale)
{
	float value;
	
	value=goldenSearch7dof(baseVol,regVol,rotTransScale,"rz",FLIRT_FINER_ANGLE_STEP / 4.0);
	value=goldenSearch7dof(baseVol,regVol,rotTransScale,"rx",FLIRT_FINER_ANGLE_STEP / 4.0);
	value=goldenSearch7dof(baseVol,regVol,rotTransScale,"ry",FLIRT_FINER_ANGLE_STEP / 4.0);

	value=goldenSearch7dof(baseVol,regVol,rotTransScale,"tx",FLIRT_TRANS_XY_STEP / 4.0);
	value=goldenSearch7dof(baseVol,regVol,rotTransScale,"ty",FLIRT_TRANS_XY_STEP / 4.0);

	value=goldenSearch7dof(baseVol,regVol,rotTransScale,"scale",FLIRT_SCALE_STEP / 4.0);

	return value;
}

float goldenSearch12dof(Volume* baseVol,Volume* regVol,Transform12 &rotTransScale,char* target,float minStep)
{
	float *p;

	if(strcmp(target,"rx")==0)
		p=&(rotTransScale.rotation.x);
	else if(strcmp(target,"ry")==0)
		p=&(rotTransScale.rotation.y);
	else if(strcmp(target,"rz")==0)
		p=&(rotTransScale.rotation.z);
	else if(strcmp(target,"tx")==0)
		p=&(rotTransScale.translation.x);
	else if(strcmp(target,"ty")==0)
		p=&(rotTransScale.translation.y);
	else if(strcmp(target,"tz")==0)
		p=&(rotTransScale.translation.z);
	else if(strcmp(target,"sx")==0)
		p=&(rotTransScale.scale.x);
	else if(strcmp(target,"sy")==0)
		p=&(rotTransScale.scale.y);
	else if(strcmp(target,"sz")==0)
		p=&(rotTransScale.scale.z);
	else if(strcmp(target,"gx")==0)
		p=&(rotTransScale.skew.x);
	else if(strcmp(target,"gy")==0)
		p=&(rotTransScale.skew.y);
	else if(strcmp(target,"gz")==0)
		p=&(rotTransScale.skew.z);
	else return 0;

	float points[4];

	if(strcmp(target,"rx")==0 || strcmp(target,"ry")==0 || strcmp(target,"rz")==0)
	{
		float tmpPoints[4]=
		{
			FLIRT_ANGLE_MIN,
			*p,
			FLIRT_ANGLE_MAX,
			0
		};

		memcpy(points,tmpPoints,sizeof(float)*4);
	}
	else if(strcmp(target,"tx")==0 || strcmp(target,"ty")==0)
	{
		float tmpPoints[4]=
		{
			FLIRT_TRANS_XY_MIN,
			*p,
			FLIRT_TRANS_XY_MAX,
			0
		};
		memcpy(points,tmpPoints,sizeof(float)*4);
	}
	else if( strcmp(target,"tz")==0)
	{
		float tmpPoints[4]=
		{
			FLIRT_TRANS_Z_MIN,
			*p,
			FLIRT_TRANS_Z_MAX,
			0
		};
		memcpy(points,tmpPoints,sizeof(float)*4);
	}
	else if(strcmp(target,"sx")==0 || strcmp(target,"sy")==0 || strcmp(target,"sz")==0)
	{
		float tmpPoints[4]=
		{
			FLIRT_SCALE_MIN,
			*p,
			FLIRT_SCALE_MAX,
			0
		};
		memcpy(points,tmpPoints,sizeof(float)*4);
	}
	else 
	{
		float tmpPoints[4]=
		{
			FLIRT_SKEW_MIN,
			*p,
			FLIRT_SKEW_MAX,
			0
		};
		memcpy(points,tmpPoints,sizeof(float)*4);
	}

	float values[4]={0,0,0,0};

	

	Matrix44* affineMat=new Matrix44();

	Matrix44 rotationMat,translationMat,scaleMat,tmpMat,skewMat;

	rotationMat.setRotateAxisAngle(1,0,0,rotTransScale.rotation.x);
	tmpMat=rotationMat;

	rotationMat.setRotateAxisAngle(0,1,0,rotTransScale.rotation.y);
	tmpMat=rotationMat*tmpMat;

	rotationMat.setRotateAxisAngle(0,0,1,rotTransScale.rotation.z);
	rotationMat=rotationMat*tmpMat;
	
	scaleMat.setScale(rotTransScale.scale.x,rotTransScale.scale.y,rotTransScale.scale.z);
	translationMat.setTranslate(rotTransScale.translation.x,rotTransScale.translation.y,rotTransScale.translation.z,baseVol->dim.x);
	skewMat.setSkew(rotTransScale.skew.x,rotTransScale.skew.y,rotTransScale.skew.z);

	tmpMat=skewMat*scaleMat;
	tmpMat=rotationMat*tmpMat;
	*affineMat=translationMat*tmpMat;
	
	//int clock1=clock();

	values[1]=calcCostFunction(baseVol,regVol,affineMat);
	
	//int clock2=clock();

	//printf("1mm cost function time:%d ms\n",clock2-clock1);

	


	float d1,d2;

	int count=0;

	while(true)
	{
		d1=points[1]-points[0];
		d2=points[2]-points[1];

		count++;

		//find the logner section, and calculate the 4th point from current 3 points.

		//d2 is longer
		if(d1 < d2)
		{
			points[3]=points[2];
			values[3]=values[2];

			//insert a point into p0...p1...<<p>>...p2
			points[2]=points[1]+d2 * 0.61803f;
			*p=points[2];
			rotationMat.setRotateAxisAngle(1,0,0,rotTransScale.rotation.x);
			tmpMat=rotationMat;

			rotationMat.setRotateAxisAngle(0,1,0,rotTransScale.rotation.y);
			tmpMat=rotationMat*tmpMat;

			rotationMat.setRotateAxisAngle(0,0,1,rotTransScale.rotation.z);
			rotationMat=rotationMat*tmpMat;
			
			scaleMat.setScale(rotTransScale.scale.x,rotTransScale.scale.y,rotTransScale.scale.z);
			translationMat.setTranslate(rotTransScale.translation.x,rotTransScale.translation.y,rotTransScale.translation.z,baseVol->dim.x);
			skewMat.setSkew(rotTransScale.skew.x,rotTransScale.skew.y,rotTransScale.skew.z);

			tmpMat=skewMat*scaleMat;
			tmpMat=rotationMat*tmpMat;

			tmpMat=rotationMat*scaleMat;
			*affineMat=translationMat*tmpMat;
			values[2]=calcCostFunction(baseVol,regVol,affineMat);
		}
		else 
		{
			points[3]=points[2];
			points[2]=points[1];
			values[3]=values[2];
			values[2]=values[1];

			//insert a point into p0...<<p>>...p1..p2
			points[1]=points[0]+ d1 * 0.38197f;

			*p=points[1];
			rotationMat.setRotateAxisAngle(1,0,0,rotTransScale.rotation.x);
			tmpMat=rotationMat;

			rotationMat.setRotateAxisAngle(0,1,0,rotTransScale.rotation.y);
			tmpMat=rotationMat*tmpMat;

			rotationMat.setRotateAxisAngle(0,0,1,rotTransScale.rotation.z);
			rotationMat=rotationMat*tmpMat;
			
			scaleMat.setScale(rotTransScale.scale.x,rotTransScale.scale.y,rotTransScale.scale.z);
			translationMat.setTranslate(rotTransScale.translation.x,rotTransScale.translation.y,rotTransScale.translation.z,baseVol->dim.x);
			skewMat.setSkew(rotTransScale.skew.x,rotTransScale.skew.y,rotTransScale.skew.z);

			tmpMat=skewMat*scaleMat;
			tmpMat=rotationMat*tmpMat;

			tmpMat=rotationMat*scaleMat;
			*affineMat=translationMat*tmpMat;
			values[1]=calcCostFunction(baseVol,regVol,affineMat);
		}

		//form new 3 points group.

		//The latest 3 points come to be the 3 poits group, shift to left <---------
		if(values[1] > values[2])
		{
			points[0]=points[1];
			points[1]=points[2];
			points[2]=points[3];

			values[0]=values[1];
			values[1]=values[2];
			values[2]=values[3];
		}
		
		

		if(d1< minStep || d2< minStep)
			break;
	}
	
	*p=points[1];

	delete affineMat;
	
	return values[1];
}


float goldenSearch9dof(Volume* baseVol,Volume* regVol,Transform12 &rotTransScale,char* target,float minStep)
{
	float *p;

	if(strcmp(target,"rx")==0)
		p=&(rotTransScale.rotation.x);
	else if(strcmp(target,"ry")==0)
		p=&(rotTransScale.rotation.y);
	else if(strcmp(target,"rz")==0)
		p=&(rotTransScale.rotation.z);
	else if(strcmp(target,"tx")==0)
		p=&(rotTransScale.translation.x);
	else if(strcmp(target,"ty")==0)
		p=&(rotTransScale.translation.y);
	else if(strcmp(target,"tz")==0)
		p=&(rotTransScale.translation.z);
	else if(strcmp(target,"sx")==0)
		p=&(rotTransScale.scale.x);
	else if(strcmp(target,"sy")==0)
		p=&(rotTransScale.scale.y);
	else if(strcmp(target,"sz")==0)
		p=&(rotTransScale.scale.z);
	else return 0;

	float points[4];

	if(strcmp(target,"rx")==0 || strcmp(target,"ry")==0 || strcmp(target,"rz")==0)
	{
		float tmpPoints[4]=
		{
			FLIRT_ANGLE_MIN,
			*p,
			FLIRT_ANGLE_MAX,
			0
		};

		memcpy(points,tmpPoints,sizeof(float)*4);
	}
	else if(strcmp(target,"tx")==0 || strcmp(target,"ty")==0)
	{
		float tmpPoints[4]=
		{
			FLIRT_TRANS_XY_MIN,
			*p,
			FLIRT_TRANS_XY_MAX,
			0
		};
		memcpy(points,tmpPoints,sizeof(float)*4);
	}
	else if( strcmp(target,"tz")==0)
	{
		float tmpPoints[4]=
		{
			FLIRT_TRANS_Z_MIN,
			*p,
			FLIRT_TRANS_Z_MAX,
			0
		};
		memcpy(points,tmpPoints,sizeof(float)*4);
	}
	else 
	{
		float tmpPoints[4]=
		{
			FLIRT_SCALE_MIN,
			*p,
			FLIRT_SCALE_MAX,
			0
		};
		memcpy(points,tmpPoints,sizeof(float)*4);
	}

	float values[4]={0,0,0,0};

	Matrix44* affineMat=new Matrix44();

	Matrix44 rotationMat,translationMat,scaleMat,tmpMat;

	rotationMat.setRotateAxisAngle(1,0,0,rotTransScale.rotation.x);
	tmpMat=rotationMat;

	rotationMat.setRotateAxisAngle(0,1,0,rotTransScale.rotation.y);
	tmpMat=rotationMat*tmpMat;

	rotationMat.setRotateAxisAngle(0,0,1,rotTransScale.rotation.z);
	rotationMat=rotationMat*tmpMat;
	
	scaleMat.setScale(rotTransScale.scale.x,rotTransScale.scale.y,rotTransScale.scale.z);
	translationMat.setTranslate(rotTransScale.translation.x,rotTransScale.translation.y,rotTransScale.translation.z,baseVol->dim.x);

	tmpMat=rotationMat*scaleMat;
	*affineMat=translationMat*tmpMat;
	values[1]=calcCostFunction(baseVol,regVol,affineMat);


	float d1,d2;

	int count=0;

	while(true)
	{
		d1=points[1]-points[0];
		d2=points[2]-points[1];

		count++;

		//find the logner section, and calculate the 4th point from current 3 points.

		//d2 is longer
		if(d1 < d2)
		{
			points[3]=points[2];
			values[3]=values[2];

			//insert a point into p0...p1...<<p>>...p2
			points[2]=points[1]+d2 * 0.61803f;
			*p=points[2];
			rotationMat.setRotateAxisAngle(1,0,0,rotTransScale.rotation.x);
			tmpMat=rotationMat;

			rotationMat.setRotateAxisAngle(0,1,0,rotTransScale.rotation.y);
			tmpMat=rotationMat*tmpMat;

			rotationMat.setRotateAxisAngle(0,0,1,rotTransScale.rotation.z);
			rotationMat=rotationMat*tmpMat;
			
			scaleMat.setScale(rotTransScale.scale.x,rotTransScale.scale.y,rotTransScale.scale.z);
			translationMat.setTranslate(rotTransScale.translation.x,rotTransScale.translation.y,rotTransScale.translation.z,baseVol->dim.x);

			tmpMat=rotationMat*scaleMat;
			*affineMat=translationMat*tmpMat;
			values[2]=calcCostFunction(baseVol,regVol,affineMat);
		}
		else 
		{
			points[3]=points[2];
			points[2]=points[1];
			values[3]=values[2];
			values[2]=values[1];

			//insert a point into p0...<<p>>...p1..p2
			points[1]=points[0]+ d1 * 0.38197f;

			*p=points[1];
			rotationMat.setRotateAxisAngle(1,0,0,rotTransScale.rotation.x);
			tmpMat=rotationMat;

			rotationMat.setRotateAxisAngle(0,1,0,rotTransScale.rotation.y);
			tmpMat=rotationMat*tmpMat;

			rotationMat.setRotateAxisAngle(0,0,1,rotTransScale.rotation.z);
			rotationMat=rotationMat*tmpMat;
			
			scaleMat.setScale(rotTransScale.scale.x,rotTransScale.scale.y,rotTransScale.scale.z);
			translationMat.setTranslate(rotTransScale.translation.x,rotTransScale.translation.y,rotTransScale.translation.z,baseVol->dim.x);

			tmpMat=rotationMat*scaleMat;
			*affineMat=translationMat*tmpMat;
			values[1]=calcCostFunction(baseVol,regVol,affineMat);
		}

		//form new 3 points group.

		//The latest 3 points come to be the 3 poits group, shift to left <---------
		if(values[1] > values[2])
		{
			points[0]=points[1];
			points[1]=points[2];
			points[2]=points[3];

			values[0]=values[1];
			values[1]=values[2];
			values[2]=values[3];
		}
		
		

		if(d1< minStep || d2< minStep)
			break;
	}
	
	*p=points[1];

	delete affineMat;
	
	return values[1];
}

float goldenSearch7dof(Volume* baseVol,Volume* regVol,Transform7 &rotTransScale,char* target,float minStep)
{
	float *p;

	if(strcmp(target,"rx")==0)
		p=&(rotTransScale.rotation.x);
	else if(strcmp(target,"ry")==0)
		p=&(rotTransScale.rotation.y);
	else if(strcmp(target,"rz")==0)
		p=&(rotTransScale.rotation.z);
	else if(strcmp(target,"tx")==0)
		p=&(rotTransScale.translateAndScale.x);
	else if(strcmp(target,"ty")==0)
		p=&(rotTransScale.translateAndScale.y);
	else if(strcmp(target,"tz")==0)
		p=&(rotTransScale.translateAndScale.z);
	else if(strcmp(target,"scale")==0)
		p=&(rotTransScale.translateAndScale.w);
	else return 0;

	float points[4];

	if(strcmp(target,"rx")==0|| strcmp(target,"ry")==0|| strcmp(target,"rz")==0)
	{
		float tmpPoints[4]=
		{
			FLIRT_ANGLE_MIN,
			*p,
			FLIRT_ANGLE_MAX,
			0
		};

		memcpy(points,tmpPoints,sizeof(float)*4);
	}
	else if(strcmp(target,"tx")==0 || strcmp(target,"ty")==0)
	{
		float tmpPoints[4]=
		{
			FLIRT_TRANS_XY_MIN,
			*p,
			FLIRT_TRANS_XY_MAX,
			0
		};
		memcpy(points,tmpPoints,sizeof(float)*4);
	}
	else if( strcmp(target,"tz")==0)
	{
		float tmpPoints[4]=
		{
			FLIRT_TRANS_Z_MIN,
			*p,
			FLIRT_TRANS_Z_MAX,
			0
		};
		memcpy(points,tmpPoints,sizeof(float)*4);
	}
	else 
	{
		float tmpPoints[4]=
		{
			FLIRT_SCALE_MIN,
			*p,
			FLIRT_SCALE_MAX,
			0
		};
		memcpy(points,tmpPoints,sizeof(float)*4);
	}

	float values[4]={0,0,0,0};

	Matrix44* affineMat=new Matrix44();

	Matrix44 rotationMat,translationMat,scaleMat,tmpMat;

	rotationMat.setRotateAxisAngle(1,0,0,rotTransScale.rotation.x);
	tmpMat=rotationMat;

	rotationMat.setRotateAxisAngle(0,1,0,rotTransScale.rotation.y);
	tmpMat=rotationMat*tmpMat;

	rotationMat.setRotateAxisAngle(0,0,1,rotTransScale.rotation.z);
	rotationMat=rotationMat*tmpMat;
	
	scaleMat.setScale(rotTransScale.translateAndScale.w,rotTransScale.translateAndScale.w,rotTransScale.translateAndScale.w);
	translationMat.setTranslate(rotTransScale.translateAndScale.x,rotTransScale.translateAndScale.y,rotTransScale.translateAndScale.z,baseVol->dim.x);

	tmpMat=rotationMat*scaleMat;
	*affineMat=translationMat*tmpMat;
	values[1]=calcCostFunction(baseVol,regVol,affineMat);


	float d1,d2;

	int count=0;

	while(true)
	{
		d1=points[1]-points[0];
		d2=points[2]-points[1];

		count++;

		//find the logner section, and calculate the 4th point from current 3 points.

		//d2 is longer
		if(d1 < d2)
		{
			points[3]=points[2];
			values[3]=values[2];

			//insert a point into p0...p1...<<p>>...p2
			points[2]=points[1]+d2 * 0.61803f;
			*p=points[2];
			rotationMat.setRotateAxisAngle(1,0,0,rotTransScale.rotation.x);
			tmpMat=rotationMat;

			rotationMat.setRotateAxisAngle(0,1,0,rotTransScale.rotation.y);
			tmpMat=rotationMat*tmpMat;

			rotationMat.setRotateAxisAngle(0,0,1,rotTransScale.rotation.z);
			rotationMat=rotationMat*tmpMat;
			
			scaleMat.setScale(rotTransScale.translateAndScale.w,rotTransScale.translateAndScale.w,rotTransScale.translateAndScale.w);
			translationMat.setTranslate(rotTransScale.translateAndScale.x,rotTransScale.translateAndScale.y,rotTransScale.translateAndScale.z,baseVol->dim.x);

			tmpMat=rotationMat*scaleMat;
			*affineMat=translationMat*tmpMat;
			values[2]=calcCostFunction(baseVol,regVol,affineMat);
		}
		else 
		{
			points[3]=points[2];
			points[2]=points[1];
			values[3]=values[2];
			values[2]=values[1];

			//insert a point into p0...<<p>>...p1..p2
			points[1]=points[0]+ d1 * 0.38197f;

			*p=points[1];
			rotationMat.setRotateAxisAngle(1,0,0,rotTransScale.rotation.x);
			tmpMat=rotationMat;

			rotationMat.setRotateAxisAngle(0,1,0,rotTransScale.rotation.y);
			tmpMat=rotationMat*tmpMat;

			rotationMat.setRotateAxisAngle(0,0,1,rotTransScale.rotation.z);
			rotationMat=rotationMat*tmpMat;
			
			scaleMat.setScale(rotTransScale.translateAndScale.w,rotTransScale.translateAndScale.w,rotTransScale.translateAndScale.w);
			translationMat.setTranslate(rotTransScale.translateAndScale.x,rotTransScale.translateAndScale.y,rotTransScale.translateAndScale.z,baseVol->dim.x);

			tmpMat=rotationMat*scaleMat;
			*affineMat=translationMat*tmpMat;
			values[1]=calcCostFunction(baseVol,regVol,affineMat);
		}

		//form new 3 points group.

		//The latest 3 points come to be the 3 poits group, shift to left <---------
		if(values[1] > values[2])
		{
			points[0]=points[1];
			points[1]=points[2];
			points[2]=points[3];

			values[0]=values[1];
			values[1]=values[2];
			values[2]=values[3];
		}
		
		

		if(d1< minStep || d2< minStep)
			break;
	}
	
	*p=points[1];

	delete affineMat;
	
	return values[1];
}

void perturbationStage(Volume* baseVol,Volume* regVol,Transform7* globalMinima,Transform7* higherResStart)
{

	Transform7 grid4mm[FLIRT_MINIMA_COUNT * 27];
	float result4mm[FLIRT_MINIMA_COUNT * 27];

	for(int i=0;i<FLIRT_MINIMA_COUNT;i++)
	{
		perturbation(baseVol,regVol,globalMinima[i],grid4mm+(i*27),result4mm+(i*27));
	}


	for(int i=0;i<FLIRT_HIGHER_RESOLUTION_START; i++)
	{
		for(int j=i;j<FLIRT_MINIMA_COUNT*27;j++)
		{
			if(result4mm[i] > result4mm[j])
			{
				Transform7 tmpTrans=grid4mm[i];
				grid4mm[i]=grid4mm[j];
				grid4mm[j]=tmpTrans;

				float tmpValue=result4mm[i];
				result4mm[i]=result4mm[j];
				result4mm[j]=tmpValue;
			}
		}

		higherResStart[i]=grid4mm[i];

		//printf("%d:rotation(%f,%f,%f),trans(%f,%f,%f),scale(%f),value is:%f\n",i,grid4mm[i].rotation.x,grid4mm[i].rotation.y,grid4mm[i].rotation.z,grid4mm[i].translateAndScale.x,grid4mm[i].translateAndScale.y,grid4mm[i].translateAndScale.z,grid4mm[i].translateAndScale.w,result4mm[i]);
	}
}

void perturbation(Volume* baseVol,Volume* regVol,Transform7 globalMinima,Transform7* grid4mm, float* result4mm )
{
	Matrix44 matXYZ,matYZ,matZ,tmpMat;

	float x,y,z;

	int pass;

	//float value;

	for(int k=-1;k<=1;k++)
	{
		z=globalMinima.rotation.z + k*(FLIRT_FINER_ANGLE_STEP)/2.0f;
		matZ.setRotateAxisAngle(0,0,1.0,z);

		for(int j=-1;j<=1;j++)
		{
			y=globalMinima.rotation.y + j*(FLIRT_FINER_ANGLE_STEP)/2.0f;
			tmpMat.setRotateAxisAngle(0,1.0,0,y);

			matYZ=matZ*tmpMat;

			for(int i=-1;i<=1;i++)
			{
				x=globalMinima.rotation.x + i*(FLIRT_FINER_ANGLE_STEP)/2.0f;

				int index=(k+1)*9 + (j+1)*3 + (i+1);

				grid4mm[index].rotation.x = x;
				grid4mm[index].rotation.y = y;
				grid4mm[index].rotation.z = z;

				tmpMat.setRotateAxisAngle(1.0,0,0,x);

				matXYZ=matYZ*tmpMat;

				//printf("%f,%f,%f:",(x*FLIRT_COARSE_ANGLE_STEP + FLIRT_ANGLE_MIN),(y*FLIRT_COARSE_ANGLE_STEP + FLIRT_ANGLE_MIN),(z*FLIRT_COARSE_ANGLE_STEP + FLIRT_ANGLE_MIN));

				pass=0;

				float4 translateAndScale=globalMinima.translateAndScale;

				while( pass < FLIRT_SEARCH_PASS )
				{
					pass++;
					result4mm[index]=goldenSearchTranslation(baseVol,regVol,&matXYZ,translateAndScale,'z',FLIRT_TRANS_Z_STEP / 2.0);
					result4mm[index]=goldenSearchTranslation(baseVol,regVol,&matXYZ,translateAndScale,'y',FLIRT_TRANS_XY_STEP / 2.0);
					result4mm[index]=goldenSearchTranslation(baseVol,regVol,&matXYZ,translateAndScale,'x',FLIRT_TRANS_XY_STEP / 2.0);
					
					result4mm[index]=goldenSearchScale(baseVol,regVol,&matXYZ,translateAndScale,FLIRT_SCALE_STEP / 2.0);
				}

				grid4mm[index].translateAndScale=translateAndScale;

			}
		}
	}
}

void sortMinima(float4* finerGrid,float* finerResults,int* minima, int count,Transform7* globalMinima)
{
	for(int i=0;i<FLIRT_MINIMA_COUNT; i++)
	{
		for(int j=i;j<count;j++)
		{
			if(finerResults[minima[i]] > finerResults[minima[j]])
			{
				int tmp=minima[i];
				minima[i]=minima[j];
				minima[j]=tmp;
			}
		}

		int3 r;
		
		r.z=minima[i]/(finerAngleCount*finerAngleCount);
		r.y=(minima[i]-r.z*(finerAngleCount*finerAngleCount))/finerAngleCount;
		r.x=minima[i] % finerAngleCount;

		globalMinima[i].translateAndScale=finerGrid[minima[i]];
		globalMinima[i].rotation.x=r.x*FLIRT_FINER_ANGLE_STEP + FLIRT_ANGLE_MIN;
		globalMinima[i].rotation.y=r.y*FLIRT_FINER_ANGLE_STEP + FLIRT_ANGLE_MIN;
		globalMinima[i].rotation.z=r.z*FLIRT_FINER_ANGLE_STEP + FLIRT_ANGLE_MIN;

		printf("%d:rotation(%f,%f,%f),trans(%f,%f,%f),scale(%f),value is:%f\n",i,globalMinima[i].rotation.x,globalMinima[i].rotation.y,globalMinima[i].rotation.z,globalMinima[i].translateAndScale.x,globalMinima[i].translateAndScale.y,globalMinima[i].translateAndScale.z,globalMinima[i].translateAndScale.w,finerResults[minima[i]]);
	}
}

bool isMinima(float* finerResults,int x,int y,int z)
{
	int index=z*finerAngleCount*finerAngleCount + y*finerAngleCount + x;

	for(int k=-1;k<=1;k+=2)
		for(int j=-1;j<=1;j+=2)
			for(int i=-1;i<=1;i+=2)
			{
				int neighborIndex = (z+k)*finerAngleCount*finerAngleCount + (y+j)*finerAngleCount + (x+i);

				if(finerResults[index] > finerResults[neighborIndex])
					return false;
			}

	return true;
}

void findMinima(float* finerResults,int* minima, int& count)
{
	count=0;
	for(int z=1;z<finerAngleCount-1;z++)
	{
		for(int y=1;y<finerAngleCount-1;y++)
		{
			for(int x=1;x<finerAngleCount-1;x++)
			{
				if(isMinima(finerResults,x,y,z))
				{
					minima[count]=z*finerAngleCount*finerAngleCount + y*finerAngleCount + x;
					count++;
				}
				
			}
		}
	}
}

void finerRegister(Volume* baseVol, Volume* regVol,float4* finerGrid,float* finerResults,float medianScale)
{
	Matrix44* affineMat=new Matrix44();

	Matrix44 translationMat,scaleMat,tmpMat;

	Matrix44 rotationXYZ,rotationYZ,rotationZ;

	for(int z=0;z<finerAngleCount;z++)
	{
		rotationZ.setRotateAxisAngle(0,0,1.0,z*FLIRT_FINER_ANGLE_STEP + FLIRT_ANGLE_MIN);

		for(int y=0;y<finerAngleCount;y++)
		{
			tmpMat.setRotateAxisAngle(0,1.0,0,y*FLIRT_FINER_ANGLE_STEP + FLIRT_ANGLE_MIN);
			rotationYZ=rotationZ*tmpMat;

			for(int x=0;x<finerAngleCount;x++)
			{
				tmpMat.setRotateAxisAngle(1.0,0,0,x*FLIRT_FINER_ANGLE_STEP + FLIRT_ANGLE_MIN);
				
				rotationXYZ=rotationYZ*tmpMat;

				int index=z*finerAngleCount*finerAngleCount + y*finerAngleCount + x;

				float scale;
				float value=0;

				float tmpValue;

				translationMat.setTranslate(finerGrid[index].x,finerGrid[index].y,finerGrid[index].z,baseVol->dim.x);

				while(true)
				{
					scaleMat.setScale(finerGrid[index].w,finerGrid[index].w,finerGrid[index].w);
					tmpMat=rotationXYZ*scaleMat;
					*affineMat=translationMat*tmpMat;
					tmpValue = calcCostFunction(baseVol,regVol,affineMat);

					if(tmpValue < value)
					{
						scale = finerGrid[index].w;
						value = tmpValue;
					}

					if(finerGrid[index].w < medianScale) 
					{
						finerGrid[index].w += FLIRT_FINER_ANGLE_STEP;
						if(finerGrid[index].w > medianScale) break;
					}
					else 
					{
						finerGrid[index].w -= FLIRT_FINER_ANGLE_STEP;
						if(finerGrid[index].w < medianScale) break;
					}
					
				}//while

				finerGrid[index].w=scale;
				finerResults[index]=value;

			}//for x
		}//for y
	}//for z
	
}

void buildFinerGrid(float4* coarseGrid,float4* finerGrid)
{
	
	int3 neighbors;
	for(int z=0;z<finerAngleCount;z++)
	{
		for(int y=0;y<finerAngleCount;y++)
		{
			for(int x=0;x<finerAngleCount;x++)
			{
				if(x % 2 ==0)
					neighbors.x=1;
				else 
					neighbors.x=2;

				if(y % 2 ==0)
					neighbors.y=1;
				else 
					neighbors.y=2;

				if(z % 2 ==0)
					neighbors.z=1;
				else 
					neighbors.z=2;

				float4 cell={0,0,0,0};

				for(int k=0;k<neighbors.z;k++)
					for(int j=0;j<neighbors.y;j++)
						for(int i=0;i<neighbors.x;i++)
						{
							int coarseIndex= (z/2 + k)*coarseAngleCount*coarseAngleCount + (y/2 + j)*coarseAngleCount + (x/2 + i);
							
							cell.x += coarseGrid[coarseIndex].x;
							cell.y += coarseGrid[coarseIndex].y;
							cell.z += coarseGrid[coarseIndex].z;
							cell.w += coarseGrid[coarseIndex].w;
						}

				cell.x /= (float)(neighbors.x * neighbors.y * neighbors.z);
				cell.y /= (float)(neighbors.x * neighbors.y * neighbors.z);
				cell.z /= (float)(neighbors.x * neighbors.y * neighbors.z);
				cell.w /= (float)(neighbors.x * neighbors.y * neighbors.z);

				finerGrid[z*finerAngleCount*finerAngleCount + y*finerAngleCount + x]=cell;
			}
		}
	}
}

void coarseRegister(Volume* baseVol, Volume* regVol,float4* coarseGrid,float& medianScale)
{
	


	Matrix44 matXYZ,matYZ,matZ,tmpMat;

	for(int z=0;z<coarseAngleCount;z++)
	{
		matZ.setRotateAxisAngle(0,0,1.0,z*FLIRT_COARSE_ANGLE_STEP + FLIRT_ANGLE_MIN);

		for(int y=0;y<coarseAngleCount;y++)
		{
			tmpMat.setRotateAxisAngle(0,1.0,0,y*FLIRT_COARSE_ANGLE_STEP + FLIRT_ANGLE_MIN);

			matYZ=matZ*tmpMat;

			for(int x=0;x<coarseAngleCount;x++)
			{
				tmpMat.setRotateAxisAngle(1.0,0,0,x*FLIRT_COARSE_ANGLE_STEP + FLIRT_ANGLE_MIN);

				matXYZ=matYZ*tmpMat;

				//printf("%f,%f,%f:",(x*FLIRT_COARSE_ANGLE_STEP + FLIRT_ANGLE_MIN),(y*FLIRT_COARSE_ANGLE_STEP + FLIRT_ANGLE_MIN),(z*FLIRT_COARSE_ANGLE_STEP + FLIRT_ANGLE_MIN));

				coarseGrid[z*coarseAngleCount*coarseAngleCount + y*coarseAngleCount + x]=localOptimize4Dof(baseVol,regVol,&matXYZ);

				medianScale+=coarseGrid[z*coarseAngleCount*coarseAngleCount + y*coarseAngleCount + x].w;

			}
		}
	}

	medianScale /= (float)(coarseAngleCount*coarseAngleCount*coarseAngleCount);

}

float4 localOptimize4Dof(Volume* baseVol, Volume* regVol, Matrix44* rotationMat)
{
	float4 coarseGrid={0,0,0,1.0};

	int pass=0;

	float v=0;

	while( pass < FLIRT_SEARCH_PASS )
	{
		pass++;

		v=goldenSearchTranslation(baseVol,regVol,rotationMat,coarseGrid,'x');
		v=goldenSearchTranslation(baseVol,regVol,rotationMat,coarseGrid,'y');
		v=goldenSearchTranslation(baseVol,regVol,rotationMat,coarseGrid,'z');
		v=goldenSearchScale(baseVol,regVol,rotationMat,coarseGrid);
	}

	////cout<<v<<endl;

	return coarseGrid;
}

float goldenSearchTranslation(Volume* baseVol,Volume* regVol,Matrix44* rotationMat,float4 &coarseGrid,char target,float minStep)
{
	float *p;

	if(target == 'x')
		p=&(coarseGrid.x);
	else if(target == 'y')
		p=&(coarseGrid.y);
	else if(target == 'z')
		p=&(coarseGrid.z);
	else return 0;

	float points[4];

	if(target == 'x' || target == 'y')
	{
		float tmpPoints[4]=
		{
			FLIRT_TRANS_XY_MIN,
			*p,
			FLIRT_TRANS_XY_MAX,
			0
		};
		memcpy(points,tmpPoints,sizeof(float)*4);
	}
	else if(target == 'z')
	{
		float tmpPoints[4]=
		{
			FLIRT_TRANS_Z_MIN,
			*p,
			FLIRT_TRANS_Z_MAX,
			0
		};
		memcpy(points,tmpPoints,sizeof(float)*4);
	}

	float values[4]={0,0,0,0};

	Matrix44* affineMat=new Matrix44();

	Matrix44 translationMat,scaleMat,tmpMat;
	
	scaleMat.setScale(coarseGrid.w,coarseGrid.w,coarseGrid.w);

	translationMat.setTranslate(coarseGrid.x,coarseGrid.y,coarseGrid.z,baseVol->dim.x);
	tmpMat=(*rotationMat)*scaleMat;
	*affineMat=translationMat*tmpMat;
	values[1]=calcCostFunction(baseVol,regVol,affineMat);


	float d1,d2;

	int count=0;

	while(true)
	{
		d1=points[1]-points[0];
		d2=points[2]-points[1];

		count++;

		//find the logner section, and calculate the 4th point from current 3 points.

		//d2 is longer
		if(d1 < d2)
		{
			points[3]=points[2];
			values[3]=values[2];

			//insert a point into p0...p1...<<p>>...p2
			points[2]=points[1]+d2 * 0.61803f;
			*p=points[2];
			translationMat.setTranslate(coarseGrid.x,coarseGrid.y,coarseGrid.z,baseVol->dim.x);
			tmpMat=(*rotationMat)*scaleMat;
			*affineMat=translationMat*tmpMat;
			values[2]=calcCostFunction(baseVol,regVol,affineMat);
		}
		else 
		{
			points[3]=points[2];
			points[2]=points[1];
			values[3]=values[2];
			values[2]=values[1];

			//insert a point into p0...<<p>>...p1..p2
			points[1]=points[0]+ d1 * 0.38197f;

			*p=points[1];
			translationMat.setTranslate(coarseGrid.x,coarseGrid.y,coarseGrid.z,baseVol->dim.x);
			tmpMat=(*rotationMat)*scaleMat;
			*affineMat=translationMat*tmpMat;
			values[1]=calcCostFunction(baseVol,regVol,affineMat);
		}

		//form new 3 points group.

		//The latest 3 points come to be the 3 poits group, shift to left <---------
		if(values[1] > values[2])
		{
			points[0]=points[1];
			points[1]=points[2];
			points[2]=points[3];

			values[0]=values[1];
			values[1]=values[2];
			values[2]=values[3];
		}
		
		

		if(d1< minStep || d2< minStep)
			break;
	}
	
	*p=points[1];

	delete affineMat;
	
	return values[1];
}

float goldenSearchScale(Volume* baseVol,Volume* regVol,Matrix44* rotationMat,float4 &coarseGrid,float minStep)
{
	float points[4]=
	{
		FLIRT_SCALE_MIN,
		coarseGrid.w,
		FLIRT_SCALE_MAX,
		0
	};
	float values[4]={0,0,0,0};

	Matrix44* affineMat=new Matrix44();

	Matrix44 translationMat,scaleMat,tmpMat;
	
	translationMat.setTranslate(coarseGrid.x,coarseGrid.y,coarseGrid.z,baseVol->dim.x);
	

	//coarseGrid.w=points[1];
	scaleMat.setScale(coarseGrid.w,coarseGrid.w,coarseGrid.w);
	tmpMat=(*rotationMat)*scaleMat;
	*affineMat=translationMat*tmpMat;
	values[1]=calcCostFunction(baseVol,regVol,affineMat);

	float d1,d2;

	while(true)
	{
		d1=points[1]-points[0];
		d2=points[2]-points[1];

		//find the logner section, and calculate the 4th point from current 3 points.

		if(d1 < d2)
		{
			points[3]=points[2];
			values[3]=values[2];

			points[2]=points[1]+d2*0.61803f;
			coarseGrid.w=points[2];
			scaleMat.setScale(coarseGrid.w,coarseGrid.w,coarseGrid.w);
			tmpMat=(*rotationMat)*scaleMat;
			*affineMat=translationMat*tmpMat;
			values[2]=calcCostFunction(baseVol,regVol,affineMat);
		}
		else 
		{
			points[3]=points[2];
			points[2]=points[1];
			values[3]=values[2];
			values[2]=values[1];

			//insert a point into p0...<<p>>...p1..p2
			points[1]=points[0]+d1 * 0.38197f;
			coarseGrid.w=points[1];
			scaleMat.setScale(coarseGrid.w,coarseGrid.w,coarseGrid.w);
			tmpMat=(*rotationMat)*scaleMat;
			*affineMat=translationMat*tmpMat;
			values[1]=calcCostFunction(baseVol,regVol,affineMat);
		}

		//form  new 3 points group.

		//The latest 3 points come to be the 3 poits group.
		if(values[1] > values[2])
		{
			points[0]=points[1];
			points[1]=points[2];
			points[2]=points[3];

			values[0]=values[1];
			values[1]=values[2];
			values[2]=values[3];
		}
		

		if(d1< minStep || d2<minStep)
			break;

	}
	
	coarseGrid.w=points[1];

	delete affineMat;
	

	return values[1];
}

float calcCostFunction(Volume* baseVol,Volume* regVol,Matrix44* affineMat)
{
	globalCount++;

	float baseEntropy=0,regEntropy=0,jointEntropy=0;

	calcEntropyGPU(baseVol,regVol,affineMat,baseEntropy,regEntropy,jointEntropy);

	/*printf("%f ",(jointEntropy-regEntropy-baseEntropy));*/

	return (jointEntropy-regEntropy-baseEntropy);
}

