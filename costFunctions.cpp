/*
 * Authored by: Chen, Shifu
 * 
 * Email: chen@gmtk.org / sf.chen@ieee.org
 *
 * The code is distributed under BSD license, you are allowed to use, modify or sell this code, but a statement is required if you used this code any where.
 * 
 */
#include "costFunctions.h"


void addToJointHist(int hist[FLIRT_BIN_COUNT][FLIRT_BIN_COUNT],float valueX, float valueY,int binCount)
{
	const unsigned int binWidth=256/binCount;

	int offsetX,offsetY;

	if(valueX > 255.0)
		offsetX = binCount -1;
	else offsetX = (int)(valueX / binWidth) ;

	if(valueY > 255.0)
		offsetY = binCount -1;
	else offsetY = (int)(valueY / binWidth);

	hist[offsetX][offsetY]++;
}

void calcEntropy(Volume* baseVol, Volume* regVol, Matrix44* dAffine, float& baseEntropy, float& regEntropy, float& jointEntropy)
{
	int binCount;
	if(baseVol->dim.x ==8.0)
		binCount=FLIRT_BIN_COUNT_8MM;
	else if(baseVol->dim.x ==4.0)
		binCount=FLIRT_BIN_COUNT_4MM;
	else 
		binCount=FLIRT_BIN_COUNT;

	const unsigned int binWidth=256/binCount;
	int x,y,z;

	//Coordination in reg volume space
	float3 regCoord;

	//Coordination after  affine transform, means the coordination in the base volume space
	float3 baseCoord;

	//fetch the value
	float baseValue, regValue;

	int baseHist[FLIRT_BIN_COUNT], regHist[FLIRT_BIN_COUNT], jointHist[FLIRT_BIN_COUNT][FLIRT_BIN_COUNT];

	memset(baseHist,0,sizeof(int)*binCount);
	memset(regHist,0,sizeof(int)*binCount);

	int i,j;
	for(i=0;i<binCount;i++)
		for(j=0;j<binCount;j++)
			jointHist[i][j]=0;


	int regIndex;

	int totalJoint =0;

	//Calculate the regHist & jointHist
	for(z=0;z<regVol->size.z;z++)
		for(y=0;y<regVol->size.y;y++)
			for(x=0;x<regVol->size.x;x++)
			{
				regCoord.x = x - regVol->center.x / regVol->dim.x;
				regCoord.y = y - regVol->center.y / regVol->dim.y;
				regCoord.z = z - regVol->center.z / regVol->dim.z;


				regIndex = getIndex(regVol ->size,x,y,z);
				
				// Fetch the value at reg Volume
				regValue = regVol->data[regIndex];

				// transform
				baseCoord=(*dAffine) * regCoord ;
				baseCoord.x += baseVol->center.x / baseVol->dim.x;
				baseCoord.y += baseVol->center.y / baseVol->dim.y;
				baseCoord.z += baseVol->center.z / baseVol->dim.z;

				//Add to the joint hist
				if(insideVolume(baseVol->size,baseCoord))
				{
					totalJoint ++;

					//Fetch the base value with interpolation.
					baseValue = interpolate(baseVol,baseCoord);

					//Add to the Joint hist
					addToJointHist(jointHist,baseValue,regValue,binCount);

				}//if
			}//for x
			/*

	for( i=0;i<binCount * binCount ; i++)
		if(jointHist[i] > 0)
			printf("%d,%d,%d\n",i% binCount,i/binCount,jointHist[i% binCount][i/binCount]);*/

	//Calculate the maginal distribution

	// i:base j:reg

	for(i=0;i<binCount;i++)
	{
		for(j=0;j<binCount;j++)
		{
			baseHist[i] += jointHist[i][j];
			regHist[j]+= jointHist[i][j];
		}
	}

	baseEntropy=0;
	regEntropy=0;
	jointEntropy=0;

	float p;

	//Calculate the base entropy
	for(i=0;i< binCount;i++)
	{
		p = (float)baseHist[i] / (float)totalJoint ;

		if(p==0)
			continue;

		baseEntropy -= p * log(p);
	}

	//Calculate the reg entropy
	for(i=0;i< binCount;i++)
	{
		p = (float)regHist[i] / (float)totalJoint ;

		if(p==0)
			continue;

		regEntropy -= p * log(p);
	}

	//Calculate the joint entropy

	for(i=0;i<binCount;i++)
	{
		for(j=0;j<binCount;j++)
		{
			if(jointHist[i][j] == 0)
				continue;

			jointEntropy -= (float)jointHist[i][j] / (float)totalJoint  * (log((float)jointHist[i][j]) - log((float)totalJoint));
		}
	}



}