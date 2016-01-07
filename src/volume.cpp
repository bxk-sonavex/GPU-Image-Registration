/*
 * Authored by: Chen, Shifu
 * 
 * Email: chen@gmtk.org / sf.chen@ieee.org
 *
 * The code is distributed under BSD license, you are allowed to use, modify or sell this code, but a statement is required if you used this code any where.
 * 
 */
#include "volume.h"

//Get the linear index by the size of Volume and the coordinate.
inline int getIndex(int3 size,int x,int y,int z)
{
	if(x>=size.x)x=size.x-1;
	if(y>=size.y)y=size.y-1;
	if(z>=size.z)z=size.z-1;

	if(x<0)x=0;
	if(y<0)y=0;
	if(z<0)z=0;

	return z*size.x*size.y + y*size.x + x;
}


//Create a Volume with a fixed size
Volume* createVolume(int3 size)
{
	Volume* vol=new Volume();
	vol->size=size;
	vol->data=new float[size.x*size.y*size.z];
	memset(vol->data,0,size.x*size.y*size.z*sizeof(float));

	return vol;
}

//Release a fixed size
void releaseVolume(Volume* vol)
{
	delete vol->data;
	delete vol;
}

//Normalize a volume, if the orginal data is 16bit based
void normalizeVolume(Volume* vol)
{
	float maxValue = vol -> maxValue;
	int i,j,k,index;
	if(maxValue > 0)
	{
		for(k=0;k<vol->size.z;k++)
			for(j=0;j<vol->size.z;j++)
				for(i=0;i<vol->size.z;i++)
				{
					index = getIndex(vol->size,i,j,k);

					vol->data[index] = vol->data[index] / maxValue * 255.0f;
				}
	}
}

//Load a Volume with a filename, a fixed size, the dataType is set to uint8 by default, change to FLIRT_UINT16 when the data is in uint16
Volume* loadVolumeRAW(char* filename, int dataType )
{
	int3 size=analyzeSizeByFilename(filename);

	//Create a Volume firstly
	Volume* vol=createVolume(size);

	//Open the file and read the data;
	FILE* fp;
	fp=fopen(filename,"rb");

	if(!fp)
	{
		//cerr<<"Error:file not exists."<<endl;
		exit(1);
	}

	int i,j,k,index;

	//read the data values
	for(k=0;k<size.z;k++)
		for(j=0;j<size.y;j++)
			for(i=0;i<size.x;i++)
			{
				index = getIndex(size,i,j,k);
				if(dataType == FLIRT_UINT8)
				{
					unsigned char value;
					fread(&value,1,1,fp);
					vol->data[index]=(float)value;
				}
				else if(dataType == FLIRT_UINT16)
				{
					unsigned short value;
					fread(&value,2,1,fp);
					vol->data[index]=(float)value;
				}
				else
				{
					//cerr<<"Error:datatype should be uint8 or uint16."<<endl;
					exit(1);
				}

				if(index == 0)
				{
					vol->maxValue = vol->data[index];
					vol->minValue = vol->data[index];
				}
				else 
				{
					if(vol->data[index] > vol->maxValue) vol->maxValue =vol->data[index];
					if(vol->data[index] < vol->minValue) vol->minValue =vol->data[index];
				}

			}
			
	if(dataType == FLIRT_UINT16)
		normalizeVolume(vol);

	fclose(fp);

	vol->size=size;

	analyzeDimensionByFilename(filename,vol);

	calcCenter(vol);

	return vol;
		
}

void writeVolumeToRaw(char* filename, Volume* vol, int dataType )
{
	FILE* fp=fopen(filename,"wb");

	for(int k=0;k<vol->size.z;k++)
		for(int j=0;j<vol->size.y;j++)
			for(int i=0;i<vol->size.x;i++)
			{
				unsigned char value=(unsigned char)vol->data[getIndex(vol->size,i,j,k)];
				fwrite(&value,sizeof(unsigned char),1,fp);
			}

	fclose(fp);


}

int3 analyzeSizeByFilename(char* filename)
{
	int num=(int)strlen(filename);

	int3 size={0,0,0};
	int xPlace[10];

	int i=0,j=0,k = 0;
	for ( i = 1; i < num-1;  i++)
	{
		if (filename[i] == 'x' && (filename[i-1]>='0' && filename[i-1]<='9') && (filename[i+1]>='0' && filename[i+1]<='9'))
		{
			xPlace[k] = i;
			k++;
		}
	}

	char str[3][10];
	int pos;

	//Width
	pos=xPlace[0]-1;
	while(pos>0 && filename[pos-1]>='0' && filename[pos-1]<='9')
		pos--;
	strncpy(str[0], filename+pos, xPlace[0]-pos);

	//Length
	pos=xPlace[1]-1;
	while(pos>0 && filename[pos-1]>='0' && filename[pos-1]<='9')
		pos--;
	strncpy(str[1], filename+pos, xPlace[1]-pos);

	//Depth
	pos=xPlace[1]+1;
	while(pos<num-1 && filename[pos+1]>='0' && filename[pos+1]<='9')
		pos++;
	strncpy(str[2], filename+xPlace[1]+1, pos-xPlace[1]);

	size.x = atoi(str[0]);
	size.y = atoi(str[1]);
	size.z = atoi(str[2]);

	return size;
}

//Find the dimension from filename
void analyzeDimensionByFilename(char* filename, Volume* vol)
{
	int num=(int)strlen(filename);

	int xPlace[10];

	int i=0,j=0,k = 0;
	for ( i = 1; i < num-1;  i++)
	{
		if (filename[i] == 'x' && (filename[i-1]>='0' && filename[i-1]<='9') && (filename[i+1]>='0' && filename[i+1]<='9'))
		{
			xPlace[k] = i;
			k++;
		}
	}

	if(k<4)
	{
		printf("ERROR: Can't extract dimension from filename.");
	}

	char str[3][10];
	int pos;

	//Width
	pos=xPlace[2]-1;
	while(pos>0 && ((filename[pos-1]>='0' && filename[pos-1]<='9') || filename[pos-1]=='.'))
		pos--;
	strncpy(str[0], filename+pos, xPlace[2]-pos);

	//Length
	pos=xPlace[3]-1;
	while(pos>0 && ((filename[pos-1]>='0' && filename[pos-1]<='9') || filename[pos-1]=='.'))
		pos--;
	strncpy(str[1], filename+pos, xPlace[3]-pos);

	//Depth
	pos=xPlace[3]+1;
	while(pos<num-1 && ((filename[pos+1]>='0' && filename[pos+1]<='9') || filename[pos+1]=='.'))
		pos++;
	strncpy(str[2], filename+xPlace[3]+1, pos-xPlace[3]);

	vol->dim.x = (float)atof(str[0]);
	vol->dim.y = (float)atof(str[1]);
	vol->dim.z = (float)atof(str[2]);
}

//Trilinear interpolation kernel,wraped from FLIRT
inline float trilinearInterpolation(float v000, float v001, float v010, 
				   float v011, float v100, float v101, 
				   float v110, float v111, 
				   float dx, float dy, float dz)

{
      float temp1, temp2, temp3, temp4, temp5, temp6;
      temp1 = (v100 - v000)*dx + v000;
      temp2 = (v101 - v001)*dx + v001;
      temp3 = (v110 - v010)*dx + v010;
      temp4 = (v111 - v011)*dx + v011;
      // second order terms
      temp5 = (temp3 - temp1)*dy + temp1;
      temp6 = (temp4 - temp2)*dy + temp2;
      // final third order term
	  float result = (temp6 - temp5)*dz + temp5;
      return (temp6 - temp5)*dz + temp5;
}

//Get the interoplated value,wraped from FLIRT
float interpolate(Volume* vol, float3 point)
{
	int ix, iy, iz;  
	ix=(int) floor(point.x);
	iy=(int) floor(point.y);
	iz=(int) floor(point.z);

	//left-top-front
	float dx=point.x-ix, dy=point.y-iy, dz=point.z-iz;

	float v000=0, v001=0, v010=0, v011=0, v100=0, v101=0, v110=0, v111=0;

	v000 = vol->data[getIndex(vol->size,ix,iy,iz)];
	v001 = vol->data[getIndex(vol->size,ix,iy,iz+1)];
	v010 = vol->data[getIndex(vol->size,ix,iy+1,iz)];
	v011 = vol->data[getIndex(vol->size,ix,iy+1,iz+1)];
	v100 = vol->data[getIndex(vol->size,ix+1,iy,iz)];
	v101 = vol->data[getIndex(vol->size,ix+1,iy,iz+1)];
	v110 = vol->data[getIndex(vol->size,ix+1,iy+1,iz)];
	v111 = vol->data[getIndex(vol->size,ix+1,iy+1,iz+1)];

	//Use trilinear interpolation
	return trilinearInterpolation(v000,v001,v010,v011,v100,v101,v110,v111,dx,dy,dz);
}



//The function for isotropic resampling
Volume* isotropicResample(Volume* baseVol,float3 targetDim)
{
	float3 step;

	//Calculate the step
	step.x = targetDim.x / baseVol->dim.x;
	step.y = targetDim.y / baseVol->dim.y;
	step.z = targetDim.z / baseVol->dim.z;

	int3 newSize;
	
	//Calculate the new size
	newSize.x = (int)(((float)baseVol->size.x /*+ 1.0*/ )/ step.x);
	newSize.y = (int)(((float)baseVol->size.y /*+ 1.0*/ )/ step.y);
	newSize.z = (int)(((float)baseVol->size.z /*+ 1.0*/ )/ step.z);

	Volume* newVol = createVolume(newSize);

	newVol->dim=targetDim;
	newVol->center=baseVol->center;

	float3 point;
	int x,y,z;

	int index;

	for (point.z=0.0, z=0; z<newSize.z; z++, point.z+=step.z) {
      for (point.y=0.0, y=0; y<newSize.y; y++, point.y+=step.y) {
		for(point.x=0.0, x=0; x<newSize.x; x++, point.x+=step.x) {
			index=getIndex(newSize,x,y,z);
			newVol->data[index] = interpolate(baseVol,point);
		}
      }
    }

	calcCenter(newVol);

	return newVol;
}

//Spacing: 1.0mm --->  2.0 mm  --->   4.0mm  --->  8.0mm
Volume* halfSample(Volume* baseVol, bool smoothed)
{
	int3 newSize;
	newSize.x = (baseVol->size.x + 1) /2 ;
	newSize.y = (baseVol->size.y + 1) /2 ;
	newSize.z = (baseVol->size.z + 1) /2 ;

	float3 newDim={baseVol->dim.x*2.0f, baseVol->dim.y*2.0f, baseVol->dim.z*2.0f};

	Volume* vol= createVolume(newSize);

	vol->dim = newDim;
	vol->center=baseVol->center;

	int x,y,z,bx,by,bz;

	int faceIndex[6],edgeIndex[12],vertexIndex[8];

	float faceValue,edgeValue,vertexValue,baseValue;

	int i,baseIndex, newIndex;

	if(smoothed)
	{
		for ( z=0, bz=0; z< vol->size.z; z++, bz+=2)
			for ( y=0, by=0; y< vol->size.y; y++, by+=2)
				for ( x=0, bx=0; x< vol->size.x; x++, bx+=2)
				{
					//Must set to 0 each iteration.
					faceValue=0.0,edgeValue=0.0,vertexValue=0.0;

					baseIndex = getIndex(baseVol->size,bx,by,bz);
					baseValue = baseVol->data[baseIndex];

					//6 faces
					faceIndex[0]=getIndex(baseVol->size,bx+1,by,bz);
					faceIndex[1]=getIndex(baseVol->size,bx-1,by,bz);
					faceIndex[2]=getIndex(baseVol->size,bx,by+1,bz);
					faceIndex[3]=getIndex(baseVol->size,bx,by-1,bz);
					faceIndex[4]=getIndex(baseVol->size,bx,by,bz+1);
					faceIndex[5]=getIndex(baseVol->size,bx,by,bz-1);

					for(i=0;i<6;i++)
						faceValue += baseVol->data[faceIndex[i]];


					//12 edges
					edgeIndex[0]=getIndex(baseVol->size,bx+1,by+1,bz);
					edgeIndex[1]=getIndex(baseVol->size,bx+1,by-1,bz);
					edgeIndex[2]=getIndex(baseVol->size,bx-1,by+1,bz);
					edgeIndex[3]=getIndex(baseVol->size,bx-1,by-1,bz);
					edgeIndex[4]=getIndex(baseVol->size,bx+1,by,bz+1);
					edgeIndex[5]=getIndex(baseVol->size,bx+1,by,bz-1);
					edgeIndex[6]=getIndex(baseVol->size,bx-1,by,bz+1);
					edgeIndex[7]=getIndex(baseVol->size,bx-1,by,bz-1);
					edgeIndex[8]=getIndex(baseVol->size,bx,by+1,bz+1);
					edgeIndex[9]=getIndex(baseVol->size,bx,by+1,bz-1);
					edgeIndex[10]=getIndex(baseVol->size,bx,by-1,bz+1);
					edgeIndex[11]=getIndex(baseVol->size,bx,by-1,bz-1);

					for(i=0;i<12;i++)
						edgeValue += baseVol->data[edgeIndex[i]];


					//8 vertexes
					vertexIndex[0]=getIndex(baseVol->size,bx+1,by+1,bz+1);
					vertexIndex[1]=getIndex(baseVol->size,bx+1,by+1,bz-1);
					vertexIndex[2]=getIndex(baseVol->size,bx+1,by-1,bz+1);
					vertexIndex[3]=getIndex(baseVol->size,bx+1,by-1,bz-1);
					vertexIndex[4]=getIndex(baseVol->size,bx-1,by+1,bz+1);
					vertexIndex[5]=getIndex(baseVol->size,bx-1,by+1,bz-1);
					vertexIndex[6]=getIndex(baseVol->size,bx-1,by-1,bz+1);
					vertexIndex[7]=getIndex(baseVol->size,bx-1,by-1,bz-1);

					for(i=0;i<8;i++)
						vertexValue += baseVol->data[vertexIndex[i]];


					newIndex = getIndex(vol->size,x,y,z);

					//Smoothed kernel
					vol->data[newIndex] = (float)(0.125 * baseValue + 0.0625 * faceValue + 0.0312 * edgeValue + (1 - 0.125 - 0.0625*6 - 0.0312 * 12)/8.0 * vertexValue);
					
				}

	}
	else // no smoothed
	{
		for ( z=0, bz=0; z< vol->size.z; z++, bz+=2)
			for ( y=0, by=0; y< vol->size.y; y++, by+=2)
				for ( x=0, bx=0; x< vol->size.x; x++, bx+=2)
				{
					vertexIndex[0]=getIndex(baseVol->size,bx,by,bz);
					vertexIndex[1]=getIndex(baseVol->size,bx+1,by,bz);
					vertexIndex[2]=getIndex(baseVol->size,bx,by+1,bz);
					vertexIndex[3]=getIndex(baseVol->size,bx,by,bz+1);
					vertexIndex[4]=getIndex(baseVol->size,bx+1,by+1,bz);
					vertexIndex[5]=getIndex(baseVol->size,bx+1,by,bz+1);
					vertexIndex[6]=getIndex(baseVol->size,bx,by+1,bz+1);
					vertexIndex[7]=getIndex(baseVol->size,bx+1,by+1,bz+1);

					vertexValue=0.0;

					for(i=0;i<8;i++)
						vertexValue += baseVol->data[vertexIndex[i]];

					newIndex = getIndex(vol->size,x,y,z);

					vol->data[newIndex] = vertexValue / 8.0f;
				}
	}

	return vol;
}

bool insideVolume(int3 size, float3 coord)
{
	if (coord.x <0 || coord.y <0 || coord.z <0 )
		return false;

	if (coord.x >size.x || coord.y > size.y || coord.z > size.z )
		return false;

	return true;
}

void calcCenter(Volume* vol)
{
	float3 center = {0,0,0};

	int x,y,z,index;

	int count=0;

	for(z=0;z<vol->size.z;z++)
		for(y=0;y<vol->size.y;y++)
			for(x=0;x<vol->size.x;x++)
			{
				index = getIndex(vol->size,x,y,z);

				if(vol->data[index] !=0)
				{
					center.x += x * vol->dim.x;
					center.y += y * vol->dim.y;
					center.z += z * vol->dim.z;
					count++;
				}

			}
	center.x /= (float)count;
	center.y /= (float)count;
	center.z /= (float)count;

	vol->center = center;
}

Volume* transformVolume(Volume* baseVol,Volume* regVol,Transform12 t12)
{
	Matrix44* affineMat=new Matrix44();

	Matrix44 rotationMat,translationMat,scaleMat,tmpMat,skewMat;

	rotationMat.setRotateAxisAngle(1,0,0,t12.rotation.x);
	tmpMat=rotationMat;

	rotationMat.setRotateAxisAngle(0,1,0,t12.rotation.y);
	tmpMat=rotationMat*tmpMat;

	rotationMat.setRotateAxisAngle(0,0,1,t12.rotation.z);
	rotationMat=rotationMat*tmpMat;

	scaleMat.setScale(t12.scale.x,t12.scale.y,t12.scale.z);
	translationMat.setTranslate(t12.translation.x,t12.translation.y,t12.translation.z,baseVol->dim.x);
	skewMat.setSkew(t12.skew.x,t12.skew.y,t12.skew.z);

	tmpMat=skewMat*scaleMat;
	tmpMat=rotationMat*tmpMat;
	*affineMat=translationMat*tmpMat;

	Matrix44 mat=*affineMat;
	mat.display();

	Volume* vol=createVolume(baseVol->size);
	vol->dim = baseVol->dim;
	vol->size=baseVol->size;

	for(int z=0;z<baseVol->size.z;z++)
	{
		for(int y=0;y<baseVol->size.y;y++)
		{
			for(int x=0;x<baseVol->size.x;x++)
			{
				float3 coord={(float)x,(float)y,(float)z};
				
				float value;

				coord.x -=baseVol->center.x;
				coord.y -=baseVol->center.y;
				coord.z -=baseVol->center.z;

				coord = mat*coord;

				//printf("%f,%f,%f\n", coord.x,coord.y,coord.z);

				coord.x +=regVol->center.x;
				coord.y +=regVol->center.y;
				coord.z +=regVol->center.z;

				value=interpolate(regVol,coord);

				vol->data[getIndex(vol->size,x,y,z)] = value;


			}
		}
	}

	calcCenter(vol);

	return vol;

}
