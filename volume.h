/*
 * Authored by: Chen, Shifu
 * 
 * Email: chen@gmtk.org / sf.chen@ieee.org
 *
 * The code is distributed under BSD license, you are allowed to use, modify or sell this code, but a statement is required if you used this code any where.
 * 
 */
#ifndef FLIRT_VOLUME_H
#define FLIRT_VOLUME_H

#include "common.h"
#include "stdio.h"
#include "stdlib.h"
#include "memory.h"
#include "string.h"

struct Volume
{
	//Dimension of this Volume
	float3 dim;

	//Size of this Volume
	int3 size;

	//Data pointer
	float* data;

	//min value
	float minValue;

	//max value
	float maxValue;

	//center of mass
	float3 center;
	
};

inline int getIndex(int3 size,int i,int j,int k);

//Create a Volume with fixed size
Volume* createVolume(int3 size);

//Release a Volume and its data
void releaseVolume(Volume* vol);

//Load a Volume with a filename, a fixed size, the dataType is set to uint8 by default, change to FLIRT_UINT16 when the data is in uint16
Volume* loadVolumeRAW(char* filename, int dataType = FLIRT_UINT8 );

void writeVolumeToRaw(char* filename, Volume* vol, int dataType = FLIRT_UINT8 );

//Find the dimension from filename
void analyzeDimensionByFilename(char* filename, Volume* vol);

//Find the dimension from filename
int3 analyzeSizeByFilename(char* filename);

Volume* isotropicResample(Volume* baseVol,float3 targetDim);

//Subsampling by 2, so the volume sequence is 1.0mm --->  2.0 mm  --->   4.0mm  --->  8.0mm
Volume* halfSample(Volume* baseVol,bool smoothed=true);

bool insideVolume(int3 size, float3 coord);

//Get the interoplated value
float interpolate(Volume* vol, float3 point);

//calculate the center of mass.

void calcCenter(Volume* vol);

//Volume* buildSimulatedData(Volume* baseVol);

Volume* transformVolume(Volume* baseVol,Volume* regVol,Transform12 t12);


#endif //#ifndef FLIRT_VOLUME_H