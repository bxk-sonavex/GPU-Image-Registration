/*
 * Authored by: Chen, Shifu
 * 
 * Email: chen@gmtk.org / sf.chen@ieee.org
 *
 * The code is distributed under BSD license, you are allowed to use, modify or sell this code, but a statement is required if you used this code any where.
 * 
 */
#include "volume.h"
#include "costFunctions.h"
#include "time.h"
#include "costFunctionsGPU.cuh"
#include "search.h"

int main()
{
	char* regFileName="data/MRT2_153359_LXZ_512x512x23_0.46875x0.468746x6_2_8bit.raw";
	char* baseFileName="data/MR_153359_LXZ_128x128x19_1.79688x1.79688x6.5_3ep_b0_1000_8bit.raw";
	

	Volume* baseVol = loadVolumeRAW(baseFileName,FLIRT_UINT8);
	Volume* regVol = loadVolumeRAW(regFileName,FLIRT_UINT8);

	float3 targetDim = {1.0,1.0,1.0};
	Volume* baseVol1 = isotropicResample(baseVol,targetDim);
	Volume* regVol1 = isotropicResample(regVol,targetDim);
	
	char name[200];

	sprintf(name,"results/Reg_inporlated_LXZ_%dx%dx%d_%fx%fx%f_3ep_b0_1000_8bit.raw",regVol1->size.x,regVol1->size.y,regVol1->size.z,regVol1->dim.x,regVol1->dim.y,regVol1->dim.z);
	writeVolumeToRaw(name,regVol1);
	sprintf(name,"results/Base_inporlated_LXZ_%dx%dx%d_%fx%fx%f_3ep_b0_1000_8bit.raw",baseVol1->size.x,baseVol1->size.y,baseVol1->size.z,baseVol1->dim.x,baseVol1->dim.y,baseVol1->dim.z);
	writeVolumeToRaw(name,baseVol1);

	Transform12 inv_t12 = registerVolume(regVol1,baseVol1);
	Volume* registeredVolume = transformVolume(baseVol1,regVol1,inv_t12);
	sprintf(name,"results/Registered_LXZ%dx%dx%d_%fx%fx%f_3ep_b0_1000_8bit.raw",registeredVolume->size.x,registeredVolume->size.y,registeredVolume->size.z,registeredVolume->dim.x,registeredVolume->dim.y,registeredVolume->dim.z);
	writeVolumeToRaw(name,registeredVolume);

	releaseVolume(baseVol);
	releaseVolume(baseVol1);
	releaseVolume(regVol1);

	return 0;
}