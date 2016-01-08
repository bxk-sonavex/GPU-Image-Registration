#ifndef FLIRT_COST_FUNCTIONS_CUH
#define FLIRT_COST_FUNCTIONS_CUH

#include "volume.h"
#include "matrix44.h"

void calcEntropyGPU(Volume* baseVol, Volume* regVol, Matrix44* dAffine, float &baseEntropy, float& regEntropy, float& jointEntropy);

void initVolumeData(Volume* baseVol, Volume* regVol);

void freeVolumeData(bool clearCpuBuffer=false);


#endif