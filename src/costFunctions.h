#ifndef FLIRT_COST_FUNCTIONS_H
#define FLIRT_COST_FUNCTIONS_H

#include "volume.h"
#include "matrix44.h"

void calcEntropy(Volume* baseVol, Volume* regVol, Matrix44* dAffine, float &baseEntropy, float& regEntropy, float& jointEntropy);


#endif

