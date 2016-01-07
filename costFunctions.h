/*
 * Authored by: Chen, Shifu
 * 
 * Email: chen@gmtk.org / sf.chen@ieee.org
 *
 * The code is distributed under BSD license, you are allowed to use, modify or sell this code, but a statement is required if you used this code any where.
 * 
 */
#ifndef FLIRT_COST_FUNCTIONS_H
#define FLIRT_COST_FUNCTIONS_H

#include "volume.h"
#include "matrix44.h"

void calcEntropy(Volume* baseVol, Volume* regVol, Matrix44* dAffine, float &baseEntropy, float& regEntropy, float& jointEntropy);


#endif

