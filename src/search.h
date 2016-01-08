#ifndef FLIRT_SEARCH_H
#define FLIRT_SEARCH_H

#include "volume.h"



Transform12 registerVolume(Volume* baseVol, Volume* regVol); 

void coarseRegister(Volume* baseVol, Volume* regVol,float4* coarseResult,float& medianScale); 

void finerRegister(Volume* baseVol, Volume* regVol,float4* finerGrid,float* finerResults,float medianScale); 

float4 localOptimize4Dof(Volume* baseVol, Volume* regVol, Matrix44* rotationMat); 

float goldenSearchTranslation(Volume* baseVol,Volume* regVol,Matrix44* rotationMat,float4 &translationAndScale,char target,float step= FLIRT_TRANS_XY_STEP);

float goldenSearchScale(Volume* baseVol,Volume* regVol,Matrix44* rotationMat,float4 &translationAndScale,float step= FLIRT_SCALE_STEP);

float calcCostFunction(Volume* baseVol,Volume* regVol,Matrix44* affineMat);

void buildFinerGrid(float4* coarseGrid,float4* finerGrid);

void findMinima(float* finerResults,int* minima, int& count);

void sortMinima(float4* finerGrid,float* finerResults,int* minima, int count,Transform7* globalMinima);

void perturbationStage(Volume* baseVol4,Volume* regVol4,Transform7* globalMinima,Transform7* higherResStart);

void perturbation(Volume* baseVol,Volume* regVol,Transform7 globalMinima,Transform7* grid4mm, float* result4mm );

Transform12 optimize2mm(Volume* baseVol2,Volume* regVol2,Transform7 higherResStart);

void perturbationZ(Volume* baseVol,Volume* regVol,Transform7& higherResStart);

float singlePerturbationZ(Volume* baseVol,Volume* regVol,Transform7& higherResStart);

float goldenSearch7dof(Volume* baseVol,Volume* regVol,Transform7 &coarseGrid,char* target,float minStep);

float goldenSearch9dof(Volume* baseVol,Volume* regVol,Transform12 &rotTransScale,char* target,float minStep);
float goldenSearch12dof(Volume* baseVol,Volume* regVol,Transform12 &rotTransScale,char* target,float minStep);

float optimize2mm9dof(Volume* baseVol,Volume* regVol,Transform12& t12);

float optimize2mm12dof(Volume* baseVol,Volume* regVol,Transform12& t12);

float optimize1mm12dof(Volume* baseVol,Volume* regVol,Transform12& t12);

float optimize1mm(Volume* baseVol,Volume* regVol,Transform12& t12);

#endif