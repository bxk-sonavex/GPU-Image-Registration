#include <stdio.h>
#include <stdlib.h>
#include <cutil.h>
#include <cstdlib>
#include <cstdio>

#include <GL/glew.h>

#if defined (__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <cuda_gl_interop.h>
#include <cutil_inline.h>
#include "bbsort.cuh"
#include "costFunctionsGPU.cuh"

float4* dRegSortedData;
float* dRegVolData;

unsigned int* dJointHist;
unsigned int* dRegHistCopy;
unsigned int* dBaseHist;
unsigned int* dRegHist;
float3* dEntropies;

float3* hEntropies;
unsigned int* hRegHist;
unsigned int* hBaseHist;
unsigned int* hJointHist;

unsigned int* hBlockPivots;
unsigned int* dBlockPivots;

cudaArray *dBaseArray;

float3 hBaseCenter,hRegCenter;

#include "costFunctionsGPU_kernel.cu"

void initRegFloat4(Volume* baseVol,int3 size)
{
	unsigned int binCount;
	if(baseVol->dim.x ==8.0)
		binCount=FLIRT_BIN_COUNT_8MM;
	else if(baseVol->dim.x ==4.0)
		binCount=FLIRT_BIN_COUNT_4MM;
	else 
		binCount=FLIRT_BIN_COUNT;

	const unsigned int binWidth=256/binCount;

	unsigned int blockSize = 512;

	unsigned int totalPoints = size.x * size.y * size.z ;

	unsigned int  blockCount;


	if(totalPoints%blockSize==0)
		blockCount=totalPoints/blockSize;
	else 
		blockCount=totalPoints/blockSize+1;

	buildRegFloat4Kernel<<<blockCount,blockSize>>>(dRegSortedData,dRegVolData,size,totalPoints);

	bbSort(dRegSortedData, totalPoints);

	float4* hRegSortedData = new float4[totalPoints];

	cudaMemcpy(hRegSortedData, dRegSortedData, sizeof(float4) * totalPoints, cudaMemcpyDeviceToHost);


	hBlockPivots[0]=0;

	unsigned int offset=0;
	
	for(unsigned int target=1;target<= binCount; target++)
	{
		while( (unsigned int)hRegSortedData[offset].w / binWidth < target)
		{
			if(offset >= totalPoints)
				break;
			else 
				offset++;
		}

		hBlockPivots[target] = offset;

	}

	cudaMemcpy(dBlockPivots,hBlockPivots , sizeof(unsigned int) * (binCount + 1), cudaMemcpyHostToDevice);

	delete hRegSortedData;


}

void initVolumeData(Volume* baseVol, Volume* regVol)
{

	int binCount;
	if(baseVol->dim.x ==8.0)
		binCount=FLIRT_BIN_COUNT_8MM;
	else if(baseVol->dim.x ==4.0)
		binCount=FLIRT_BIN_COUNT_4MM;
	else 
		binCount=FLIRT_BIN_COUNT;

	//const unsigned int binWidth=256/binCount;


	hEntropies = new float3[binCount];
	hRegHist = new unsigned int[binCount];
	hBaseHist=new unsigned int[binCount];
	hJointHist=new unsigned int[binCount*binCount];
	hBlockPivots=new unsigned int[binCount + 1];


	float3 baseCenter={baseVol->center.x / baseVol->dim.x, baseVol->center.y / baseVol->dim.y, baseVol->center.z / baseVol->dim.z};
	float3 regCenter={regVol->center.x / regVol->dim.x, regVol->center.y / regVol->dim.y, regVol->center.z / regVol->dim.z};

	cudaMemcpyToSymbol( cBaseCenter, &baseCenter, sizeof(float3)) ;

	cutilCheckMsg("cudaMemcpyToSymbol  baseVol failed");

	cudaMemcpyToSymbol( cRegCenter, &regCenter, sizeof(float3)) ;

	cutilCheckMsg("cudaMemcpyToSymbol  regVol failed");

	cudaExtent baseVolumeSize = make_cudaExtent(baseVol -> size.x, baseVol -> size.y, baseVol -> size.z);


	// create 3D array
    cudaChannelFormatDesc baseDesc = cudaCreateChannelDesc<float>();
    cutilSafeCall( cudaMalloc3DArray(&dBaseArray, &baseDesc, baseVolumeSize) );

    // copy data to 3D array
    cudaMemcpy3DParms baseCopyParams = {0};
    baseCopyParams.srcPtr   = make_cudaPitchedPtr((void*)baseVol->data, baseVolumeSize.width*sizeof(float), baseVolumeSize.width, baseVolumeSize.height);
    baseCopyParams.dstArray = dBaseArray;
    baseCopyParams.extent   = baseVolumeSize;
    baseCopyParams.kind     = cudaMemcpyHostToDevice;
    cutilSafeCall( cudaMemcpy3D(&baseCopyParams) );

    // set texture parameters
    tBaseVolData.normalized = false;                      // access with normalized texture coordinates
    if(!FLIRT_HARDWARE_INTERPOLATION || baseVol->dim.x >=4.0)
		tBaseVolData.filterMode = cudaFilterModePoint;      
	else 
		tBaseVolData.filterMode = cudaFilterModeLinear;   
    tBaseVolData.addressMode[0] = cudaAddressModeWrap;   // wrap texture coordinates
    tBaseVolData.addressMode[1] = cudaAddressModeWrap;
    tBaseVolData.addressMode[2] = cudaAddressModeWrap;

    // bind array to 3D texture
    cutilSafeCall(cudaBindTextureToArray(tBaseVolData, dBaseArray, baseDesc));

	const unsigned int regTotalSize=(regVol -> size.x * regVol ->size.y * regVol -> size. z);

	CUDA_SAFE_CALL(cudaMalloc((void**)&dRegVolData, sizeof(float) * regTotalSize));
	CUDA_SAFE_CALL(cudaMemcpy(dRegVolData, regVol->data,sizeof(float) * regTotalSize, cudaMemcpyHostToDevice));



	//Alloc memory
	CUDA_SAFE_CALL(cudaMalloc((void**)&dJointHist, sizeof(unsigned int) * binCount * binCount));
	CUDA_SAFE_CALL(cudaMalloc((void**)&dRegHistCopy, sizeof(unsigned int) * binCount ));
	CUDA_SAFE_CALL(cudaMalloc((void**)&dBaseHist, sizeof(unsigned int) * binCount));
	CUDA_SAFE_CALL(cudaMalloc((void**)&dRegHist, sizeof(unsigned int) * binCount));


	CUDA_SAFE_CALL(cudaMalloc((void**)&dEntropies, sizeof(float3) * binCount));


	CUDA_SAFE_CALL(cudaMalloc((void**)&dBlockPivots, sizeof(unsigned int) * (binCount + 1)));


	CUDA_SAFE_CALL(cudaMalloc((void**)&dRegSortedData, sizeof(float4) * regTotalSize));

	


	initRegFloat4(baseVol,regVol -> size);


}

void freeVolumeData(bool clearCpuBuffer)
{
	//Unbind textures
	cudaUnbindTexture( tRegVolData );

	//Free arrays
	cudaFreeArray( dBaseArray );


	CUDA_SAFE_CALL(cudaFree(dJointHist));
	CUDA_SAFE_CALL(cudaFree(dRegHistCopy));
	CUDA_SAFE_CALL(cudaFree(dBaseHist));
	CUDA_SAFE_CALL(cudaFree(dRegHist));
	CUDA_SAFE_CALL(cudaFree(dEntropies));

	CUDA_SAFE_CALL(cudaFree(dRegVolData));
	CUDA_SAFE_CALL(cudaFree(dRegSortedData));
	CUDA_SAFE_CALL(cudaFree(dBlockPivots));

	/*CUDA_SAFE_CALL(cudaFree(cBaseCenter));
	CUDA_SAFE_CALL(cudaFree(cRegCenter));*/


	delete hEntropies;
	delete hRegHist;
	delete hJointHist;
	delete hBaseHist;
	delete hBlockPivots;
}

void calcEntropyGPU(Volume* baseVol, Volume* regVol, Matrix44* dAffine, float &baseEntropy, float& regEntropy, float& jointEntropy)
{
	//unsigned int timer;

	baseEntropy=0;
	regEntropy=0;
	jointEntropy=0;

	unsigned int binCount;

	if(baseVol->dim.x ==8.0)
		binCount=FLIRT_BIN_COUNT_8MM;
	else if(baseVol->dim.x ==4.0)
		binCount=FLIRT_BIN_COUNT_4MM;
	else 
		binCount=FLIRT_BIN_COUNT;

	//printf("%d,%d,%d\n",regVol->size.x,regVol->size.y,regVol->size.z);
	
	
	//cudaMemset(dBaseHist, 0,sizeof(unsigned int) * binCount);

	/*float mat[16];
	for(int i=0;i<4;i++)
		memcpy(mat+(i*4),dAffine->data[i],sizeof(float)*4);*/

	cutilCheckMsg("cudaMemcpyToSymbol  before cAffine failed");

	cudaMemcpyToSymbol( cAffine,dAffine->data , sizeof(float)*16) ;

	cutilCheckMsg("cudaMemcpyToSymbol  after cAffine failed");


	unsigned int totalPoints = regVol->size.x * regVol->size.y * regVol->size.z ;

	unsigned int blockCount = binCount;

	unsigned int binShift = 0;

	

	////if the first sub bin is too big, calculate it invidually

	if(hBlockPivots[1] > binCount * binCount)
	{

		binShift=1;

		if(!FLIRT_HARDWARE_INTERPOLATION || baseVol->dim.x >=4.0)
			calcFirstBinKernel<<< binCount,WARP_SIZE * WARP_COUNT, binCount * WARP_COUNT * sizeof(unsigned int) >>>(dRegSortedData,hBlockPivots[1],baseVol->size,regVol->size,dJointHist);
		else
			calcFirstBinHardwareInterpolatedKernel<<< binCount,WARP_SIZE * WARP_COUNT, binCount * WARP_COUNT * sizeof(unsigned int) >>>(dRegSortedData,hBlockPivots[1],baseVol->size,regVol->size,dJointHist);

		//CUDA_SAFE_CALL( cudaThreadSynchronize() );

		cutilCheckMsg("calcFirstBinKernel failed");
		

		unsigned int step = binCount * binCount /2;


		//Reduce (binCount * binCount) to binCount
		while(step > binCount/2)
		{
			if( step > 1024 )
				blockCount = binCount/2;
			else blockCount = 16;

			int blockSize = step / blockCount;
			reduceFirstBinKernel<<<blockCount,blockSize>>>(dJointHist,step);
			//CUDA_SAFE_CALL( cudaThreadSynchronize() );
			step /=2 ;
		}

		cudaMemcpy(hRegHist, dJointHist, sizeof(unsigned int)*binCount, cudaMemcpyDeviceToHost);

		unsigned int firstBinSum=0;

		for(unsigned int i=0;i<binCount;i++)
			firstBinSum += hRegHist[i];

		cudaMemcpy(dRegHist, &firstBinSum, sizeof(unsigned int), cudaMemcpyHostToDevice);

		blockCount = binCount -1;


	}

	//CUDA_SAFE_CALL( cudaThreadSynchronize() );

	cutilCheckMsg("before calcJointHistKernel failed");
	
	// calculate the joint hist in parallel
	if(!FLIRT_HARDWARE_INTERPOLATION || baseVol->dim.x >=4.0)
		calcJointHistKernel<<< blockCount,WARP_SIZE * WARP_COUNT, binCount * WARP_COUNT * sizeof(unsigned int) >>>(dRegSortedData,dBlockPivots,baseVol->size,regVol->size,dJointHist,dRegHist,totalPoints,binShift);
	else
		calcJointHistHardwareInterpolatedKernel<<< blockCount,WARP_SIZE * WARP_COUNT, binCount * WARP_COUNT * sizeof(unsigned int) >>>(dRegSortedData,dBlockPivots,baseVol->size,regVol->size,dJointHist,dRegHist,totalPoints,binShift);

	cutilCheckMsg("after calcJointHistKernel failed");


	//CUT_SAFE_CALL(cutCreateTimer(&timer));
	//CUT_SAFE_CALL(cutStartTimer(timer));

	cudaMemcpy(hRegHist, dRegHist, sizeof(unsigned int)*binCount, cudaMemcpyDeviceToHost);

	unsigned int jointSum=0;

	for(unsigned int i=0;i<binCount;i++)
		jointSum += hRegHist[i];

	float logJointSum = log((float)jointSum);
	
	//CUDA_SAFE_CALL( cudaThreadSynchronize() );

	//CUT_SAFE_CALL(cutStopTimer(timer));

	//float t = cutGetAverageTimerValue(timer);

	// calculate the base distribution in parallel

	//calcBaseHistKernel<<<binCount,binCount>>>(dJointHist,dBaseHist);

	//calcJointEntropies<<<binCount,binCount,sizeof(float)*binCount>>>(dJointHist,(float)jointSum,dEntropies,logJointSum);
	//calcBaseRegEntropies<<<8,binCount/8>>>(dBaseHist,dRegHist,(float)jointSum,dEntropies,logJointSum);

	
	//cudaMemcpy(hEntropies, dEntropies, sizeof(float3)*binCount, cudaMemcpyDeviceToHost);
	
	

	cudaMemcpy(hJointHist, dJointHist, sizeof(unsigned int)*binCount*binCount, cudaMemcpyDeviceToHost);
	cudaMemcpy(hBaseHist, dBaseHist, sizeof(unsigned int)*binCount, cudaMemcpyDeviceToHost);

	memset(hBaseHist,0,sizeof(unsigned int)*binCount);

	//float3 result={0,0,0};

	for(unsigned int i=0;i<binCount;i++)
	{	
		for(unsigned int j=0;j<binCount;j++)
		{
			if(hJointHist[i*binCount + j])jointEntropy+= -(float)hJointHist[i*binCount + j]/jointSum * (log((float)hJointHist[i*binCount + j])-logJointSum);
			hBaseHist[i]+=hJointHist[j*binCount + i];
		}
		
		if(hBaseHist[i])baseEntropy+= -(float)hBaseHist[i]/jointSum * (log((float)hBaseHist[i])-logJointSum);
		if(hRegHist[i])regEntropy+= -(float)hRegHist[i]/jointSum * (log((float)hRegHist[i])-logJointSum);
	}

	//printf("%f ",(jointEntropy-regEntropy-baseEntropy));

	//printf("JointSum:%d,MI:%f,%f,%f\n ",jointSum,jointEntropy,regEntropy,baseEntropy);

	


	
	//printf("%f ",t);

}
