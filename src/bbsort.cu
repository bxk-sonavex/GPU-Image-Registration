#include <stdio.h>
#include <stdlib.h>

#include <cutil.h>

#include "vector_types.h"

#include "bbsort.cuh"

#include "bbsort_kernel.cu"


float cmpKernel(float4 v){
	return v.w;
}



bool assignSliceToBuckets(unsigned int* sliceCount,int sliceSize,unsigned int* bucketOffset,unsigned int* bucketOfSlice,unsigned int* bucketSizes,unsigned int* sliceOffsetInBucket,int& bucketsCount,float step)
{
	int i=0;

	bool overflow=false;

	int tmpSum=0;

	bucketOffset[0]=0;

	for(i=0;i<sliceSize; i++){
		if(sliceCount[i] >BLOCK_SIZE)
		{
			overflow=true;
		}

		tmpSum += sliceCount[i];
		bucketOfSlice[i]=bucketsCount;
		bucketSizes[bucketsCount] = tmpSum;
		sliceOffsetInBucket[i]=tmpSum -sliceCount[i];
		if(tmpSum > BLOCK_SIZE )
		{	
			if(i != 0)
			{
				bucketOfSlice[i]=bucketsCount+1;
				bucketSizes[bucketsCount] -= sliceCount[i];
				sliceOffsetInBucket[i]=0;
				bucketOffset[bucketsCount+1]=bucketOffset[bucketsCount] + tmpSum -  sliceCount[i];

				bucketsCount++;
				tmpSum=sliceCount[i];
				bucketSizes[bucketsCount] = tmpSum;
			}
			else 
			{
				bucketOffset[bucketsCount+1]=bucketOffset[bucketsCount] + tmpSum ;
				sliceOffsetInBucket[i]=0;
				tmpSum=0;
				bucketsCount++;
			}
		}

	}
	bucketsCount++;

	return overflow;

}

void reduceMinMax(float4* dData,int size,float& result,bool isMax)
{

	int step;
	step=(size%2==0)?
		(size/2):(size/2 +1);

	int blockSize=BLOCK_SIZE;

	int blockCount;

	int length=size;

	float4 float4result;

	while(step > 0)
	{
		if(step%BLOCK_SIZE==0)
			blockCount=step/BLOCK_SIZE;
		else 
			blockCount=step/BLOCK_SIZE+1;

		if(isMax)
			reduceMaxD<<<blockCount,blockSize>>>(dData,step,length);
		else 
			reduceMinD<<<blockCount,blockSize>>>(dData,step,length);

		length=step;

		step=(step%2==0 || step==1)?(step/2):(step/2 +1);
	}

	CUDA_SAFE_CALL(cudaMemcpy(&float4result, dData, sizeof(float4), cudaMemcpyDeviceToHost));

	result=cmpKernel(float4result);
}

void evaluateDisorder(float4* dData,int size,float maxValue, float minValue, int& listOrder)
{
	int blockCount;

	if((size-1) % BLOCK_SIZE ==0)blockCount=size/BLOCK_SIZE;
	else blockCount=size/BLOCK_SIZE+1;

	float* dDiffData;
	CUDA_SAFE_CALL(cudaMalloc((void**)&dDiffData, sizeof(float) * size));

	calDifferenceD<<<blockCount,BLOCK_SIZE,(BLOCK_SIZE)*sizeof(float4)>>>(dData,dDiffData,size);


	float sum=0;

	int step;
	step=(size%2==0)?
		(size/2):(size/2 +1);

	int blockSize=BLOCK_SIZE;

	int length=size;

	while(step > 0)
	{

		if(step%BLOCK_SIZE==0)
			blockCount=step/BLOCK_SIZE;
		else 
			blockCount=step/BLOCK_SIZE+1;

		reduceSumD<<<blockCount,blockSize>>>(dDiffData,step,length);

		length=step;

		step=(step%2==0 || step==1)?(step/2):(step/2 +1);
	}

	CUDA_SAFE_CALL(cudaMemcpy(&sum, dDiffData, sizeof(float), cudaMemcpyDeviceToHost));

	if( sum < (maxValue - minValue) * size / 10)
		listOrder=NEARLY_SORTED;
	else 
		listOrder=DISORDERLY;

	CUDA_SAFE_CALL(cudaFree(dDiffData));
}


void bbSort(float4* dData,int size,int listOrder)
{
	unsigned int timer;

	CUT_SAFE_CALL(cutCreateTimer(&timer));

	//CUDA_SAFE_CALL( cudaThreadSynchronize() );
	
	//CUT_SAFE_CALL(cutStartTimer(timer));

	float minValue,maxValue;

	float4*  dTmpData;
	
	
	CUDA_SAFE_CALL(cudaMalloc((void**)&dTmpData, sizeof(float4) * size));
	
	CUDA_SAFE_CALL(cudaMemcpy(dTmpData, dData, sizeof(float4) * size, cudaMemcpyDeviceToDevice));


	reduceMinMax(dTmpData,size,maxValue,true);

	
	

	CUDA_SAFE_CALL(cudaMemcpy(dTmpData, dData, sizeof(float4) * size, cudaMemcpyDeviceToDevice));

	

	reduceMinMax(dTmpData,size,minValue,false);

	//CUDA_SAFE_CALL( cudaThreadSynchronize() );

	//CUT_SAFE_CALL(cutStopTimer(timer));  

	//float t = cutGetAverageTimerValue(timer);

	//printf( "Max min used: %fms\n",t);

	


	if(minValue == maxValue)
	{
		CUDA_SAFE_CALL(cudaFree(dTmpData));
		//printf(" flat, size is:%d\n",size);
		return ;
	}

	if(listOrder == AUTO_EVALUATE )
	{
		evaluateDisorder(dData,size,maxValue,minValue,listOrder);
	}

	

	int blockCount;

	if(size%BLOCK_SIZE==0)blockCount=size/BLOCK_SIZE;
	else blockCount=size/BLOCK_SIZE+1;

	float sliceStep =(float)(50.0*((double)(maxValue-minValue)/(double)size));

	int sliceSize = (int)((maxValue-minValue)/sliceStep) + 10;

	unsigned int* dSliceCounts;
	unsigned int* dOffsetInSlice;

	CUDA_SAFE_CALL(cudaMalloc((void**)&dOffsetInSlice, sizeof(unsigned int) * size));
	CUDA_SAFE_CALL(cudaMalloc((void**)&dSliceCounts, sizeof(unsigned int) * sliceSize));
	CUDA_SAFE_CALL(cudaMemset(dSliceCounts,0, sizeof(int) * sliceSize));

	//CUDA_SAFE_CALL( cudaThreadSynchronize() );
	
	//CUT_SAFE_CALL(cutResetTimer(timer));

	if(listOrder == NEARLY_SORTED)
	{
		assignElementToSlicesNearlySortedD<<<blockCount, BLOCK_SIZE>>>(dData,size,dSliceCounts,dOffsetInSlice,minValue,sliceStep,sliceSize,blockCount);
		//printf("NEARLY SORTED\n");
	}
	else 
		assignElementToSlicesD<<<blockCount, BLOCK_SIZE>>>(dData,size,dSliceCounts,dOffsetInSlice,minValue,sliceStep,sliceSize);

	//CUDA_SAFE_CALL( cudaThreadSynchronize() );

	//CUT_SAFE_CALL(cutStopTimer(timer));  

	//t = cutGetAverageTimerValue(timer);

	//printf( "slice dvision used: %fms\n",t);

	unsigned int* hSliceCounts=new unsigned int[sliceSize];
	CUDA_SAFE_CALL(cudaMemcpy(hSliceCounts, dSliceCounts, sizeof(unsigned int) * sliceSize, cudaMemcpyDeviceToHost));

	int looseBucketSize=size/100;

	unsigned int* hBucketOffsets=new unsigned int[looseBucketSize];
	unsigned int* hBucketSizes=new unsigned int[looseBucketSize];
	unsigned int* hBucketOfSlices=new unsigned int[sliceSize];
	unsigned int* hSliceOffsetInBucket=new unsigned int[sliceSize];
	int bucketsCount=0;

	memset(hBucketSizes,0,sizeof(int) * looseBucketSize);
	memset(hSliceOffsetInBucket,0,sizeof(unsigned int) * sliceSize);

	bool overflow;

	

	overflow = assignSliceToBuckets(hSliceCounts,sliceSize,hBucketOffsets,hBucketOfSlices,hBucketSizes,hSliceOffsetInBucket,bucketsCount,sliceStep);

	
	
	unsigned int* dBucketOffsets;
	unsigned int* dBucketSizes;

	unsigned int* dBucketOfSlices;
	unsigned int* dSliceOffsetInBucket;

	CUDA_SAFE_CALL(cudaMalloc((void**)&dBucketOfSlices, sizeof(unsigned int) * sliceSize));
	CUDA_SAFE_CALL(cudaMalloc((void**)&dSliceOffsetInBucket, sizeof(unsigned int) * sliceSize));
	CUDA_SAFE_CALL(cudaMalloc((void**)&dBucketOffsets, sizeof(unsigned int) * bucketsCount));
	CUDA_SAFE_CALL(cudaMalloc((void**)&dBucketSizes, sizeof(unsigned int) * bucketsCount));


	CUDA_SAFE_CALL(cudaMemcpy(dBucketOfSlices, hBucketOfSlices, sizeof(unsigned int) * sliceSize, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(dSliceOffsetInBucket, hSliceOffsetInBucket, sizeof(unsigned int) * sliceSize, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(dBucketOffsets, hBucketOffsets, sizeof(unsigned int) * bucketsCount, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(dBucketSizes, hBucketSizes, sizeof(unsigned int) * bucketsCount, cudaMemcpyHostToDevice));

	cudaBindTexture(0,tBucketOffsets,dBucketOffsets);
	cudaBindTexture(0,tBucketSizes,dBucketSizes);
	cudaBindTexture(0,tBucketOfSlices,dBucketOfSlices);
	cudaBindTexture(0,tSliceOffsetInBucket,dSliceOffsetInBucket);

	

	//CUDA_SAFE_CALL( cudaThreadSynchronize() );
	
	//CUT_SAFE_CALL(cutResetTimer(timer));

	assignElementToBucketD<<<blockCount, BLOCK_SIZE>>>(dData,dTmpData,size,dOffsetInSlice,minValue,sliceStep);

	//CUDA_SAFE_CALL( cudaThreadSynchronize() );

	//CUT_SAFE_CALL(cutStopTimer(timer));  

	//t = cutGetAverageTimerValue(timer);

	//printf( "bucket dvision used: %fms\n",t);

	


	CUDA_SAFE_CALL( cudaThreadSynchronize() );

	bitonicSortD<<<bucketsCount, BLOCK_SIZE, sizeof(float4) * BLOCK_SIZE>>>(dTmpData);

    CUDA_SAFE_CALL(cudaMemcpy(dData, dTmpData, sizeof(float4) * size, cudaMemcpyDeviceToDevice));

	
	
	if(overflow){
		for(int i=0;i<bucketsCount;i++)
		{
			if(hBucketSizes[i] > BLOCK_SIZE)
			{
				bbSort(dData + hBucketOffsets[i],hBucketSizes[i],listOrder);
			}
		}
	}

	
	
	delete hBucketOffsets;
	delete hBucketOfSlices;
	delete hSliceCounts;
	delete hBucketSizes;

	CUDA_SAFE_CALL(cudaFree(dOffsetInSlice));
	CUDA_SAFE_CALL(cudaFree(dSliceCounts));
	CUDA_SAFE_CALL(cudaFree(dTmpData));

	cudaUnbindTexture( tBucketSizes );
	CUDA_SAFE_CALL(cudaFree(dBucketSizes));

	cudaUnbindTexture( tBucketOffsets );
	CUDA_SAFE_CALL(cudaFree(dBucketOffsets));

	cudaUnbindTexture( tBucketOfSlices );
	CUDA_SAFE_CALL(cudaFree(dBucketOfSlices));

	cudaUnbindTexture( tSliceOffsetInBucket );
	CUDA_SAFE_CALL(cudaFree(dSliceOffsetInBucket));
}
