/*
 * Authored by: Chen, Shifu
 * 
 * Email: chen@gmtk.org / sf.chen@ieee.org
 *
 * The code is distributed under BSD license, you are allowed to use, modify or sell this code, but a statement is required if you used this code any where.
 * 
 */
#ifndef _BBSORT_KERNEL_H_
#define _BBSORT_KERNEL_H_

#include "math_constants.h"

texture<unsigned int, 1, cudaReadModeElementType> tBucketSizes;
texture<unsigned int, 1, cudaReadModeElementType> tBucketOffsets;
texture<unsigned int, 1, cudaReadModeElementType> tBucketOfSlices;
texture<unsigned int, 1, cudaReadModeElementType> tSliceOffsetInBucket;

__device__ float dCmpKernel(float4 v){
	return v.w;
}


__global__ static void reduceMaxD(float4 * dData,int step,int length)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

	if(index + step >=length)
		return ;
	dData[index] = dCmpKernel(dData[index])>dCmpKernel(dData[index+step])?dData[index]:dData[index+step];
}

__global__ static void reduceMinD(float4 * dData,int step,int length)
{
	
    int index = blockIdx.x * blockDim.x + threadIdx.x;

	if(index + step >=length)
		return ;

	dData[index] = dCmpKernel(dData[index])<dCmpKernel(dData[index+step])?dData[index]:dData[index+step];
}

__global__ static void reduceSumD(float * dDiffData,int step,int length)
{

    int index = blockIdx.x * blockDim.x + threadIdx.x;

	if(index + step >=length)
		return ;

	dDiffData[index] += dDiffData[index+step];
}

__global__ static void calDifferenceD(float4 * dData,float * dDiffData,int size)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

	if(index > size-1)
		return ;

    const unsigned int tid = threadIdx.x;
	
	extern __shared__ float4 sData[];

	sData[tid]=dData[index];

	__syncthreads();
	
	if(tid < blockDim.x -1)
		dDiffData[index] = abs(dCmpKernel(sData[tid+1]) - dCmpKernel(sData[tid]));
	else 
		dDiffData[index] =0;
	
}

__device__ inline void dSwap(float4 & a, float4 & b)
{
    float4 tmp = a;
    a = b;
    b = tmp;
}

__global__ static void bitonicSortD(float4 * datas)
{
    extern __shared__ float4 shared[];

	const unsigned int bid=blockIdx.x;

    const unsigned int tid = threadIdx.x;

	__shared__ unsigned int count;
	__shared__ unsigned int offset;

	if(tid == 0)
	{
		count=tex1Dfetch(tBucketSizes,bid);
		offset=tex1Dfetch(tBucketOffsets,bid);
	}

	__syncthreads();

    if(tid < count)
		shared[tid] = datas[tid+offset];
	else 
	{
		shared[tid].x=0x1f800000;
		shared[tid].y=0x1f800000;
		shared[tid].z=0x1f800000;
		shared[tid].w=0x1f800000;
	}

    __syncthreads();

    for (unsigned int k = 2; k <= BLOCK_SIZE; k *= 2)
    {
        for (unsigned int j = k / 2; j>0; j /= 2)
        {
            unsigned int ixj = tid ^ j;
            

            if (ixj > tid)
            {
                if ((tid & k) == 0)
                {
                    if (dCmpKernel(shared[tid]) > dCmpKernel(shared[ixj]))
                    {
                        dSwap(shared[tid], shared[ixj]);
                    }
                }
                else
                {
                    if (dCmpKernel(shared[tid]) < dCmpKernel(shared[ixj]))
                    {
                        dSwap(shared[tid], shared[ixj]);
                    }
                }
            }
            
            __syncthreads();
        }
    }
    if(tid < count)
		datas[tid+offset] = shared[tid];
}

__global__ void assignElementToSlicesD(float4* dDatas,int number,unsigned int* dSliceCounts,unsigned int* dOffsetInSlice,float minValue,float step,int sliceSize)
{
	unsigned int index= __mul24(blockIdx.x,blockDim.x) + threadIdx.x;

	if(index > number-1)
		return ;

	unsigned int s=(int)((dCmpKernel(dDatas[index]) - minValue)/ step);

	unsigned int offset=atomicInc(dSliceCounts + s,0xFFFFFFF);

	dOffsetInSlice[index] = offset;

}

__global__ void assignElementToSlicesNearlySortedD(float4* dDatas,int number,unsigned int* dSliceCounts,unsigned int* dOffsetInSlice,float minValue,float step,int sliceSize,int blockCount)
{
	unsigned int index= blockIdx.x + blockCount * threadIdx.x;

	if(index > number-1)
		return ;

	unsigned int s=(int)((dCmpKernel(dDatas[index]) - minValue)/ step);

	unsigned int offset=atomicInc(dSliceCounts + s,0xFFFFFFF);

	dOffsetInSlice[index] = offset;

}

__global__ void assignElementToBucketD(float4* dDatas,float4*  dNewDatas,int number,unsigned int* dOffsetInSlice,float minValue,float step)
{

	unsigned int index= __mul24(blockIdx.x,blockDim.x) + threadIdx.x;

	if(index > number-1)
		return ;

	unsigned int s=(int)((dCmpKernel(dDatas[index]) - minValue)/ step);

	unsigned int b=tex1Dfetch(tBucketOfSlices,s);

	unsigned int offset =tex1Dfetch(tBucketOffsets,b) + tex1Dfetch(tSliceOffsetInBucket,s) + dOffsetInSlice[index];

	dNewDatas[offset] =dDatas[index];

}

#endif // _BBSORT_KERNEL_H_
