#include "math_constants.h"
#include <cutil_math.h>


texture<float, 3, cudaReadModeElementType> tBaseVolData;
texture<float, 3, cudaReadModeElementType> tRegVolData;
__constant__ __device__ float cAffine[16];
__constant__ __device__ float3 cBaseCenter[1];
__constant__ __device__ float3 cRegCenter[1];

#ifndef __DEVICE_EMULATION__
    #define WARP_LOG2SIZE 5
#else
    #define WARP_LOG2SIZE 0
#endif

__device__ int dGetIndex(int3 size, int x, int y ,int z)
{
	if(x>=size.x)x=size.x-1;
	if(y>=size.y)y=size.y-1;
	if(z>=size.z)z=size.z-1;

	if(x<0)x=0;
	if(y<0)y=0;
	if(z<0)z=0;

	return z*size.x*size.y + y*size.x + x;
}


//Peformance the  affine transform
__device__ float3 dTransform(float x, float y, float z)
{
	float3 result;

	result.x = cAffine[0]*x + cAffine[1]*y + cAffine[2]*z + cAffine[3];

	result.y = cAffine[4]*x + cAffine[5]*y + cAffine[6]*z + cAffine[7];

	result.z = cAffine[8]*x + cAffine[9]*y + cAffine[10]*z + cAffine[11];

	return result;

}

//Trilinear interpolation kernel
__device__ float dTrilinearInterpolation(float v000, float v001, float v010, 
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
      return (temp6 - temp5)*dz + temp5;
}

//Get the interoplated value
__device__ float dInterpolate(int3 size,  float3 point)
{
	int ix, iy, iz;  
	ix=(int) floor(point.x);
	iy=(int) floor(point.y);
	iz=(int) floor(point.z);

	//left-top-front
	float dx=point.x-ix, dy=point.y-iy, dz=point.z-iz;

	float v000=0, v001=0, v010=0, v011=0, v100=0, v101=0, v110=0, v111=0;

	v000 = tex3D(tBaseVolData,ix,iy,iz);
	v001 = tex3D(tBaseVolData,ix,iy,iz+1);
	v010 = tex3D(tBaseVolData,ix,iy+1,iz);
	v011 = tex3D(tBaseVolData,ix,iy+1,iz+1);
	v100 = tex3D(tBaseVolData,ix+1,iy,iz);
	v101 = tex3D(tBaseVolData,ix+1,iy,iz+1);
	v110 = tex3D(tBaseVolData,ix+1,iy+1,iz);
	v111 = tex3D(tBaseVolData,ix+1,iy+1,iz+1);

	//Use trilinear interpolation
	return dTrilinearInterpolation(v000,v001,v010,v011,v100,v101,v110,v111,dx,dy,dz);
}


__global__ static void buildRegFloat4Kernel(float4* dRegSortedData,float* dRegVolData, int3 regSize,unsigned int totalPoints)
{

	unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;

	if(index >= totalPoints) 
		return ;

	int x,y,z;
		
	z = index / (regSize.x * regSize.y);
	y = (index - z * (regSize.x * regSize.y)) / regSize.x;
	x = index - z * (regSize.x * regSize.y) - y * regSize.x ; 

	float4 p=make_float4((float)x,(float)y,(float)z,dRegVolData[index]);

	dRegSortedData[index]=p;
}
__global__ static void calcFirstBinKernel(float4* dRegSortedData,unsigned int firstBinPoints, int3 baseSize, int3 regSize,unsigned int* dJointHist)
{
	const unsigned int binCount=gridDim.x;
	
	const unsigned int binWidth=256/binCount;

	unsigned int index;

	const unsigned int tid = threadIdx.x;

	const unsigned int globalTid = blockIdx.x * blockDim.x + tid;
	
	const unsigned int gridSize=binCount * blockDim.x; 

	const unsigned int warpId = tid / WARP_SIZE;

	


	#ifndef __DEVICE_EMULATION__
        const unsigned int threadTag = tid << (32 - WARP_LOG2SIZE);
    #else
        const unsigned int threadTag = 0;
    #endif

	volatile extern __shared__ unsigned int sJointHist[];

	const unsigned int iteration = binCount * WARP_COUNT / blockDim.x;


	//Init the shared memory, no bank conflict
	for(int i=0;i< iteration;i++)
		sJointHist[i * blockDim.x + tid] =0;
	__syncthreads();


	for( index = globalTid; index < firstBinPoints; index += gridSize)
	{

		//unsigned int regDataIndex = index + startIndex;
		
		//Fetch the coordination and shift to the mass center
		float x = dRegSortedData[index].x - cRegCenter[0].x;
		float y = dRegSortedData[index].y - cRegCenter[0].y;
		float z = dRegSortedData[index].z - cRegCenter[0].z;
		

		float baseValue;

		int baseOffset;

		//Transform and shift to base volume's coordinate
		float3 baseCoord = dTransform(x, y, z) + cBaseCenter[0];


		//Inside the base volume
		if( baseCoord.x >=0 && baseCoord.x <= baseSize.x && baseCoord.y >=0 && baseCoord.y <= baseSize.y && baseCoord.z >=0 && baseCoord.z <= baseSize.z)
		{
			//baseValue = tex3D(tBaseVolData,baseCoord.x,baseCoord.y,baseCoord.z);
			baseValue = dInterpolate(baseSize,baseCoord);

			baseOffset = (int)baseValue / binWidth ;

			unsigned int taggedValue, offset, val;
			do
			{
				offset = warpId * binCount + baseOffset;

				val = sJointHist[offset] & 0x07FFFFFF;

				taggedValue = threadTag | (val + 1);

				sJointHist[offset] = taggedValue;

			} while (sJointHist[offset] != taggedValue);

		}

	}

	
	__syncthreads();

	for(int i=0;i< iteration;i++)
		sJointHist[i * blockDim.x + tid] = sJointHist[i * binCount + tid] & 0x07FFFFFF;

	__syncthreads();
	
	//reduce 8 warps to one
	if( tid  < binCount )
	{
		for(int i=0;i< 4;i++)
		{
			sJointHist[tid + i * binCount] += sJointHist[ 4 * binCount + tid + i * binCount];
		}
		__syncthreads();

		for(int i=0;i< 2;i++)
		{
			sJointHist[tid + i * binCount] += sJointHist[ 2 * binCount + tid + i * binCount];
		}
		__syncthreads();

		for(int i=0;i< 1;i++)
		{
			sJointHist[tid + i * binCount] += sJointHist[1 * binCount + tid + i * binCount];
		}
		__syncthreads();

		dJointHist[blockIdx.x * binCount + tid ] = sJointHist[ tid ];
	}
}



__global__ static void calcJointHistKernel(float4* dRegSortedData,unsigned int* dBlockPivots, int3 baseSize, int3 regSize,unsigned int* dJointHist,unsigned int* dRegHist, int totalPoints, unsigned int binShift)
{
	const unsigned int binCount=gridDim.x + binShift;
	
	const unsigned int binWidth=256/binCount;

	const unsigned int blockSize=blockDim.x; 

	unsigned int index;

	const unsigned int tid = threadIdx.x;

	const unsigned int warpId = tid / WARP_SIZE;

	const unsigned int currentBin = blockIdx.x + binShift;

	const unsigned int startIndex = dBlockPivots[currentBin];

	const unsigned int endIndex = dBlockPivots[currentBin+1];

	/*if(startIndex == endIndex)
		return ;*/

	#ifndef __DEVICE_EMULATION__
        const unsigned int threadTag = tid << (32 - WARP_LOG2SIZE);
    #else
        const unsigned int threadTag = 0;
    #endif

	volatile extern __shared__ unsigned int sJointHist[];

	// iteration equals to totalBinCounts / blockSize
	const unsigned int iteration = binCount * WARP_COUNT / blockSize;

	//Init the shared memory, no bank conflict
	for(int i=0;i< iteration;i++)
		sJointHist[i * blockSize + tid] =0;

	__syncthreads();


	for( index = startIndex + tid; index < endIndex; index += blockSize)
	{

		//unsigned int regDataIndex = index + startIndex;
		

		//Fetch the coordination and shift to the mass center
		float x = dRegSortedData[index].x - cRegCenter[0].x;
		float y = dRegSortedData[index].y - cRegCenter[0].y;
		float z = dRegSortedData[index].z - cRegCenter[0].z;
		

		float baseValue;

		int baseOffset;

		//Transform and shift to base volume's coordinate
		float3 baseCoord = dTransform(x, y, z) + cBaseCenter[0];

		//Inside the base volume
		if( baseCoord.x >=0 && baseCoord.x <= (float)baseSize.x && baseCoord.y >=0 && baseCoord.y <= (float)baseSize.y && baseCoord.z >=0 && baseCoord.z <= (float)baseSize.z)
		{
			//baseValue = tex3D(tBaseVolData,baseCoord.x,baseCoord.y,baseCoord.z);
			baseValue = dInterpolate(baseSize,baseCoord);

			baseOffset = (int)baseValue / binWidth ;

			unsigned int taggedValue, offset, val;
			do
			{
				offset = warpId * binCount + baseOffset;

				val = sJointHist[offset] & 0x07FFFFFF;

				taggedValue = threadTag | (val + 1);

				sJointHist[offset] = taggedValue;

			} while (sJointHist[offset] != taggedValue);

		}

	}

	
	__syncthreads();

	for(int i=0;i< iteration;i++)
		sJointHist[i * blockSize + tid] = sJointHist[i * blockSize + tid] & 0x07FFFFFF;

	__syncthreads();
	
	
	if( tid  < binCount )
	{
		for(int i=0;i< 4;i++)
		{
			sJointHist[tid + i * binCount] += sJointHist[ 4 * binCount + tid + i * binCount];
		}
		__syncthreads();

		for(int i=0;i< 2;i++)
		{
			sJointHist[tid + i * binCount] += sJointHist[ 2 * binCount + tid + i * binCount];
		}
		__syncthreads();

		for(int i=0;i< 1;i++)
		{
			sJointHist[tid + i * binCount] += sJointHist[1 * binCount + tid + i * binCount];
		}
		__syncthreads();

		dJointHist[currentBin * binCount + tid ] = sJointHist[ tid ];
	}

	unsigned int step = binCount/2;

	while(step > 0)
	{
		if(tid < step)sJointHist[tid] += sJointHist[tid + step];
		__syncthreads();
		step /= 2;
	}

	if(tid ==0)
		dRegHist[currentBin] = sJointHist[0];

	//Marginal Histogram of Reg

}


__global__ static void calcFirstBinHardwareInterpolatedKernel(float4* dRegSortedData,unsigned int firstBinPoints, int3 baseSize, int3 regSize,unsigned int* dJointHist)
{
	const unsigned int binCount=gridDim.x;
	
	const unsigned int binWidth=256/binCount;

	unsigned int index;

	const unsigned int tid = threadIdx.x;

	const unsigned int globalTid = blockIdx.x * blockDim.x + tid;
	
	const unsigned int gridSize=binCount * blockDim.x; 

	const unsigned int warpId = tid / WARP_SIZE;

	


	#ifndef __DEVICE_EMULATION__
        const unsigned int threadTag = tid << (32 - WARP_LOG2SIZE);
    #else
        const unsigned int threadTag = 0;
    #endif

	volatile extern __shared__ unsigned int sJointHist[];

	const unsigned int iteration = binCount * WARP_COUNT / blockDim.x;


	//Init the shared memory, no bank conflict
	for(int i=0;i< iteration;i++)
		sJointHist[i * blockDim.x + tid] =0;
	__syncthreads();


	for( index = globalTid; index < firstBinPoints; index += gridSize)
	{

		//unsigned int regDataIndex = index + startIndex;
		
		//Fetch the coordination and shift to the mass center
		float x = dRegSortedData[index].x - cRegCenter[0].x;
		float y = dRegSortedData[index].y - cRegCenter[0].y;
		float z = dRegSortedData[index].z - cRegCenter[0].z;
		

		float baseValue;

		int baseOffset;

		//Transform and shift to base volume's coordinate
		float3 baseCoord = dTransform(x, y, z) + cBaseCenter[0];


		//Inside the base volume
		if( baseCoord.x >=0 && baseCoord.x <= baseSize.x && baseCoord.y >=0 && baseCoord.y <= baseSize.y && baseCoord.z >=0 && baseCoord.z <= baseSize.z)
		{
			baseValue = tex3D(tBaseVolData,baseCoord.x,baseCoord.y,baseCoord.z);
			//baseValue = dInterpolate(baseSize,baseCoord);

			baseOffset = (int)baseValue / binWidth ;

			unsigned int taggedValue, offset, val;
			do
			{
				offset = warpId * binCount + baseOffset;

				val = sJointHist[offset] & 0x07FFFFFF;

				taggedValue = threadTag | (val + 1);

				sJointHist[offset] = taggedValue;

			} while (sJointHist[offset] != taggedValue);

		}

	}

	
	__syncthreads();

	for(int i=0;i< iteration;i++)
		sJointHist[i * blockDim.x + tid] = sJointHist[i * binCount + tid] & 0x07FFFFFF;

	__syncthreads();
	
	//reduce 8 warps to one
	if( tid  < binCount )
	{
		for(int i=0;i< 4;i++)
		{
			sJointHist[tid + i * binCount] += sJointHist[ 4 * binCount + tid + i * binCount];
		}
		__syncthreads();

		for(int i=0;i< 2;i++)
		{
			sJointHist[tid + i * binCount] += sJointHist[ 2 * binCount + tid + i * binCount];
		}
		__syncthreads();

		for(int i=0;i< 1;i++)
		{
			sJointHist[tid + i * binCount] += sJointHist[1 * binCount + tid + i * binCount];
		}
		__syncthreads();

		dJointHist[blockIdx.x * binCount + tid ] = sJointHist[ tid ];
	}
}



__global__ static void calcJointHistHardwareInterpolatedKernel(float4* dRegSortedData,unsigned int* dBlockPivots, int3 baseSize, int3 regSize,unsigned int* dJointHist,unsigned int* dRegHist, int totalPoints, unsigned int binShift)
{
	const unsigned int binCount=gridDim.x + binShift;
	
	const unsigned int binWidth=256/binCount;

	const unsigned int blockSize=blockDim.x; 

	unsigned int index;

	const unsigned int tid = threadIdx.x;

	const unsigned int warpId = tid / WARP_SIZE;

	const unsigned int currentBin = blockIdx.x + binShift;

	const unsigned int startIndex = dBlockPivots[currentBin];

	const unsigned int endIndex = dBlockPivots[currentBin+1];

	


	/*if(startIndex == endIndex)
		return ;*/

	#ifndef __DEVICE_EMULATION__
        const unsigned int threadTag = tid << (32 - WARP_LOG2SIZE);
    #else
        const unsigned int threadTag = 0;
    #endif

	volatile extern __shared__ unsigned int sJointHist[];

	// iteration equals to totalBinCounts / blockSize
	const unsigned int iteration = binCount * WARP_COUNT / blockSize;

	//Init the shared memory, no bank conflict
	for(int i=0;i< iteration;i++)
		sJointHist[i * blockSize + tid] =0;

	__syncthreads();


	for( index = startIndex + tid; index < endIndex; index += blockSize)
	{

		//unsigned int regDataIndex = index + startIndex;
		

		//Fetch the coordination and shift to the mass center
		float x = dRegSortedData[index].x - cRegCenter[0].x;
		float y = dRegSortedData[index].y - cRegCenter[0].y;
		float z = dRegSortedData[index].z - cRegCenter[0].z;
		

		float baseValue;

		int baseOffset;

		//Transform and shift to base volume's coordinate
		float3 baseCoord = dTransform(x, y, z) + cBaseCenter[0];



		//Inside the base volume
		if( baseCoord.x >=0 && baseCoord.x <= (float)baseSize.x && baseCoord.y >=0 && baseCoord.y <= (float)baseSize.y && baseCoord.z >=0 && baseCoord.z <= (float)baseSize.z)
		{
			baseValue = tex3D(tBaseVolData,baseCoord.x,baseCoord.y,baseCoord.z);
			//baseValue = dInterpolate(baseSize,baseCoord);

			baseOffset = (int)baseValue / binWidth ;

			unsigned int taggedValue, offset, val;
			do
			{
				offset = warpId * binCount + baseOffset;

				val = sJointHist[offset] & 0x07FFFFFF;

				taggedValue = threadTag | (val + 1);

				sJointHist[offset] = taggedValue;

			} while (sJointHist[offset] != taggedValue);

		}

	}

	
	__syncthreads();

	for(int i=0;i< iteration;i++)
		sJointHist[i * blockSize + tid] = sJointHist[i * blockSize + tid] & 0x07FFFFFF;

	__syncthreads();
	
	
	if( tid  < binCount )
	{
		for(int i=0;i< 4;i++)
		{
			sJointHist[tid + i * binCount] += sJointHist[ 4 * binCount + tid + i * binCount];
		}
		__syncthreads();

		for(int i=0;i< 2;i++)
		{
			sJointHist[tid + i * binCount] += sJointHist[ 2 * binCount + tid + i * binCount];
		}
		__syncthreads();

		for(int i=0;i< 1;i++)
		{
			sJointHist[tid + i * binCount] += sJointHist[1 * binCount + tid + i * binCount];
		}
		__syncthreads();

		dJointHist[currentBin * binCount + tid ] = sJointHist[ tid ];
	}

	unsigned int step = binCount/2;

	while(step > 0)
	{
		if(tid < step)sJointHist[tid] += sJointHist[tid + step];
		__syncthreads();
		step /= 2;
	}

	if(tid ==0)
		dRegHist[currentBin] = sJointHist[0];

	//Marginal Histogram of Reg

}

//__global__ static void calcBaseHistKernel(unsigned int * dJointHist,unsigned int * dBaseHist)
//{
//	//histIndex = reg * bin_counts + base;
//	int x = blockIdx.x;
//	int y = threadIdx.x;
//
//	unsigned int count = dJointHist[ x * blockDim.x + y];
//
//	atomicAdd(dBaseHist + y, count);  // y: base
//}

//Reductions
__global__ static void reduceFirstBinKernel(unsigned int * data,int step)
{

    unsigned int index= __mul24(blockIdx.x,blockDim.x) + threadIdx.x;

	data[index] += data[index + step];
}



// Entropy lists
// One block one bin
__global__ static void calcJointEntropies(unsigned int * dJointHist,float jointSum,float3* dEntropies,float logJointSum)
{
	const unsigned int binCount=gridDim.x;

	unsigned int index= __mul24(blockIdx.x,blockDim.x) + threadIdx.x;

	const unsigned int tid =  threadIdx.x;

	extern __shared__ float sBinEntropies[];

	sBinEntropies[tid] = (float)dJointHist[blockIdx.x * binCount + tid];

	if(sBinEntropies[tid] !=0)
		sBinEntropies[tid] = -sBinEntropies[tid]/jointSum * (log(sBinEntropies[tid]) - logJointSum);

	__syncthreads();
	
	
	//Reduction

	unsigned int step = binCount/2;

	while(step > 0)
	{
		if(tid < step)sBinEntropies[tid] += sBinEntropies[tid + step];
		
		__syncthreads();

		step /= 2;
	}

	if(tid ==0)
		dEntropies[blockIdx.x].z = sBinEntropies[0];
}

__global__ static void calcBaseRegEntropies(unsigned int * dBaseHist,unsigned int *  dRegHist,float jointSum,float3* dEntropies,float logJointSum)
{
	unsigned int index= __mul24(blockIdx.x,blockDim.x) + threadIdx.x;



	if(dRegHist[index] !=0)
	{
		dEntropies[index].x = -(float)dRegHist[index]/jointSum * (log((float)dRegHist[index])-logJointSum);
	}

	if(dBaseHist[index] !=0)
	{
		dEntropies[index].y = -(float)dBaseHist[index]/jointSum * (log((float)dBaseHist[index])-logJointSum);
	}
	
}
