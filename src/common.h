#ifndef FLIRT_COMMON_H
#define FLIRT_COMMON_H

#include "vector_types.h"
#include "stdio.h"
#include "math.h"
#include "matrix44.h"
#include "options.h"

//Data type of orginal data, uint8 or uint16
#define FLIRT_UINT8 0
#define FLIRT_UINT16 1

template< typename T>
T Max(T x, T y)
{
	return x > y ? x : y;
}

template< typename T>
T Min(T x, T y)
{
	return x < y ? x : y;
}

typedef struct
{
	float3 rotation;
	float4 translateAndScale;
}Transform7;

typedef struct
{
	float3 rotation;
	float3 translation;
	float3 scale;
}Transform9;

typedef struct
{
	float3 rotation;
	float3 translation;
	float3 scale;
	float3 skew;
}Transform12;


#endif