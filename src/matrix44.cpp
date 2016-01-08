#include "matrix44.h"
#include "memory.h"
#include "math.h"
#include "stdio.h"

//Get a identity matrix by default
Matrix44::Matrix44()
{
	data[0][0] = 1;
	data[0][1] = 0;
	data[0][2] = 0;
	data[0][3] = 0;

	data[1][0] = 0;
	data[1][1] = 1;
	data[1][2] = 0;
	data[1][3] = 0;

	data[2][0] = 0;
	data[2][1] = 0;
	data[2][2] = 1;
	data[2][3] = 0;

	data[3][0] = 0;
	data[3][1] = 0;
	data[3][2] = 0;
	data[3][3] = 1;
}

Matrix44::Matrix44(float m[4][4])
{
	int i,j;
	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
			data[i][j]=m[i][j];
	calcDet();
}

Matrix44::Matrix44(float a11, float a12, float a13, float a14,
				   float a21, float a22, float a23, float a24,
				   float a31, float a32, float a33, float a34,
				   float a41, float a42, float a43, float a44)
{
	data[0][0] = a11;
	data[0][1] = a12;
	data[0][2] = a13;
	data[0][3] = a14;
	
	data[1][0] = a21;
	data[1][1] = a22;
	data[1][2] = a23;
	data[1][3] = a24;
	
	data[2][0] = a31;
	data[2][1] = a32;
	data[2][2] = a33;
	data[2][3] = a34;
	
	data[3][0] = a41;
	data[3][1] = a42;
	data[3][2] = a43;
	data[3][3] = a44;

	calcDet();
}

void Matrix44::setValues(float a11, float a12, float a13, float a14,
					     float a21, float a22, float a23, float a24,
					     float a31, float a32, float a33, float a34,
					     float a41, float a42, float a43, float a44)
{
	data[0][0] = a11;
	data[0][1] = a12;
	data[0][2] = a13;
	data[0][3] = a14;
	
	data[1][0] = a21;
	data[1][1] = a22;
	data[1][2] = a23;
	data[1][3] = a24;
	
	data[2][0] = a31;
	data[2][1] = a32;
	data[2][2] = a33;
	data[2][3] = a34;
	
	data[3][0] = a41;
	data[3][1] = a42;
	data[3][2] = a43;
	data[3][3] = a44;

	calcDet();
}

void Matrix44::setValues(float m[4][4])
{
	int i,j;
	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
			data[i][j]=m[i][j];
	calcDet();
}

///////////////////////// operator + - * with Matrix44 or value ///////////////////////////////////

const Matrix44 operator+ (const Matrix44& A, const Matrix44& B)
{
	Matrix44 mat;

	int i,j;

	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		{
			mat.data[i][j]=A.data[i][j] + B.data[i][j];
		}

	return mat;
}

const Matrix44 operator- (const Matrix44& A, const Matrix44& B)
{
	Matrix44 mat;

	int i,j;

	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		{
			mat.data[i][j]=A.data[i][j] - B.data[i][j];
			
		}

	return mat;
}

const Matrix44 operator* (const Matrix44& A, const Matrix44& B)
{
	Matrix44 mat;

	int i,j,k;

	for(i=0;i<4;i++)
	{
		for(j=0;j<4;j++)
		{
			mat.data[i][j]=0.0;
			for(k=0;k<4;k++)
				mat.data[i][j] += A.data[i][k]* B.data[k][j];
			
			//printf("%f ",mat.data[i][j]);
		}

		//printf("\n");
	}

	return mat;
}

const Matrix44 operator* (const Matrix44& A, const float& d)
{
	Matrix44 mat;

	int i,j;

	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		{
			mat.data[i][j] = A.data[i][j]*d;
		}

	return mat;
}

/////////////////////////// operator * with float3 and float4 ///////////////////////////////////////

const float4 operator* (const Matrix44& A, const float4& vec)
{
	
	float4 result={0.0,0.0,0.0,0.0};

	result.x += A.data[0][0]*vec.x + A.data[0][1]*vec.y + A.data[0][2]*vec.z + A.data[0][3]*vec.w;
	result.y += A.data[1][0]*vec.x + A.data[1][1]*vec.y + A.data[1][2]*vec.z + A.data[1][3]*vec.w;
	result.z += A.data[2][0]*vec.x + A.data[2][1]*vec.y + A.data[2][2]*vec.z + A.data[2][3]*vec.w;
	result.w += A.data[3][0]*vec.x + A.data[3][1]*vec.y + A.data[3][2]*vec.z + A.data[3][3]*vec.w;

	return result;
}

const float3 operator* (const Matrix44& A, const float3& vec)
{
	
	//Make a vec4 with this  vec3, the 4th element is set to 1.0
	float4 vec4;
	vec4.x=vec.x; vec4.y=vec.y; vec4.z=vec.z; vec4.w=1.0;

	float3 result;

	vec4= A * vec4;

	//Write back to vec3
	result.x = vec4.x;result.y = vec4.y;result.z = vec4.z;

	return result;
}

/////////////////////////// operator += -= *= ///////////////////////////////////////

const Matrix44 Matrix44::operator+=(const Matrix44& mat)
{
	int i,j;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			data[i][j] += mat.data[i][j];
	
	return *this;
}

const Matrix44  Matrix44::operator-=(const Matrix44& mat)
{
	int i,j;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			data[i][j] -= mat.data[i][j];
	
	return *this;
}


const Matrix44 Matrix44::operator*=(const float d)
{
	int i,j;

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			data[i][j] *= d;
	
	return *this;
}


///////////////// Determinant ////////////////////////////

void Matrix44::calcDet(){
	det=(data[0][0] * det3(1, 2, 3, 1, 2, 3)
	      - data[0][1] * det3(1, 2, 3, 0, 2, 3)
		  + data[0][2] * det3(1, 2, 3, 0, 1, 3)
		  - data[0][3] * det3(1, 2, 3, 0, 1, 2));
}

float Matrix44::det3(int r1, int r2, int r3, int c1, int c2, int c3) const
{
	return (data[r1][c1] * det2(r2, r3, c2, c3)
		  - data[r1][c2] * det2(r2, r3, c1, c3)
		  + data[r1][c3] * det2(r2, r3, c1, c2));
}
float Matrix44::getDet(){

	calcDet();
	return det;
}

/////////////////////  //////////////////////////////

Matrix44 Matrix44::transpose() 
{
	return Matrix44(data[0][0], data[1][0], data[2][0], data[3][0],
		              data[0][1], data[1][1], data[2][1], data[3][1],
		              data[0][2], data[1][2], data[2][2], data[3][2],
		              data[0][3], data[1][3], data[2][3], data[3][3]);
}

Matrix44 Matrix44::adjoint() 
{
	Matrix44 A;
	
	A.data[0][0] =  det3(1, 2, 3, 1, 2, 3);
	A.data[0][1] = -det3(0, 2, 3, 1, 2, 3);
	A.data[0][2] =  det3(0, 1, 3, 1, 2, 3);
	A.data[0][3] = -det3(0, 1, 2, 1, 2, 3);
	
	A.data[1][0] = -det3(1, 2, 3, 0, 2, 3);
	A.data[1][1] =  det3(0, 2, 3, 0, 2, 3);
	A.data[1][2] = -det3(0, 1, 3, 0, 2, 3);
	A.data[1][3] =  det3(0, 1, 2, 0, 2, 3);
	
	A.data[2][0] =  det3(1, 2, 3, 0, 1, 3);
	A.data[2][1] = -det3(0, 2, 3, 0, 1, 3);
	A.data[2][2] =  det3(0, 1, 3, 0, 1, 3);
	A.data[2][3] = -det3(0, 1, 2, 0, 1, 3);
	
	A.data[3][0] = -det3(1, 2, 3, 0, 1, 2);
	A.data[3][1] =  det3(0, 2, 3, 0, 1, 2);
	A.data[3][2] = -det3(0, 1, 3, 0, 1, 2);
	A.data[3][3] =  det3(0, 1, 2, 0, 1, 2);
	
	return A;
}

Matrix44 Matrix44::inverse()
{
	if(det == 0)
	{
		//cerr<<"Warning: devided by zero!"<<endl;
		return Matrix44();
	}
	calcDet();

	Matrix44 mat = adjoint();
	mat =mat* (1.0f/det);
	return mat;
}

void Matrix44::setIdentity()
{
	setValues(1.f, 0.f, 0.f, 0.f,
			  0.f, 1.f, 0.f, 0.f,
			  0.f, 0.f, 1.f, 0.f,
		      0.f, 0.f, 0.f, 1.f);
	det=1.0;
}

MatrixStruct Matrix44::makeStruct()
{
	MatrixStruct matStruct;

	memcpy(&matStruct.data,&data,sizeof(float)*16);

	return matStruct;
	
}

void Matrix44::display()
{
	int i,j;

	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			printf("%f ",data[i][j]);
			
		}
		printf("\n");
	}
}

Matrix44 makeClass(MatrixStruct matStruct)
{
	Matrix44 mat;

	mat.setValues(matStruct.data);

	return mat;
}
void Matrix44::setTranslate(float tx, float ty, float tz,float dim)
{
	data[0][0] = 1.f;
	data[0][1] = 0.f;
	data[0][2] = 0.f;
	data[0][3] = tx/dim;
	
	data[1][0] = 0.f;
	data[1][1] = 1.f;
	data[1][2] = 0.f;
	data[1][3] = ty/dim;
	
	data[2][0] = 0.f;
	data[2][1] = 0.f;
	data[2][2] = 1.f;
	data[2][3] = tz/dim;
	
	data[3][0] = 0.f;
	data[3][1] = 0.f;
	data[3][2] = 0.f;
	data[3][3] = 1.f;
}


void Matrix44::setRotateAxisAngle(float ax, float ay, float az, float angle)
{
	float q[4];
	float radAngle = (float)(angle * PI / 180.0);
	float c = cos(radAngle/2);
	float s = sin(radAngle/2);
	
	q[0] = c;
	q[1] = s * ax;
	q[2] = s * ay;
	q[3] = s * az;
	
	data[0][0] = 1.f - 2.f * (q[2] * q[2] + q[3] * q[3]);
	data[0][1] =       2.f * (q[1] * q[2] + q[3] * q[0]);
	data[0][2] =       2.f * (q[3] * q[1] - q[2] * q[0]);
	data[0][3] = 0.f;
	
	data[1][0] =       2.f * (q[1] * q[2] - q[3] * q[0]);
	data[1][1] = 1.f - 2.f * (q[3] * q[3] + q[1] * q[1]);
	data[1][2] =       2.f * (q[2] * q[3] + q[1] * q[0]);
	data[1][3] = 0.f;
	
	data[2][0] =       2.f * (q[3] * q[1] + q[2] * q[0]);
	data[2][1] =       2.f * (q[2] * q[3] - q[1] * q[0]);
	data[2][2] = 1.f - 2.f * (q[2] * q[2] + q[1] * q[1]);
	data[2][3] = 0.f;
	
	data[3][0] = 0.f;
	data[3][1] = 0.f;
	data[3][2] = 0.f;
	data[3][3] = 1.f;
}


void Matrix44::setScale(float sx, float sy, float sz)
{
	data[0][0] = sx;
	data[0][1] = 0.f;
	data[0][2] = 0.f;
	data[0][3] = 0.f;
	
	data[1][0] = 0.f;
	data[1][1] = sy;
	data[1][2] = 0.f;
	data[1][3] = 0.f;
	
	data[2][0] = 0.f;
	data[2][1] = 0.f;
	data[2][2] = sz;
	data[2][3] = 0.f;
	
	data[3][0] = 0.f;
	data[3][1] = 0.f;
	data[3][2] = 0.f;
	data[3][3] = 1.f;
}

void Matrix44::setSkew(float gx, float gy, float gz)
{
	data[0][0] = 1;
	data[0][1] = gx*gz;
	data[0][2] = gx;
	data[0][3] = 0.f;
	
	data[1][0] = gx;
	data[1][1] = 1;
	data[1][2] = 0.f;
	data[1][3] = 0.f;
	
	data[2][0] = 0.f;
	data[2][1] = gz;
	data[2][2] = 1;
	data[2][3] = 0.f;
	
	data[3][0] = 0.f;
	data[3][1] = 0.f;
	data[3][2] = 0.f;
	data[3][3] = 1.f;
}