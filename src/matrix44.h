/*
 * Authored by: Chen, Shifu
 * 
 * Email: chen@gmtk.org / sf.chen@ieee.org
 *
 * The code is distributed under BSD license, you are allowed to use, modify or sell this code, but a statement is required if you used this code any where.
 * 
 */
#ifndef FLIRT_MATRIX44_H
#define FLIRT_MATRIX44_H

#include "vector_types.h"


#define PI 3.1415926

//Matrix44 in struct
struct MatrixStruct
{
	float data[4][4];
};

//Matrix44 in class
class Matrix44
{
public:
	float data[4][4];
	float det;

public:
	Matrix44();
	Matrix44(float m[4][4]);
	Matrix44(float a11, float a12, float a13, float a14,
		       float a21, float a22, float a23, float a24,
		       float a31, float a32, float a33, float a34,
		       float a41, float a42, float a43, float a44);

	//Operators: + - * with Matrix44
	friend const Matrix44 operator+ (const Matrix44& A, const Matrix44& B);
	friend const Matrix44 operator- (const Matrix44& A, const Matrix44& B);
	friend const Matrix44 operator* (const Matrix44& A, const Matrix44& B);
	friend const Matrix44 operator* (const Matrix44& A, const float& d);

	//Operator * with float3 and float4
	friend const  float4 operator* (const Matrix44& A, const float4& vec);
	friend const  float3 operator* (const Matrix44& A, const float3& vec);

	//Operators: += -= *=
	const Matrix44 operator+=(const Matrix44& mat);
	const Matrix44 operator-=(const Matrix44& mat);
	const Matrix44 operator*=(const float d);

	//Determinant
	void calcDet();
	float getDet();
	float det3(int r1, int r2, int r3, int c1, int c2, int c3) const;
	float det2(int r1, int r2, int c1, int c2) const
	{ 
		return (data[r1][c1]*data[r2][c2] - data[r2][c1]*data[r1][c2]); 
	}

	// & transpose & adjoint & identity
	 Matrix44 transpose() ;
	 Matrix44 adjoint() ;
	 Matrix44 inverse() ;
	void setIdentity();


	void setValues(float m[4][4]);
	void setValues(float a11, float a12, float a13, float a14,
		       float a21, float a22, float a23, float a24,
		       float a31, float a32, float a33, float a34,
		       float a41, float a42, float a43, float a44);

	MatrixStruct makeStruct();

	void setTranslate(float tx, float ty, float tz,float dim);
	void setRotateAxisAngle(float ax, float ay, float az, float angle);
	void setScale(float sx, float sy, float sz);
	void setSkew(float gx, float gy, float gz);

	void display();
};

Matrix44 makeClass(MatrixStruct matStruct);


#endif