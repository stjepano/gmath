#include "gmath.h"

#include <stdlib.h>
#include <string.h>

#include <math.h>

Vec3::Vec3() {
	EL[0] = EL[1] = EL[2] = 0.0f;
	EL[3] = 1.0f;
}

Vec3::Vec3(float x, float y, float z, float w) {
	EL[0] = x;
	EL[1] = y;
	EL[2] = z;
	EL[3] = w;
}

Vec3::Vec3(float elements[]) {
	memcpy(EL,elements,sizeof(float)*4);
}

Vec3::Vec3(const Vec3& o) {
	memcpy(EL,o.EL,sizeof(float)*4);
}

Vec3& Vec3::operator =(const Vec3& o) {
	memcpy(EL,o.EL,sizeof(float)*4);
	return (*this);
}

Vec3 Vec3::operator +(const Vec3& o) const {
	return Vec3(EL[0] + o.EL[0], EL[1] + o.EL[1], EL[2] + o.EL[2]);
}

Vec3 Vec3::operator -(const Vec3& o) const {
	return Vec3(EL[0] - o.EL[0], EL[1] - o.EL[1], EL[2] - o.EL[2]);
}

Vec3& Vec3::operator +=(const Vec3& o ) {
	EL[0] += o.EL[0];
	EL[1] += o.EL[1];
	EL[2] += o.EL[2];
	return (*this);
}

Vec3& Vec3::operator -=(const Vec3& o ) {
	EL[0] -= o.EL[0];
	EL[1] -= o.EL[1];
	EL[2] -= o.EL[2];
	return (*this);
}

float Vec3::length() const {
	return sqrtf( SQR(EL[0]) + SQR(EL[1]) + SQR(EL[2]) );
}

float Vec3::lengthSqr() const {
	return SQR(EL[0]) + SQR(EL[1]) + SQR(EL[2]);
}

Vec3& Vec3::scale(float s) {
	EL[0] *= s;
	EL[1] *= s;
	EL[2] *= s;
	return *this;
}

Vec3& Vec3::normalize() {
	float s = 1.0f/length();
	return scale(s);
}

float Dot(const Vec3& a, const Vec3& b) {
	return a.EL[0]*b.EL[0] + a.EL[1]*b.EL[1] + a.EL[2]*b.EL[2];
}

Vec3 Cross(const Vec3& a, const Vec3& b) {
	return Vec3(
		a.EL[1]*b.EL[2] - a.EL[2]*b.EL[1],
		a.EL[2]*b.EL[0] - a.EL[0]*b.EL[2],
		a.EL[0]*b.EL[1] - a.EL[1]*b.EL[0]
	);
}


Vec3 Scale(const Vec3& a, float s) {
	return Vec3(a.EL[0]*s,a.EL[1]*s,a.EL[2]*s);
}

Vec3 Normalize(const Vec3& a) {
	float s = 1.0f/a.length();
	return Scale(a,s);
}



/************************************
	Mat44
*************************************/


Mat44::Mat44() {
	E00 = 1.0f; E01 = 0.0f; E02 = 0.0f; E03 = 0.0f;
	E10 = 0.0f; E11 = 1.0f; E12 = 0.0f; E13 = 0.0f;
	E20 = 0.0f; E21 = 0.0f; E22 = 1.0f; E23 = 0.0f;
	E30 = 0.0f; E31 = 0.0f; E32 = 0.0f; E33 = 1.0f;
}

Mat44::Mat44(float elements[]) {
	memcpy(EL, elements, sizeof(float)*16);
}

Mat44::Mat44(const Mat44& o) {
	memcpy(EL, o.EL, sizeof(float)*16);
}

Mat44& Mat44::operator =(const Mat44& o) {
	memcpy(EL,o.EL,sizeof(float)*16);
	return *this;
}

Mat44 Mat44::operator+(const Mat44& o) const {
	float el[16];
	for ( int i = 0; i < 16; i++) {
		el[i] = this->EL[i] + o.EL[i];
	}
	return Mat44(el);
}

Mat44 Mat44::operator -(const Mat44& o) const {
	float el[16];
	for ( int i = 0; i < 16; i++) {
		el[i] = this->EL[i] - o.EL[i];
	}
	return Mat44(el);
}

Mat44& Mat44::operator +=(const Mat44& o) {
	for ( int i = 0; i < 16; i++ ) {
		EL[i] += o.EL[i];
	}
	return *this;
}

Mat44& Mat44::operator -=(const Mat44& o) {
	for ( int i = 0; i < 16; i++ ) {
		EL[i] -= o.EL[i];
	}
	return *this;
}

float& Mat44::el(int r, int c) {
	return EL[c*4+r];
}

const float Mat44::el(int r, int c) const {
	return EL[c*4+r];
}

Mat44 Mat44::operator *(const Mat44& o) const {
	Mat44 result;
	
	for ( int col = 0; col < 4; col++) {
		for (int row = 0; row < 4; row++) {
			result.el(row,col) = 	this->el(row,col+0) * o.el(row+0,col) + 
									this->el(row,col+1) * o.el(row+1,col) + 
									this->el(row,col+2) * o.el(row+2,col) + 
									this->el(row,col+3) * o.el(row+3,col);
		}
	}
	
	return result;
}

Mat44& Mat44::operator *=(const Mat44& o) {
	float results[16];
	for ( int col = 0; col < 4; col++) {
		for (int row = 0; row < 4; row++ ) {
			results[col*4+row] = 	this->el(row,col+0) * o.el(row+0,col) + 
									this->el(row,col+1) * o.el(row+1,col) + 
									this->el(row,col+2) * o.el(row+2,col) + 
									this->el(row,col+3) * o.el(row+3,col);
		}
	}
	memcpy(EL,results,sizeof(float)*16);
	return *this;
}


Mat44& Mat44::transpose() {
	float results[16];
	for (int col = 0; col < 4; col++) {
		for (int row = 0; row < 4; row++) {
			results[col*4+row] = this->EL[row*4+col];
		}
	}
	memcpy(EL,results,sizeof(float)*16);
	return *this;
}


float Mat44::determinant() const {
	/* a1 a11 a14 a4 - a1 a10 a15 a4 - a11 a13 a2 a4 + a10 a13 a3 a4 - a0 a11 a14 a5 + a0 a10 a15 a5 + a11 a12 a2 a5 - a10 a12 a3 a5 - a1 a11 a12 a6 + a0 a11 a13 a6 + a1 a10 a12 a7 - a0 a10 a13 a7 - a15 a2 a5 a8 + a14 a3 a5 a8 + a1 a15 a6 a8 - a13 a3 a6 a8 - a1 a14 a7 a8 + a13 a2 a7 a8 + a15 a2 a4 a9 - a14 a3 a4 a9 - a0 a15 a6 a9 + a12 a3 a6 a9 + a0 a14 a7 a9 - a12 a2 a7 a9 */
	float det = EL[1] * EL[11] * EL[14] * EL[4] - EL[1] * EL[10] * EL[15] * EL[4] - EL[11] * EL[13] * EL[2] * EL[4] +
				EL[10] * EL[13] * EL[3] * EL[4] - EL[0] * EL[11] * EL[14] * EL[5] + EL[0] * EL[10] * EL[15] * EL[5] + 
				EL[11] * EL[12] * EL[2] * EL[5] - EL[10] * EL[12] * EL[3] * EL[5] - EL[1] * EL[11] * EL[12] * EL[6] + 
				EL[0] * EL[11] * EL[13] * EL[6] + EL[1] * EL[10] * EL[12] * EL[7] - EL[0] * EL[10] * EL[13] * EL[7] - 
				EL[15] * EL[2] * EL[5] * EL[8] + EL[14] * EL[3] * EL[5] * EL[8] + EL[1] * EL[15] * EL[6] * EL[8] - 
				EL[13] * EL[3] * EL[6] * EL[8] - EL[1] * EL[14] * EL[7] * EL[8] + EL[13] * EL[2] * EL[7] * EL[8] + 
				EL[15] * EL[2] * EL[4] * EL[9] - EL[14] * EL[3] * EL[4] * EL[9] - EL[0] * EL[15] * EL[6] * EL[9] + 
				EL[12] * EL[3] * EL[6] * EL[9] + EL[0] * EL[14] * EL[7] * EL[9] - EL[12] * EL[2] * EL[7] * EL[9];
	return det;
}

Mat44& Mat44::invert() {
	float oneOverDet = 1.0f/determinant();
	float inv[16];
	inv[0] =   EL[5]*EL[10]*EL[15] - EL[5]*EL[11]*EL[14] - EL[9]*EL[6]*EL[15] + EL[9]*EL[7]*EL[14] + EL[13]*EL[6]*EL[11] - EL[13]*EL[7]*EL[10];
	inv[4] =  -EL[4]*EL[10]*EL[15] + EL[4]*EL[11]*EL[14] + EL[8]*EL[6]*EL[15] - EL[8]*EL[7]*EL[14] - EL[12]*EL[6]*EL[11] + EL[12]*EL[7]*EL[10];
	inv[8] =   EL[4]*EL[9]*EL[15] - EL[4]*EL[11]*EL[13] - EL[8]*EL[5]*EL[15] + EL[8]*EL[7]*EL[13] + EL[12]*EL[5]*EL[11] - EL[12]*EL[7]*EL[9];
	inv[12] = -EL[4]*EL[9]*EL[14] + EL[4]*EL[10]*EL[13] + EL[8]*EL[5]*EL[14] - EL[8]*EL[6]*EL[13] - EL[12]*EL[5]*EL[10] + EL[12]*EL[6]*EL[9];
	inv[1] =  -EL[1]*EL[10]*EL[15] + EL[1]*EL[11]*EL[14] + EL[9]*EL[2]*EL[15] - EL[9]*EL[3]*EL[14] - EL[13]*EL[2]*EL[11] + EL[13]*EL[3]*EL[10];
	inv[5] =   EL[0]*EL[10]*EL[15] - EL[0]*EL[11]*EL[14] - EL[8]*EL[2]*EL[15] + EL[8]*EL[3]*EL[14] + EL[12]*EL[2]*EL[11] - EL[12]*EL[3]*EL[10];
	inv[9] =  -EL[0]*EL[9]*EL[15] + EL[0]*EL[11]*EL[13] + EL[8]*EL[1]*EL[15] - EL[8]*EL[3]*EL[13] - EL[12]*EL[1]*EL[11] + EL[12]*EL[3]*EL[9];
	inv[13] =  EL[0]*EL[9]*EL[14] - EL[0]*EL[10]*EL[13] - EL[8]*EL[1]*EL[14] + EL[8]*EL[2]*EL[13] + EL[12]*EL[1]*EL[10] - EL[12]*EL[2]*EL[9];
	inv[2] =   EL[1]*EL[6]*EL[15] - EL[1]*EL[7]*EL[14] - EL[5]*EL[2]*EL[15] + EL[5]*EL[3]*EL[14] + EL[13]*EL[2]*EL[7] - EL[13]*EL[3]*EL[6];
	inv[6] =  -EL[0]*EL[6]*EL[15] + EL[0]*EL[7]*EL[14] + EL[4]*EL[2]*EL[15] - EL[4]*EL[3]*EL[14] - EL[12]*EL[2]*EL[7] + EL[12]*EL[3]*EL[6];
	inv[10] =  EL[0]*EL[5]*EL[15] - EL[0]*EL[7]*EL[13] - EL[4]*EL[1]*EL[15] + EL[4]*EL[3]*EL[13] + EL[12]*EL[1]*EL[7] - EL[12]*EL[3]*EL[5];
	inv[14] = -EL[0]*EL[5]*EL[14] + EL[0]*EL[6]*EL[13] + EL[4]*EL[1]*EL[14] - EL[4]*EL[2]*EL[13] - EL[12]*EL[1]*EL[6] + EL[12]*EL[2]*EL[5];
	inv[3] =  -EL[1]*EL[6]*EL[11] + EL[1]*EL[7]*EL[10] + EL[5]*EL[2]*EL[11] - EL[5]*EL[3]*EL[10] - EL[9]*EL[2]*EL[7] + EL[9]*EL[3]*EL[6];
	inv[7] =   EL[0]*EL[6]*EL[11] - EL[0]*EL[7]*EL[10] - EL[4]*EL[2]*EL[11] + EL[4]*EL[3]*EL[10] + EL[8]*EL[2]*EL[7] - EL[8]*EL[3]*EL[6];
	inv[11] = -EL[0]*EL[5]*EL[11] + EL[0]*EL[7]*EL[9] + EL[4]*EL[1]*EL[11] - EL[4]*EL[3]*EL[9] - EL[8]*EL[1]*EL[7] + EL[8]*EL[3]*EL[5];
	inv[15] =  EL[0]*EL[5]*EL[10] - EL[0]*EL[6]*EL[9] - EL[4]*EL[1]*EL[10] + EL[4]*EL[2]*EL[9] + EL[8]*EL[1]*EL[6] - EL[8]*EL[2]*EL[5];
	
	for ( int i = 0; i < 16; i++) 
		inv[i] *= oneOverDet;
		
	memcpy(EL,inv,sizeof(16));
	return *this;
}


void TransformVector(const Mat44& m, const Vec3& vIn, Vec3& vOut) {
	vOut.X = m.E00*vIn.X + m.E01*vIn.Y + m.E02*vIn.Z + m.E03*vIn.W;
	vOut.Y = m.E10*vIn.X + m.E11*vIn.Y + m.E12*vIn.Z + m.E13*vIn.W;
	vOut.Z = m.E20*vIn.X + m.E21*vIn.Y + m.E22*vIn.Z + m.E23*vIn.W;
}

void TransformVectors(const Mat44& m, int n, Vec3* in, Vec3* out) {
	for ( int i = 0; i < n; i++ ) {
		const Vec3& vIn = in[i];
		Vec3& vOut = out[i];
		vOut.X = m.E00*vIn.X + m.E01*vIn.Y + m.E02*vIn.Z + m.E03*vIn.W;
		vOut.Y = m.E10*vIn.X + m.E11*vIn.Y + m.E12*vIn.Z + m.E13*vIn.W;
		vOut.Z = m.E20*vIn.X + m.E21*vIn.Y + m.E22*vIn.Z + m.E23*vIn.W;
	}
}


void Transpose( const Mat44& m, Mat44& out ) {
	out = m;
	out.transpose();
}

void Inverse( const Mat44& m, Mat44& out ) {
	out = m;
	out.invert();
}

Mat44 ConstructRotXMatrix(float angleDeg) {
	float angle = DEG2RAD(angleDeg);
	
	Mat44 m;
	float c = cosf(angle);
	float s = sinf(angle);
	
	m.E11 = c; m.E12 = -s;
	m.E21 = s; m.E22 = c;
	
	return m;
}

Mat44 ConstructRotYMatrix(float angleDeg) {
	float angle = DEG2RAD(angleDeg);
	
	Mat44 m;
	float c = cosf(angle);
	float s = sinf(angle);
	
	m.E00 = c; m.E02 = s;
	m.E20 = -s; m.E22 = c;
	
	return m;
}

Mat44 ConstructRotZMatrix(float angleDeg) {
	float angle = DEG2RAD(angleDeg);
	
	Mat44 m;
	float c = cosf(angle);
	float s = sinf(angle);
	
	m.E00 = c; m.E01 = -s;
	m.E10 = s; m.E11 = c;
	return m;
}

Mat44 ConstructRotMatrix(float angleDeg, const Vec3& axis) {
	float angle = DEG2RAD(angleDeg);
	
	Mat44 m;
	float c = cosf(angle);
	float s = sinf(angle);
	
	float ax = axis.X;
	float ay = axis.Y;
	float az = axis.Z;
	
	m.E00 = c+(1-c)*SQR(ax);
	m.E01 = (1-c)*ax*ay - s*az;
	m.E02 = (1-c)*ax*az + s*ay;
	
	m.E10 = (1-c)*ax*ay + s*az;
	m.E11 = c+(1-c)*SQR(ay);
	m.E12 = (1-c)*ay*az - s*ax;
	
	m.E20 = (1-c)*ax*az - s*ay;
	m.E21 = (1-c)*ay*az + s*ax;
	m.E22 = c+(1-c)*SQR(az);
	
	return m;
}


Mat44 ConstructScaleMatrix(float sx, float sy, float sz) {
	Mat44 m;
	m.E00 = sx;
	m.E11 = sy;
	m.E22 = sz;
	return m;
}

Mat44 ConstructTranslationMatrix(float tx, float ty, float tz) {
	Mat44 m;
	m.E03 = tx;
	m.E13 = ty;
	m.E23 = tz;
	return m;
}