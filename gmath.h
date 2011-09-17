#ifndef GMATH_H_
#define GMATH_H_

const float PI = 3.14159265f;

#define DEG2RAD( deg ) ((deg)*0.0174532925f)
#define RAD2DEG( rad ) ((rad)/0.0174532925f)
#define SQR(x) ((x)*(x))


struct Vec3 {
public:
	union {
		struct {
			float X,Y,Z,W;
		};
		float EL[4];
	};
	
	Vec3();
	Vec3(float x, float y, float z, float w = 1.0f);
	Vec3(const Vec3& other);
	Vec3(float elements[]);
	
	Vec3& operator=(const Vec3& o);
	
	Vec3 operator+(const Vec3& o) const;
	Vec3 operator-(const Vec3& o) const;
	
	Vec3& operator+=(const Vec3& o);
	Vec3& operator-=(const Vec3& o);
	
	float length() const;
	float lengthSqr() const;
	Vec3& normalize();
	Vec3& scale(float s);
};

Vec3 Cross(const Vec3& a, const Vec3& b);
float Dot(const Vec3& a, const Vec3& b);

Vec3 Normalize(const Vec3& a);
Vec3 Scale(const Vec3& a, float s);



struct Mat44 {
public:
	// opengl column major ordering, first 4 elements - first column
	union {
		struct {
			float E00, E10, E20, E30;
			float E01, E11, E21, E31;
			float E02, E12, E22, E32;
			float E03, E13, E23, E33;
		};
		float EL[16];	
	};
	
	Mat44();
	Mat44(float elements[]);
	Mat44(const Mat44& other);
	
	Mat44& operator=(const Mat44& other);
	
	Mat44 operator+(const Mat44& other) const;
	Mat44 operator-(const Mat44& other) const;
	
	Mat44& operator+=(const Mat44& other);
	Mat44& operator-=(const Mat44& other);
	
	Mat44 operator*(const Mat44& o) const;
	
	Mat44& operator*=(const Mat44& o);
	
	float& el(int r, int c);
	const float el(int r, int c) const;
	
	float determinant() const;
	
	Mat44& transpose();
	Mat44& invert();
};

void TransformVector( const Mat44& m, const Vec3& vec, Vec3* out);
void TransformVectors( const Mat44& m, int n, const Vec3* in, Vec3* out);

void Transpose( const Mat44& m, Mat44* out );
void Inverse( const Mat44& m, Mat44* out );


Mat44 ConstructRotXMatrix(float angleDeg);
Mat44 ConstructRotYMatrix(float angleDeg);
Mat44 ConstructRotZMatrix(float angleDeg);
Mat44 ConstructRotMatrix(float angleDeg, const Vec3& axis);
Mat44 ConstructScaleMatrix(float sx, float sy, float sz);
Mat44 ConstructTranslationMatrix(float tx, float ty, float tz);

#endif