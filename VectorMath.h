#ifndef __VECTOR_MATH_H__
#define __VECTOR_MATH_H__

#include <stdio.h>
#include <vector>
#include <string>
#include <cmath>
#include <array>
#include <algorithm>

#define VM_PI 3.14159265359

#ifdef ENABLE_HP
#define C_TYPE long double
#else
#define C_TYPE float
#endif

namespace VectorMath {

	inline C_TYPE maxCompFromVector(std::vector<C_TYPE>& value) {
		C_TYPE res = !value.empty() ? value[0] : 0;
		for (unsigned int i = 1; i < value.size(); ++i)
			if (value[i] > res) res = value[i];
		return res;
	}

	class VectorX {
	public:
		unsigned short n;
		std::vector<C_TYPE> v;

		VectorX() : n(0) {}
	};

	class Vector2 : public VectorX {
	public:
		Vector2() { n = 2; v = { 0, 0 }; }
		Vector2(C_TYPE x) { n = 2; v = { x, x }; }
		Vector2(C_TYPE x, C_TYPE y) { n = 2; v = { x, y }; }

		inline Vector2 operator+(Vector2& b) {
			return Vector2(	v[0] + b.v[0],
							v[1] + b.v[1]	);
		}

		inline Vector2 operator-(Vector2& b) {
			return Vector2(	v[0] - b.v[0],
							v[1] - b.v[1]	);
		}

		inline Vector2 operator*(Vector2& b) {
			return Vector2(	v[0] * b.v[0],
							v[1] * b.v[1]	);
		}

		inline Vector2 operator/(Vector2& b) {
			return Vector2(	v[0] / b.v[0],
							v[1] / b.v[1]	);
		}

		inline Vector2 normalize() {
			C_TYPE scalar = maxCompFromVector(v);
			return Vector2(v[0] / scalar, v[1] / scalar);
		}

		inline C_TYPE dot(Vector2& b) {
			return (v[0] * b.v[0] + v[1] * b.v[1]);
		}

		inline C_TYPE length() {
			return sqrt(v[0] * v[0] + v[1] * v[1]);
		}

		friend std::ostream& operator<<(std::ostream& os, const Vector2& value);
	};

	

	class Vector3 : public VectorX {
	public:
		Vector3() { n = 3; v = { 0, 0, 0 }; }
		Vector3(C_TYPE x) { n = 3; v = { x, x, x }; }
		Vector3(C_TYPE x, C_TYPE y, C_TYPE z) { n = 3; v = { x, y, z }; }

		inline Vector3 operator+(Vector3& b) {
			return Vector3(	v[0] + b.v[0],
							v[1] + b.v[1],
							v[2] + b.v[2]	);
		}

		inline Vector3 operator-(Vector3& b) {
			return Vector3(	v[0] - b.v[0],
							v[1] - b.v[1],
							v[2] - b.v[2]	);
		}

		inline Vector3 operator*(Vector3& b) {
			return Vector3(	v[0] * b.v[0],
							v[1] * b.v[1],
							v[2] * b.v[2]	);
		}

		inline Vector3 operator/(Vector3& b) {
			return Vector3(	v[0] / b.v[0],
							v[1] / b.v[1],
							v[2] / b.v[2]	);
		}

		inline Vector3 normalize() {
			C_TYPE scalar = maxCompFromVector(v);
			return Vector3(v[0] / scalar, v[1] / scalar, v[2] / scalar);
		}

		inline C_TYPE dot(Vector3& b) {
			return (v[0] * b.v[0] + v[1] * b.v[1] + v[2] * b.v[2]);
		}

		inline C_TYPE length() {
			return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
		}

		inline Vector3 cross(Vector3& b) {
			return Vector3(	v[1] * b.v[2] - v[2] * b.v[1],
							v[0] * b.v[2] - v[2] * b.v[0],
							v[0] * b.v[1] - v[1] * b.v[0]	);
		}

		friend std::ostream& operator<<(std::ostream& os, const Vector3& value);
	};

	class Vector4 : public VectorX {
	public:
		Vector4() { n = 4; v = { 0, 0, 0, 0 }; }
		Vector4(C_TYPE x) { n = 4; v = { x, x, x, x }; }
		Vector4(C_TYPE x, C_TYPE y, C_TYPE z, C_TYPE w) { n = 4; v = { x, y, z, w }; }

		inline Vector4 operator+(Vector4& b) {
			return Vector4(v[0] + b.v[0],
				v[1] + b.v[1],
				v[2] + b.v[2],
				v[3] + b.v[3]	);
		}

		inline Vector4 operator-(Vector4& b) {
			return Vector4(v[0] - b.v[0],
				v[1] - b.v[1],
				v[2] - b.v[2],
				v[3] - b.v[3]	);
		}

		inline Vector4 operator*(Vector4& b) {
			return Vector4(v[0] * b.v[0],
				v[1] * b.v[1],
				v[2] * b.v[2],
				v[3] * b.v[3]	);
		}

		inline Vector4 operator/(Vector4& b) {
			return Vector4(v[0] / b.v[0],
				v[1] / b.v[1],
				v[2] / b.v[2],
				v[3] / b.v[3]	);
		}

		inline Vector4 normalize() {
			C_TYPE scalar = maxCompFromVector(v);
			return Vector4(v[0] / scalar, v[1] / scalar, v[2] / scalar, v[3] / scalar);
		}

		inline C_TYPE dot(Vector4& b) {
			return (v[0] * b.v[0] + v[1] * b.v[1] + v[2] * b.v[2] + v[3] * b.v[3]);
		}

		inline C_TYPE length() {
			return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3]);
		}

		friend std::ostream& operator<<(std::ostream& os, const Vector4& value);
	};

	inline std::ostream& output(std::ostream& os, std::vector<C_TYPE>& v) {
		for (const auto& value : v) os << value << " ";
		os << "\n";
		return os;
	}

	inline void normalize(std::vector<C_TYPE>& v) {
		C_TYPE scalar = maxCompFromVector(v);
		for (auto& c : v) c /= scalar;
	}

	inline std::ostream& operator<<(std::ostream& os, Vector2& value) { return output(os, value.v); }
	inline std::ostream& operator<<(std::ostream& os, Vector3& value) { return output(os, value.v); }
	inline std::ostream& operator<<(std::ostream& os, Vector4& value) { return output(os, value.v); }

	inline void normalize(Vector2& value) { normalize(value.v); }
	inline void normalize(Vector4& value) { normalize(value.v); }

	inline Vector3 normalize(Vector3 value) {
		Vector3 r = Vector3(value);
		normalize(r.v);
		return r;
	}

	inline Vector3 operator-(Vector3& a, const Vector3& b) { return a - b; }

	inline C_TYPE dot(Vector2& a, Vector2& b) { return a.v[0] * b.v[0] + a.v[1] * b.v[1]; }
	inline C_TYPE dot(Vector3& a, Vector3& b) { return a.v[0] * b.v[0] + a.v[1] * b.v[1] + a.v[2] * b.v[2]; }
	inline C_TYPE dot(Vector4& a, Vector4& b) { return a.v[0] * b.v[0] + a.v[1] * b.v[1] + a.v[2] * b.v[2] + a.v[3] * b.v[3]; }

	inline C_TYPE length(Vector2& a) { return sqrt(a.v[0] * a.v[0] + a.v[1] * a.v[1]); }
	inline C_TYPE length(Vector3& a) { return sqrt(a.v[0] * a.v[0] + a.v[1] * a.v[1] + a.v[2] * a.v[2]); }
	inline C_TYPE length(Vector4& a) { return sqrt(a.v[0] * a.v[0] + a.v[1] * a.v[1] + a.v[2] * a.v[2] + a.v[3] * a.v[3]); }

	inline Vector3 cross(const Vector3& a, const Vector3& b) {
		return Vector3(		a.v[1] * b.v[2] - a.v[2] * b.v[1],
							a.v[0] * b.v[2] - a.v[2] * b.v[0],
							a.v[0] * b.v[1] - a.v[1] * b.v[0]	);
	}

	class MatrixX : public VectorX {
	protected:
		virtual void init() = 0;

	public:
		MatrixX() {}
	};

	class Matrix2x2 : public MatrixX {
	protected:
		inline virtual void init() {
			n = 4;
			v = {
				1, 0,
				0, 1
			};
		}

	public:
		Matrix2x2() { init(); }
		Matrix2x2(std::array<C_TYPE, 4> value) { init(); for (unsigned int i = 0; i < n; ++i) v[i] = value[i]; }

		inline Matrix2x2 operator+(Matrix2x2 b) {
			Matrix2x2 r = Matrix2x2();
			for (unsigned int i = 0; i < n; ++i) r.v[i] = v[i] + b.v[i];
			return r;
		}

		inline Matrix2x2 operator-(Matrix2x2 b) {
			Matrix2x2 r = Matrix2x2();
			for (unsigned int i = 0; i < n; ++i) r.v[i] = v[i] - b.v[i];
			return r;
		}

		inline Matrix2x2 operator*(Matrix2x2 b) {
			Matrix2x2 r = Matrix2x2();
			for (unsigned int i = 0; i < n; ++i) r.v[i] = v[i] * b.v[i];
			return r;
		}

		inline Matrix2x2& operator*=(Matrix2x2 b) {
			for (unsigned int i = 0; i < n; ++i) v[i] = v[i] * b.v[i];
			return *this;
		}

		inline Matrix2x2 identity() {
			return Matrix2x2({
				1, 0,
				0, 1
				});
		}

		inline Matrix2x2 transpose() {
			return Matrix2x2({
					v[0],	v[2],
					v[1],	v[3]
				});
		}
	};

	class Matrix3x3 : public MatrixX {
	protected:
		inline virtual void init() {
			n = 9;
			v = {
				1, 0, 0,
				0, 1, 0,
				0, 0, 1
			};
		}

	public:
		Matrix3x3() { init(); }
		Matrix3x3(std::array<C_TYPE, 9> value) { init(); for (unsigned int i = 0; i < n; ++i) v[i] = value[i]; }

		inline Matrix3x3 operator+(Matrix3x3 b) {
			Matrix3x3 r = Matrix3x3();
			for (unsigned int i = 0; i < n; ++i) r.v[i] = v[i] + b.v[i];
			return r;
		}

		inline Matrix3x3 operator-(Matrix3x3 b) {
			Matrix3x3 r = Matrix3x3();
			for (unsigned int i = 0; i < n; ++i) r.v[i] = v[i] - b.v[i];
			return r;
		}

		inline Matrix3x3 operator*(Matrix3x3 b) {
			Matrix3x3 r = Matrix3x3();
			for (unsigned int i = 0; i < n; ++i) r.v[i] = v[i] * b.v[i];
			return r;
		}

		inline Matrix3x3& operator*=(Matrix3x3 b) {
			for (unsigned int i = 0; i < n; ++i) v[i] = v[i] * b.v[i];
			return *this;
		}

		inline Matrix3x3 identity() {
			return Matrix3x3({
				1, 0, 0,
				0, 1, 0,
				0, 0, 1
				});
		}

		inline Matrix3x3 transpose() {
			return Matrix3x3({
					v[0],	v[3],	v[6],
					v[1],	v[4],	v[7],
					v[2],	v[5],	v[8]
				});
		}
	};

	/*
	* Homogeneous matrix class 
	*/
	class Matrix4x4 : public MatrixX {
	protected:
		inline virtual void init() {
			n = 16;
			v = {
				1, 0, 0, 0,
				0, 1, 0, 0,
				0, 0, 1, 0,
				0, 0, 0, 1
			};
		}

	public:
		Matrix4x4() { init(); }
		Matrix4x4(std::array<C_TYPE, 16> value) { init(); for (unsigned int i = 0; i < n; ++i) v[i] = value[i]; }
		Matrix4x4(const Matrix4x4& m) {
			v = m.v;
			n = m.n;
		}

		inline Matrix4x4 operator+(Matrix4x4 b) {
			Matrix4x4 r = Matrix4x4();
			for (unsigned int i = 0; i < n; ++i) r.v[i] = v[i] + b.v[i];
			return r;
		}

		inline Matrix4x4 operator-(Matrix4x4 b) {
			Matrix4x4 r = Matrix4x4();
			for (unsigned int i = 0; i < n; ++i) r.v[i] = v[i] - b.v[i];
			return r;
		}

		inline Matrix4x4 operator*(Matrix4x4 b) {
			Matrix4x4 r = Matrix4x4();
			for (unsigned int i = 0; i < n; ++i) r.v[i] = v[i] * b.v[i];
			return r;
		}

		inline Matrix4x4& operator*=(Matrix4x4 b) {
			for (unsigned int i = 0; i < n; ++i) v[i] = v[i] * b.v[i];
			return *this;
		}

		inline Matrix4x4 identity() {
			return Matrix4x4({
				1, 0, 0, 0,
				0, 1, 0, 0,
				0, 0, 1, 0,
				0, 0, 0, 1
			});
		}

		inline Matrix4x4 transpose() {
			Matrix4x4 res = Matrix4x4(*this);
			Vector4 tempData[4][4];

			unsigned int col = 0;
			const unsigned int stride = sqrt(n);
			for (unsigned int i = 0; i < n; ++i) {
				unsigned int row = i % stride;
				tempData[col][row] = v[i];
				res.v[i] = tempData[row][col].v[i];
				if (row == 0) ++col;
			}

			return res;
		}
	};

	inline Matrix4x4 translate(Matrix4x4& m, C_TYPE x, C_TYPE y, C_TYPE z) {
			return Matrix4x4({
					m.v[0],	m.v[1],	m.v[2],	x,
					m.v[4],	m.v[5],	m.v[6],	y,
					m.v[8],	m.v[9],	m.v[10],z,
					m.v[12],m.v[13],m.v[14],m.v[15]
				});
		}

		inline Matrix4x4 translate(Matrix4x4& m, Vector3 vec) {
			return Matrix4x4({
					m.v[0],	m.v[1],	m.v[2],	vec.v[0],
					m.v[4],	m.v[5],	m.v[6],	vec.v[1],
					m.v[8],	m.v[9],	m.v[10],vec.v[2],
					m.v[12],m.v[13],m.v[14],m.v[15]
				});
		}

		inline Matrix4x4 scale(Matrix4x4& m, C_TYPE x) {
			return Matrix4x4({
					x,		m.v[1],	m.v[2],	m.v[3],
					m.v[4],	x,		m.v[6],	m.v[7],
					m.v[8],	m.v[9],	x,		m.v[11],
					m.v[12],m.v[13],m.v[14],m.v[15]
				});
		}

		inline Matrix4x4 scale(Matrix4x4& m, C_TYPE x, C_TYPE y, C_TYPE z) {
			return Matrix4x4({
					x,		m.v[1],	m.v[2],	m.v[3],
					m.v[4],	y,		m.v[6],	m.v[7],
					m.v[8],	m.v[9],	z,		m.v[11],
					m.v[12],m.v[13],m.v[14],m.v[15]
				});
		}

		inline Matrix4x4 scale(Matrix4x4& m, Vector3 vec) {
			return Matrix4x4({
					vec.v[0],	m.v[1],		m.v[2],		m.v[3],
					m.v[4],		vec.v[1],	m.v[6],		m.v[7],
					m.v[8],		m.v[9],		vec.v[2],	m.v[11],
					m.v[12],	m.v[13],	m.v[14],	m.v[15]
				});
		}

		inline Matrix4x4 rotateX(Matrix4x4& m, C_TYPE rad) {
			return Matrix4x4({
					m.v[0],	m.v[1],		m.v[2],		m.v[3],
					m.v[4],	cos(rad),	-sin(rad),	m.v[7],
					m.v[8],	sin(rad),	cos(rad),	m.v[11],
					m.v[12],m.v[13],	m.v[14],	m.v[15]
				});
		}

		inline Matrix4x4 rotateY(Matrix4x4& m, C_TYPE rad) {
			return Matrix4x4({
					cos(rad),	m.v[1],		sin(rad),	m.v[3],
					m.v[4],		m.v[5],		m.v[6],		m.v[7],
					-sin(rad),	m.v[9],		cos(rad),	m.v[11],
					m.v[12],	m.v[13],	m.v[14],	m.v[15]
				});
		}

		inline Matrix4x4 rotateZ(Matrix4x4& m, C_TYPE rad) {
			return Matrix4x4({
					cos(rad),	-sin(rad),	m.v[2],		m.v[3],
					sin(rad),	cos(rad),	m.v[6],		m.v[7],
					m.v[8],		m.v[9],		m.v[10],	m.v[11],
					m.v[12],	m.v[13],	m.v[14],	m.v[15]
				});
		}

		inline Matrix4x4& rotate(Matrix4x4& m, Vector3 axis, C_TYPE rad) {
			Matrix4x4 r = Matrix4x4();

			if (axis.v[0]) r *= rotateX(m, rad);
			if (axis.v[1]) r *= rotateY(m, rad);
			if (axis.v[2]) r *= rotateZ(m, rad);

			return r;
		}

		inline Matrix4x4 lookAt(Vector3& eye, Vector3& center, const Vector3& up, const Matrix4x4& m) {
			Vector3 zaxis = normalize(center - eye);
			Vector3 xaxis = normalize(cross(up, zaxis));
			Vector3 yaxis = cross(zaxis, xaxis);

			return Matrix4x4({
					xaxis.v[0],	yaxis.v[0], zaxis.v[0], -dot(xaxis, eye), 
					xaxis.v[1], yaxis.v[1], zaxis.v[1], -dot(yaxis, eye),
					xaxis.v[2],	yaxis.v[2], zaxis.v[2], -dot(zaxis, eye),
					0, 0, 0, 1
				});

			return m;
		}

		inline Matrix4x4& ortho(const C_TYPE &b, const C_TYPE &t, const C_TYPE &l, const C_TYPE &r, const C_TYPE &f, const C_TYPE &n, Matrix4x4& m) {
			m.v[0] = 2 / (r - l);
			m.v[1] = 0;
			m.v[2] = 0;
			m.v[3] = -(r + l) / (r - l);

			m.v[4] = 0;
			m.v[5] = 2 / (t - b);
			m.v[6] = 0;
			m.v[7] = -(t + b) / (t - b);

			m.v[8] = 0;
			m.v[9] = 0;
			m.v[10] = -2 / (f - n);
			m.v[11] = -(f + n) / (f - n);

			m.v[12] = 0;
			m.v[13] = 0;
			m.v[14] = 0;
			m.v[15] = 1;

			return m;
		}

		inline Matrix4x4& perspective(const C_TYPE &angleOfView, const C_TYPE &n, const C_TYPE& f, Matrix4x4& m) {
			
			float scale = 1 / tan(angleOfView * .5f * VM_PI / 180);

			m.v[0] = scale;
			m.v[5] = scale;
			m.v[10] = -f / (f- n);
			m.v[11] = -f * n / (f - n);
			m.v[14] = -1;
			m.v[15] = 0;

			return m;
		}
};

#endif