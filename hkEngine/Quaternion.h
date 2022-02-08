#pragma once

#include "Matrix.hpp"
#include "Vector.hpp"

namespace hkEngine {
	
	class Transform;
/*
* Quaternion struct
*/
class Quaternion {
public:
	// Quaternion Public Methods
	static Quaternion FromTransform(const Transform& t);
	Transform ToTransform() const;
	Quaternion(const Transform& t);
	Quaternion() : v(0, 0, 0), w(1) {}
	Quaternion& operator+=(const Quaternion& q) {
		v += q.v;
		w += q.w;
		return *this;
	}
	friend Quaternion operator+(const Quaternion& q1, const Quaternion& q2) {
		Quaternion ret = q1;
		return ret += q2;
	}
	Quaternion& operator-=(const Quaternion& q) {
		v -= q.v;
		w -= q.w;
		return *this;
	}
	Quaternion operator-() const {
		Quaternion ret;
		ret.v = -v;
		ret.w = -w;
		return ret;
	}
	friend Quaternion operator-(const Quaternion& q1, const Quaternion& q2) {
		Quaternion ret = q1;
		return ret -= q2;
	}
	Quaternion& operator*=(float f) {
		v *= f;
		w *= f;
		return *this;
	}
	Quaternion operator*(float f) const {
		Quaternion ret = *this;
		ret.v *= f;
		ret.w *= f;
		return ret;
	}
	Quaternion& operator/=(float f) {
		v /= f;
		w /= f;
		return *this;
	}
	Quaternion operator/(float f) const {
		Quaternion ret = *this;
		ret.v /= f;
		ret.w /= f;
		return ret;
	}

	friend std::ostream& operator<<(std::ostream& os, const Quaternion& q) {
		os <<  q.v.x, q.v.y, q.v.z, q.w;
		return os;
	}


	// Quaternion Public Data
	float3 v;
	float w;
};


#pragma region Inline Functions

inline float Dot(const Quaternion& q1, const Quaternion& q2) {
	return Dot(q1.v, q2.v) + q1.w * q2.w;
}

inline Quaternion Normalize(const Quaternion& q) {
	return q / std::sqrt(Dot(q, q));
}

// Quaternion Inline Functions
inline Quaternion operator*(float f, const Quaternion& q) { return q * f; }
Quaternion Slerp(float t, const Quaternion& q1, const Quaternion& q2);


#pragma endregion

}	//hkEngine

