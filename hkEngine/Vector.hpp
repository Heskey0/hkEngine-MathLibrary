#pragma once



#include "logging.hpp"
#include "Scalar.hpp"

namespace hkEngine{




/*
* Vector2 class
*/
template <typename T>
class Vector2 {
public:
	// Vector2 Public Methods
	Vector2() { x = y = 0; }
	Vector2(T xx, T yy) : x(xx), y(yy) { HK_CHECK(!HasNaNs()); }
	bool HasNaNs() const { return isNaN(x) || isNaN(y); }

	Vector2<T> operator+(const Vector2<T>& v) const {
		HK_CHECK(!v.HasNaNs());
		return Vector2(x + v.x, y + v.y);
	}

	Vector2<T>& operator+=(const Vector2<T>& v) {
		HK_CHECK(!v.HasNaNs());
		x += v.x;
		y += v.y;
		return *this;
	}
	Vector2<T> operator-(const Vector2<T>& v) const {
		HK_CHECK(!v.HasNaNs());
		return Vector2(x - v.x, y - v.y);
	}

	Vector2<T>& operator-=(const Vector2<T>& v) {
		HK_CHECK(!v.HasNaNs());
		x -= v.x;
		y -= v.y;
		return *this;
	}
	bool operator==(const Vector2<T>& v) const { return x == v.x && y == v.y; }
	bool operator!=(const Vector2<T>& v) const { return x != v.x || y != v.y; }
	template <typename U>
	Vector2<T> operator*(U f) const {
		return Vector2<T>(f * x, f * y);
	}

	template <typename U>
	Vector2<T>& operator*=(U f) {
		HK_CHECK(!isNaN(f));
		x *= f;
		y *= f;
		return *this;
	}
	template <typename U>
	Vector2<T> operator/(U f) const {
		HK_CHECK(f != 0);
		float inv = (float)1 / f;
		return Vector2<T>(x * inv, y * inv);
	}

	template <typename U>
	Vector2<T>& operator/=(U f) {
		HK_CHECK(f != 0);
		float inv = (float)1 / f;
		x *= inv;
		y *= inv;
		return *this;
	}
	Vector2<T> operator-() const { return Vector2<T>(-x, -y); }
	T operator[](int i) const {
		HK_CHECK(i >= 0 && i <= 1);
		if (i == 0) return x;
		return y;
	}

	T& operator[](int i) {
		HK_CHECK(i >= 0 && i <= 1);
		if (i == 0) return x;
		return y;
	}
	float LengthSquared() const { return x * x + y * y; }
	float Length() const { return std::sqrt(LengthSquared()); }

	// Vector2 Public Data
	T x, y;
};




/*
* Vector3 class
*/
template <typename T>
class Vector3 {
public:
	// Vector3 Public Methods
	T operator[](int i) const {
		HK_CHECK(i >= 0 && i <= 2);
		if (i == 0) return x;
		if (i == 1) return y;
		return z;
	}
	T& operator[](int i) {
		HK_CHECK(i >= 0 && i <= 2);
		if (i == 0) return x;
		if (i == 1) return y;
		return z;
	}
	Vector3() { x = y = z = 0; }
	Vector3(T x, T y, T z) : x(x), y(y), z(z) { HK_CHECK(!HasNaNs()); }
	bool HasNaNs() const { return isNaN(x) || isNaN(y) || isNaN(z); }

	Vector3<T> operator+(const Vector3<T>& v) const {
		HK_CHECK(!v.HasNaNs());
		return Vector3(x + v.x, y + v.y, z + v.z);
	}
	Vector3<T>& operator+=(const Vector3<T>& v) {
		HK_CHECK(!v.HasNaNs());
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}
	Vector3<T> operator-(const Vector3<T>& v) const {
		HK_CHECK(!v.HasNaNs());
		return Vector3(x - v.x, y - v.y, z - v.z);
	}
	Vector3<T>& operator-=(const Vector3<T>& v) {
		HK_CHECK(!v.HasNaNs());
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	}
	bool operator==(const Vector3<T>& v) const {
		return x == v.x && y == v.y && z == v.z;
	}
	bool operator!=(const Vector3<T>& v) const {
		return x != v.x || y != v.y || z != v.z;
	}
	template <typename U>
	Vector3<T> operator*(U s) const {
		return Vector3<T>(s * x, s * y, s * z);
	}
	template <typename U>
	Vector3<T>& operator*=(U s) {
		HK_CHECK(!isNaN(s));
		x *= s;
		y *= s;
		z *= s;
		return *this;
	}
	template <typename U>
	Vector3<T> operator/(U f) const {
		HK_CHECK(f != 0);
		float inv = (float)1 / f;
		return Vector3<T>(x * inv, y * inv, z * inv);
	}

	template <typename U>
	Vector3<T>& operator/=(U f) {
		HK_CHECK(f != 0);
		float inv = (float)1 / f;
		x *= inv;
		y *= inv;
		z *= inv;
		return *this;
	}
	Vector3<T> operator-() const { return Vector3<T>(-x, -y, -z); }
	float LengthSquared() const { return x * x + y * y + z * z; }
	float Length() const { return std::sqrt(LengthSquared()); }

	// Vector3 Public Data
	T x, y, z;
};



/*
* Inline functions
*/
#pragma region Inline Functions


template <typename T>
inline std::ostream& operator<<(std::ostream& os, const Vector3<T>& v) {
	os << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
	return os;
}

template <typename T, typename U>
inline Vector3<T> operator*(U s, const Vector3<T>& v) {
	return v * s;
}
template <typename T>
Vector3<T> Abs(const Vector3<T>& v) {
	return Vector3<T>(std::abs(v.x), std::abs(v.y), std::abs(v.z));
}

template <typename T>
inline T Dot(const Vector3<T>& v1, const Vector3<T>& v2) {
	HK_CHECK(!v1.HasNaNs() && !v2.HasNaNs());
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

template <typename T>
inline T AbsDot(const Vector3<T>& v1, const Vector3<T>& v2) {
	HK_CHECK(!v1.HasNaNs() && !v2.HasNaNs());
	return std::abs(Dot(v1, v2));
}

template <typename T>
inline Vector3<T> Cross(const Vector3<T>& v1, const Vector3<T>& v2) {
	HK_CHECK(!v1.HasNaNs() && !v2.HasNaNs());
	double v1x = v1.x, v1y = v1.y, v1z = v1.z;
	double v2x = v2.x, v2y = v2.y, v2z = v2.z;
	return Vector3<T>((v1y * v2z) - (v1z * v2y), (v1z * v2x) - (v1x * v2z),
		(v1x * v2y) - (v1y * v2x));
}

template <typename T>
inline Vector3<T> Normalize(const Vector3<T>& v) {
	return v / v.Length();
}
template <typename T>
T MinComponent(const Vector3<T>& v) {
	return std::min(v.x, std::min(v.y, v.z));
}

template <typename T>
T MaxComponent(const Vector3<T>& v) {
	return std::max(v.x, std::max(v.y, v.z));
}

template <typename T>
int MaxDimension(const Vector3<T>& v) {
	return (v.x > v.y) ? ((v.x > v.z) ? 0 : 2) : ((v.y > v.z) ? 1 : 2);
}

template <typename T>
Vector3<T> Min(const Vector3<T>& p1, const Vector3<T>& p2) {
	return Vector3<T>(std::min(p1.x, p2.x), std::min(p1.y, p2.y),
		std::min(p1.z, p2.z));
}

template <typename T>
Vector3<T> Max(const Vector3<T>& p1, const Vector3<T>& p2) {
	return Vector3<T>(std::max(p1.x, p2.x), std::max(p1.y, p2.y),
		std::max(p1.z, p2.z));
}

template <typename T>
Vector3<T> Permute(const Vector3<T>& v, int x, int y, int z) {
	return Vector3<T>(v[x], v[y], v[z]);
}

template <typename T>
inline void CoordinateSystem(const Vector3<T>& v1, Vector3<T>* v2,
	Vector3<T>* v3) {
	if (std::abs(v1.x) > std::abs(v1.y))
		*v2 = Vector3<T>(-v1.z, 0, v1.x) / std::sqrt(v1.x * v1.x + v1.z * v1.z);
	else
		*v2 = Vector3<T>(0, v1.z, -v1.y) / std::sqrt(v1.y * v1.y + v1.z * v1.z);
	*v3 = Cross(v1, *v2);
}

template <typename T>
inline Vector3<T> Faceforward(const Vector3<T>& v, const Vector3<T>& v2) {
	return (Dot(v, v2) < 0.f) ? -v : v;
}



template <typename T, typename U>
inline Vector2<T> operator*(U f, const Vector2<T>& v) {
	return v * f;
}
template <typename T>
inline float Dot(const Vector2<T>& v1, const Vector2<T>& v2) {
	HK_CHECK(!v1.HasNaNs() && !v2.HasNaNs());
	return v1.x * v2.x + v1.y * v2.y;
}

template <typename T>
inline float AbsDot(const Vector2<T>& v1, const Vector2<T>& v2) {
	HK_CHECK(!v1.HasNaNs() && !v2.HasNaNs());
	return std::abs(Dot(v1, v2));
}

template <typename T>
inline Vector2<T> Normalize(const Vector2<T>& v) {
	return v / v.Length();
}
template <typename T>
Vector2<T> Abs(const Vector2<T>& v) {
	return Vector2<T>(std::abs(v.x), std::abs(v.y));
}



template <typename T>
inline float Distance(const Vector3<T>& p1, const Vector3<T>& p2) {
	return (p1 - p2).Length();
}

template <typename T>
inline float DistanceSquared(const Vector3<T>& p1, const Vector3<T>& p2) {
	return (p1 - p2).LengthSquared();
}

template <typename T>
Vector3<T> Lerp(float t, const Vector3<T>& p0, const Vector3<T>& p1) {
	return (1 - t) * p0 + t * p1;
}



template <typename T>
inline float Distance(const Vector2<T>& p1, const Vector2<T>& p2) {
	return (p1 - p2).Length();
}

template <typename T>
inline float DistanceSquared(const Vector2<T>& p1, const Vector2<T>& p2) {
	return (p1 - p2).LengthSquared();
}


#pragma endregion

typedef Vector2<float> float2;
typedef Vector2<int> int2;
typedef Vector3<float> float3;
typedef Vector3<int> int3;

}