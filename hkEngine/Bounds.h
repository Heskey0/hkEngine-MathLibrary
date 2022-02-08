#pragma once


#include "Ray.hpp"
#include "Vector.hpp"



//Bounding Box

namespace hkEngine {


/*
* Bounds2 class
*/
template <typename T>
class Bounds2 {
public:
	// Bounds2 Public Methods
	Bounds2() {
		T minNum = std::numeric_limits<T>::lowest();
		T maxNum = std::numeric_limits<T>::max();
		pMin = Vector2<T>(maxNum, maxNum);
		pMax = Vector2<T>(minNum, minNum);
	}
	explicit Bounds2(const Vector2<T>& p) : pMin(p), pMax(p) {}
	Bounds2(const Vector2<T>& p1, const Vector2<T>& p2) {
		pMin = Vector2<T>(std::min(p1.x, p2.x), std::min(p1.y, p2.y));
		pMax = Vector2<T>(std::max(p1.x, p2.x), std::max(p1.y, p2.y));
	}
	template <typename U>
	explicit operator Bounds2<U>() const {
		return Bounds2<U>((Vector2<U>)pMin, (Vector2<U>)pMax);
	}

	Vector2<T> Diagonal() const { return pMax - pMin; }
	T Area() const {
		Vector2<T> d = pMax - pMin;
		return (d.x * d.y);
	}
	int MaximumExtent() const {
		Vector2<T> diag = Diagonal();
		if (diag.x > diag.y)
			return 0;
		else
			return 1;
	}
	inline const Vector2<T>& operator[](int i) const {
		HK_CHECK(i == 0 || i == 1);
		return (i == 0) ? pMin : pMax;
	}
	inline Vector2<T>& operator[](int i) {
		HK_CHECK(i == 0 || i == 1);
		return (i == 0) ? pMin : pMax;
	}
	bool operator==(const Bounds2<T>& b) const {
		return b.pMin == pMin && b.pMax == pMax;
	}
	bool operator!=(const Bounds2<T>& b) const {
		return b.pMin != pMin || b.pMax != pMax;
	}
	Vector2<T> Lerp(const float2& t) const {
		return Vector2<T>(Lerp(t.x, pMin.x, pMax.x),
			Lerp(t.y, pMin.y, pMax.y));
	}
	Vector2<T> Offset(const Vector2<T>& p) const {
		Vector2<T> o = p - pMin;
		if (pMax.x > pMin.x) o.x /= pMax.x - pMin.x;
		if (pMax.y > pMin.y) o.y /= pMax.y - pMin.y;
		return o;
	}
	void BoundingSphere(Vector2<T>* c, float* rad) const {
		*c = (pMin + pMax) / 2;
		*rad = Inside(*c, *this) ? Distance(*c, pMax) : 0;
	}
	friend std::ostream& operator<<(std::ostream& os, const Bounds2<T>& b) {
		os << "[ " << b.pMin << " - " << b.pMax << " ]";
		return os;
	}

	// Bounds2 Public Data
	Vector2<T> pMin, pMax;
};



/*
* Bounds3 class
*/
template <typename T>
class Bounds3 {
public:
	// Bounds3 Public Methods
	Bounds3() {
		T minNum = std::numeric_limits<T>::lowest();
		T maxNum = std::numeric_limits<T>::max();
		pMin = Vector3<T>(maxNum, maxNum, maxNum);
		pMax = Vector3<T>(minNum, minNum, minNum);
	}
	explicit Bounds3(const Vector3<T>& p) : pMin(p), pMax(p) {}
	Bounds3(const Vector3<T>& p1, const Vector3<T>& p2)
		: pMin(std::min(p1.x, p2.x), std::min(p1.y, p2.y),
			std::min(p1.z, p2.z)),
		pMax(std::max(p1.x, p2.x), std::max(p1.y, p2.y),
			std::max(p1.z, p2.z)) {}
	const Vector3<T>& operator[](int i) const;
	Vector3<T>& operator[](int i);
	bool operator==(const Bounds3<T>& b) const {
		return b.pMin == pMin && b.pMax == pMax;
	}
	bool operator!=(const Bounds3<T>& b) const {
		return b.pMin != pMin || b.pMax != pMax;
	}
	Vector3<T> Corner(int corner) const {
		HK_CHECK(corner >= 0 && corner < 8);
		return Vector3<T>(
			(*this)[(corner & 1)].x,
			(*this)[(corner & 2) ? 1 : 0].y,
			(*this)[(corner & 4) ? 1 : 0].z
			);
	}
	Vector3<T> Diagonal() const { return pMax - pMin; }
	T SurfaceArea() const {
		Vector3<T> d = Diagonal();
		return 2 * (d.x * d.y + d.x * d.z + d.y * d.z);
	}
	T Volume() const {
		Vector3<T> d = Diagonal();
		return d.x * d.y * d.z;
	}
	int MaximumExtent() const {
		Vector3<T> d = Diagonal();
		if (d.x > d.y && d.x > d.z)
			return 0;
		else if (d.y > d.z)
			return 1;
		else
			return 2;
	}
	Vector3<T> Lerp(const float3& t) const {
		return Vector3<T>(Lerp(t.x, pMin.x, pMax.x),
			Lerp(t.y, pMin.y, pMax.y),
			Lerp(t.z, pMin.z, pMax.z));
	}
	Vector3<T> Offset(const Vector3<T>& p) const {
		Vector3<T> o = p - pMin;
		if (pMax.x > pMin.x) o.x /= pMax.x - pMin.x;
		if (pMax.y > pMin.y) o.y /= pMax.y - pMin.y;
		if (pMax.z > pMin.z) o.z /= pMax.z - pMin.z;
		return o;
	}
	void BoundingSphere(Vector3<T>* center, float* radius) const {
		*center = (pMin + pMax) / 2;
		*radius = Inside(*center, *this) ? Distance(*center, pMax) : 0;
	}
	template <typename U>
	explicit operator Bounds3<U>() const {
		return Bounds3<U>((Vector3<U>)pMin, (Vector3<U>)pMax);
	}
	bool IntersectP(const Ray& ray, float* hitt0 = nullptr,
		float* hitt1 = nullptr) const;
	bool IntersectP(const Ray& ray, const float3& invDir, const int dirIsNeg[3]) const;
	friend std::ostream& operator<<(std::ostream& os, const Bounds3<T>& b) {
		os << "[ " << b.pMin << " - " << b.pMax << " ]";
		return os;
	}

	// Bounds3 Public Data
	Vector3<T> pMin, pMax;
};

typedef Bounds2<float> Bounds2f;
typedef Bounds2<int> Bounds2i;
typedef Bounds3<float> Bounds3f;
typedef Bounds3<int> Bounds3i;

/*
* Bounds2iIterator class
*/
class Bounds2iIterator : public std::forward_iterator_tag {
public:
	Bounds2iIterator(const Bounds2i& b, const int2& pt)
		: p(pt), bounds(&b) {}
	Bounds2iIterator operator++() {
		advance();
		return *this;
	}
	Bounds2iIterator operator++(int) {
		Bounds2iIterator old = *this;
		advance();
		return old;
	}
	bool operator==(const Bounds2iIterator& bi) const {
		return p == bi.p && bounds == bi.bounds;
	}
	bool operator!=(const Bounds2iIterator& bi) const {
		return p != bi.p || bounds != bi.bounds;
	}

	int2 operator*() const { return p; }

private:
	void advance() {
		++p.x;
		if (p.x == bounds->pMax.x) {
			p.x = bounds->pMin.x;
			++p.y;
		}
	}
	int2 p;
	const Bounds2i* bounds;
};


/*
* inline functions
*/
#pragma region InlineFunctions




template <typename T>
Bounds3<T> Union(const Bounds3<T>& b, const Vector3<T>& p) {
	Bounds3<T> ret;
	ret.pMin = Min(b.pMin, p);
	ret.pMax = Max(b.pMax, p);
	return ret;
}

template <typename T>
Bounds3<T> Union(const Bounds3<T>& b1, const Bounds3<T>& b2) {
	Bounds3<T> ret;
	ret.pMin = Min(b1.pMin, b2.pMin);
	ret.pMax = Max(b1.pMax, b2.pMax);
	return ret;
}

template <typename T>
Bounds3<T> Intersect(const Bounds3<T>& b1, const Bounds3<T>& b2) {
	// Important: assign to pMin/pMax directly and don't run the Bounds2()
	// constructor, since it takes min/max of the points passed to it.  In
	// turn, that breaks returning an invalid bound for the case where we
	// intersect non-overlapping bounds (as we'd like to happen).
	Bounds3<T> ret;
	ret.pMin = Max(b1.pMin, b2.pMin);
	ret.pMax = Min(b1.pMax, b2.pMax);
	return ret;
}

template <typename T>
bool Overlaps(const Bounds3<T>& b1, const Bounds3<T>& b2) {
	bool x = (b1.pMax.x >= b2.pMin.x) && (b1.pMin.x <= b2.pMax.x);
	bool y = (b1.pMax.y >= b2.pMin.y) && (b1.pMin.y <= b2.pMax.y);
	bool z = (b1.pMax.z >= b2.pMin.z) && (b1.pMin.z <= b2.pMax.z);
	return (x && y && z);
}

template <typename T>
bool Inside(const Vector3<T>& p, const Bounds3<T>& b) {
	return (p.x >= b.pMin.x && p.x <= b.pMax.x && p.y >= b.pMin.y &&
		p.y <= b.pMax.y && p.z >= b.pMin.z && p.z <= b.pMax.z);
}

template <typename T>
bool InsideExclusive(const Vector3<T>& p, const Bounds3<T>& b) {
	return (p.x >= b.pMin.x && p.x < b.pMax.x&& p.y >= b.pMin.y &&
		p.y < b.pMax.y&& p.z >= b.pMin.z && p.z < b.pMax.z);
}

template <typename T, typename U>
inline Bounds3<T> Expand(const Bounds3<T>& b, U delta) {
	return Bounds3<T>(b.pMin - Vector3<T>(delta, delta, delta),
		b.pMax + Vector3<T>(delta, delta, delta));
}

// Minimum squared distance from point to box; returns zero if point is
// inside.
template <typename T, typename U>
inline float DistanceSquared(const Vector3<T>& p, const Bounds3<U>& b) {
	float dx = std::max({ float(0), b.pMin.x - p.x, p.x - b.pMax.x });
	float dy = std::max({ float(0), b.pMin.y - p.y, p.y - b.pMax.y });
	float dz = std::max({ float(0), b.pMin.z - p.z, p.z - b.pMax.z });
	return dx * dx + dy * dy + dz * dz;
}

template <typename T, typename U>
inline float Distance(const Vector3<T>& p, const Bounds3<U>& b) {
	return std::sqrt(DistanceSquared(p, b));
}

inline Bounds2iIterator begin(const Bounds2i& b) {
	return Bounds2iIterator(b, b.pMin);
}

inline Bounds2iIterator end(const Bounds2i& b) {
	// Normally, the ending point is at the minimum x value and one past
	// the last valid y value.
	int2 pEnd(b.pMin.x, b.pMax.y);
	// However, if the bounds are degenerate, override the end point to
	// equal the start point so that any attempt to iterate over the bounds
	// exits out immediately.
	if (b.pMin.x >= b.pMax.x || b.pMin.y >= b.pMax.y)
		pEnd = b.pMin;
	return Bounds2iIterator(b, pEnd);
}

template <typename T>
inline const Vector3<T>& Bounds3<T>::operator[](int i) const {
	HK_CHECK(i == 0 || i == 1);
	return (i == 0) ? pMin : pMax;
}

template <typename T>
inline Vector3<T>& Bounds3<T>::operator[](int i) {
	HK_CHECK(i == 0 || i == 1);
	return (i == 0) ? pMin : pMax;
}

template <typename T>
Bounds2<T> Union(const Bounds2<T>& b, const Vector2<T>& p) {
	Bounds2<T> ret;
	ret.pMin = Min(b.pMin, p);
	ret.pMax = Max(b.pMax, p);
	return ret;
}

template <typename T>
Bounds2<T> Union(const Bounds2<T>& b, const Bounds2<T>& b2) {
	Bounds2<T> ret;
	ret.pMin = Min(b.pMin, b2.pMin);
	ret.pMax = Max(b.pMax, b2.pMax);
	return ret;
}

template <typename T>
Bounds2<T> Intersect(const Bounds2<T>& b1, const Bounds2<T>& b2) {
	// Important: assign to pMin/pMax directly and don't run the Bounds2()
	// constructor, since it takes min/max of the points passed to it.  In
	// turn, that breaks returning an invalid bound for the case where we
	// intersect non-overlapping bounds (as we'd like to happen).
	Bounds2<T> ret;
	ret.pMin = Max(b1.pMin, b2.pMin);
	ret.pMax = Min(b1.pMax, b2.pMax);
	return ret;
}

template <typename T>
bool Overlaps(const Bounds2<T>& ba, const Bounds2<T>& bb) {
	bool x = (ba.pMax.x >= bb.pMin.x) && (ba.pMin.x <= bb.pMax.x);
	bool y = (ba.pMax.y >= bb.pMin.y) && (ba.pMin.y <= bb.pMax.y);
	return (x && y);
}

template <typename T>
bool Inside(const Vector2<T>& pt, const Bounds2<T>& b) {
	return (pt.x >= b.pMin.x && pt.x <= b.pMax.x && pt.y >= b.pMin.y &&
		pt.y <= b.pMax.y);
}

template <typename T>
bool InsideExclusive(const Vector2<T>& pt, const Bounds2<T>& b) {
	return (pt.x >= b.pMin.x && pt.x < b.pMax.x&& pt.y >= b.pMin.y &&
		pt.y < b.pMax.y);
}

template <typename T, typename U>
Bounds2<T> Expand(const Bounds2<T>& b, U delta) {
	return Bounds2<T>(b.pMin - Vector2<T>(delta, delta),
		b.pMax + Vector2<T>(delta, delta));
}



inline float3 OffsetRayOrigin(const float3& p, const float3& pError,
	const float3& n, const float3& w) {
	float d = Dot(Abs(n), pError);
	float3 offset = d * float3(n);
	if (Dot(w, n) < 0) offset = -offset;
	float3 po = p + offset;
	// Round offset point _po_ away from _p_
	for (int i = 0; i < 3; ++i) {
		if (offset[i] > 0)
			po[i] = NextfloatUp(po[i]);
		else if (offset[i] < 0)
			po[i] = NextfloatDown(po[i]);
	}
	return po;
}

inline float3 SphericalDirection(float sinTheta, float cosTheta, float phi) {
	return float3(sinTheta * std::cos(phi), sinTheta * std::sin(phi),
		cosTheta);
}

inline float3 SphericalDirection(float sinTheta, float cosTheta, float phi,
	const float3& x, const float3& y,
	const float3& z) {
	return sinTheta * std::cos(phi) * x + sinTheta * std::sin(phi) * y +
		cosTheta * z;
}

inline float SphericalTheta(const float3& v) {
	return std::acos(Clamp(v.z, -1, 1));
}

inline float SphericalPhi(const float3& v) {
	float p = std::atan2(v.y, v.x);
	return (p < 0) ? (p + 2 * Pi) : p;
}

#pragma endregion

}
