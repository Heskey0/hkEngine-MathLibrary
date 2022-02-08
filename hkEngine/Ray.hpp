#pragma once


#include "Scalar.hpp"
#include "Vector.hpp"
#include "Medium.h"

namespace hkEngine {

class Ray {
public:
	// Ray Public Methods
	Ray() : tMax(Infinity), time(0.f), medium(nullptr) {}
	Ray(const float3& o, const float3& d, float tMax = Infinity,
		float time = 0.f, const Medium* medium = nullptr)
		: o(o), d(d), tMax(tMax), time(time), medium(medium) {}
	float3 operator()(float t) const { return o + d * t; }
	bool HasNaNs() const { return (o.HasNaNs() || d.HasNaNs() || isNaN(tMax)); }

	// Ray Public Data
	float3 o;
	float3 d;
	mutable float tMax;
	float time;
	const Medium* medium;
};

class RayDifferential : public Ray {
public:
	// RayDifferential Public Methods
	RayDifferential() { hasDifferentials = false; }
	RayDifferential(const float3& o, const float3& d, float tMax = Infinity,
		float time = 0.f, const Medium* medium = nullptr)
		: Ray(o, d, tMax, time, medium) {
		hasDifferentials = false;
	}
	RayDifferential(const Ray& ray) : Ray(ray) { hasDifferentials = false; }
	bool HasNaNs() const {
		return Ray::HasNaNs() ||
			(hasDifferentials &&
				(rxOrigin.HasNaNs() || ryOrigin.HasNaNs() ||
					rxDirection.HasNaNs() || ryDirection.HasNaNs()));
	}
	void ScaleDifferentials(float s) {
		rxOrigin = o + (rxOrigin - o) * s;
		ryOrigin = o + (ryOrigin - o) * s;
		rxDirection = d + (rxDirection - d) * s;
		ryDirection = d + (ryDirection - d) * s;
	}


	// RayDifferential Public Data
	bool hasDifferentials;
	float3 rxOrigin, ryOrigin;
	float3 rxDirection, ryDirection;
};


}