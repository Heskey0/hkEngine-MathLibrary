#pragma once

#include "Bounds.h"

namespace hkEngine {

template <typename T>
inline bool Bounds3<T>::IntersectP(const Ray& ray, float* hitt0,
	float* hitt1) const {
	float t0 = 0, t1 = ray.tMax;
	for (int i = 0; i < 3; ++i) {
		// Update interval for _i_th bounding box slab
		float invRayDir = 1 / ray.d[i];
		float tNear = (pMin[i] - ray.o[i]) * invRayDir;
		float tFar = (pMax[i] - ray.o[i]) * invRayDir;

		// Update parametric interval from slab intersection $t$ values
		if (tNear > tFar) std::swap(tNear, tFar);

		// Update _tFar_ to ensure robust ray--bounds intersection
		tFar *= 1 + 2 * gamma(3);
		t0 = tNear > t0 ? tNear : t0;
		t1 = tFar < t1 ? tFar : t1;
		if (t0 > t1) return false;
	}
	if (hitt0) *hitt0 = t0;
	if (hitt1) *hitt1 = t1;
	return true;
}

template <typename T>
inline bool Bounds3<T>::IntersectP(const Ray& ray, const float3& invDir,
	const int dirIsNeg[3]) const {
	const Bounds3f& bounds = *this;
	// Check for ray intersection against $x$ and $y$ slabs
	float tMin = (bounds[dirIsNeg[0]].x - ray.o.x) * invDir.x;
	float tMax = (bounds[1 - dirIsNeg[0]].x - ray.o.x) * invDir.x;
	float tyMin = (bounds[dirIsNeg[1]].y - ray.o.y) * invDir.y;
	float tyMax = (bounds[1 - dirIsNeg[1]].y - ray.o.y) * invDir.y;

	// Update _tMax_ and _tyMax_ to ensure robust bounds intersection
	tMax *= 1 + 2 * gamma(3);
	tyMax *= 1 + 2 * gamma(3);
	if (tMin > tyMax || tyMin > tMax) return false;
	if (tyMin > tMin) tMin = tyMin;
	if (tyMax < tMax) tMax = tyMax;

	// Check for ray intersection against $z$ slab
	float tzMin = (bounds[dirIsNeg[2]].z - ray.o.z) * invDir.z;
	float tzMax = (bounds[1 - dirIsNeg[2]].z - ray.o.z) * invDir.z;

	// Update _tzMax_ to ensure robust bounds intersection
	tzMax *= 1 + 2 * gamma(3);
	if (tMin > tzMax || tzMin > tMax) return false;
	if (tzMin > tMin) tMin = tzMin;
	if (tzMax < tMax) tMax = tzMax;
	return (tMin < ray.tMax) && (tMax > 0);
}

}	//hkEngine