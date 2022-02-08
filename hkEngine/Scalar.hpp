#pragma once

#include <type_traits>
#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <assert.h>
#include <string.h>
#include <iomanip>

#include <intrin.h>


namespace hkEngine {

	// Global Constants

#define Maxfloat std::numeric_limits<float>::max()
#define Infinity std::numeric_limits<float>::infinity()


static float MachineEpsilon = std::numeric_limits<float>::epsilon() * 0.5;
static float ShadowEpsilon = 0.0001f;
static double Pi = 3.14159265358979323846;
static double InvPi = 0.31830988618379067154;
static double Inv2Pi = 0.15915494309189533577;
static double Inv4Pi = 0.07957747154594766788;
static double PiOver2 = 1.57079632679489661923;
static double PiOver4 = 0.78539816339744830961;
static double Sqrt2 = 1.41421356237309504880;


template <typename T>
inline bool isNaN(const T x) {
	return std::isnan(x);
}
template <>
inline bool isNaN(const int x) {
	return false;
}


// Global Inline Functions
inline uint32_t floatToBits(float f) {
	uint32_t ui;
	memcpy(&ui, &f, sizeof(float));
	return ui;
}

inline float BitsTofloat(uint32_t ui) {
	float f;
	memcpy(&f, &ui, sizeof(uint32_t));
	return f;
}

inline uint64_t floatToBits(double f) {
	uint64_t ui;
	memcpy(&ui, &f, sizeof(double));
	return ui;
}

inline double BitsTofloat(uint64_t ui) {
	double f;
	memcpy(&f, &ui, sizeof(uint64_t));
	return f;
}

inline float NextfloatUp(float v) {
	// Handle infinity and negative zero for _NextfloatUp()_
	if (std::isinf(v) && v > 0.) return v;
	if (v == -0.f) v = 0.f;

	// Advance _v_ to next higher float
	uint32_t ui = floatToBits(v);
	if (v >= 0)
		++ui;
	else
		--ui;
	return BitsTofloat(ui);
}

inline float NextfloatDown(float v) {
	// Handle infinity and positive zero for _NextfloatDown()_
	if (std::isinf(v) && v < 0.) return v;
	if (v == 0.f) v = -0.f;
	uint32_t ui = floatToBits(v);
	if (v > 0)
		--ui;
	else
		++ui;
	return BitsTofloat(ui);
}

inline double NextfloatUp(double v, int delta = 1) {
	if (std::isinf(v) && v > 0.) return v;
	if (v == -0.f) v = 0.f;
	uint64_t ui = floatToBits(v);
	if (v >= 0.)
		ui += delta;
	else
		ui -= delta;
	return BitsTofloat(ui);
}

inline double NextfloatDown(double v, int delta = 1) {
	if (std::isinf(v) && v < 0.) return v;
	if (v == 0.f) v = -0.f;
	uint64_t ui = floatToBits(v);
	if (v > 0.)
		ui -= delta;
	else
		ui += delta;
	return BitsTofloat(ui);
}

inline float gamma(int n) {
	return (n * MachineEpsilon) / (1 - n * MachineEpsilon);
}

inline float GammaCorrect(float value) {
	if (value <= 0.0031308f) return 12.92f * value;
	return 1.055f * std::pow(value, (float)(1.f / 2.4f)) - 0.055f;
}

inline float InverseGammaCorrect(float value) {
	if (value <= 0.04045f) return value * 1.f / 12.92f;
	return std::pow((value + 0.055f) * 1.f / 1.055f, (float)2.4f);
}

template <typename T, typename U, typename V>
inline T Clamp(T val, U low, V high) {
	if (val < low)
		return low;
	else if (val > high)
		return high;
	else
		return val;
}

template <typename T>
inline T Mod(T a, T b) {
	T result = a - (a / b) * b;
	return (T)((result < 0) ? result + b : result);
}

template <>
inline float Mod(float a, float b) {
	return std::fmod(a, b);
}

inline float Radians(float deg) { return (Pi / 180) * deg; }

inline float Degrees(float rad) { return (180 / Pi) * rad; }

inline float Log2(float x) {
	const double invLog2 = 1.442695040888963387004650940071;
	return std::log(x) * invLog2;
}

inline int Log2Int(uint32_t v) {
	unsigned long lz = 0;
	if (_BitScanReverse(&lz, v)) return lz;
	return 0;
}

inline int Log2Int(int32_t v) { return Log2Int((uint32_t)v); }

inline int Log2Int(uint64_t v) {
	unsigned long lz = 0;
#if defined(_WIN64)
	_BitScanReverse64(&lz, v);
#else
	if (_BitScanReverse(&lz, v >> 32))
		lz += 32;
	else
		_BitScanReverse(&lz, v & 0xffffffff);
#endif // _WIN64
	return lz;
}

inline int Log2Int(int64_t v) { return Log2Int((uint64_t)v); }

template <typename T>
inline bool IsPowerOf2(T v) {
	return v && !(v & (v - 1));
}

inline int32_t RoundUpPow2(int32_t v) {
	v--;
	v |= v >> 1;
	v |= v >> 2;
	v |= v >> 4;
	v |= v >> 8;
	v |= v >> 16;
	return v + 1;
}

inline int64_t RoundUpPow2(int64_t v) {
	v--;
	v |= v >> 1;
	v |= v >> 2;
	v |= v >> 4;
	v |= v >> 8;
	v |= v >> 16;
	v |= v >> 32;
	return v + 1;
}

inline int CountTrailingZeros(uint32_t v) {
	unsigned long index;
	if (_BitScanForward(&index, v))
		return index;
	else
		return 32;
}

template <typename Predicate>
int FindInterval(int size, const Predicate& pred) {
	int first = 0, len = size;
	while (len > 0) {
		int half = len >> 1, middle = first + half;
		// Bisect range based on value of _pred_ at _middle_
		if (pred(middle)) {
			first = middle + 1;
			len -= half + 1;
		}
		else
			len = half;
	}
	return Clamp(first - 1, 0, size - 2);
}

inline float Lerp(float t, float v1, float v2) { return (1 - t) * v1 + t * v2; }

inline bool Quadratic(float a, float b, float c, float* t0, float* t1) {
	// Find quadratic discriminant
	double discrim = (double)b * (double)b - 4 * (double)a * (double)c;
	if (discrim < 0) return false;
	double rootDiscrim = std::sqrt(discrim);

	// Compute quadratic _t_ values
	double q;
	if (b < 0)
		q = -.5 * (b - rootDiscrim);
	else
		q = -.5 * (b + rootDiscrim);
	*t0 = q / a;
	*t1 = c / q;
	if (*t0 > *t1) std::swap(*t0, *t1);
	return true;
}

inline float ErfInv(float x) {
	float w, p;
	x = Clamp(x, -.99999f, .99999f);
	w = -std::log((1 - x) * (1 + x));
	if (w < 5) {
		w = w - 2.5f;
		p = 2.81022636e-08f;
		p = 3.43273939e-07f + p * w;
		p = -3.5233877e-06f + p * w;
		p = -4.39150654e-06f + p * w;
		p = 0.00021858087f + p * w;
		p = -0.00125372503f + p * w;
		p = -0.00417768164f + p * w;
		p = 0.246640727f + p * w;
		p = 1.50140941f + p * w;
	}
	else {
		w = std::sqrt(w) - 3;
		p = -0.000200214257f;
		p = 0.000100950558f + p * w;
		p = 0.00134934322f + p * w;
		p = -0.00367342844f + p * w;
		p = 0.00573950773f + p * w;
		p = -0.0076224613f + p * w;
		p = 0.00943887047f + p * w;
		p = 1.00167406f + p * w;
		p = 2.83297682f + p * w;
	}
	return p * x;
}

inline float Erf(float x) {
	// constants
	float a1 = 0.254829592f;
	float a2 = -0.284496736f;
	float a3 = 1.421413741f;
	float a4 = -1.453152027f;
	float a5 = 1.061405429f;
	float p = 0.3275911f;

	// Save the sign of x
	int sign = 1;
	if (x < 0) sign = -1;
	x = std::abs(x);

	// A&S formula 7.1.26
	float t = 1 / (1 + p * x);
	float y =
		1 -
		(((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * std::exp(-x * x);

	return sign * y;
}


}	//hkEngine