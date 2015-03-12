#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <math.h>
//#include <cstring>
#include <android/bitmap.h>
#include <android/log.h>
//#include <RayTracerContext.h>

#define  LOG_TAG    "LIBRAY"
#define	 CL_TAG		"OPENCL"
#define  LOGI(...)  __android_log_print(ANDROID_LOG_DEBUG,LOG_TAG,__VA_ARGS__)
#define  LOGE(...)  __android_log_print(ANDROID_LOG_DEBUG,LOG_TAG,__VA_ARGS__)
#define  LOGC(...)  __android_log_print(ANDROID_LOG_DEBUG,CL_TAG,__VA_ARGS__)



#define IOR_AIR 1.0f

#define DEF_R 1.0f
#define DEF_G 0.0f
#define DEF_B 0.0f
#define DEF_KD 1.0f
#define DEF_KS 0.5f
#define DEF_SHINE 100.0f
#define DEF_T 0.0f;
#define DEF_IOR 1.0f

inline float det(const float &ax, const float &ay, const float &az,
				 const float &bx, const float &by, const float &bz,
				 const float &cx, const float &cy, const float &cz) {
	return ax*(by*cz - bz*cy) + ay*(bz*cx - bx*cz) + az*(bx*cy - by*cx);
}

inline uint32_t argb_to_int(const uint32_t &a, const uint32_t &r, const uint32_t g, const uint32_t &b) {
	return (a<<24) | (b<<16) | (g<<8) | r;
}

inline float fast_rsqrt(float number) {
	long i;
	float x2, y;
	const float threehalfs = 1.5f;

	x2 = number * 0.5F;
	y  = number;
	i  = *(long*)&y;
	i  = 0x5f3759df - (i >> 1);
	y  = *(float*)&i;
	y  = y * (threehalfs - (x2 * y * y));   // 1st iteration
//      y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed

	return y;
}

struct Vertex {
	float x;
	float y;
	float z;
};

inline void operator *= (const float &a, Vertex &_v) {
	_v.x *= a;
	_v.y *= a;
	_v.z *= a;
}

inline void operator *= (Vertex &_v, const float &a) {
	_v.x *= a;
	_v.y *= a;
	_v.z *= a;
}

inline void operator += (Vertex &v1, const Vertex &v2) {
	v1.x = v1.x + v2.x;
	v1.y = v1.y + v2.y;
	v1.z = v1.z + v2.z;
}

inline void operator -= (Vertex &v1, const Vertex &v2) {
	v1.x = v1.x - v2.x;
	v1.y = v1.y - v2.y;
	v1.z = v1.z - v2.z;
}

inline void cross(const Vertex &v1, const Vertex &v2, Vertex &result) {
	result.x = v1.y*v2.z - v1.z*v2.y;
	result.y = v1.z*v2.x - v1.x*v2.z;
	result.z = v1.x*v2.y - v1.y*v2.x;
}

inline Vertex operator * (const float &a, const Vertex &v) {
	Vertex va;
	va.x = a * v.x;
	va.y = a * v.y;
	va.z = a * v.z;
	return va;
}

inline Vertex operator * (const Vertex &v, const float &a) {
	Vertex va;
	va.x = a * v.x;
	va.y = a * v.y;
	va.z = a * v.z;
	return va;
}

inline Vertex operator - (const Vertex &v) {
	Vertex vm;
	vm.x = -v.x;
	vm.y = -v.y;
	vm.z = -v.z;
	return vm;
}

inline Vertex operator - (const Vertex &v1, const Vertex &v2) {
	Vertex vm;
	vm.x = v1.x - v2.x;
	vm.y = v1.y - v2.y;
	vm.z = v1.z - v2.z;
	return vm;
}

inline Vertex operator + (const Vertex &v1, const Vertex &v2) {
	Vertex vm;
	vm.x = v1.x + v2.x;
	vm.y = v1.y + v2.y;
	vm.z = v1.z + v2.z;
	return vm;
}

inline float inv_norm(const Vertex &_v) {
	return fast_rsqrt((_v.x*_v.x + _v.y*_v.y) + _v.z*_v.z);
	//return 1.0f / sqrt((_v.x*_v.x + _v.y*_v.y) + _v.z*_v.z);
}

inline Vertex inverse(const Vertex &_v) {
    Vertex retv;
    retv.x = 1.0f / _v.x;
    retv.y = 1.0f / _v.y;
    retv.z = 1.0f / _v.z;
    return retv;
}

inline float norm(const Vertex &_v) {
	return sqrt((_v.x*_v.x + _v.y*_v.y) + _v.z*_v.z);
}

inline void normalize(Vertex &_v) {
	float norm = inv_norm(_v);
	_v.x *= norm;
	_v.y *= norm;
	_v.z *= norm;
}

inline float cosine1(const Vertex &v1, const Vertex &v2) { // Only normalize v1.
	return (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z) / norm(v1);
}

inline float cosine(const Vertex &v1, const Vertex &v2) { // Normalize both.
	return (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z) / (norm(v1) * norm(v2));
}

inline float dot(const Vertex &v1, const Vertex &v2) {
	return (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
}

struct PointLight {
	Vertex pos;
	float r;
	float g;
	float b;
};

struct Material {
	float r;
	float g;
	float b;
	float kd;
	float ks;
	float shine;
	float T;
	float ior;
};

struct Triangle {
	int v1_idx;
	int v2_idx;
	int v3_idx;
	int mat;
	float c[3];
	Vertex n;
	float n_inv_norm;
	int n_calculated; // 0 or 1
};

// For the ray stack
struct RayItem {
	Vertex dir;
	Vertex orig;
	int recursion_number;
	float intensity_coef;
};

// For passing of camera data to the CL device.
struct CamData_CL {
    Vertex pos;
    Vertex proj_orig;
    Vertex x_step, y_step;
};

struct thread_range {
	int start_x;
	int start_y;
	int end_x;
	int end_y;
};

struct Ray {
    Vertex orig;
    Vertex dir;
    Vertex inv_dir;
};

#endif
