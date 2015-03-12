#include "raytracer.h"

extern "C" {

JNIEXPORT void JNICALL Java_com_janlangr_nativeraytracerlib_RayTracerLibrary_provideKernelFile(JNIEnv* env, jclass clazz, jstring fileStr) {
	const char* fileChar = env->GetStringUTFChars(fileStr, NULL);
	context.provideKernelChar(fileChar);
}

JNIEXPORT jint JNICALL Java_com_janlangr_nativeraytracerlib_RayTracerLibrary_provideBitmap(JNIEnv* env, jobject thiz, jobject bitmap) {
	int result = AndroidBitmap_getInfo(env, bitmap, &(context.bitmapInfo));
	if(result < 0) {
		//LOGE("AndroidBitmap_getInfo error");
		return 1;
	}

	context.setProjectionResolution(context.bitmapInfo.width, context.bitmapInfo.height);
	return 0;
}

JNIEXPORT jint JNICALL Java_com_janlangr_nativeraytracerlib_RayTracerLibrary_addTriangle(JNIEnv* env, jobject thiz, jint v1_idx, jint v2_idx, jint v3_idx, jint mat_idx) {
	return context.addTriangle(v1_idx, v2_idx, v3_idx, mat_idx);
}

JNIEXPORT jint JNICALL Java_com_janlangr_nativeraytracerlib_RayTracerLibrary_addVertex(JNIEnv* env, jobject thiz, jfloat x, jfloat y, jfloat z) {
	return context.addVertex(x, y, z);
}

JNIEXPORT jint JNICALL Java_com_janlangr_nativeraytracerlib_RayTracerLibrary_addMaterial(JNIEnv* env, jobject thiz,
		jfloat r, jfloat g, jfloat b, jfloat kd, jfloat ks, jfloat shine, jfloat _T, jfloat ior) {
	return context.addMaterial(r, g, b, kd, ks, shine, _T, ior);
}

JNIEXPORT jint JNICALL Java_com_janlangr_nativeraytracerlib_RayTracerLibrary_addPointLight(JNIEnv* env, jobject thiz,
		jfloat x, jfloat y, jfloat z, jfloat r, jfloat g, jfloat b) {
	return context.addPointLight(x, y, z, r, g, b);
}

// DYNAMIC OBJECT METHODS //////////////
/// Create a new dynamic object within the context.
JNIEXPORT jint JNICALL Java_com_janlangr_nativeraytracerlib_RayTracerLibrary_createDynamicObject(JNIEnv* env, jclass clazz, jint numVert, jint numTri) {
	return context.createDynamicObject(numVert, numTri);
}

/// Add a vertex to the specified dynamic object.
JNIEXPORT jint JNICALL Java_com_janlangr_nativeraytracerlib_RayTracerLibrary_addDynamicVertex(JNIEnv* env, jclass clazz, jint dynObjIdx, jfloat x, jfloat y, jfloat z) {
	return context.addDynamicVertex(dynObjIdx, x, y, z);
}

/// Add a triangle to the specified dynamic object.
JNIEXPORT jint JNICALL Java_com_janlangr_nativeraytracerlib_RayTracerLibrary_addDynamicTriangle(JNIEnv* env, jclass clazz, jint dynObjIdx,
                                                                                                jint v1_idx, jint v2_idx, jint v3_idx, jint mat_idx) {
	return context.addDynamicTriangle(dynObjIdx, v1_idx, v2_idx, v3_idx, mat_idx);
}

JNIEXPORT void JNICALL Java_com_janlangr_nativeraytracerlib_RayTracerLibrary_initialize(JNIEnv* env, jobject thiz) {
	context.initialize();
}

JNIEXPORT void JNICALL Java_com_janlangr_nativeraytracerlib_RayTracerLibrary_finish(JNIEnv* env, jclass clazz) {
	context.finish();
}

JNIEXPORT jint JNICALL Java_com_janlangr_nativeraytracerlib_RayTracerLibrary_render(JNIEnv* env, jobject thiz, jobject bitmap, jfloat delta) {
	if(AndroidBitmap_lockPixels(env, bitmap, &context.pixels) < 0) {
		//LOGE("AndroidBitmap_lockPixels error");
		return -1;
	}

	context.rayTrace_entry(delta);

	AndroidBitmap_unlockPixels(env, bitmap);

	return 0;
}

}
