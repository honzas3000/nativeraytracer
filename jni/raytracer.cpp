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

JNIEXPORT jint JNICALL Java_com_janlangr_nativeraytracerlib_RayTracerLibrary_addTriangle(JNIEnv* env, jclass clazz, jint meshID,
		jfloat v1x, jfloat v1y, jfloat v1z, jfloat v2x, jfloat v2y, jfloat v2z, jfloat v3x, jfloat v3y, jfloat v3z, jint mat_idx) {
	return context.addTriangle(meshID, v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z, mat_idx);
}

//JNIEXPORT jint JNICALL Java_com_janlangr_nativeraytracerlib_RayTracerLibrary_addVertex(JNIEnv* env, jclass clazz, jint meshID, jfloat x, jfloat y, jfloat z) {
//	return context.addVertex(meshID, x, y, z);
//}

JNIEXPORT jint JNICALL Java_com_janlangr_nativeraytracerlib_RayTracerLibrary_addMaterial(JNIEnv* env, jclass clazz,
		jfloat r, jfloat g, jfloat b, jfloat kd, jfloat ks, jfloat shine, jfloat _T, jfloat ior) {
	return context.addMaterial(r, g, b, kd, ks, shine, _T, ior);
}

JNIEXPORT jint JNICALL Java_com_janlangr_nativeraytracerlib_RayTracerLibrary_addPointLight(JNIEnv* env, jclass clazz,
		jfloat x, jfloat y, jfloat z, jfloat r, jfloat g, jfloat b) {
	return context.addPointLight(x, y, z, r, g, b);
}

// OBJECT METHODS //////////////
/// Create a new mesh template.
JNIEXPORT jint JNICALL Java_com_janlangr_nativeraytracerlib_RayTracerLibrary_createMeshTemplate(JNIEnv* env, jclass clazz, jint numVert, jint numTri) {
	return context.createMeshTemplate(numVert, numTri);
}
/// Create a new object based on a mesh template.
JNIEXPORT jint JNICALL Java_com_janlangr_nativeraytracerlib_RayTracerLibrary_createObject(JNIEnv* env, jclass clazz, jint meshID) {
	return context.createObject(meshID);
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
