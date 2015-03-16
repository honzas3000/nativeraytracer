package com.janlangr.nativeraytracerlib;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import android.content.Context;
import android.graphics.Bitmap;
import android.util.Log;

public class RayTracerLibrary {
	public static Bitmap bitmap = null;
	
	static {
		int res_x = 512, res_y = 256;
        System.loadLibrary("raytracer");
        bitmap = Bitmap.createBitmap(res_x, res_y, Bitmap.Config.ARGB_8888);
        changeBitmap(bitmap);
        
    }
	
	/**
	 * Provide a bitmap for the ray tracer to draw to. The number
	 * of rays (width * height) is determined by the dimensions
	 * of this bitmap.
	 */
	public static boolean changeBitmap(Bitmap bmp) {
		if(bmp.getConfig() == Bitmap.Config.ARGB_8888) {
			bitmap = bmp;
			if(provideBitmap(bmp) == 0) { // Pass the bitmap's reference to the C++ code.
				return true;
			} else {
				return false;
			}
		} else {
			return false;
		}
	}
	
	
	
	public static int rayTraceScene(float delta) {
		return render(bitmap, delta);
	}
	
	/**
	 * Specify a new triangle using indices of already existing vertices.
	 * @param v1_idx
	 * @param v2_idx
	 * @param v3_idx
	 * @return
	 */
	public static native int addTriangle(int objIndex, float v1x, float v1y, float v1z, float v2x, float v2y, float v2z, float v3x, float v3y, float v3z, int mat_idx);
	//public static native int addSphere(int c_idx, int r_idx, int mat_idx);
	//public static native int addRadius(float r);
//	public static native int addVertex(int objIndex, float x, float y, float z);
	public static native int addMaterial(float r, float g, float b, float kd, float ks, float shine, float T, float ior);
	public static native int addPointLight(float x, float y, float z, float r, float g, float b);
	
	// Meshes and instances of meshes (objects).
	public static native int createMeshTemplate(int numVert, int numTri);
	public static native int createObject(int meshID);
//	public static native int addDynamicVertex(int dynObjIdx, float x, float y, float z);
//	public static native int addDynamicTriangle(int dynObjIdx, int v1_idx, int v2_idx, int v3_idx, int mat_idx);
	
//	public static native void enqueueApplyTransform(int dynObjIdx, float rot_x, float rot_y, float rot_z, float tr_x, float tr_y, float tr_z);
	
	// Native methods.
	private static native int provideBitmap(Bitmap bmp);
	public static native void provideKernelFile(String fileStr);
	
	public static void init(Context act) {
		InputStream is = null;
		try {
			is = act.getResources().getAssets().open("raytrace_kernel.cl");
		} catch (IOException e) {
			Log.d("FILE","File not loaded.");
		}
		
		String fileStr = "";
		if(is != null) {
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
			// Parse data.
			String line;
			
			try {
				while((line = br.readLine()) != null) {
					fileStr += line+"\n";
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		//Log.d("OPENCL", "Java: "+fileStr);
		
		provideKernelFile(fileStr);
		initialize();
	}
	
	private static native void initialize();
	public static native void finish();
	private static native int render(Bitmap bmp, float delta);
}
