#ifndef _RAY_TRACER_CONTEXT
#define _RAY_TRACER_CONTEXT

//#define ONE_THREAD
//#define CL
#define RAY_THREADS

//#include <vector>
//#define __CL_ENABLE_EXCEPTIONS
#define CL_USE_DEPRECATED_OPENCL_1_1_APIS
#define __NO_STD_VECTOR
#define __NO_STD_STRING

#include <android/bitmap.h>
/*#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>*/
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <CL/cl.hpp>
#include <pthread.h>


#include <Structures.h>
#include <Camera.h>
#include <BVH.h>
#include <ObjectBVH.h>

#define VERT_ARRAY_INIT_SIZE 100000
#define TRI_ARRAY_INIT_SIZE 100000

struct ThreadCellStartCoord {
	int start_x, end_x;
	int start_y, end_y;
};

class RayTracerContext {
public:
	AndroidBitmapInfo bitmapInfo;
	void* pixels;
	uint32_t* pixels_uint;
private:
	// camera data
	Camera cam;

	// Data structures
	BVH bvh;
	RayItem* rayStack;
	int numRayStack;
	int rayStack_size;

	// MODEL DATA //
	// Vertices
	int numVert;
	int vertArray_size;
	Vertex* vertArray;

	// Triangles
	int numTri;
	int triArray_size;
	Triangle* triArray;

	// Materials
	int numMat;
	int matArray_size;
	Material* matArray;

	// Point lights
	int numPointLights;
	int ptLightArray_size;
	PointLight* ptLightArray;

	// DYNAMIC OBJECT DATA
	DynamicObject* dynObjArray;
	int numDynObj;
	int dynObjArray_size;

	ObjectBVH objbvh;

	// parameters
	float ambient;
	int recursion_depth;
	const float rec_ray_offset;

	// test
	const float deltaAngle;
	float angle;

	// OpenCL
	CamData_CL camData_CL;
	cl_mem camData_cl_mem;
	//
    cl_mem vertArray_cl_mem;
    cl_mem triArray_cl_mem;
    cl_mem matArray_cl_mem;
    cl_mem ptLightArray_cl_mem;

	uint32_t* pixels_cl;

	const char* kernelFile;

	// OPENCL VARIABLES
	cl_int errorCode; ///< Error code variables used to catch errors.
	cl_platform_id* platformIDs; ///< Platform id.
	cl_device_id* deviceIDs; ///< Device id.
	cl_context gpuContext; ///< Context.
	cl_command_queue commandQueue; ///< Command queue for the used context.
	cl_program programObject; ///< Command queue for the used context.
	cl_kernel kernel_trace; ///< The ray tracing kernel.
	cl_mem pixelBuffer_mem; // Memory object for the pixel buffer

	// Work domain sizes.
	size_t global_work_size[2];
	size_t local_work_size[2];

	int packets_x, packets_y;
	int raypacket_x, raypacket_y;

	// Pthreads
#define NUM_THREADS 4
#define CELL_SIZE 8
#define PACKET_SIZE 8
#define NUM_RAYS_PACKET PACKET_SIZE*PACKET_SIZE

	pthread_t threads[NUM_THREADS];
	//pthread_args thread_args[NUM_THREADS];
	thread_range thread_ranges[NUM_THREADS];

	// Cell grid
	ThreadCellStartCoord* threadCellStartCoords[NUM_THREADS];
	int numThreadCells[NUM_THREADS];

    // Scene info
    AABB sceneBox;
    Vertex camPos;
    float camDist;
public:
	// OpenCL
	void printError(cl_int errorCode);
	void getDeviceInfo(const cl_device_id &deviceID);

    // Assign the source code of the OpenCL kernel file to the "kernelFile" variable.
	void provideKernelChar(const char* fileChar) {
		kernelFile = fileChar;
	}

	void initCL() {
		// Initialize OpenCL.
		cl_uint num_platforms;
		errorCode = clGetPlatformIDs(3, platformIDs, &num_platforms);
		printError(errorCode);
		LOGC("Platforms: %d", (int)num_platforms);

		cl_uint num_devices;
		errorCode = clGetDeviceIDs(platformIDs[0], CL_DEVICE_TYPE_GPU, 3, deviceIDs, &num_devices);
		printError(errorCode);
		LOGC("Devices: %d", (int)num_devices);
		// Get device info.
		for(int i=0; i<(int)num_devices; i++) {
			LOGC("----------------------");
			getDeviceInfo(deviceIDs[i]);
		}

		gpuContext = clCreateContext(NULL, 1, &deviceIDs[0], NULL, NULL, &errorCode);
		printError(errorCode);
		commandQueue = clCreateCommandQueue(gpuContext, deviceIDs[0], 0, &errorCode);
		printError(errorCode);



		size_t pixels_cl_size = sizeof(uint32_t) * cam.res_x * cam.res_y;
		LOGC("size of pixels_cl: %d", (int)pixels_cl_size);

		// Memory objects.

		// Pixel buffer memory object, the rendering kernel will write computed
		// colors here.
		pixelBuffer_mem = clCreateBuffer(gpuContext, CL_MEM_WRITE_ONLY, pixels_cl_size, NULL, &errorCode); printError(errorCode); printError(errorCode);
		// Camera data memory object.
		camData_cl_mem = clCreateBuffer(gpuContext, CL_MEM_READ_ONLY, sizeof(CamData_CL), NULL, &errorCode); printError(errorCode);
		// Model data memory objects.
		vertArray_cl_mem = clCreateBuffer(gpuContext, CL_MEM_READ_ONLY, numVert * sizeof(Vertex), NULL, &errorCode); printError(errorCode);
		triArray_cl_mem = clCreateBuffer(gpuContext, CL_MEM_READ_ONLY, numTri * sizeof(Triangle), NULL, &errorCode); printError(errorCode);
		matArray_cl_mem = clCreateBuffer(gpuContext, CL_MEM_READ_ONLY, numMat * sizeof(Material), NULL, &errorCode); printError(errorCode);
		ptLightArray_cl_mem = clCreateBuffer(gpuContext, CL_MEM_READ_ONLY, numPointLights * sizeof(PointLight), NULL, &errorCode); printError(errorCode);

		//LOGC("kernel file: %s", kernelFile);
		const size_t source_size[] = {(strlen(kernelFile) + 1) * sizeof(char)};
		LOGC("kernel char* size: %d", (strlen(kernelFile) + 1) * sizeof(char));


		programObject = clCreateProgramWithSource(gpuContext, 1, &kernelFile, source_size, &errorCode);
		LOGC("programObject:");
		printError(errorCode);
		errorCode = clBuildProgram(programObject, 1, deviceIDs, "-cl-denorms-are-zero -cl-no-signed-zeros -cl-fast-relaxed-math", NULL, NULL); // -I ./
		LOGC("buildProgram: %d", errorCode);
		printError(errorCode);

		kernel_trace = clCreateKernel(programObject, "traceRay", &errorCode);
		printError(errorCode);

		// Model data.
		errorCode = clEnqueueWriteBuffer(commandQueue, vertArray_cl_mem, CL_TRUE, 0, numVert * sizeof(Vertex), (void*)&vertArray, 0, NULL, NULL); printError(errorCode);
		errorCode = clEnqueueWriteBuffer(commandQueue, triArray_cl_mem, CL_TRUE, 0, numTri * sizeof(Triangle), (void*)&triArray, 0, NULL, NULL); printError(errorCode);
		errorCode = clEnqueueWriteBuffer(commandQueue, matArray_cl_mem, CL_TRUE, 0, numMat * sizeof(Material), (void*)&matArray, 0, NULL, NULL); printError(errorCode);
		errorCode = clEnqueueWriteBuffer(commandQueue, ptLightArray_cl_mem, CL_TRUE, 0, numPointLights * sizeof(PointLight), (void*)&ptLightArray, 0, NULL, NULL); printError(errorCode);





		size_t ker_pref_mult = 0;
		errorCode = clGetKernelWorkGroupInfo(kernel_trace, deviceIDs[0], CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(size_t), (void*)&ker_pref_mult, NULL);
		LOGC("CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE: %d", ker_pref_mult);
		printError(errorCode);
		size_t ker_wg_size = 0;
		errorCode = clGetKernelWorkGroupInfo(kernel_trace, deviceIDs[0], CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), (void*)&ker_wg_size, NULL);
		LOGC("CL_KERNEL_WORK_GROUP_SIZE: %d", ker_wg_size);
		printError(errorCode);

		raypacket_x = 8;
		raypacket_y = 8;

		// Set the global and local work size.
		int res_x = (int)cam.res_x;
		int res_y = (int)cam.res_y;

		packets_x = res_x / raypacket_x;
		packets_y = res_y / raypacket_y;

		LOGC("total packets: %d * %d", packets_x, packets_y);

		int global_x = packets_x, global_y = packets_y;

		int wg_x = ker_wg_size / 8;
		int wg_y = 8;

		LOGC("local_work_size: %d, %d", wg_x, wg_y);

		if(global_x % wg_x > 0) {
			global_x = (global_x / wg_x + 1) * wg_x;
		}
		if(global_y % wg_y > 0) {
			global_y = (global_y / wg_y + 1) * wg_y;
		}

		global_work_size[0] = global_x;
		global_work_size[1] = global_y;

		LOGC("global_work_size: %d, %d", global_x, global_y);

		local_work_size[0] = wg_x;
		local_work_size[1] = wg_y;
	}

	void initThreadsCells() {
		int res_x = cam.res_x;
		int res_y = cam.res_y;

		int block = NUM_THREADS * CELL_SIZE;

		int numCells_x;
		if(res_x % block > 0) {
			numCells_x = res_x / block + 1;
		} else {
			numCells_x = res_x / block;
		}
		int numCells_y;
		if(res_y % CELL_SIZE > 0) {
			numCells_y = res_y / CELL_SIZE + 1;
		} else {
			numCells_y = res_y / CELL_SIZE;
		}

		int numCells_total = numCells_x * numCells_y;

		// Prepare all the cell coordinates for all the threads.
		// The pattern for 4 threads would be like this:
		// 123412351234...
		// 341234123412...
		// 123412341234...
		// where each number represents a square of CELL_SIZE * CELL_SIZE rays.
		for(int k=0; k<NUM_THREADS; k++) {
			int numCells = 0;
			threadCellStartCoords[k] = new ThreadCellStartCoord[numCells_total];

			int block_2 = block / 2;
			for(int i=-1; i<numCells_x; i++) {
				for(int j=0; j<numCells_y; j++) {
					int start_x = i * block + k * CELL_SIZE + block_2 * (j % 2); // offset the odd lines by half block.
					int start_y = j * CELL_SIZE;

					if(start_x >= res_x || start_y >= res_y || start_x < 0) {
						continue;
					} else {
						int end_x = fmin(start_x + CELL_SIZE, res_x);
						int end_y = fmin(start_y + CELL_SIZE, res_y);

						threadCellStartCoords[k][numCells].start_x = start_x;
						threadCellStartCoords[k][numCells].start_y = start_y;
						threadCellStartCoords[k][numCells].end_x = end_x;
						threadCellStartCoords[k][numCells].end_y = end_y;
						numCells++;

						//LOGI("New cell: %d, %d x %d, %d", start_x, end_x, start_y, end_y);
					}
				}
				//LOGI("----------------------------");
			}
			numThreadCells[k] = numCells;
			//LOGI("----------------------------");
			//LOGI("----------------------------");
		}
	}

	void initThreads() {
		int res_x = cam.res_x;
		int res_y = cam.res_y;

		int x_parts, y_parts;

		switch(NUM_THREADS) {
			case 1: {
				x_parts = 1;
				y_parts = 1;
			} break;
			case 2: {
				x_parts = 2;
				y_parts = 1;
			} break;
			case 4: {
				x_parts = 2;
				y_parts = 2;
			} break;
			case 8: {
				x_parts = 4;
				y_parts = 2;
			} break;
			default: break;
		}

		int res_x_part = res_x / x_parts;
		int res_y_part = res_y / y_parts;

		for(int i=0; i<x_parts; i++) {
			for(int j=0; j<y_parts; j++) {
				int idx = i + j * x_parts;

				LOGI("Init thread args: %d", idx);

				thread_ranges[idx].start_x = i * res_x_part;
				thread_ranges[idx].end_x = (i + 1) * res_x_part;
				thread_ranges[idx].start_y = j * res_y_part;
				thread_ranges[idx].end_y = (j + 1) * res_y_part;

				LOGI("Init thread args: %d, %d x %d, %d", thread_ranges[idx].start_x, thread_ranges[idx].end_x, thread_ranges[idx].start_y, thread_ranges[idx].end_y);
			}
		}
	}

	void initialize() {
	    LOGI("Vertices: %d", numVert);
	    LOGI("Triangles: %d", numTri);

	    // Create the BVH for the current scene.
		bvh.createBVH(vertArray, numVert, triArray, numTri);

		for(int i=0; i<numDynObj; i++) {
            dynObjArray[i].initBuildTree();
		}

		objbvh.createBVH(dynObjArray, numDynObj);

		// Init camera parameters.
        camPos = sceneBox.getCenter();
        camDist = fmax(sceneBox.getXsize(), fmax(sceneBox.getYsize(), sceneBox.getZsize())) * 1.2f;

#ifdef CL
		initCL();
#endif

#ifdef RAY_THREADS
		initThreads();
		initThreadsCells();
		int x = 10;
#endif
	}



	RayTracerContext() : deltaAngle(0.8f), angle(0.0f), numPointLights(2), ambient(0.0f), recursion_depth(2), rec_ray_offset(1e-5f) {
		rayStack_size = 0x7FFFFFFF >> (31 - recursion_depth);
		rayStack = new RayItem[rayStack_size];

		// Vertex array
		vertArray = new Vertex[VERT_ARRAY_INIT_SIZE];
		vertArray_size = VERT_ARRAY_INIT_SIZE;
		numVert = 0;

		// Triangle array
		triArray = new Triangle[TRI_ARRAY_INIT_SIZE];
		triArray_size = TRI_ARRAY_INIT_SIZE;
		numTri = 0;

		matArray = new Material[10];
		matArray_size = 10;
		numMat = 1;
		// Fill in the default material.
		matArray[0].r = DEF_R;
		matArray[0].g = DEF_G;
		matArray[0].b = DEF_B;
		matArray[0].kd = DEF_KD;
		matArray[0].ks = DEF_KS;
		matArray[0].shine = DEF_SHINE;
		matArray[0].T = DEF_T;
		matArray[0].ior = DEF_IOR;

		ptLightArray = new PointLight[1];
		ptLightArray_size = 1;
		numPointLights = 0;

		// Dynamic object array.
        dynObjArray = new DynamicObject[2];
        dynObjArray_size = 2;
        numDynObj = 0;

		// Setup the camera.
		Vertex cam_pos = {0.0f, 0.0f, -5.0f};
		Vertex cam_dir = {0.0f, 0.0f, 1.0f};
		Vertex cam_up = {0.0f, 1.0f, 0.0f};
		cam.setup(960, 540, 1.85f, cam_pos, cam_dir, cam_up);

		// OpenCL
		// OPENCL VARIABLES
		platformIDs = new cl_platform_id[3]; ///< Platform id.
		deviceIDs = new cl_device_id[3]; ///< Device id.
		gpuContext = NULL; ///< Context.
		commandQueue = NULL; ///< Command queue for the used context.
		programObject = NULL; ///< Command queue for the used context.
		kernel_trace = NULL; ///< The ray tracing kernel.
		pixelBuffer_mem = NULL; // Memory object for the pixel buffer
	}

	~RayTracerContext();

	void setProjectionResolution(const int &_res_x, const int &_res_y);

	int addPointLight(const float &_x, const float &_y, const float &_z,
					  const float &_r, const float &_g, const float &_b);
	int addVertex(const float &_x, const float &_y, const float &_z);
	int addTriangle(const int &v1, const int &v2, const int &v3, const int &mat_idx);
	int addMaterial(const float &r, const float &g, const float &b,
					const float &kd, const float &ks, const float &shine,
					const float &_T, const float &ior);

    int createDynamicObject(const int &_numVert, const int &_numTri) {
        if(dynObjArray_size - numDynObj == 0) { // Reallocate the array.
            dynObjArray_size <<= 1;
            DynamicObject* new_dynObjArray = new DynamicObject[dynObjArray_size];
            memcpy(new_dynObjArray, dynObjArray, numDynObj * sizeof(DynamicObject));
            delete[] dynObjArray;
            dynObjArray = new_dynObjArray;
        }

        dynObjArray[numDynObj].initialize(_numVert, _numTri);

        numDynObj++;

        return numDynObj - 1;
    }

    /// Add a vertex to the dynamic object with index dynObjIndex.
    int addDynamicVertex(const int &dynObjIndex, const float &_x, const float &_y, const float &_z) {
        if(dynObjIndex >= 0 && dynObjIndex < numDynObj) {
            return dynObjArray[dynObjIndex].addVertex(_x, _y, _z);
        }
        return -1;
    }

    /// Add a triangle to the dynamic object with index dynObjIndex.
    /**
        The vertex indices point to the private vertex array of the object.
        The material index points to the global material array of the context.
    */
    int addDynamicTriangle(const int &dynObjIndex, const int &v1, const int &v2, const int &v3, const int &mat_idx) {
        if(dynObjIndex >= 0 && dynObjIndex < numDynObj) {
            return dynObjArray[dynObjIndex].addTriangle(v1, v2, v3, mat_idx);
        }
        return -1;
    }

	void updateTriangleCenter(const int &tri_idx);

    /// Render using one CPU thread
	void rayTraceScene(const float &delta);
	///

	/// Render using OpenCL
	void rayTraceSceneCL(const float &delta);
    ///

    /// Render using multiple CPU threads
    // slave thread functions
	void* rayTraceSceneThread(void* screenRange);
	void* rayTraceScenePacketsThread(void* thread_id_context);

    // master function
	void rayTraceSceneThreads(const float &delta);
	///

    /// Main function called when a scene is to be rendered
	void rayTrace_entry(const float &delta) {
#ifdef CL
		rayTraceSceneCL(delta);
#endif
#ifdef RAY_THREADS
		rayTraceSceneThreads(delta);
#endif
#ifdef ONE_THREAD
		rayTraceScene(delta);
#endif
	}

	void finish() {
#ifdef CL
	clReleaseKernel(kernel_trace);
	clReleaseProgram(programObject);
	clReleaseCommandQueue(commandQueue);
	clReleaseContext(gpuContext);
	clReleaseMemObject(camData_cl_mem);
	clReleaseMemObject(vertArray_cl_mem);
	clReleaseMemObject(triArray_cl_mem);
	clReleaseMemObject(matArray_cl_mem);
	clReleaseMemObject(ptLightArray_cl_mem);
	clReleaseMemObject(pixelBuffer_mem);
#endif

		//pthread_exit(NULL);
	}
};

struct pthread_args {
	int id;
	RayTracerContext* context;
};

#endif
