#include <RayTracerContext.h>

//#define USE_RAY_PACKETS
//#define SECONDARY_RAYS
//#define LIGHTS

RayTracerContext::~RayTracerContext() {
//	delete[] rayStack;

//	delete[] vertArray;
//	delete[] triArray;
	delete[] matArray;

	delete[] pixels_cl;

	for(int i=0; i<NUM_THREADS; i++) {
		delete[] threadCellStartCoords[i];
	}
}


// Initialization
void RayTracerContext::setProjectionResolution(const int &_res_x, const int &_res_y) {
	cam.changeResolution(_res_x, _res_y);

	LOGI("Changing resolution to: %d, %d.", _res_x, _res_y);

	pixels_cl = new uint32_t[_res_x * _res_y];
}

int RayTracerContext::addPointLight(const float &_x, const float &_y, const float &_z,
					   const float &_r, const float &_g, const float &_b) {
	if(ptLightArray_size - numPointLights == 0) {
		//LOGI("REALLOC: PointLight array.");

		ptLightArray_size <<= 1; // double the size
		PointLight* new_ptLightArray = new PointLight[ptLightArray_size];
		memcpy(new_ptLightArray, ptLightArray, numPointLights * sizeof(PointLight));

		delete[] ptLightArray;
		ptLightArray = new_ptLightArray;
	}
	ptLightArray[numPointLights].pos.x = _x;
	ptLightArray[numPointLights].pos.y = _y;
	ptLightArray[numPointLights].pos.z = _z;
	ptLightArray[numPointLights].r = _r;
	ptLightArray[numPointLights].g = _g;
	ptLightArray[numPointLights].b = _b;
	numPointLights++;

	//LOGI("PointLight added.");

	return numPointLights - 1;
}

int RayTracerContext::addTriangle(const int &meshID, const float &v1x, const float &v1y, const float &v1z,
                    const float &v2x, const float &v2y, const float &v2z,
                    const float &v3x, const float &v3y, const float &v3z, const int &mat_idx) {
    if(mat_idx >= numMat || meshID >= numMesh) return -1;
    meshArray[meshID].addTriangle(v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z, mat_idx);
}

int RayTracerContext::addMaterial(const float &r, const float &g, const float &b,
					const float &kd, const float &ks, const float &shine,
					const float &_T, const float &ior) {
	if(matArray_size - numMat == 0) { // Reallocate the array.
		//LOGI("REALLOC: Material array.");

		matArray_size <<= 1;
		Material* new_matArray = new Material[matArray_size];
		memcpy(new_matArray, matArray, numMat * sizeof(Material));
		/*for(int i=0; i<MAT_COMPONENTS*numMat; i++) {
			new_matArray[i] = matArray[i];
		}*/
		delete[] matArray;
		matArray = new_matArray;
	}
	//int index = MAT_COMPONENTS*numMat;
	matArray[numMat].r = r;
	matArray[numMat].g = g;
	matArray[numMat].b = b;
	matArray[numMat].kd = kd;
	matArray[numMat].ks = ks;
	matArray[numMat].shine = shine;
	matArray[numMat].T = _T;
	matArray[numMat].ior = ior;
	numMat++;

	//LOGI("Material added.");

	return numMat - 1;
}

// CL
void RayTracerContext::printError(cl_int errorCode) {
	switch (errorCode) {
	case CL_SUCCESS:                            LOGC("Success!"); break;
	/*case CL_DEVICE_NOT_FOUND:                   cout << "Device not found." << endl; break;
	case CL_DEVICE_NOT_AVAILABLE:               cout << "Device not available" << endl; break;
	case CL_COMPILER_NOT_AVAILABLE:             cout << "Compiler not available" << endl; break;
	case CL_MEM_OBJECT_ALLOCATION_FAILURE:      cout << "Memory object allocation failure" << endl; break;
	case CL_OUT_OF_RESOURCES:                   cout << "Out of resources" << endl; break;
	case CL_OUT_OF_HOST_MEMORY:                 cout << "Out of host memory" << endl; break;
	case CL_PROFILING_INFO_NOT_AVAILABLE:       cout << "Profiling information not available" << endl; break;
	case CL_MEM_COPY_OVERLAP:                   cout << "Memory copy overlap" << endl; break;
	case CL_IMAGE_FORMAT_MISMATCH:              cout << "Image format mismatch" << endl; break;
	case CL_IMAGE_FORMAT_NOT_SUPPORTED:         cout << "Image format not supported" << endl; break;
	case CL_BUILD_PROGRAM_FAILURE:              cout << "Program build failure" << endl; break;
	case CL_MAP_FAILURE:                        cout << "Map failure" << endl; break;
	case CL_INVALID_VALUE:                      cout << "Invalid value" << endl; break;
	case CL_INVALID_DEVICE_TYPE:                cout << "Invalid device type" << endl; break;
	case CL_INVALID_PLATFORM:                   cout << "Invalid platform" << endl; break;
	case CL_INVALID_DEVICE:                     cout << "Invalid device" << endl; break;
	case CL_INVALID_CONTEXT:                    cout << "Invalid context" << endl; break;
	case CL_INVALID_QUEUE_PROPERTIES:           cout << "Invalid queue properties" << endl; break;
	case CL_INVALID_COMMAND_QUEUE:              cout << "Invalid command queue" << endl; break;
	case CL_INVALID_HOST_PTR:                   cout << "Invalid host pointer" << endl; break;
	case CL_INVALID_MEM_OBJECT:                 cout << "Invalid memory object" << endl; break;
	case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:    cout << "Invalid image format descriptor" << endl; break;
	case CL_INVALID_IMAGE_SIZE:                 cout << "Invalid image size" << endl; break;
	case CL_INVALID_SAMPLER:                    cout << "Invalid sampler" << endl; break;
	case CL_INVALID_BINARY:                     cout << "Invalid binary" << endl; break;
	case CL_INVALID_BUILD_OPTIONS:              cout << "Invalid build options" << endl; break;
	case CL_INVALID_PROGRAM:                    cout << "Invalid program" << endl; break;
	case CL_INVALID_PROGRAM_EXECUTABLE:         cout << "Invalid program executable" << endl; break;
	case CL_INVALID_KERNEL_NAME:                cout << "Invalid kernel name" << endl; break;
	case CL_INVALID_KERNEL_DEFINITION:          cout << "Invalid kernel definition" << endl; break;
	case CL_INVALID_KERNEL:                     cout << "Invalid kernel" << endl; break;
	case CL_INVALID_ARG_INDEX:                  cout << "Invalid argument index" << endl; break;
	case CL_INVALID_ARG_VALUE:                  cout << "Invalid argument value" << endl; break;
	case CL_INVALID_ARG_SIZE:                   cout << "Invalid argument size" << endl; break;
	case CL_INVALID_KERNEL_ARGS:                cout << "Invalid kernel arguments" << endl; break;
	case CL_INVALID_WORK_DIMENSION:             cout << "Invalid work dimension" << endl; break;
	case CL_INVALID_WORK_GROUP_SIZE:            cout << "Invalid work group size" << endl; break;
	case CL_INVALID_WORK_ITEM_SIZE:             cout << "Invalid work item size" << endl; break;
	case CL_INVALID_GLOBAL_OFFSET:              cout << "Invalid global offset" << endl; break;
	case CL_INVALID_EVENT_WAIT_LIST:            cout << "Invalid event wait list" << endl; break;
	case CL_INVALID_EVENT:                      cout << "Invalid event" << endl; break;
	case CL_INVALID_OPERATION:                  cout << "Invalid operation" << endl; break;
	case CL_INVALID_GL_OBJECT:                  cout << "Invalid OpenGL object" << endl; break;
	case CL_INVALID_BUFFER_SIZE:                cout << "Invalid buffer size" << endl; break;
	case CL_INVALID_MIP_LEVEL:                  cout << "Invalid mip-map level" << endl; break;
	default: cout << "Unknown"; break;*/
	default: LOGC("Unknown error type."); break;
	}
}

void RayTracerContext::getDeviceInfo(const cl_device_id &deviceID) {
	size_t size;

	// Vendor.
	char* vendor;
	clGetDeviceInfo(deviceID, CL_DEVICE_VENDOR, 0, NULL, &size);
	vendor = (char*)malloc(sizeof(char)*size);
	clGetDeviceInfo(deviceID, CL_DEVICE_VENDOR, size, vendor, NULL);
	LOGC("Device vendor: %s", vendor);
	free(vendor);

	// Device type.
	cl_device_type* device_type;
	clGetDeviceInfo(deviceID, CL_DEVICE_TYPE, 0, NULL, &size);
	device_type = (cl_device_type*)malloc(sizeof(cl_device_type)*size);
	clGetDeviceInfo(deviceID, CL_DEVICE_TYPE, size, device_type, NULL);
	if (*device_type & CL_DEVICE_TYPE_GPU) {
		LOGC("Device is a GPU.");
	}
	else if (*device_type & CL_DEVICE_TYPE_CPU) {
		LOGC("Device is a CPU.");
	}
	free(device_type);

	// Supported OpenCL version.
	char* device_version;
	clGetDeviceInfo(deviceID, CL_DEVICE_VERSION, 0, NULL, &size);
	device_version = (char*)malloc(sizeof(char)*size);
	clGetDeviceInfo(deviceID, CL_DEVICE_VERSION, size, device_version, NULL);
	LOGC("Device OpenCL version: %s", device_version);
	free(device_version);

	// Command queue properties.
	cl_command_queue_properties* c_c_properties;
	clGetDeviceInfo(deviceID, CL_DEVICE_QUEUE_PROPERTIES, 0, NULL, &size);
	c_c_properties = (cl_command_queue_properties*)malloc(sizeof(cl_command_queue_properties)*size);
	clGetDeviceInfo(deviceID, CL_DEVICE_QUEUE_PROPERTIES, size, c_c_properties, NULL);
	if (*c_c_properties & CL_QUEUE_PROFILING_ENABLE) {
		LOGC("Profiling enable.");
	}
	if (*c_c_properties & CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE) {
		LOGC("Out of order exec mode enable.");
	}
	free(c_c_properties);

	// Timer resolution.
	size_t* timer_res;
	clGetDeviceInfo(deviceID, CL_DEVICE_PROFILING_TIMER_RESOLUTION, 0, NULL, &size);
	timer_res = (size_t*)malloc(sizeof(size_t)*size);
	clGetDeviceInfo(deviceID, CL_DEVICE_PROFILING_TIMER_RESOLUTION, size, timer_res, NULL);
	LOGC("Timer resolution [ns]: %d", *timer_res);
	free(timer_res);

	// Device profile.
	char* device_profile;
	clGetDeviceInfo(deviceID, CL_DEVICE_PROFILE, 0, NULL, &size);
	device_profile = (char*)malloc(sizeof(char)*size);
	clGetDeviceInfo(deviceID, CL_DEVICE_PROFILE, size, device_profile, NULL);
	LOGC("Device OpenCL profile: %s", device_profile);
	free(device_profile);

	// Work-group dimensions.
	cl_uint* dim;
	clGetDeviceInfo(deviceID, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, 0, NULL, &size);
	dim = (cl_uint*)malloc(sizeof(cl_uint)*size);
	clGetDeviceInfo(deviceID, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, size, dim, NULL);
	LOGC("Work-group dimensions: %d", *dim);

	// Sizes of work-group dimensions.
	size_t* sizes;
	clGetDeviceInfo(deviceID, CL_DEVICE_MAX_WORK_ITEM_SIZES, 0, NULL, &size);
	sizes = (size_t*)malloc(sizeof(size_t)*size);
	clGetDeviceInfo(deviceID, CL_DEVICE_MAX_WORK_ITEM_SIZES, size, sizes, NULL);
	for (int i = 0; i < *dim; i++) {
		LOGC("Dimension %d size: %d", i, sizes[i]);
		//cout << "Dimension " << i << " size: " << sizes[i] << endl;
	}
	free(sizes);
	free(dim);

	// Work-group size.
	size_t* w_size;
	clGetDeviceInfo(deviceID, CL_DEVICE_MAX_WORK_GROUP_SIZE, 0, NULL, &size);
	w_size = (size_t*)malloc(sizeof(size_t)*size);
	clGetDeviceInfo(deviceID, CL_DEVICE_MAX_WORK_GROUP_SIZE, size, w_size, NULL);
	LOGC("Work group size: %d", *w_size);
	free(w_size);

	// Memory object max size.
	cl_ulong* mem_obj_size;
	clGetDeviceInfo(deviceID, CL_DEVICE_MAX_MEM_ALLOC_SIZE, 0, NULL, &size);
	mem_obj_size = (cl_ulong*)malloc(sizeof(cl_ulong)*size);
	clGetDeviceInfo(deviceID, CL_DEVICE_MAX_MEM_ALLOC_SIZE, size, mem_obj_size, NULL);
	LOGC("Memory object size [B]: %d", (int)(*mem_obj_size / 1000000));
	free(mem_obj_size);

	// Compute units.
	cl_uint* compute_units;
	clGetDeviceInfo(deviceID, CL_DEVICE_MAX_COMPUTE_UNITS, 0, NULL, &size);
	compute_units = (cl_uint*)malloc(sizeof(cl_uint)*size);
	clGetDeviceInfo(deviceID, CL_DEVICE_MAX_COMPUTE_UNITS, size, compute_units, NULL);
	LOGC("Max compute units: %d", *compute_units);
	free(compute_units);

	// Device frequency.
	cl_uint* freq;
	clGetDeviceInfo(deviceID, CL_DEVICE_MAX_CLOCK_FREQUENCY, 0, NULL, &size);
	freq = (cl_uint*)malloc(sizeof(cl_uint)*size);
	clGetDeviceInfo(deviceID, CL_DEVICE_MAX_CLOCK_FREQUENCY, size, freq, NULL);
	LOGC("Frequency [Mhz]: %d", *freq);
	free(freq);

	// Global memory size.
	cl_ulong* global_mem;
	clGetDeviceInfo(deviceID, CL_DEVICE_GLOBAL_MEM_SIZE, 0, NULL, &size);
	global_mem = (cl_ulong*)malloc(sizeof(cl_ulong)*size);
	clGetDeviceInfo(deviceID, CL_DEVICE_GLOBAL_MEM_SIZE, size, global_mem, NULL);
	LOGC("Global memory [B]: %d", (int)(*global_mem / 1000000));
	free(global_mem);
}

void RayTracerContext::rayTrace_thread_grid(void* thread_id_context) {
//	pthread_args* thread_args = (pthread_args*)thread_id_context;
//	int id = thread_args->id;
//
//	//thread_range* screen_range = &thread_ranges[id];
//
//    int nearestObj_cam, nearestObj_light;
//	int nearestStaticTri_cam, nearestDynamicTri_cam;
//	int nearestStaticTri_light, nearestDynamicTri_light;
//	float t_cam, t_light, t_cam_dyn;
//
//	int nearestToLightIdx, nearestToCamIdx;
//
//	Vertex pixel_w;
//	Vertex col_step_back = -((float)(CELL_SIZE) * cam.y_step);
//	Vertex contactPoint, ray, inv_ray, lightRay, inv_lightRay, surfaceNormal, proj_surfaceNormal, invertedNormal;
//
//	Vertex reflectedRay, refractedRay;
//
//	int numRayStack_thread;
//	RayItem rayStack_thread[rayStack_size]; // Each thread needs its own stack.
//	RayItem current_rayItem;
//
//	float r, g, b, r_cur, g_cur, b_cur;
//
//	int bmp_rows;
//
//	// Acquire the pixel array;
//	uint32_t* pixels_uint_thread = (uint32_t*) pixels;
//
//	for(int k=0; k<numThreadCells[id]; k++) {
//		int start_x = threadCellStartCoords[id][k].start_x;
//		int start_y = threadCellStartCoords[id][k].start_y;
//		int end_x = threadCellStartCoords[id][k].end_x;
//		int end_y = threadCellStartCoords[id][k].end_y;
//
//		pixel_w = (cam.proj_origin - cam.pos) + (cam.x_step * start_x) + (cam.y_step * start_y);
//		// Trace rays for all pixels in <start_x, end_x) x <start_y, end_y).
//		for(int i=start_x; i<end_x; i++) {
//			bmp_rows = i + (start_y * bitmapInfo.width);
//			for(int j=start_y; j<end_y; j++) {
//				r = 0.0f;
//				g = 0.0f;
//				b = 0.0f;
//
//				numRayStack_thread = 0;
//				// The first primary ray.
//				//pixel_w = cam.proj_origin + i*cam.x_step + j*cam.y_step - cam.pos;
//
//				memcpy(&ray, &pixel_w, sizeof(Vertex));
//				//ray -= cam.pos;
//
//				rayStack_thread[numRayStack_thread].dir.x = ray.x;
//				rayStack_thread[numRayStack_thread].dir.y = ray.y;
//				rayStack_thread[numRayStack_thread].dir.z = ray.z;
//				rayStack_thread[numRayStack_thread].orig.x = cam.pos.x;
//				rayStack_thread[numRayStack_thread].orig.y = cam.pos.y;
//				rayStack_thread[numRayStack_thread].orig.z = cam.pos.z;
//				rayStack_thread[numRayStack_thread].intensity_coef = 1.0f;
//				rayStack_thread[numRayStack_thread].recursion_number = 0;
//
//				numRayStack_thread++;
//
//				while(numRayStack_thread > 0) {
//					current_rayItem = rayStack_thread[numRayStack_thread-1];
//					numRayStack_thread--;
//
//					memcpy(&ray, &(current_rayItem.dir), sizeof(Vertex));
//
//					normalize(ray);
//
//					inv_ray.x = 1.0f / ray.x;
//					inv_ray.y = 1.0f / ray.y;
//					inv_ray.z = 1.0f / ray.z;
//
//					nearestStaticTri_cam = -1;
//					nearestToLightIdx = -2;
//
////                    nearestToCamObjIdx = -1;
//
////					t_cam = bvh.findNearestPrimitive_new(ray.x, ray.y, ray.z,
////														 inv_ray.x, inv_ray.y, inv_ray.z,
////														 current_rayItem.orig.x, current_rayItem.orig.y, current_rayItem.orig.z,
////														 &nearestStaticTri_cam);
//
////                    t_cam_dyn = obvh.findNearestObject(ray.x, ray.y, ray.z,
////                                                       inv_ray.x, inv_ray.y, inv_ray.z,
////                                                       current_rayItem.orig.x, current_rayItem.orig.y, current_rayItem.orig.z,
////                                                       &nearestObj_cam, &nearestDynamicTri_cam);
//
//
//					if(nearestStaticTri_cam >= 0) {
//						//LOGI("I see a triangle!");
//
//						contactPoint.x = current_rayItem.orig.x + t_cam * ray.x;
//						contactPoint.y = current_rayItem.orig.y + t_cam * ray.y;
//						contactPoint.z = current_rayItem.orig.z + t_cam * ray.z;
//
//						Material* mat = NULL;
////						if(triArray[nearestStaticTri_cam].n_calculated) {
////							memcpy(&surfaceNormal, &(triArray[nearestStaticTri_cam].n), sizeof(Vertex));
////						} else {
////							// update triangle normal, don't forget to divide by inv_norm
////						}
//
//						// Ambient
//						r_cur = ambient * mat->r;
//						g_cur = ambient * mat->g;
//						b_cur = ambient * mat->b;
//	#ifdef LIGHTS
//						// Find all the lights that light up the nearest primitive.
//						for(int k=0; k<numPointLights; k++) {
//							memcpy(&lightRay, &contactPoint, sizeof(Vertex));
//							lightRay -= ptLightArray[k].pos;
//							//lightRay = contactPoint - Vec3(pointLightArray[k*6], pointLightArray[k*6 + 1], pointLightArray[k*6 + 2]);
//
////							t_light = bvh.findNearestPrimitive(lightRay.x, lightRay.y, lightRay.z, ptLightArray[k].pos.x, ptLightArray[k].pos.y, ptLightArray[k].pos.z, nearestToLightIdx);
//							//t_light = findNearestPrimitive(lightRay.x, lightRay.y, lightRay.z, pointLightArray[k*6], pointLightArray[k*6 + 1], pointLightArray[k*6 + 2], nearestToLightIdx);
//
//							if(nearestStaticTri_cam == nearestToLightIdx) { // && type_cam == type_light
//								float lightRay_inv_norm = inv_norm(lightRay);
//								float proj = dot(surfaceNormal, lightRay);
//
//								// Diffuse
//								float diffuse_cos = -proj * lightRay_inv_norm;
//								if(diffuse_cos > 0.0f) { // light_cos > 0.0f
//									float coef = mat->kd * diffuse_cos;
//									r_cur += coef * ptLightArray[k].r * mat->r;
//									g_cur += coef * ptLightArray[k].g * mat->g;
//									b_cur += coef * ptLightArray[k].b * mat->b;
//								}
//
//								// Specular (reflected dot toObserver)
//								memcpy(&proj_surfaceNormal, &surfaceNormal, sizeof(Vertex));
//								proj_surfaceNormal *= (proj + proj);
//
//								lightRay -= proj_surfaceNormal;
//
//								float spec_cos = cosine(lightRay, ray);
//								if(spec_cos <= 0.0f) { // This might be stupid and useless.
//									float coef = mat->ks * pow(spec_cos, mat->shine);
//									r_cur += coef * ptLightArray[k].r;
//									g_cur += coef * ptLightArray[k].g;
//									b_cur += coef * ptLightArray[k].b;
//								}
//							}
//						}
//	#endif
//						r_cur *= current_rayItem.intensity_coef;
//						g_cur *= current_rayItem.intensity_coef;
//						b_cur *= current_rayItem.intensity_coef;
//
//	#ifdef SECONDARY_RAYS
//						if(current_rayItem.recursion_number < recursion_depth) {
//							//float ray_inv_norm = inv_norm(ray);
//							float proj = dot(surfaceNormal, ray);
//
//							// Reflected ray
//							if(proj < 0.0f && mat->ks > 0.0f) { // can determine earlier
//								memcpy(&reflectedRay, &ray, sizeof(Vertex));
//								memcpy(&proj_surfaceNormal, &surfaceNormal, sizeof(Vertex));
//
//								proj_surfaceNormal *= (proj + proj);
//								reflectedRay -= proj_surfaceNormal;
//
//								normalize(reflectedRay);
//								//reflectedRay = Vec3::normalize(ray + 2.0f * Vec3::dot(-1.0f*ray, surfaceNormal) * surfaceNormal);
//
//								rayStack_thread[numRayStack_thread].dir = reflectedRay;
//								rayStack_thread[numRayStack_thread].orig.x = reflectedRay.x * rec_ray_offset + contactPoint.x;
//								rayStack_thread[numRayStack_thread].orig.y = reflectedRay.y * rec_ray_offset + contactPoint.y;
//								rayStack_thread[numRayStack_thread].orig.z = reflectedRay.z * rec_ray_offset + contactPoint.z;
//								rayStack_thread[numRayStack_thread].intensity_coef = mat->ks * current_rayItem.intensity_coef;
//								rayStack_thread[numRayStack_thread].recursion_number = current_rayItem.recursion_number + 1;
//
//								numRayStack_thread++;
//							}
//
//							// Refracted ray
//							if(mat->T > 0.0f) {
//
//								memcpy(&invertedNormal, &surfaceNormal, sizeof(Vertex));
//
//								//Vec3 normalizedRay = Vec3::normalize(ray), refractedRay, invertedNormal;
//
//								float cos_theta_1 = cosine1(ray, surfaceNormal), ior, sqrterm;
//
//								if(cos_theta_1 < 0.0f) {
//									ior = IOR_AIR / mat->ior;
//									//invertedNormal = surfaceNormal;
//								} else {
//									ior = mat->ior / IOR_AIR;
//									cos_theta_1 *= -1.0f;
//									//invertedNormal = -1.0f*surfaceNormal;
//									invertedNormal *= -1.0f;
//								}
//								sqrterm = 1.0f - ior*ior*(1.0f - cos_theta_1*cos_theta_1);
//
//								if(sqrterm > 0.0f) {
//									sqrterm = cos_theta_1*ior + sqrt(sqrterm);
//
//									memcpy(&refractedRay, &ray, sizeof(Vertex));
//									refractedRay *= ior;
//									refractedRay += invertedNormal * (-sqrterm);
//
//									//refractedRay = invertedNormal * (-sqrterm) + normalizedRay * ior;
//								}
//
//								//Vec3 refractedRay = Vec3::normalize(ray);
//								rayStack_thread[numRayStack_thread].dir = refractedRay;
//								rayStack_thread[numRayStack_thread].orig.x = refractedRay.x * rec_ray_offset + contactPoint.x;
//								rayStack_thread[numRayStack_thread].orig.y = refractedRay.y * rec_ray_offset + contactPoint.y;
//								rayStack_thread[numRayStack_thread].orig.z = refractedRay.z * rec_ray_offset + contactPoint.z;
//								rayStack_thread[numRayStack_thread].intensity_coef = mat->T * current_rayItem.intensity_coef;
//								rayStack_thread[numRayStack_thread].recursion_number = current_rayItem.recursion_number + 1;
//
//								numRayStack_thread++;
//							}
//						}
//	#endif
//						r += r_cur;
//						g += g_cur;
//						b += b_cur;
//					} else {
//						r += 0.2f * current_rayItem.intensity_coef;
//						g += 0.2f * current_rayItem.intensity_coef;
//						b += 0.2f * current_rayItem.intensity_coef;
//					}
//				}
//
//				// Limit the color components to 1.0f.
//				r = r > 1.0f ? 1.0f : r;
//				g = g > 1.0f ? 1.0f : g;
//				b = b > 1.0f ? 1.0f : b;
//
//                // Set the color in the color buffer.
//                pixels_uint_thread[i + (j)*bitmapInfo.width] = argb_to_int(255, (int)(255.0f * r), (int)(255.0f * g), (int)(255.0f * b));
//
//				bmp_rows += bitmapInfo.width;
//
//				pixel_w += cam.y_step;
//			} // End of j
//			pixel_w += col_step_back;
//			pixel_w += cam.x_step;
//		} // End of i
//	}
}

void RayTracerContext::rayTrace_thread_horSlabs(void* thread_id_context) {
	pthread_args* thread_args = (pthread_args*)thread_id_context;
	int id = thread_args->id;

	//thread_range* screen_range = &thread_ranges[id];

    int nearestObj_cam, nearestObj_light;
    int nearestTri_cam, nearestTri_light;
	float t_cam, t_light;

	int nearestToLightIdx, nearestToCamIdx;

	Vertex pixel_w;
	Vertex row_step_back = -((float)(bitmapInfo.width) * cam.x_step);
	Vertex pixel_skip_y = ((float)((NUM_THREADS - 1) * LINES_PER_SLAB) * cam.y_step);
	Vertex contactPoint, ray, inv_ray, lightRay, inv_lightRay, surfaceNormal, proj_surfaceNormal, invertedNormal;

	Ray primaryRay, curRay;
	primaryRay.orig = cam.pos;

	Vertex reflectedRay, refractedRay;

	int numRayStack_thread;
	RayItem rayStack_thread[rayStack_size]; // Each thread needs its own stack.
	RayItem current_rayItem;

	float r, g, b, r_cur, g_cur, b_cur;

	int color_buffer_index;

	// Acquire the pixel array;
	uint32_t* pixels_uint_thread = (uint32_t*) pixels;

	// Find the first pixel for this thread.
	pixel_w = (cam.proj_origin - cam.pos) + (cam.y_step * (float)(id * LINES_PER_SLAB));

    // All the horizontal slabs this thread has to render.
	for(unsigned int k = id*LINES_PER_SLAB; k<bitmapInfo.height; k += NUM_THREADS*LINES_PER_SLAB) {
        color_buffer_index = k * bitmapInfo.width;
        // All the image rows in this slab.
        unsigned int j_end = k + LINES_PER_SLAB;
		for(unsigned int j=k; j<j_end; j++) {
			for(unsigned int i=0; i<bitmapInfo.width; i++) {
				r = 0.0f;
				g = 0.0f;
				b = 0.0f;

				numRayStack_thread = 0;
				// The first primary ray.
				//pixel_w = cam.proj_origin + i*cam.x_step + j*cam.y_step - cam.pos;

				memcpy(&(primaryRay.dir), &pixel_w, sizeof(Vertex));
				//ray -= cam.pos;

				rayStack_thread[numRayStack_thread].ray = primaryRay;
				rayStack_thread[numRayStack_thread].intensity_coef = 1.0f;
				rayStack_thread[numRayStack_thread].recursion_number = 0;

				numRayStack_thread++;

				while(numRayStack_thread > 0) {
					current_rayItem = rayStack_thread[--numRayStack_thread];

					memcpy(&curRay, &(current_rayItem.ray), sizeof(Ray));

//					normalize(curRay.dir);

					curRay.inv_dir = inverse(curRay.dir);

					nearestObj_cam = -1;
					nearestTri_cam = -1;
					nearestToLightIdx = -2;

//                    nearestToCamObjIdx = -1;

                    // We need t_cam, the intersected Triangle and the index of the intersected object.
                    t_cam = objbvh.traverse(curRay, &nearestObj_cam, &nearestTri_cam);

//					t_cam = bvh.findNearestPrimitive_new(ray.x, ray.y, ray.z,
//														 inv_ray.x, inv_ray.y, inv_ray.z,
//														 current_rayItem.orig.x, current_rayItem.orig.y, current_rayItem.orig.z,
//														 &nearestStaticTri_cam);

//                    t_cam_dyn = obvh.findNearestObject(ray.x, ray.y, ray.z,
//                                                       inv_ray.x, inv_ray.y, inv_ray.z,
//                                                       current_rayItem.orig.x, current_rayItem.orig.y, current_rayItem.orig.z,
//                                                       &nearestObj_cam, &nearestDynamicTri_cam);


					if(nearestObj_cam >= 0) {
						//LOGI("I see a triangle!");

//						contactPoint.x = current_rayItem.ray.orig.x + t_cam * ray.x;
//						contactPoint.y = current_rayItem.ray.orig.y + t_cam * ray.y;
//						contactPoint.z = current_rayItem.ray.orig.z + t_cam * ray.z;

						Material* mat = NULL;
//						if(triArray[nearestStaticTri_cam].n_calculated) {
//							memcpy(&surfaceNormal, &(triArray[nearestStaticTri_cam].n), sizeof(Vertex));
//						} else {
//							// update triangle normal, don't forget to divide by inv_norm
//						}

						// Ambient
//						float _r = 0.0f, _g = 0.0f, _b = 0.0f;
//						switch(id) {
//                            case 0: {
//                                _r = 0; _g = 1; _b = 0;
//                            } break;
//                            case 1: {
//                                _r = 1; _g = 0; _b = 0;
//                            } break;
//                            case 2: {
//                                _r = 0; _g = 0; _b = 1;
//                            } break;
//                            case 3: {
//                                _r = 1; _g = 1; _b = 1;
//                            } break;
//                        }
//                        r_cur = _r;
//						g_cur = _g;
//						b_cur = _b;

//						r_cur = ambient * mat->r;
//						g_cur = ambient * mat->g;
//						b_cur = ambient * mat->b;

                        // test
                        switch(nearestObj_cam) {
                            case 0: {
                                r_cur = 0.0f;
                                g_cur = 1.0f;
                                b_cur = 0.0f;
                            } break;
                            case 1: {
                                r_cur = 0.0f;
                                g_cur = 0.0f;
                                b_cur = 1.0f;
                            } break;
                        }

	#ifdef LIGHTS
						// Find all the lights that light up the nearest primitive.
						for(int k=0; k<numPointLights; k++) {
							memcpy(&lightRay, &contactPoint, sizeof(Vertex));
							lightRay -= ptLightArray[k].pos;
							//lightRay = contactPoint - Vec3(pointLightArray[k*6], pointLightArray[k*6 + 1], pointLightArray[k*6 + 2]);

//							t_light = bvh.findNearestPrimitive(lightRay.x, lightRay.y, lightRay.z, ptLightArray[k].pos.x, ptLightArray[k].pos.y, ptLightArray[k].pos.z, nearestToLightIdx);
							//t_light = findNearestPrimitive(lightRay.x, lightRay.y, lightRay.z, pointLightArray[k*6], pointLightArray[k*6 + 1], pointLightArray[k*6 + 2], nearestToLightIdx);

							if(nearestStaticTri_cam == nearestToLightIdx) { // && type_cam == type_light
								float lightRay_inv_norm = inv_norm(lightRay);
								float proj = dot(surfaceNormal, lightRay);

								// Diffuse
								float diffuse_cos = -proj * lightRay_inv_norm;
								if(diffuse_cos > 0.0f) { // light_cos > 0.0f
									float coef = mat->kd * diffuse_cos;
									r_cur += coef * ptLightArray[k].r * mat->r;
									g_cur += coef * ptLightArray[k].g * mat->g;
									b_cur += coef * ptLightArray[k].b * mat->b;
								}

								// Specular (reflected dot toObserver)
								memcpy(&proj_surfaceNormal, &surfaceNormal, sizeof(Vertex));
								proj_surfaceNormal *= (proj + proj);

								lightRay -= proj_surfaceNormal;

								float spec_cos = cosine(lightRay, ray);
								if(spec_cos <= 0.0f) { // This might be stupid and useless.
									float coef = mat->ks * pow(spec_cos, mat->shine);
									r_cur += coef * ptLightArray[k].r;
									g_cur += coef * ptLightArray[k].g;
									b_cur += coef * ptLightArray[k].b;
								}
							}
						}
	#endif
						r_cur *= current_rayItem.intensity_coef;
						g_cur *= current_rayItem.intensity_coef;
						b_cur *= current_rayItem.intensity_coef;

    #ifdef SECONDARY_RAYS
						if(current_rayItem.recursion_number < recursion_depth) {

							//float ray_inv_norm = inv_norm(ray);
							float proj = dot(surfaceNormal, ray);

							// Reflected ray
							if(proj < 0.0f && mat->ks > 0.0f) { // can determine earlier
								memcpy(&reflectedRay, &ray, sizeof(Vertex));
								memcpy(&proj_surfaceNormal, &surfaceNormal, sizeof(Vertex));

								proj_surfaceNormal *= (proj + proj);
								reflectedRay -= proj_surfaceNormal;

								normalize(reflectedRay);
								//reflectedRay = Vec3::normalize(ray + 2.0f * Vec3::dot(-1.0f*ray, surfaceNormal) * surfaceNormal);

								rayStack_thread[numRayStack_thread].dir = reflectedRay;
								rayStack_thread[numRayStack_thread].orig.x = reflectedRay.x * rec_ray_offset + contactPoint.x;
								rayStack_thread[numRayStack_thread].orig.y = reflectedRay.y * rec_ray_offset + contactPoint.y;
								rayStack_thread[numRayStack_thread].orig.z = reflectedRay.z * rec_ray_offset + contactPoint.z;
								rayStack_thread[numRayStack_thread].intensity_coef = mat->ks * current_rayItem.intensity_coef;
								rayStack_thread[numRayStack_thread].recursion_number = current_rayItem.recursion_number + 1;

								numRayStack_thread++;
							}

							// Refracted ray
							if(mat->T > 0.0f) {

								memcpy(&invertedNormal, &surfaceNormal, sizeof(Vertex));

								//Vec3 normalizedRay = Vec3::normalize(ray), refractedRay, invertedNormal;

								float cos_theta_1 = cosine1(ray, surfaceNormal), ior, sqrterm;

								if(cos_theta_1 < 0.0f) {
									ior = IOR_AIR / mat->ior;
									//invertedNormal = surfaceNormal;
								} else {
									ior = mat->ior / IOR_AIR;
									cos_theta_1 *= -1.0f;
									//invertedNormal = -1.0f*surfaceNormal;
									invertedNormal *= -1.0f;
								}
								sqrterm = 1.0f - ior*ior*(1.0f - cos_theta_1*cos_theta_1);

								if(sqrterm > 0.0f) {
									sqrterm = cos_theta_1*ior + sqrt(sqrterm);

									memcpy(&refractedRay, &ray, sizeof(Vertex));
									refractedRay *= ior;
									refractedRay += invertedNormal * (-sqrterm);

									//refractedRay = invertedNormal * (-sqrterm) + normalizedRay * ior;
								}

								//Vec3 refractedRay = Vec3::normalize(ray);
								rayStack_thread[numRayStack_thread].dir = refractedRay;
								rayStack_thread[numRayStack_thread].orig.x = refractedRay.x * rec_ray_offset + contactPoint.x;
								rayStack_thread[numRayStack_thread].orig.y = refractedRay.y * rec_ray_offset + contactPoint.y;
								rayStack_thread[numRayStack_thread].orig.z = refractedRay.z * rec_ray_offset + contactPoint.z;
								rayStack_thread[numRayStack_thread].intensity_coef = mat->T * current_rayItem.intensity_coef;
								rayStack_thread[numRayStack_thread].recursion_number = current_rayItem.recursion_number + 1;

								numRayStack_thread++;
							}
						}
	#endif
						r += r_cur;
						g += g_cur;
						b += b_cur;
					} else {
					    float _r = 0, _g = 0, _b = 0;
					    switch(id) {
                            case 0: {
                                _r = 0; _g = 1; _b = 0;
                            } break;
                            case 1: {
                                _r = 1; _g = 0; _b = 0;
                            } break;
                            case 2: {
                                _r = 0; _g = 0; _b = 1;
                            } break;
                            case 3: {
                                _r = 1; _g = 1; _b = 1;
                            } break;
                        }

						r += 0.1f * current_rayItem.intensity_coef;
						g += 0.1f * current_rayItem.intensity_coef;
						b += 0.1f * current_rayItem.intensity_coef;

//                        r += _r;
//                        g += _g;
//                        b += _b;
					}
				}

				// Limit the color components to 1.0f.
				r = r > 1.0f ? 1.0f : r;
				g = g > 1.0f ? 1.0f : g;
				b = b > 1.0f ? 1.0f : b;

                // Set the color in the color buffer.
//                if(color_buffer_index < bitmapInfo.width * bitmapInfo.height) {
//                    pixels_uint_thread[color_buffer_index] = argb_to_int(255, (int)(255.0f * r), (int)(255.0f * g), (int)(255.0f * b));
//                }
                pixels_uint_thread[color_buffer_index] = argb_to_int(255, (int)(255.0f * r), (int)(255.0f * g), (int)(255.0f * b));

				pixel_w += cam.x_step; // Step one pixel to the right.
				color_buffer_index++;
			} // End of i
			pixel_w += row_step_back; // Step back to the first column.
			pixel_w += cam.y_step; // Advance to the next row.
		} // End of j
		pixel_w += pixel_skip_y;
	} // End of k
}

void* rayTraceSceneThread_wrapper(void* thread_id_context) {
	RayTracerContext* cntxt = (RayTracerContext*)(((pthread_args*)thread_id_context)->context);

#ifdef USE_RAY_PACKETS
    cntxt->rayTraceScenePacketsThread(thread_id_context);
#else
    #ifdef GRID
        cntxt->rayTrace_thread_grid(thread_id_context);
    #else
        cntxt->rayTrace_thread_horSlabs(thread_id_context);
    #endif
#endif

    return NULL;
}

void RayTracerContext::rayTraceSceneThreads(const float &delta) {
	// Camera movement
	angle += deltaAngle * delta;
	Vertex new_dir = {-camDist * cos(angle), -camHeight, -camDist * sin(angle)};
	Vertex new_pos = {camPos.x + camDist * cos(angle), camPos.y + camHeight, camPos.z + camDist * sin(angle)};
//	Vertex new_dir = {-19725.9f + 19726.8f, 10746.7f - 10747.1f, 5461.97f - 5462.21f};
//	Vertex new_pos = {-19726.8f, 10747.1f, 5462.21f};
	cam.changeDir(new_dir);
	cam.changePos(new_pos);

	pthread_args* thread_args = new pthread_args[NUM_THREADS];

	for(int i=0; i<NUM_THREADS; i++) {
		thread_args[i].id = i;
		thread_args[i].context = this;
		int rc = pthread_create(&threads[i], NULL, rayTraceSceneThread_wrapper, (void*)&(thread_args[i]));
	}
	for(int i=0; i<NUM_THREADS; i++) {
		pthread_join(threads[i], NULL);
	}

	delete[] thread_args;
}

//// Single-threaded ray tracing of the scene, no packets.
//void RayTracerContext::rayTraceScene(const float &delta) {
//	// Camera movement
//	angle += deltaAngle * delta;
//	Vertex new_dir = {-5.0f * cos(angle), -1.0f, 8.0f};
//	Vertex new_pos = {5.0f * cos(angle), 1.0f, -8.0f};
//	cam.changeDir(new_dir);
//	cam.changePos(new_pos);
//
//	int nearestToCamIdx, nearestToLightIdx;
//	float t_cam, t_light;
//
//	Vertex pixel_w = cam.proj_origin - cam.pos;
//	Vertex col_step_back = -((float)(bitmapInfo.height) * cam.y_step);
//	Vertex contactPoint, ray, inv_ray, lightRay, inv_lightRay, surfaceNormal, proj_surfaceNormal, invertedNormal;
//
//	Vertex reflectedRay, refractedRay;
//	RayItem current_rayItem;
//
//	float r, g, b, r_cur, g_cur, b_cur;
//
//	// Acquire the pixel array;
//	pixels_uint = (uint32_t*) pixels;
//
//	int bmp_rows;
//
//	// Cast a ray for each pixel on the projection plane.
//	for(uint32_t i=0; i<bitmapInfo.width; i++) {
//		bmp_rows = i;
//		for(uint32_t j=0; j<bitmapInfo.height; j++) {
//			r = 0.0f;
//			g = 0.0f;
//			b = 0.0f;
//
//			numRayStack = 0;
//			// The first primary ray.
//			//pixel_w = cam.proj_origin + i*cam.x_step + j*cam.y_step;
//			//Vec3 ray = pixel_w - cam.pos;
//
//			memcpy(&ray, &pixel_w, sizeof(Vertex));
//			//ray -= cam.pos;
//
//			rayStack[numRayStack].dir.x = ray.x;
//			rayStack[numRayStack].dir.y = ray.y;
//			rayStack[numRayStack].dir.z = ray.z;
//			rayStack[numRayStack].orig.x = cam.pos.x;
//			rayStack[numRayStack].orig.y = cam.pos.y;
//			rayStack[numRayStack].orig.z = cam.pos.z;
//			rayStack[numRayStack].intensity_coef = 1.0f;
//			rayStack[numRayStack].recursion_number = 0;
//
//			numRayStack++;
//
//			while(numRayStack > 0) {
//				current_rayItem = rayStack[numRayStack-1];
//				numRayStack--;
//
//				memcpy(&ray, &(current_rayItem.dir), sizeof(Vertex));
//
//				normalize(ray);
//
//				inv_ray.x = 1.0f / ray.x;
//				inv_ray.y = 1.0f / ray.y;
//				inv_ray.z = 1.0f / ray.z;
//
//				nearestToCamIdx = -1;
//				nearestToLightIdx = -2;
//
//				//t_cam = bvh.findNearestPrimitive(ray.x, ray.y, ray.z, current_rayItem.orig.x, current_rayItem.orig.y, current_rayItem.orig.z, nearestToCamIdx);
//				t_cam = bvh.findNearestPrimitive_new(ray.x, ray.y, ray.z,
//													 inv_ray.x, inv_ray.y, inv_ray.z,
//													 current_rayItem.orig.x, current_rayItem.orig.y, current_rayItem.orig.z,
//													 &nearestToCamIdx);
//				if(nearestToCamIdx >= 0) {
//					//LOGI("I see a triangle!");
//
//					contactPoint.x = current_rayItem.orig.x + t_cam * ray.x;
//					contactPoint.y = current_rayItem.orig.y + t_cam * ray.y;
//					contactPoint.z = current_rayItem.orig.z + t_cam * ray.z;
//
//					Material* mat = &(matArray[triArray[nearestToCamIdx].mat]);
//					if(triArray[nearestToCamIdx].n_calculated) {
//						memcpy(&surfaceNormal, &(triArray[nearestToCamIdx].n), sizeof(Vertex));
//					} else {
//						// update triangle normal, don't forget to divide by inv_norm
//					}
//
//					// Ambient
//					r_cur = ambient * mat->r;
//					g_cur = ambient * mat->g;
//					b_cur = ambient * mat->b;
//
//					// Find all the lights that light up the nearest primitive.
//					for(int k=0; k<numPointLights; k++) {
//						memcpy(&lightRay, &contactPoint, sizeof(Vertex));
//						lightRay -= ptLightArray[k].pos;
//						//lightRay = contactPoint - Vec3(pointLightArray[k*6], pointLightArray[k*6 + 1], pointLightArray[k*6 + 2]);
//
//						t_light = bvh.findNearestPrimitive(lightRay.x, lightRay.y, lightRay.z, ptLightArray[k].pos.x, ptLightArray[k].pos.y, ptLightArray[k].pos.z, nearestToLightIdx);
//						//t_light = findNearestPrimitive(lightRay.x, lightRay.y, lightRay.z, pointLightArray[k*6], pointLightArray[k*6 + 1], pointLightArray[k*6 + 2], nearestToLightIdx);
//
//						if(nearestToCamIdx == nearestToLightIdx) { // && type_cam == type_light
//							float lightRay_inv_norm = inv_norm(lightRay);
//							float proj = dot(surfaceNormal, lightRay);
//
//							// Diffuse
//							float diffuse_cos = -proj * lightRay_inv_norm;
//							if(diffuse_cos > 0.0f) { // light_cos > 0.0f
//								float coef = mat->kd * diffuse_cos;
//								r_cur += coef * ptLightArray[k].r * mat->r;
//								g_cur += coef * ptLightArray[k].g * mat->g;
//								b_cur += coef * ptLightArray[k].b * mat->b;
//							}
//
//							// Specular (reflected dot toObserver)
//							memcpy(&proj_surfaceNormal, &surfaceNormal, sizeof(Vertex));
//							proj_surfaceNormal *= (proj + proj);
//
//							lightRay -= proj_surfaceNormal;
//
//							float spec_cos = cosine(lightRay, ray);
//							if(spec_cos <= 0.0f) { // This might be stupid and useless.
//								float coef = mat->ks * pow(spec_cos, mat->shine);
//								r_cur += coef * ptLightArray[k].r;
//								g_cur += coef * ptLightArray[k].g;
//								b_cur += coef * ptLightArray[k].b;
//							}
//						}
//					}
//					r_cur *= current_rayItem.intensity_coef;
//					g_cur *= current_rayItem.intensity_coef;
//					b_cur *= current_rayItem.intensity_coef;
//
//					if(current_rayItem.recursion_number < recursion_depth) {
//						//float ray_inv_norm = inv_norm(ray);
//						float proj = dot(surfaceNormal, ray);
//
//						// Reflected ray
//						if(proj < 0.0f && mat->ks > 0.0f) { // can determine earlier
//							memcpy(&reflectedRay, &ray, sizeof(Vertex));
//							memcpy(&proj_surfaceNormal, &surfaceNormal, sizeof(Vertex));
//
//							proj_surfaceNormal *= (proj + proj);
//							reflectedRay -= proj_surfaceNormal;
//
//							normalize(reflectedRay);
//							//reflectedRay = Vec3::normalize(ray + 2.0f * Vec3::dot(-1.0f*ray, surfaceNormal) * surfaceNormal);
//
//							rayStack[numRayStack].dir = reflectedRay;
//							rayStack[numRayStack].orig.x = reflectedRay.x * rec_ray_offset + contactPoint.x;
//							rayStack[numRayStack].orig.y = reflectedRay.y * rec_ray_offset + contactPoint.y;
//							rayStack[numRayStack].orig.z = reflectedRay.z * rec_ray_offset + contactPoint.z;
//							rayStack[numRayStack].intensity_coef = mat->ks * current_rayItem.intensity_coef;
//							rayStack[numRayStack].recursion_number = current_rayItem.recursion_number + 1;
//
//							numRayStack++;
//						}
//
//						// Refracted ray
//						if(mat->T > 0.0f) {
//							memcpy(&invertedNormal, &surfaceNormal, sizeof(Vertex));
//
//							//Vec3 normalizedRay = Vec3::normalize(ray), refractedRay, invertedNormal;
//
//							float cos_theta_1 = cosine1(ray, surfaceNormal), ior, sqrterm;
//
//							if(cos_theta_1 < 0.0f) {
//								ior = IOR_AIR / mat->ior;
//								//invertedNormal = surfaceNormal;
//							} else {
//								ior = mat->ior / IOR_AIR;
//								cos_theta_1 *= -1.0f;
//								//invertedNormal = -1.0f*surfaceNormal;
//								invertedNormal *= -1.0f;
//							}
//							sqrterm = 1.0f - ior*ior*(1.0f - cos_theta_1*cos_theta_1);
//
//							if(sqrterm > 0.0f) {
//								sqrterm = cos_theta_1*ior + sqrt(sqrterm);
//
//								memcpy(&refractedRay, &ray, sizeof(Vertex));
//								refractedRay *= ior;
//								refractedRay += invertedNormal * (-sqrterm);
//
//								//refractedRay = invertedNormal * (-sqrterm) + normalizedRay * ior;
//							}
//
//							//Vec3 refractedRay = Vec3::normalize(ray);
//							rayStack[numRayStack].dir = refractedRay;
//							rayStack[numRayStack].orig.x = refractedRay.x * rec_ray_offset + contactPoint.x;
//							rayStack[numRayStack].orig.y = refractedRay.y * rec_ray_offset + contactPoint.y;
//							rayStack[numRayStack].orig.z = refractedRay.z * rec_ray_offset + contactPoint.z;
//							rayStack[numRayStack].intensity_coef = mat->T * current_rayItem.intensity_coef;
//							rayStack[numRayStack].recursion_number = current_rayItem.recursion_number + 1;
//
//							numRayStack++;
//						}
//					}
//
//					r += r_cur;
//					g += g_cur;
//					b += b_cur;
//				} else {
//					r += 0.6f * current_rayItem.intensity_coef;
//					g += 0.6f * current_rayItem.intensity_coef;
//					b += 0.6f * current_rayItem.intensity_coef;
//				}
//			}
//
//			//} else { // Background.
//			//	pixels_uint[i + (j)*bitmapInfo.width] = argb_to_int(255, 125, 125, 125);
//			//}
//			// Limit the color components to 1.0f.
//			r = r > 1.0f ? 1.0f : r;
//			g = g > 1.0f ? 1.0f : g;
//			b = b > 1.0f ? 1.0f : b;
//			pixels_uint[bmp_rows] = argb_to_int(255, (int)(255.0f * r), (int)(255.0f * g), (int)(255.0f * b));
//
//			bmp_rows += bitmapInfo.width;
//
//			pixel_w += cam.y_step;
//		} // End of j
//		pixel_w += col_step_back;
//		pixel_w += cam.x_step;
//	} // End of i
//}

//void RayTracerContext::rayTraceSceneCL(const float &delta) {
//    // Camera movement
//	angle += deltaAngle * delta;
//	Vertex new_dir = {-5.0f * cos(angle), -1.0f, 8.0f};
//	Vertex new_pos = {5.0f * cos(angle), 1.0f, -8.0f};
////	Vertex new_dir = {-19725.9f + 19726.8f, 10746.7f - 10747.1f, 5461.97f - 5462.21f};
////	Vertex new_pos = {-19726.8f, 10747.1f, 5462.21f};
//	cam.changeDir(new_dir);
//	cam.changePos(new_pos);
//
//    // Set the camera data for transfer to GPU memory.
//	camData_CL.pos = cam.pos;
//	camData_CL.proj_orig = cam.proj_origin;
//	camData_CL.x_step = cam.x_step;
//	camData_CL.y_step = cam.y_step;
//
//    // PARALLEL PART //////////////////////////////////////////////
//	cl_int pixel_count = cam.res_x * cam.res_y;
//	cl_int res_x = (int)(cam.res_x);
//	cl_int res_y = (int)(cam.res_y);
//
//    //size_t global_work_size[] = { (size_t)res_x, (size_t)res_y }; ///< The total number of required work items.
//	//size_t global_work_size[] = { (size_t)(cam.res_x), (size_t)(cam.res_y) }; ///< The total number of required work items.
//	//size_t local_work_size[] = { (size_t)128 }; ///< The number of work items in a work group.
//
//	// Load new data into GPU memory.
//	// Camera
//    errorCode = clEnqueueWriteBuffer(commandQueue, camData_cl_mem, CL_TRUE, 0, sizeof(CamData_CL), (void*)&camData_CL, 0, NULL, NULL); //printError(errorCode);
//
//
//	// Set the kernel call arguments. TODO: For some data, this can only be done once during initialization?
//	errorCode = clSetKernelArg(kernel_trace, 0, sizeof(cl_mem), (void*)&pixelBuffer_mem); //printError(errorCode);
//	errorCode = clSetKernelArg(kernel_trace, 1, sizeof(cl_int), (void*)&(packets_x)); //printError(errorCode);
//	errorCode = clSetKernelArg(kernel_trace, 2, sizeof(cl_int), (void*)&(packets_y)); //printError(errorCode);
//	errorCode = clSetKernelArg(kernel_trace, 3, sizeof(cl_mem), (void*)&camData_cl_mem); //printError(errorCode);
//	errorCode = clSetKernelArg(kernel_trace, 4, sizeof(cl_mem), (void*)&vertArray_cl_mem); //printError(errorCode);
//	errorCode = clSetKernelArg(kernel_trace, 5, sizeof(cl_int), (void*)&numVert); //printError(errorCode);
//	errorCode = clSetKernelArg(kernel_trace, 6, sizeof(cl_mem), (void*)&triArray_cl_mem); //printError(errorCode);
//	errorCode = clSetKernelArg(kernel_trace, 7, sizeof(cl_int), (void*)&numTri); //printError(errorCode);
//	errorCode = clSetKernelArg(kernel_trace, 8, sizeof(cl_mem), (void*)&matArray_cl_mem); //printError(errorCode);
//	errorCode = clSetKernelArg(kernel_trace, 9, sizeof(cl_int), (void*)&numMat); //printError(errorCode);
//	errorCode = clSetKernelArg(kernel_trace, 10, sizeof(cl_mem), (void*)&ptLightArray_cl_mem); //printError(errorCode);
//	errorCode = clSetKernelArg(kernel_trace, 11, sizeof(cl_int), (void*)&numPointLights); //printError(errorCode);
//	errorCode = clSetKernelArg(kernel_trace, 12, sizeof(cl_int), (void*)&pixel_count); //printError(errorCode);
//	errorCode = clSetKernelArg(kernel_trace, 12, sizeof(cl_int), (void*)&pixel_count); //printError(errorCode);
//	errorCode = clSetKernelArg(kernel_trace, 13, sizeof(cl_int), (void*)&raypacket_x); //printError(errorCode);
//	errorCode = clSetKernelArg(kernel_trace, 14, sizeof(cl_int), (void*)&raypacket_y); //printError(errorCode);
//	errorCode = clSetKernelArg(kernel_trace, 15, numTri * sizeof(Triangle), NULL); //printError(errorCode);
//	errorCode = clSetKernelArg(kernel_trace, 16, numVert * sizeof(Vertex), NULL); //printError(errorCode);
//
//    //LOGC("clEnqueueNDRangeKernel:");
//	errorCode = clEnqueueNDRangeKernel(commandQueue, kernel_trace, 2, NULL, global_work_size, local_work_size, 0, NULL, NULL); //printError(errorCode);
//	errorCode = clFinish(commandQueue); //printError(errorCode);
//
//	// Read the new pixels.
//	errorCode = clEnqueueReadBuffer(commandQueue, pixelBuffer_mem, CL_TRUE, 0, sizeof(uint32_t) * cam.res_x * cam.res_y, pixels, 0, NULL, NULL); //printError(errorCode);
//
//}







