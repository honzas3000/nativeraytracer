#include <RayTracerContext.h>

void* RayTracerContext::rayTraceScenePacketsThread(void* thread_id_context) {
	pthread_args* thread_args = (pthread_args*)thread_id_context;
	int id = thread_args->id;

	//thread_range* screen_range = &thread_ranges[id];

	int nearestToCamIdx, nearestToLightIdx;
	float t_cam, t_light;

	Vertex pixel_w;
	Vertex col_step_back = -((float)(CELL_SIZE) * cam.y_step);
	Vertex contactPoint, ray, inv_ray, lightRay, inv_lightRay, surfaceNormal, proj_surfaceNormal, invertedNormal;

	Vertex reflectedRay, refractedRay;

	int numRayStack_thread;
	RayItem rayStack_thread[rayStack_size]; // Each thread needs its own stack.
	RayItem current_rayItem;

	float r_cur, g_cur, b_cur;
	float r[NUM_RAYS_PACKET], g[NUM_RAYS_PACKET], b[NUM_RAYS_PACKET];
	Ray rays[NUM_RAYS_PACKET];
	// All primary rays start at the camera center:
	for(int i=0; i<NUM_RAYS_PACKET; i++) {
        rays[i].orig = cam.pos;
	}
	int nearestPrimitives[NUM_RAYS_PACKET];

	int bmp_rows;

	// Acquire the pixel array;
	uint32_t* pixels_uint_thread = (uint32_t*) pixels;

	for(int k=0; k<numThreadCells[id]; k++) {
		int start_x = threadCellStartCoords[id][k].start_x;
		int start_y = threadCellStartCoords[id][k].start_y;
		int end_x = threadCellStartCoords[id][k].end_x;
		int end_y = threadCellStartCoords[id][k].end_y;

		pixel_w = (cam.proj_origin - cam.pos) + (cam.x_step * start_x) + (cam.y_step * start_y);
		// Trace ray packets for all pixels in <start_x, end_x) x <start_y, end_y).

		// Instead of individual rays, send packets
		for(int i=start_x; i<end_x; i+=PACKET_SIZE) {
			bmp_rows = i + (start_y * bitmapInfo.width);
			for(int j=start_y; j<end_y; j+=PACKET_SIZE) {
				memset(&r, 0, NUM_RAYS_PACKET*sizeof(float)); //r = 0.0f;
				memset(&g, 0, NUM_RAYS_PACKET*sizeof(float)); //g = 0.0f;
				memset(&b, 0, NUM_RAYS_PACKET*sizeof(float)); //b = 0.0f;

                // Send primary rays as packets of PACKET_SIZE x PACKET_SIZE.
                // Create the ray packet.
                Ray* _ray;
                Vertex packet_w;
                for(int x=0; x<PACKET_SIZE; x++) {
                    packet_w = pixel_w + x * cam.x_step;
                    for(int y=0; y<PACKET_SIZE; y++) {
                        _ray = &(rays[x + y*PACKET_SIZE]);
                        memcpy(&(_ray->dir), &packet_w, sizeof(Vertex));
                        _ray->inv_dir = inverse(_ray->dir);
                        // r->orig is already set
                        packet_w += cam.y_step;
                    }
                }

                memset(&nearestPrimitives, -1, NUM_RAYS_PACKET*sizeof(int));

                // Find the nearest intersected primitives for all the rays in the packet.
                //bvh.findNearestPrimitive_Packets(&rays[0], &nearestPrimitives[0], PACKET_SIZE);

                // Solve each of the rays individually for shadow and secondary rays.
                for(int m=0; m<PACKET_SIZE; m++) {
                    for(int n=0; n<PACKET_SIZE; n++) {
                        int pixel_index = bmp_rows + m + n * bitmapInfo.width;
                        int tri_idx = nearestPrimitives[m + n * PACKET_SIZE];
                        if(tri_idx >= 0) { // nearestPrimitives[m + n * PACKET_SIZE] >= 0
                            int re = (matArray[triArray[tri_idx].mat].r * 255.0f);
                            int gr = (matArray[triArray[tri_idx].mat].g * 255.0f);
                            int bl = (matArray[triArray[tri_idx].mat].b * 255.0f);
                            if(pixel_index < bitmapInfo.width * bitmapInfo.height) { // Just to make sure, remove later.
                                pixels_uint_thread[pixel_index] = argb_to_int(255, re, gr, bl);
                            }
                        } else {
                            if(pixel_index < bitmapInfo.width * bitmapInfo.height) { // Just to make sure, remove later.
                                pixels_uint_thread[pixel_index] = argb_to_int(255, 50, 50, 50);
                            }
                        }
                    }
                }

                // Test drawing of cells and packets.
//                for(int m=0; m<PACKET_SIZE; m++) {
//                    for(int n=0; n<PACKET_SIZE; n++) {
//                        int pixel_index = bmp_rows + m + n * bitmapInfo.width;
//                        if((m == 0 || m == PACKET_SIZE - 1) && (n == 0 || n == PACKET_SIZE - 1)) { // nearestPrimitives[m + n * PACKET_SIZE] >= 0
//                            if(pixel_index < bitmapInfo.width * bitmapInfo.height) { // Just to make sure, remove later.
//                                pixels_uint_thread[pixel_index] = argb_to_int(255, (int)(255.0f), (int)(255.0f), (int)(255.0f));
//                            }
//                        } else {
//                            int re, gr, bl;
//                            switch(id) {
//                                case 0: {re = 255; gr = 0; bl = 0;}; break;
//                                case 1: {re = 255; gr = 255; bl = 0;}; break;
//                                case 2: {re = 255; gr = 0; bl = 255;}; break;
//                                case 3: {re = 50; gr = 50; bl = 50;}; break;
//                                default: {re = 255; gr = 255; bl = 255;}; break;
//                            }
//                            if(pixel_index < bitmapInfo.width * bitmapInfo.height) { // Just to make sure, remove later.
//                                pixels_uint_thread[pixel_index] = argb_to_int(255, re, gr, bl);
//                            }
//                        }
//                    }
//                }



				bmp_rows += PACKET_SIZE * bitmapInfo.width;

				pixel_w += PACKET_SIZE * cam.y_step;
			} // End of j
			pixel_w += col_step_back;
			pixel_w += PACKET_SIZE * cam.x_step;
		} // End of i
	}

	//LOGI("Completed thread exec.: %d, %d x %d, %d", screen_range->start_x, screen_range->end_x, screen_range->start_y, screen_range->end_y);

	return NULL;
}
