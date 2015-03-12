#ifndef _BVH
#define _BVH

#include <AABB.h>
#include <Structures.h>


struct BVH_node {
    AABB box; // 24B
	int* primArray;
	int numPrim;
	int leftChild;
	int rightChild;

	bool isLeaf; // 1B

	BVH_node() : leftChild(-1), rightChild(-1), isLeaf(false) {}
};

struct nodeBucket {
	int* primArray;
	int numPrim;

	int leftPrimitiveCount;
	int rightPrimitiveCount;
	AABB leftBox, rightBox;

	nodeBucket() : numPrim(0), leftPrimitiveCount(0), rightPrimitiveCount(0) {}
};

class BVH {
private:
	float C_t, C_i;
	int divBuckets, numPrimInLeaf;
	nodeBucket* bucketList;

	// Scene data
	// Vertices
	int numVert;
	Vertex* vertArray;
	// Triangles
	int numTri; // The actual number of triangles.
	Triangle* triArray;

	int* primArray;

	// BVH node list
	BVH_node* nodeList;
	int nodeList_size;
	int numBvhNodes;

	// For queries
	//BVH_node** nodeStack;
public:
	BVH() : C_t(2.0f), C_i(1.0f), divBuckets(10), numPrimInLeaf(1) {
		nodeList = new BVH_node[100000];
		nodeList_size = 100000;

		bucketList = new nodeBucket[divBuckets];
	}

	~BVH() {
		delete[] primArray;
		//delete[] nodeStack;
		delete[] bucketList;
		delete[] nodeList;
	}
private:
	void clearBucketList() {
		for(int i=0; i<divBuckets; i++) {
			bucketList[i].leftBox.reset();
			bucketList[i].rightBox.reset();
			bucketList[i].leftPrimitiveCount = 0;
			bucketList[i].rightPrimitiveCount = 0;
			bucketList[i].primArray = NULL;
			bucketList[i].numPrim = 0;
		}
	}

	void containPrimitives(AABB &box, int* &_primArray, const int &_numPrim) {
		for(int i=0; i<_numPrim; i++) {
			int triIndex = _primArray[i];
			box.contain(vertArray[triArray[triIndex].v1_idx]);
			box.contain(vertArray[triArray[triIndex].v2_idx]);
			box.contain(vertArray[triArray[triIndex].v3_idx]);
		}
	}

	int getDividingIndexBucketing(int* &_primArray, const int &_primCount,
								  const float &_parentArea, const int &_divAxis) {
		float bestSahValue = FLT_MAX;
		int dividingIndex = 1;

		int curDivBuckets = divBuckets > _primCount ? _primCount : divBuckets;

		clearBucketList();

		float left, right;
		left = triArray[_primArray[0]].c[_divAxis];
		//left = cenArray[3*_primArray[0] + _divAxis];
		right = triArray[_primArray[_primCount - 1]].c[_divAxis];

		float bucketSize = (right - left) / curDivBuckets;
		// Distribute primitives into buckets
		if(bucketSize > 0.0f) {
			int bucketIndex = -1;
			int new_index;

			for(int i=0; i<_primCount; i++) {
				new_index = (int)((triArray[_primArray[i]].c[_divAxis] - left) / bucketSize);
				new_index = new_index >= curDivBuckets ? curDivBuckets - 1 : new_index < 0 ? 0 : new_index;

				if(new_index > bucketIndex) { // Move on to the next bucket.
					bucketIndex = new_index;
					bucketList[bucketIndex].primArray = _primArray + i; // The primitive list of this bucket starts here.
					bucketList[bucketIndex].numPrim++; // And has one more primitive.
				} else {
					bucketList[bucketIndex].numPrim++; // And has one more primitive.
				}
			}
		} else {
			bucketList[0].primArray = _primArray;
			bucketList[0].numPrim = _primCount;
		}

		// Contain from left to right.
		containPrimitives(bucketList[0].leftBox, bucketList[0].primArray, bucketList[0].numPrim);
		bucketList[0].leftPrimitiveCount = bucketList[0].numPrim;
		for(int i=1; i<curDivBuckets; i++) {
			bucketList[i].leftBox.contain(bucketList[i-1].leftBox);
			containPrimitives(bucketList[i].leftBox, bucketList[i-1].primArray, bucketList[i-1].numPrim);
			bucketList[i].leftPrimitiveCount = bucketList[i].numPrim + bucketList[i-1].leftPrimitiveCount;
		}
		// Contain from right to left.
		containPrimitives(bucketList[curDivBuckets-1].leftBox, bucketList[curDivBuckets-1].primArray, bucketList[curDivBuckets-1].numPrim);
		bucketList[curDivBuckets-1].rightPrimitiveCount = bucketList[curDivBuckets-1].numPrim;
		for(int i=curDivBuckets - 2; i>=0; i--) {
			bucketList[i].rightBox.contain(bucketList[i+1].rightBox);
			containPrimitives(bucketList[i].rightBox, bucketList[i+1].primArray, bucketList[i+1].numPrim);
			bucketList[i].rightPrimitiveCount = bucketList[i].numPrim + bucketList[i+1].rightPrimitiveCount;
		}

		// Find the best dividing plane.
		for(int i=0; i<curDivBuckets-1; i++) {
			float sahVal = C_t +
					(C_i * bucketList[i].leftBox.getArea() * bucketList[i].leftPrimitiveCount +
					C_i * bucketList[i+1].rightBox.getArea() * bucketList[i+1].rightPrimitiveCount) / _parentArea;
			if(sahVal < bestSahValue) {
				bestSahValue = sahVal;
				dividingIndex = bucketList[i].leftPrimitiveCount;
			}
		}
		return dividingIndex;
	}

	float rayTriangleIntersection(const int &tri_idx,
								  const float &rayX, const float &rayY, const float &rayZ,
								  const float &originX, const float &originY, const float &originZ) {
		Vertex v1 = {vertArray[triArray[tri_idx].v1_idx].x, vertArray[triArray[tri_idx].v1_idx].y, vertArray[triArray[tri_idx].v1_idx].z};
		Vertex v2 = {vertArray[triArray[tri_idx].v2_idx].x, vertArray[triArray[tri_idx].v2_idx].y, vertArray[triArray[tri_idx].v2_idx].z};
		Vertex v3 = {vertArray[triArray[tri_idx].v3_idx].x, vertArray[triArray[tri_idx].v3_idx].y, vertArray[triArray[tri_idx].v3_idx].z};

		//Vertex b = {originX - v1.x, originY - v1.y, originZ - v1.z};
		Vertex vec1 = {v2.x - v1.x, v2.y - v1.y, v2.z - v1.z};
		Vertex vec2 = {v3.x - v1.x, v3.y - v1.y, v3.z - v1.z};

		Vertex D = {rayX, rayY, rayZ};
		Vertex P;
		cross(D, vec2, P);

		float det = dot(vec1, P);
		if(det > -0.000001f && det < 0.000001f) return -1.0f;
		float inv_det = 1.0f / det;

		Vertex O = {originX, originY, originZ};
		Vertex T = O - v1;

		float u = dot(T, P) * inv_det;
		if(u < 0.0f || u > 1.0f) return -1.0f;

		Vertex Q;
		cross(T, vec1, Q);

		float v = dot(D, Q) * inv_det;
		if(v < 0.0f || u + v > 1.0f) return -1.0f;

		return dot(vec2, Q) * inv_det; // t



/*
		float detA = det(vec1.x, vec1.y, vec1.z, vec2.x, vec2.y, vec2.z, -rayX, -rayY, -rayZ);

		if (detA == 0.0f) { // < 0.00000001
			return -1.0f;
		}

		float detA_inv = 1.0f / detA;

		float detA_o = det(vec1.x, vec1.y, vec1.z, vec2.x, vec2.y, vec2.z, b.x, b.y, b.z);

		float o = detA_o * detA_inv;
		if (o < 0) { // The object is too close or behind the camera.
			return -1.0f;
		}

		float detA_m = det(b.x, b.y, b.z, vec2.x, vec2.y, vec2.z, -rayX, -rayY, -rayZ);
		float detA_n = det(vec1.x, vec1.y, vec1.z, b.x, b.y, b.z, -rayX, -rayY, -rayZ);
		float m = detA_m * detA_inv;
		float n = detA_n * detA_inv;

		if(m < 0.0f || n < 0.0f || m + n > 1.0f) {
			return -1.0f;
		} else {
			return o;
		}
*/
	}

	bool comparePrimitives(const int &tri1_idx, const int &tri2_idx, const int &axis) {
		return triArray[tri1_idx].c[axis] <= triArray[tri2_idx].c[axis];
		//return cenArray[3*tri1_idx + axis] <= cenArray[3*tri2_idx + axis];
	}

	void sortPrimitiveArray(int* &_primArray, const int &_numPrim, const int &axis) {
		int numSwaps = 1;
		while(numSwaps > 0) {
			numSwaps = 0;
			for(int i=0; i<_numPrim-1; i++) {
				if(!comparePrimitives(_primArray[i], _primArray[i + 1], axis)) {
					int temp = _primArray[i];
					_primArray[i] = _primArray[i + 1];
					_primArray[i + 1] = temp;
					numSwaps++;
				}
			}
		}
	}
public:
	void createBVH(Vertex* &_vertArray, const int &_numVert, Triangle* &_triArray, const int &_numTri) {
		//int numNodes = 0;

		vertArray = _vertArray;
		numVert = _numVert;
		triArray = _triArray;
		numTri = _numTri;

		// An array of triangle indices, to be sorted and partitioned into nodes.
		primArray = new int[numTri];
		//memset(primArray, 0, numTri * sizeof(int));
		for(int i=0; i<numTri; i++) {
			primArray[i] = i;
		}

		BVH_node* root = &nodeList[0];
		numBvhNodes = 1;

		root->primArray = primArray;
		root->numPrim = numTri;
		containPrimitives(root->box, root->primArray, root->numPrim);
		//numNodes++;

		// If the number of primitives is lower than the number of primitives for a leaf.
		if(numTri <= numPrimInLeaf) {
			root->isLeaf = true;
			//nodeStack = new BVH_node*[numBvhNodes];
			return;
		}

		root->isLeaf = false;

		int* toSubdivideStack = new int[numTri];
		int toSubdivideStack_size = numTri;
		int numToSubdivide = 0;
		BVH_node* curNode;
		int* primSubList;
		int numPrimSubList;

		toSubdivideStack[numToSubdivide] = 0;
		numToSubdivide++;

		while(numToSubdivide > 0) {
			curNode = &nodeList[toSubdivideStack[numToSubdivide-1]];
			numToSubdivide--;

			if(curNode->numPrim <= numPrimInLeaf) {
				curNode->isLeaf = true; // Do not subdivide further.
			} else {
				if(nodeList_size - numBvhNodes < 2) {
					LOGI("Reallocating nodeList.");
					// Reallocate the node array.
					nodeList_size *= 2;
					BVH_node* new_nodeList = new BVH_node[nodeList_size];
					memcpy(new_nodeList, nodeList, numBvhNodes * sizeof(BVH_node));
					delete[] nodeList;
					nodeList = new_nodeList;

					curNode = &nodeList[toSubdivideStack[numToSubdivide]];
				}

				BVH_node* leftChild = &nodeList[numBvhNodes];
				numBvhNodes++;
				BVH_node* rightChild = &nodeList[numBvhNodes];
				numBvhNodes++;

				/*BVH_node* leftChild = new BVH_node();
				BVH_node* rightChild = new BVH_node();
				numBvhNodes += 2;*/

				//LOGI("Two children created.");

				int divAxis = curNode->box.getWidestAxis();
				//LOGI("divAxis obtained.");

				// The primitive sub list to work with.
				primSubList = curNode->primArray;
				numPrimSubList = curNode->numPrim;
				//LOGI("primSubList obtained.");

				sortPrimitiveArray(primSubList, numPrimSubList, divAxis);

				//LOGI("Primitive sublist sorted.");

				int divIndex = getDividingIndexBucketing(primSubList, numPrimSubList, curNode->box.getArea(), divAxis);

				//LOGI("divIndex: %d, numPrim: %d", divIndex, numPrimSubList);
				//int divIndex = 1;

				leftChild->primArray = primSubList;
				leftChild->numPrim = divIndex;
				rightChild->primArray = primSubList + divIndex;
				rightChild->numPrim = numPrimSubList - divIndex;

				// left child, contain triangles in the box
				containPrimitives(leftChild->box, leftChild->primArray, leftChild->numPrim);
				// right child, contain triangles in the box
				containPrimitives(rightChild->box, rightChild->primArray, rightChild->numPrim);

				curNode->primArray = NULL;
				curNode->numPrim = 0;
				/*curNode->leftChild = leftChild;
				curNode->rightChild = rightChild;*/
				curNode->leftChild = numBvhNodes - 2;
				curNode->rightChild = numBvhNodes - 1;

				if(toSubdivideStack_size - numToSubdivide < 2) {
					//LOGI("Reallocating the toSubdivide stack.");
					// Reallocate the stack.
					toSubdivideStack_size *= 2;
					int* new_toSubdivideStack = new int[toSubdivideStack_size];
					memcpy(new_toSubdivideStack, toSubdivideStack, numToSubdivide * sizeof(int));
					/*for(int i=0; i<numToSubdivide; i++) {
						new_toSubdivideStack[i] = toSubdivideStack[i];
					}*/
					delete[] toSubdivideStack;
					toSubdivideStack = new_toSubdivideStack;
				}
				toSubdivideStack[numToSubdivide] = numBvhNodes - 2;
				numToSubdivide++;
				toSubdivideStack[numToSubdivide] = numBvhNodes - 1;
				numToSubdivide++;

				//LOGI("Added children to stack.");
			}
		}
		//nodeStack = new BVH_node*[numBvhNodes];

        delete[] toSubdivideStack;

		LOGI("BVH DONE, %d nodes.", numBvhNodes);
	}

	float findNearestPrimitive(const float &x, const float &y, const float &z,
			const float &origin_x, const float &origin_y, const float &origin_z, int &nearest) {
		BVH_node* nodeStack[numBvhNodes];
		int nodeStackSize = 0;

		BVH_node* currentNode;
		float min_dist = FLT_MAX;
		float dist;

		nodeStack[nodeStackSize] = &nodeList[0]; // push back
		nodeStackSize++;
		while(nodeStackSize > 0) {
			currentNode = nodeStack[nodeStackSize-1]; // pop back
			nodeStackSize--;
			// If the ray intersects this box, check its left and right child.
			if(currentNode->box.checkRayIntersect(origin_x, origin_y, origin_z, x, y, z)) {
				if(!currentNode->isLeaf) {
					nodeStack[nodeStackSize] = &nodeList[currentNode->leftChild];
					nodeStackSize++;
					nodeStack[nodeStackSize] = &nodeList[currentNode->rightChild];
					nodeStackSize++;
				} else {
					int* leafPrimArray = currentNode->primArray;
					int numLeafPrim = currentNode->numPrim;

					for(int i=0; i<numLeafPrim; i++) {
						int pr = currentNode->primArray[i];
						/*int v1_idx = triArray[4*pr];
						int v2_idx = triArray[4*pr + 1];
						int v3_idx = triArray[4*pr + 2];*/

						dist = rayTriangleIntersection(pr, x, y, z, origin_x, origin_y, origin_z);

						// kontroluje se tady a predtim v rayTriangleIntersection
						if(dist >= 0.0f && dist < min_dist) {
							min_dist = dist;
							nearest = pr;
						}
					}
				}
			}
		}
		return min_dist;
	}

	float findNearestPrimitive_new(const float &x, const float &y, const float &z,
								   const float &inv_x, const float &inv_y, const float &inv_z,
								   const float &origin_x, const float &origin_y, const float &origin_z, int* nearest) {
		BVH_node* nodeStack[numBvhNodes];
		int nodeStackSize = 1;

		BVH_node* currentNode;
		float min_dist = FLT_MAX;
		float dist;
		float box_t;

		int* leafPrimArray;
		int numLeafPrim;

		nodeStack[0] = &nodeList[0]; // push back
		//nodeStackSize++;
		while(nodeStackSize > 0) {
			currentNode = nodeStack[--nodeStackSize]; // pop back
			// If the ray intersects this box, check its left and right child.
			//box_t = currentNode->box.checkRayIntersect_new(origin_x, origin_y, origin_z, inv_x, inv_y, inv_z);
			bool box_t_is = currentNode->box.checkRayIntersect_new(origin_x, origin_y, origin_z, inv_x, inv_y, inv_z);
			if(box_t_is) {
				if(!currentNode->isLeaf) {
					nodeStack[nodeStackSize] = &nodeList[currentNode->leftChild];
					nodeStackSize++;
					nodeStack[nodeStackSize] = &nodeList[currentNode->rightChild];
					nodeStackSize++;
				} else {
					leafPrimArray = currentNode->primArray;
					numLeafPrim = currentNode->numPrim;

					for(int i=0; i<numLeafPrim; i++) {
						int pr = currentNode->primArray[i];
						dist = rayTriangleIntersection(currentNode->primArray[i], x, y, z, origin_x, origin_y, origin_z);

						// kontroluje se tady a predtim v rayTriangleIntersection
						if(dist >= 0.0f && dist < min_dist) {
							min_dist = dist;
							*nearest = pr;
						}
					}
				}
			}
		}
		return min_dist;
	}

	#define TYPE_SQUARE 0b00001111
	#define TYPE_TRI_BOTTOM_LEFT 0b00001101
	#define TYPE_TRI_BOTTOM_RIGHT 0b00001110
	#define TYPE_TRI_TOP_RIGHT 0b00000111
	#define TYPE_TRI_TOP_LEFT 0b00001011
	#define TYPE_LINE_LEFT 0b00001001
	#define TYPE_LINE_RIGHT 0b00000110
	#define TYPE_LINE_BOTTOM 0b00001100
	#define TYPE_LINE_TOP 0b00000011
	#define TYPE_SINGLE 0b00010000
	#define TYPE_MISS 0b00000000

    // A packet of rays
	struct RayPacketJob {
        BVH_node* nodeToCheck;
        char type;
	};

	struct RayJob {
        BVH_node* nodeToCheck;
        int rayIdx;
	};

	float findNearestPrimitive_Packets(Ray* rays, int* nearest, int packet_size) {
	    int numRays = packet_size * packet_size;

        RayPacketJob rayPacketJobStack[numBvhNodes];
        int stackSize_packet = 0;
        RayPacketJob curJob;

        RayJob rayJobStack[numBvhNodes * numRays];
        int stackSize_ray = 0;
        RayJob curRayJob;

        Ray* ray;

		float min_dist[numRays];
		for(int i=0; i<numRays; i++) min_dist[i] = FLT_MAX;
		float dist;
		float box_t;

		int* leafPrimArray;
		int numLeafPrim;

		rayPacketJobStack[0].nodeToCheck = &nodeList[0]; // BVH root
		rayPacketJobStack[0].type = TYPE_SQUARE;
		stackSize_packet++;

        while(stackSize_packet > 0) {
            memcpy(&curJob, &rayPacketJobStack[--stackSize_packet], sizeof(RayPacketJob));
            //curJob = rayPacketJobStack[--stackSize];

            char newType = 0b00000000;
            AABB* box = &(curJob.nodeToCheck->box);

            ray = &rays[0];
            if(box->checkRayIntersect_new(ray->orig.x, ray->orig.y, ray->orig.z, ray->inv_dir.x, ray->inv_dir.y, ray->inv_dir.z)) {
                newType |= 0b00001000;
            }
            ray = &rays[packet_size];
            if(box->checkRayIntersect_new(ray->orig.x, ray->orig.y, ray->orig.z, ray->inv_dir.x, ray->inv_dir.y, ray->inv_dir.z)) {
                newType |= 0b00000100;
            }
            ray = &rays[numRays - 1];
            if(box->checkRayIntersect_new(ray->orig.x, ray->orig.y, ray->orig.z, ray->inv_dir.x, ray->inv_dir.y, ray->inv_dir.z)) {
                newType |= 0b00000010;
            }
            ray = &rays[numRays - packet_size];
            if(box->checkRayIntersect_new(ray->orig.x, ray->orig.y, ray->orig.z, ray->inv_dir.x, ray->inv_dir.y, ray->inv_dir.z)) {
                newType |= 0b00000001;
            }

            if(newType == TYPE_SQUARE) {
                if(curJob.nodeToCheck->isLeaf) {
                    leafPrimArray = curJob.nodeToCheck->primArray;
                    numLeafPrim = curJob.nodeToCheck->numPrim;
                    // All rays with all triangles in the AABB.
                    for(int r=0; r<numRays; r++) {
                        ray = &rays[r];
                        for(int i=0; i<numLeafPrim; i++) {
                            int pr = leafPrimArray[i];
                            dist = rayTriangleIntersection(pr, ray->dir.x, ray->dir.y, ray->dir.z,
                                                           ray->orig.x, ray->orig.y, ray->orig.z);
                            // kontroluje se tady a predtim v rayTriangleIntersection FIXED
                            if(dist >= 0.0f && dist < min_dist[r]) {
                                min_dist[r] = dist;
                                nearest[r] = pr;
                            }
                        }
                    }
                } else {
                    rayPacketJobStack[stackSize_packet].nodeToCheck = &(nodeList[curJob.nodeToCheck->leftChild]);
                    stackSize_packet++;
                    rayPacketJobStack[stackSize_packet].nodeToCheck = &(nodeList[curJob.nodeToCheck->rightChild]);
                    stackSize_packet++;
                }
            } else {
//                for(int i=0; i<numRays; i++) {
//                    rayJobStack[stackSize_ray].nodeToCheck = curJob.nodeToCheck;
//                    rayJobStack[stackSize_ray].rayIdx = i;
//                    stackSize_ray++;
//                }
            }
        }

//        // individual rays.
//        while(stackSize_ray > 0) {
//            memcpy(&curRayJob, &rayJobStack[--stackSize_ray], sizeof(RayJob));
//
//            ray = &rays[curRayJob.rayIdx];
//            if(curRayJob.nodeToCheck->box.checkRayIntersect_new(ray->orig.x, ray->orig.y, ray->orig.z, ray->inv_dir.x, ray->inv_dir.y, ray->inv_dir.z)) {
//                if(curRayJob.nodeToCheck->isLeaf) {
//                    leafPrimArray = curRayJob.nodeToCheck->primArray;
//                    numLeafPrim = curRayJob.nodeToCheck->numPrim;
//
//                    ray = &rays[curRayJob.rayIdx];
//                    for(int i=0; i<numLeafPrim; i++) {
//                        int pr = leafPrimArray[i];
//                        dist = rayTriangleIntersection(pr, ray->dir.x, ray->dir.y, ray->dir.z,
//                                                       ray->orig.x, ray->orig.y, ray->orig.z);
//                        // kontroluje se tady a predtim v rayTriangleIntersection FIXED
//                        if(dist >= 0.0f && dist < min_dist[curRayJob.rayIdx]) {
//                            min_dist[curRayJob.rayIdx] = dist;
//                            nearest[curRayJob.rayIdx] = pr;
//                        }
//                    }
//                } else {
//                    rayJobStack[stackSize_ray].nodeToCheck = &(nodeList[curRayJob.nodeToCheck->leftChild]);
//                    stackSize_ray++;
//                    rayJobStack[stackSize_ray].nodeToCheck = &(nodeList[curRayJob.nodeToCheck->rightChild]);
//                    stackSize_ray++;
//                }
//            }
//        }
	}

    AABB* getRootAABB() {
        return &(nodeList[0].box);
    }
};

#endif



















