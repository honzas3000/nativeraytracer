#ifndef OBJECTBVH_H
#define OBJECTBVH_H

#include <AABB.h>
#include <Structures.h>
#include <DynObj.h>

class ObjectBVH {
private:
    struct ObjectBVH_node {
        AABB box; // 24B
        int* primArray;
        int numPrim;
        int leftChild;
        int rightChild;

        bool isLeaf; // 1B

        ObjectBVH_node() : leftChild(-1), rightChild(-1), isLeaf(false) {}
    };

    struct nodeBucket {
        int* primArray;
        int numPrim;

        int leftPrimitiveCount;
        int rightPrimitiveCount;
        AABB leftBox, rightBox;

        nodeBucket() : numPrim(0), leftPrimitiveCount(0), rightPrimitiveCount(0) {}
    };

	float C_t, C_i;
	int divBuckets, numPrimInLeaf;
	nodeBucket* bucketList;

	// Scene data
	int numDynObj;
	DynamicObject* dynObjArray;

	int* primArray;

	// BVH node list
	ObjectBVH_node* nodeList;
	int nodeList_size;
	int numBvhNodes;

	// For queries
	//BVH_node** nodeStack;
public:
	ObjectBVH() : C_t(1.0f), C_i(1.0f), divBuckets(10), numPrimInLeaf(1) {
		nodeList = new ObjectBVH_node[1000];
		nodeList_size = 1000;

		bucketList = new nodeBucket[divBuckets];
	}

	~ObjectBVH() {
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
			box.contain(*(dynObjArray[_primArray[i]].objectBox));
		}
	}

	int getDividingIndexBucketing(int* &_primArray, const int &_primCount,
								  const float &_parentArea, const int &_divAxis) {
		float bestSahValue = FLT_MAX;
		int dividingIndex = 1;

		int curDivBuckets = divBuckets > _primCount ? _primCount : divBuckets;

		clearBucketList();

		float left, right;
		left = dynObjArray[_primArray[0]].center[_divAxis];
		//left = cenArray[3*_primArray[0] + _divAxis];
		right = dynObjArray[_primArray[_primCount - 1]].center[_divAxis];

		float bucketSize = (right - left) / curDivBuckets;
		// Distribute primitives into buckets
		if(bucketSize > 0.0f) {
			int bucketIndex = -1;
			int new_index;

			for(int i=0; i<_primCount; i++) {
				new_index = (int)((dynObjArray[_primArray[i]].center[_divAxis] - left) / bucketSize);
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

	bool comparePrimitives(const int &tri1_idx, const int &tri2_idx, const int &axis) {
		return dynObjArray[tri1_idx].center[axis] <= dynObjArray[tri2_idx].center[axis];
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
	void createBVH(DynamicObject* _dynObjArray, const int &_numDynObj) {
		//int numNodes = 0;

		dynObjArray = _dynObjArray;
		numDynObj = _numDynObj;

		// An array of object indices, to be sorted and partitioned into nodes.
		primArray = new int[numDynObj];
		//memset(primArray, 0, numTri * sizeof(int));
		for(int i=0; i<numDynObj; i++) {
			primArray[i] = i;
		}

		ObjectBVH_node* root = &nodeList[0];
		numBvhNodes = 1;

		root->primArray = primArray;
		root->numPrim = numDynObj;
		containPrimitives(root->box, root->primArray, root->numPrim);
		//numNodes++;

		// If the number of primitives is lower than the number of primitives for a leaf.
		if(numDynObj <= numPrimInLeaf) {
			root->isLeaf = true;
			//nodeStack = new BVH_node*[numBvhNodes];
			return;
		}

		root->isLeaf = false;

		int* toSubdivideStack = new int[numDynObj];
		int toSubdivideStack_size = numDynObj;
		int numToSubdivide = 0;
		ObjectBVH_node* curNode;
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
					ObjectBVH_node* new_nodeList = new ObjectBVH_node[nodeList_size];
					memcpy(new_nodeList, nodeList, numBvhNodes * sizeof(ObjectBVH_node));
					delete[] nodeList;
					nodeList = new_nodeList;

					curNode = &nodeList[toSubdivideStack[numToSubdivide]];
				}

				ObjectBVH_node* leftChild = &nodeList[numBvhNodes];
				numBvhNodes++;
				ObjectBVH_node* rightChild = &nodeList[numBvhNodes];
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

	float findNearestObject(const float &x, const float &y, const float &z,
                            const float &inv_x, const float &inv_y, const float &inv_z,
                            const float &origin_x, const float &origin_y, const float &origin_z,
                            int* nearestObject, int* nearestTriangle) {
		ObjectBVH_node* nodeStack[numBvhNodes];
		int nodeStackSize = 1;

		ObjectBVH_node* currentNode;
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
			float t_box;
			bool box_t_is = currentNode->box.checkRayIntersect_new_t(origin_x, origin_y, origin_z, inv_x, inv_y, inv_z, &t_box);
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
						int obj = currentNode->primArray[i];
						int nearest;
						dist = dynObjArray[obj].findNearestPrimitive(x, y, z, inv_x, inv_y, inv_z, origin_x, origin_y, origin_z, &nearest);

						// kontroluje se tady a predtim v rayTriangleIntersection
						if(dist >= 0.0f && dist < min_dist) {
							min_dist = dist;
							*nearestObject = obj;
							*nearestTriangle = nearest;
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
//	struct RayPacketJob {
//        ObjectBVH_node* nodeToCheck;
//        char type;
//	};
//
//	struct RayJob {
//        ObjectBVH_node* nodeToCheck;
//        int rayIdx;
//	};

//	float findNearestPrimitive_Packets(Ray* rays, int* nearest, int packet_size) {
//	    int numRays = packet_size * packet_size;
//
//        RayPacketJob rayPacketJobStack[numBvhNodes];
//        int stackSize_packet = 0;
//        RayPacketJob curJob;
//
//        RayJob rayJobStack[numBvhNodes * numRays];
//        int stackSize_ray = 0;
//        RayJob curRayJob;
//
//        Ray* ray;
//
//		float min_dist[numRays];
//		for(int i=0; i<numRays; i++) min_dist[i] = FLT_MAX;
//		float dist;
//		float box_t;
//
//		int* leafPrimArray;
//		int numLeafPrim;
//
//		rayPacketJobStack[0].nodeToCheck = &nodeList[0]; // BVH root
//		rayPacketJobStack[0].type = TYPE_SQUARE;
//		stackSize_packet++;
//
//        while(stackSize_packet > 0) {
//            memcpy(&curJob, &rayPacketJobStack[--stackSize_packet], sizeof(RayPacketJob));
//            //curJob = rayPacketJobStack[--stackSize];
//
//            char newType = 0b00000000;
//            AABB* box = &(curJob.nodeToCheck->box);
//
//            ray = &rays[0];
//            if(box->checkRayIntersect_new(ray->orig.x, ray->orig.y, ray->orig.z, ray->inv_dir.x, ray->inv_dir.y, ray->inv_dir.z)) {
//                newType |= 0b00001000;
//            }
//            ray = &rays[packet_size];
//            if(box->checkRayIntersect_new(ray->orig.x, ray->orig.y, ray->orig.z, ray->inv_dir.x, ray->inv_dir.y, ray->inv_dir.z)) {
//                newType |= 0b00000100;
//            }
//            ray = &rays[numRays - 1];
//            if(box->checkRayIntersect_new(ray->orig.x, ray->orig.y, ray->orig.z, ray->inv_dir.x, ray->inv_dir.y, ray->inv_dir.z)) {
//                newType |= 0b00000010;
//            }
//            ray = &rays[numRays - packet_size];
//            if(box->checkRayIntersect_new(ray->orig.x, ray->orig.y, ray->orig.z, ray->inv_dir.x, ray->inv_dir.y, ray->inv_dir.z)) {
//                newType |= 0b00000001;
//            }
//
//            if(newType == TYPE_SQUARE) {
//                if(curJob.nodeToCheck->isLeaf) {
//                    leafPrimArray = curJob.nodeToCheck->primArray;
//                    numLeafPrim = curJob.nodeToCheck->numPrim;
//                    // All rays with all triangles in the AABB.
//                    for(int r=0; r<numRays; r++) {
//                        ray = &rays[r];
//                        for(int i=0; i<numLeafPrim; i++) {
//                            int pr = leafPrimArray[i];
//                            dist = rayTriangleIntersection(pr, ray->dir.x, ray->dir.y, ray->dir.z,
//                                                           ray->orig.x, ray->orig.y, ray->orig.z);
//                            // kontroluje se tady a predtim v rayTriangleIntersection FIXED
//                            if(dist >= 0.0f && dist < min_dist[r]) {
//                                min_dist[r] = dist;
//                                nearest[r] = pr;
//                            }
//                        }
//                    }
//                } else {
//                    rayPacketJobStack[stackSize_packet].nodeToCheck = &(nodeList[curJob.nodeToCheck->leftChild]);
//                    stackSize_packet++;
//                    rayPacketJobStack[stackSize_packet].nodeToCheck = &(nodeList[curJob.nodeToCheck->rightChild]);
//                    stackSize_packet++;
//                }
//            } else {
////                for(int i=0; i<numRays; i++) {
////                    rayJobStack[stackSize_ray].nodeToCheck = curJob.nodeToCheck;
////                    rayJobStack[stackSize_ray].rayIdx = i;
////                    stackSize_ray++;
////                }
//            }
//        }
//
////        // individual rays.
////        while(stackSize_ray > 0) {
////            memcpy(&curRayJob, &rayJobStack[--stackSize_ray], sizeof(RayJob));
////
////            ray = &rays[curRayJob.rayIdx];
////            if(curRayJob.nodeToCheck->box.checkRayIntersect_new(ray->orig.x, ray->orig.y, ray->orig.z, ray->inv_dir.x, ray->inv_dir.y, ray->inv_dir.z)) {
////                if(curRayJob.nodeToCheck->isLeaf) {
////                    leafPrimArray = curRayJob.nodeToCheck->primArray;
////                    numLeafPrim = curRayJob.nodeToCheck->numPrim;
////
////                    ray = &rays[curRayJob.rayIdx];
////                    for(int i=0; i<numLeafPrim; i++) {
////                        int pr = leafPrimArray[i];
////                        dist = rayTriangleIntersection(pr, ray->dir.x, ray->dir.y, ray->dir.z,
////                                                       ray->orig.x, ray->orig.y, ray->orig.z);
////                        // kontroluje se tady a predtim v rayTriangleIntersection FIXED
////                        if(dist >= 0.0f && dist < min_dist[curRayJob.rayIdx]) {
////                            min_dist[curRayJob.rayIdx] = dist;
////                            nearest[curRayJob.rayIdx] = pr;
////                        }
////                    }
////                } else {
////                    rayJobStack[stackSize_ray].nodeToCheck = &(nodeList[curRayJob.nodeToCheck->leftChild]);
////                    stackSize_ray++;
////                    rayJobStack[stackSize_ray].nodeToCheck = &(nodeList[curRayJob.nodeToCheck->rightChild]);
////                    stackSize_ray++;
////                }
////            }
////        }
//	}

};

#endif // OBJECTBVH_H
