#ifndef _TRIANGLEBVH_H
#define _TRIANGLEBVH_H

#include <stdlib.h>
//#include <arm_neon.h>

#include <AABB.h>
#include <Structures.h>
#include <Triangle.h>
#include <BVHStructures.h>
//#include <RigidObject.h>

// A BVH for triangle primitives.
class TriangleBVH {
    // Construction data.
	float C_t, C_i;
	int divBuckets, numPrimInLeaf;
	nodeBucket* bucketList;

	// Scene data
	Triangle* objArray;
	int objArray_size;
	int numObj;

	int* primIndexArray; // An array of primitive indices.

	// BVH node list
	BVH_node* nodeList;
	int nodeList_size;
	int numBvhNodes;
public:
	TriangleBVH() : C_t(1.0f), C_i(1.0f), divBuckets(10), numPrimInLeaf(3), nodeList(NULL), nodeList_size(0), numBvhNodes(0) {
		bucketList = new nodeBucket[divBuckets];
	}

	~TriangleBVH() {
		delete[] primIndexArray;
		//delete[] nodeStack;
		delete[] bucketList;
		delete[] nodeList;
	}

	void setParams(const float &_C_t, const float &_C_i, const int &_numPrimInLeaf) {
        numPrimInLeaf = _numPrimInLeaf;
        C_t = _C_t;
        C_i = _C_i;
	}
private:
	void clearBucketList() {
		for(unsigned int i=0; i<divBuckets; i++) {
			bucketList[i].leftBox.reset();
			bucketList[i].rightBox.reset();
			bucketList[i].leftPrimitiveCount = 0;
			bucketList[i].rightPrimitiveCount = 0;
			bucketList[i].primIndexArray = NULL;
			bucketList[i].numPrim = 0;
		}
	}

	void containPrimitives(AABB &box, int* &_primIndexArray, const int &_numPrim) {
		for(unsigned int i=0; i<_numPrim; i++) {
            box.contain(objArray[_primIndexArray[i]].v1);
            box.contain(objArray[_primIndexArray[i]].v2);
            box.contain(objArray[_primIndexArray[i]].v3);
		}
	}

	int getDividingIndexBucketing(int* &_primIndexArray, const int &_primCount,
								  const float &_parentArea, const int &_divAxis) {
		float bestSahValue = FLT_MAX;
		int dividingIndex = 1;

		int curDivBuckets = divBuckets > _primCount ? _primCount : divBuckets;

//		LOGI("curDivBuckets: %d.", curDivBuckets);

		clearBucketList();

		float left, right;
		left = objArray[_primIndexArray[0]].c[_divAxis];
		right = objArray[_primIndexArray[_primCount - 1]].c[_divAxis];

//		LOGI("Got left and right.");

		float bucketSize = (right - left) / (float)curDivBuckets;
		// Distribute primitives into buckets
		if(bucketSize > 0.0f) {
			int bucketIndex = -1;
			int new_index;

			for(unsigned int i=0; i<_primCount; i++) {
				new_index = (int)((objArray[_primIndexArray[i]].c[_divAxis] - left) / bucketSize);
				new_index = new_index >= curDivBuckets ? curDivBuckets - 1 : new_index < 0 ? 0 : new_index;

				if(new_index > bucketIndex) { // Move on to the next bucket.
					bucketIndex = new_index;
					bucketList[bucketIndex].primIndexArray = _primIndexArray + i; // The primitive list of this bucket starts here.
					bucketList[bucketIndex].numPrim++; // And has one more primitive.
				} else {
					bucketList[bucketIndex].numPrim++; // And has one more primitive.
				}
			}
		} else {
			bucketList[0].primIndexArray = _primIndexArray;
			bucketList[0].numPrim = _primCount;
		}

//		LOGI("Distributed.");

		// Contain from left to right.
		containPrimitives(bucketList[0].leftBox, bucketList[0].primIndexArray, bucketList[0].numPrim);
		bucketList[0].leftPrimitiveCount = bucketList[0].numPrim;
//		LOGI("Contained left 0.");
		for(unsigned int i=1; i<curDivBuckets; i++) {
			bucketList[i].leftBox.contain(bucketList[i-1].leftBox);
			containPrimitives(bucketList[i].leftBox, bucketList[i].primIndexArray, bucketList[i].numPrim);
			bucketList[i].leftPrimitiveCount = bucketList[i].numPrim + bucketList[i-1].leftPrimitiveCount;
//			LOGI("Contained left %d.", i);
		}
		// Contain from right to left.
		containPrimitives(bucketList[curDivBuckets-1].leftBox, bucketList[curDivBuckets-1].primIndexArray, bucketList[curDivBuckets-1].numPrim);
		bucketList[curDivBuckets-1].rightPrimitiveCount = bucketList[curDivBuckets-1].numPrim;
//		LOGI("Contained right 0.");
		for(int i=curDivBuckets - 2; i>=0; i--) {
			bucketList[i].rightBox.contain(bucketList[i+1].rightBox);
			containPrimitives(bucketList[i].rightBox, bucketList[i].primIndexArray, bucketList[i].numPrim);
			bucketList[i].rightPrimitiveCount = bucketList[i].numPrim + bucketList[i+1].rightPrimitiveCount;
//			LOGI("Contained right %d.", i);
		}

//		LOGI("Contained both ways.");

		// Find the best dividing plane.
		for(unsigned int i=0; i<curDivBuckets-1; i++) {
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

	bool comparePrimitives(const int &prim_1, const int &prim_2, const int &_axis) {
	    return objArray[prim_1].c[_axis] <= objArray[prim_2].c[_axis];
	}

	void sortPrimitiveArray(int* &_primIndexArray, const int &_numPrim, const int &_axis) {
		int numSwaps = 1;
		while(numSwaps > 0) {
			numSwaps = 0;
			for(unsigned int i=0; i<_numPrim-1; i++) {
				if(!comparePrimitives(_primIndexArray[i], _primIndexArray[i + 1], _axis)) {
					int temp = _primIndexArray[i];
					_primIndexArray[i] = _primIndexArray[i + 1];
					_primIndexArray[i + 1] = temp;
					numSwaps++;
				}
			}
		}
	}
public:
	void createBVH(Triangle* _objArray, const int &_numObj) {
//        nodeList = (BVH_node*)(aligned_alloc(32, sizeof(BVH_node) * (_numTri << 1)));
        if(_numObj == 0) {
            nodeList_size = 1;
            nodeList = new BVH_node[nodeList_size];
            numBvhNodes = 1;
            return;
        }

        objArray = _objArray;
		numObj = _numObj;

		// An array of object indices, to be sorted and partitioned into nodes.
		primIndexArray = new int[numObj];
		for(unsigned int i=0; i<numObj; i++) {
			primIndexArray[i] = i;
		}

        // If the number of primitives is lower than the number of primitives for a leaf.
		if(numObj <= numPrimInLeaf) {
            nodeList_size = 1;
            nodeList = new BVH_node[nodeList_size];
            numBvhNodes = 1;

            BVH_node* root = &nodeList[0];

            root->primIndexArray = primIndexArray;
            root->numPrim_isLeaf = numObj;
            containPrimitives(root->box, root->primIndexArray, root->numPrim_isLeaf);

			return;
		}

        if(nodeList == NULL) {
//            LOGI("Creating new nodeList of size %d.", numObj << 1);
            nodeList_size = numObj << 1;
            nodeList = new BVH_node[nodeList_size];
        }

		BVH_node* root = &nodeList[0];
		numBvhNodes = 2; // Skip one, Jacco said.

		root->primIndexArray = primIndexArray;
		root->numPrim_isLeaf = numObj;
		containPrimitives(root->box, root->primIndexArray, root->numPrim_isLeaf);
//		LOGI("Root created and primitives contained.");

		int* toSubdivideStack = new int[nodeList_size];
		int toSubdivideStack_size = nodeList_size;
		int numToSubdivide = 0;
		BVH_node* curNode;
		int* primSubList;
		int numPrimSubList;

		toSubdivideStack[numToSubdivide] = 0;
		numToSubdivide++;

		while(numToSubdivide > 0) {
//            LOGI("Checking node %d.", toSubdivideStack[numToSubdivide-1]);
			curNode = &nodeList[toSubdivideStack[--numToSubdivide]];

			if(curNode->numPrim_isLeaf > numPrimInLeaf) { // Subdivide further.
				if(nodeList_size - numBvhNodes < 2) {
//					LOGI("Reallocating nodeList.");
					// Reallocate the node array.
					nodeList_size <<= 1;
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

//				LOGI("Two children created.");

				int divAxis = curNode->box.getWidestAxis();
//				LOGI("divAxis obtained.");

				// The primitive sub list to work with.
				primSubList = curNode->primIndexArray;
				numPrimSubList = curNode->numPrim_isLeaf;
//				LOGI("primSubList obtained.");

				sortPrimitiveArray(primSubList, numPrimSubList, divAxis);

//				LOGI("Primitive sublist sorted.");

				int divIndex = getDividingIndexBucketing(primSubList, numPrimSubList, curNode->box.getArea(), divAxis);
//                int divIndex = 1;

//				LOGI("divIndex: %d, numPrim: %d", divIndex, numPrimSubList);

				leftChild->primIndexArray = primSubList;
				leftChild->numPrim_isLeaf = divIndex;
				rightChild->primIndexArray = primSubList + divIndex;
				rightChild->numPrim_isLeaf = numPrimSubList - divIndex;

				// left child, contain primitives in the box
				containPrimitives(leftChild->box, leftChild->primIndexArray, leftChild->numPrim_isLeaf);
				// right child, contain primitives in the box
				containPrimitives(rightChild->box, rightChild->primIndexArray, rightChild->numPrim_isLeaf);

				curNode->primIndexArray = NULL;
				curNode->numPrim_isLeaf = 0;
				/*curNode->leftChild = leftChild;
				curNode->rightChild = rightChild;*/
				curNode->leftChild = numBvhNodes - 2;
				//curNode->rightChild = numBvhNodes - 1;

				if(toSubdivideStack_size - numToSubdivide < 2) {
//					LOGI("Reallocating the toSubdivide stack.");
					// Reallocate the stack.
					toSubdivideStack_size <<= 1;
					int* new_toSubdivideStack = new int[toSubdivideStack_size];
					memcpy(new_toSubdivideStack, toSubdivideStack, numToSubdivide * sizeof(int));
					delete[] toSubdivideStack;
					toSubdivideStack = new_toSubdivideStack;
				}
				toSubdivideStack[numToSubdivide] = numBvhNodes - 2;
				numToSubdivide++;
				toSubdivideStack[numToSubdivide] = numBvhNodes - 1;
				numToSubdivide++;

//				LOGI("Added children to stack.");
			}
		}
		//nodeStack = new BVH_node*[numBvhNodes];

        delete[] toSubdivideStack;

		LOGI("BVH DONE, %d nodes.", numBvhNodes);
	}

	float traverse(const Ray &ray, int* nearestTri) {
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
			bool box_t_is = currentNode->box.intersect(ray);
			if(box_t_is) {
				if(!currentNode->numPrim_isLeaf) {
					nodeStack[nodeStackSize] = &nodeList[currentNode->leftChild];
					nodeStackSize++;
					nodeStack[nodeStackSize] = &nodeList[currentNode->leftChild + 1];
					nodeStackSize++;
				} else {
					leafPrimArray = currentNode->primIndexArray;
					numLeafPrim = currentNode->numPrim_isLeaf;

					for(unsigned int i=0; i<numLeafPrim; i++) {
						int pr = currentNode->primIndexArray[i];
						//dist = rayTriangleIntersection(pr, ray);
						dist = objArray[pr].intersect(ray);

						// kontroluje se tady a predtim v rayTriangleIntersection
						if(dist >= 0.0f && dist < min_dist) { //
							min_dist = dist;
//							*nearest = pr;
                            *nearestTri = pr;
						}
					}
				}
			}
		}
		return min_dist;
	}

//	#define TYPE_SQUARE 0b00001111
//	#define TYPE_TRI_BOTTOM_LEFT 0b00001101
//	#define TYPE_TRI_BOTTOM_RIGHT 0b00001110
//	#define TYPE_TRI_TOP_RIGHT 0b00000111
//	#define TYPE_TRI_TOP_LEFT 0b00001011
//	#define TYPE_LINE_LEFT 0b00001001
//	#define TYPE_LINE_RIGHT 0b00000110
//	#define TYPE_LINE_BOTTOM 0b00001100
//	#define TYPE_LINE_TOP 0b00000011
//	#define TYPE_SINGLE 0b00010000
//	#define TYPE_MISS 0b00000000
//
//    // A packet of rays
//	struct RayPacketJob {
//        BVH_node* nodeToCheck;
//        char type;
//	};
//
//	struct RayJob {
//        BVH_node* nodeToCheck;
//        int rayIdx;
//	};
//
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
//                if(curJob.nodeToCheck->numPrim_isLeaf) {
//                    leafPrimArray = curJob.nodeToCheck->primArray;
//                    numLeafPrim = curJob.nodeToCheck->numPrim_isLeaf;
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
//                    rayPacketJobStack[stackSize_packet].nodeToCheck = &(nodeList[curJob.nodeToCheck->leftChild + 1]);
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
//	}

    AABB* getRootAABB() {
        return &(nodeList[0].box);
    }
};

#endif



















