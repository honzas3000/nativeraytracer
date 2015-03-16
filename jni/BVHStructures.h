#ifndef _BVHSTRUCTURES_H
#define _BVHSTRUCTURES_H

#include <AABB.h>

// BVH structures /////////////////////
struct BVH_node {
    AABB box; // 24B
	int* primIndexArray; // Pointer into the primitive index array.
	int numPrim_isLeaf; // If the number of primitives is 0, the node is interior.
	int leftChild; // Right child is leftChild + 1.

	BVH_node() : leftChild(-1), numPrim_isLeaf(0) {}
};

struct nodeBucket {
	int* primIndexArray;
	int numPrim;

	int leftPrimitiveCount;
	int rightPrimitiveCount;
	AABB leftBox, rightBox;

	nodeBucket() : numPrim(0), leftPrimitiveCount(0), rightPrimitiveCount(0) {}
};

#endif // _BVHSTRUCTURES_H
