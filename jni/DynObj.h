#ifndef DYNOBJ_H
#define DYNOBJ_H

#include <Structures.h>
#include <Matrix.h>

class DynamicObject {
private:
    // Data structures
	BVH bvh;

	// MODEL DATA //
	// Vertices
	int numVert;
	int vertArray_size;
	Vertex* vertArray;

	// Triangles
	int numTri;
	int triArray_size;
	Triangle* triArray;

	// One model transform for the entire object.
	Matrix4 modelMatrix;


public:
    AABB* objectBox;
    float center[3];

	DynamicObject() {

	}

	void initialize(const int &_numVert, const int &_numTri) {
        numVert = 0;
        vertArray_size = _numVert;
        vertArray = new Vertex[vertArray_size];

        numTri = 0;
        triArray_size = _numTri;
        triArray = new Triangle[triArray_size];
	}

	void updateTriangleCenter(const int &tri_idx) {
	Vertex v1, v2, v3;

	int v1_idx = triArray[tri_idx].v1_idx;
	int v2_idx = triArray[tri_idx].v2_idx;
	int v3_idx = triArray[tri_idx].v3_idx;

	v1.x = vertArray[v1_idx].x;
	v1.y = vertArray[v1_idx].y;
	v1.z = vertArray[v1_idx].z;
	v2.x = vertArray[v2_idx].x;
	v2.y = vertArray[v2_idx].y;
	v2.z = vertArray[v2_idx].z;
	v3.x = vertArray[v3_idx].x;
	v3.y = vertArray[v3_idx].y;
	v3.z = vertArray[v3_idx].z;

	float min_x = v1.x < v2.x ? (v3.x < v1.x ? v3.x : v1.x) : (v3.x < v2.x ? v3.x : v2.x);
	float min_y = v1.y < v2.y ? (v3.y < v1.y ? v3.y : v1.y) : (v3.y < v2.y ? v3.y : v2.y);
	float min_z = v1.z < v2.z ? (v3.z < v1.z ? v3.z : v1.z) : (v3.z < v2.z ? v3.z : v2.z);
	float max_x = v1.x > v2.x ? (v3.x > v1.x ? v3.x : v1.x) : (v3.x > v2.x ? v3.x : v2.x);
	float max_y = v1.y > v2.y ? (v3.y > v1.y ? v3.y : v1.y) : (v3.y > v2.y ? v3.y : v2.y);
	float max_z = v1.z > v2.z ? (v3.z > v1.z ? v3.z : v1.z) : (v3.z > v2.z ? v3.z : v2.z);

	triArray[tri_idx].c[0] = (max_x + min_x) * 0.5f;
	triArray[tri_idx].c[1] = (max_y + min_y) * 0.5f;
	triArray[tri_idx].c[2] = (max_z + min_z) * 0.5f;
}

	int addVertex(const float &_x, const float &_y, const float &_z) {
        if(vertArray_size - numVert == 0) { // Reallocate the array.
            //LOGI("REALLOC: Vertex array.");

            vertArray_size <<= 1;
            Vertex* new_vertArray = new Vertex[vertArray_size];
            memcpy(new_vertArray, vertArray, numVert * sizeof(Vertex));
            /*for(int i=0; i<3*numVert; i++) {
                new_vertArray[i] = vertArray[i];
            }*/
            delete[] vertArray;
            vertArray = new_vertArray;
        }
        vertArray[numVert].x = _x;
        vertArray[numVert].y = _y;
        vertArray[numVert].z = _z;
        numVert++;

        return numVert - 1;
    }

    int addTriangle(const int &v1, const int &v2, const int &v3, const int &mat_idx) {
        // Check if the vertices exist.
        if(v1 >= numVert || v2 >= numVert || v3 >= numVert) {
            return -1; // TODO: Change to a predefined error constant
        }
        // Check if the material exists.
        int matToUse = mat_idx;
        if(triArray_size - numTri == 0) { // Reallocate the array.
            //LOGI("REALLOC: Triangle array.");

            triArray_size <<= 1;
            Triangle* new_triArray = new Triangle[triArray_size];
            memcpy(new_triArray, triArray, numTri * sizeof(Triangle));

            /*for(int i=0; i<4*numTri; i++) {
                new_triArray[i] = triArray[i];
            }*/
            delete[] triArray;
            triArray = new_triArray;
        }
        triArray[numTri].v1_idx = v1;
        triArray[numTri].v2_idx = v2;
        triArray[numTri].v3_idx = v3;
        triArray[numTri].mat = matToUse;

        // Calculate and store the inverse norm of the triangle normal calculated as cross(v3-v1, v2-v1).
        Vertex normal;
        Vertex v3_1, v2_1;
        v3_1.x = vertArray[v3].x - vertArray[v1].x;
        v3_1.y = vertArray[v3].y - vertArray[v1].y;
        v3_1.z = vertArray[v3].z - vertArray[v1].z;

        v2_1.x = vertArray[v2].x - vertArray[v1].x;
        v2_1.y = vertArray[v2].y - vertArray[v1].y;
        v2_1.z = vertArray[v2].z - vertArray[v1].z;

        cross(v3_1, v2_1, normal);
    //	cross(v2_1, v3_1, normal);

        memcpy(&(triArray[numTri].n), &normal, sizeof(Vertex));
        triArray[numTri].n_inv_norm = inv_norm(normal);
        triArray[numTri].n_calculated = 1;

        triArray[numTri].n.x *= triArray[numTri].n_inv_norm;
        triArray[numTri].n.y *= triArray[numTri].n_inv_norm;
        triArray[numTri].n.z *= triArray[numTri].n_inv_norm;

        updateTriangleCenter(numTri);
        numTri++;

        return numTri - 1;
    }

    /// Load a new model transform for this object.
	/**
        Load a new model transform matrix for this object
        that will be used for this object until a new call to
        loadNewTransform.

	*/
	void loadNewTransform(const Matrix4 &newModelMatrix) {
        memcpy(&modelMatrix, &newModelMatrix, sizeof(Matrix4));
	}

    /// Apply the current model transform.
    /**
        Apply the model transform on this object represented
        by the current modelMatrix. Subsequently update all triangle
        centers and rebuild the BVH of this object.
    */
	void applyModelTransform() {
        // Transform all the vertices.
        for(int i=0; i<numVert; i++) {
            vertArray[i] = modelMatrix * vertArray[i];
        }
        // Update triangle centers.
        for(int i=0; i<numTri; i++) {
            updateTriangleCenter(i);
        }
        // Rebuild the BVH for this object.
        bvh.createBVH(vertArray, numVert, triArray, numTri);
        objectBox = bvh.getRootAABB();

        // Update the center (used for binning in the object BVH construction)
        center[0] = objectBox->getCenter().x;
        center[1] = objectBox->getCenter().y;
        center[2] = objectBox->getCenter().z;
	}

	void initBuildTree() {
        bvh.createBVH(vertArray, numVert, triArray, numTri);
        objectBox = bvh.getRootAABB();

        // Update the center (used for binning in the object BVH construction)
        center[0] = objectBox->getCenter().x;
        center[1] = objectBox->getCenter().y;
        center[2] = objectBox->getCenter().z;
	}

    float findNearestPrimitive(const float &x, const float &y, const float &z,
                               const float &inv_x, const float &inv_y, const float &inv_z,
                               const float &origin_x, const float &origin_y, const float &origin_z,
                               int* nearestTriangle) {
        return bvh.findNearestPrimitive_new(x, y, z, inv_x, inv_y, inv_z, origin_x, origin_y, origin_z, nearestTriangle);
    }
};

#endif
