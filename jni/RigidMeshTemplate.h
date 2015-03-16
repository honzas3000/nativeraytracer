#ifndef _RIGIDMESHTEMPLATE_H
#define _RIGIDMESHTEMPLATE_H

#include <Structures.h>
#include <Triangle.h>
#include <AABB.h>
#include <TriangleBVH.h>

class RigidMeshTemplate {
//    // Vertices
//	int numVert;
//	int vertArray_size;
//	Vertex* vertArray;
	// Triangles
	int numTri;
	int triArray_size;
	Triangle* triArray;

	// The BVH of this mesh, only needs to be built once without any
	// model transform applied.
	TriangleBVH bvh;
public:
    RigidMeshTemplate() : numTri(0) {
        bvh.setParams(2.0f, 1.0f, 3);
    }

    ~RigidMeshTemplate() {
//        delete[] vertArray;
        delete[] triArray;
    }

    void init(const int _numVert, const int _numTri) {
//        vertArray_size = _numVert;
//        vertArray = new Vertex[vertArray_size];

        triArray_size = _numTri;
        triArray = new Triangle[triArray_size];
    }

    float traverse(const Ray &ray, int* nearestTri) {
        return bvh.traverse(ray, nearestTri);
    }

    void buildBVH() {
        LOGI("Mesh has %d triangles, starting BVH build.", numTri);
        bvh.createBVH(triArray, numTri);

        LOGI("MeshTemplate AABB:");
        bvh.getRootAABB()->print();
    }

    AABB* getRootAABB() {
        return bvh.getRootAABB();
    }

    void updateTriangleCenter(const int &tri_idx) {
        Vertex v1, v2, v3;

        v1 = triArray[tri_idx].v1;
        v2 = triArray[tri_idx].v2;
        v3 = triArray[tri_idx].v3;

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

//	int addVertex(const float &_x, const float &_y, const float &_z) {
//        if(vertArray_size - numVert == 0) { // Reallocate the array.
//            vertArray_size <<= 1;
//            Vertex* new_vertArray = new Vertex[vertArray_size];
//            memcpy(new_vertArray, vertArray, numVert * sizeof(Vertex));
//            delete[] vertArray;
//            vertArray = new_vertArray;
//        }
//        vertArray[numVert].x = _x;
//        vertArray[numVert].y = _y;
//        vertArray[numVert].z = _z;
//        numVert++;
//
//        return numVert - 1;
//    }

    int addTriangle(const float &v1x, const float &v1y, const float &v1z,
                    const float &v2x, const float &v2y, const float &v2z,
                    const float &v3x, const float &v3y, const float &v3z, const int &mat_idx) {
//        // Check if the vertices exist.
//        if(v1 >= numVert || v2 >= numVert || v3 >= numVert) {
//            return -1; // TODO: Change to a predefined error constant
//        }
        if(triArray_size - numTri == 0) { // Reallocate the array.
            triArray_size <<= 1;
            Triangle* new_triArray = new Triangle[triArray_size];
            memcpy(new_triArray, triArray, numTri * sizeof(Triangle));
            delete[] triArray;
            triArray = new_triArray;
        }

        Vertex v1 = {v1x, v1y, v1z};
        Vertex v2 = {v2x, v2y, v2z};
        Vertex v3 = {v3x, v3y, v3z};

        triArray[numTri].v1 = v1;
        triArray[numTri].v2 = v2;
        triArray[numTri].v3 = v3;
        triArray[numTri].mat = mat_idx;

        // Calculate and store the inverse norm of the triangle normal calculated as cross(v3-v1, v2-v1).
        Vertex normal;
        Vertex v3_1, v2_1;
        v3_1 = triArray[numTri].v3 - triArray[numTri].v1;
        v2_1 = triArray[numTri].v2 - triArray[numTri].v1;

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
};

#endif
