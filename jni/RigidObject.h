#ifndef _RIGIDOBJECT_H
#define _RIGIDOBJECT_H

#include <Structures.h>
#include <Matrix.h>
#include <RigidMeshTemplate.h>
#include <AABB.h>

class RigidObject {
    // Mesh geometry template
    RigidMeshTemplate* meshArray;
    int meshTmpl;
	// One model transform for the entire object. Used to transform rays.
	Matrix4 mMat, inv_mMat;

	inline void refitWorldAABB() {
	    box_world.reset();
        // Transform the vertices of the mesh template's root AABB.
        AABB* box_root = meshArray[meshTmpl].getRootAABB();
        Vertex corners[8];

        Vertex v0 = {box_root->data[0], box_root->data[1], box_root->data[2]};
        Vertex v1 = {box_root->data[0], box_root->data[1], box_root->data[5]};
        Vertex v2 = {box_root->data[0], box_root->data[4], box_root->data[2]};
        Vertex v3 = {box_root->data[0], box_root->data[4], box_root->data[5]};
        Vertex v4 = {box_root->data[3], box_root->data[1], box_root->data[2]};
        Vertex v5 = {box_root->data[3], box_root->data[1], box_root->data[5]};
        Vertex v6 = {box_root->data[3], box_root->data[4], box_root->data[2]};
        Vertex v7 = {box_root->data[3], box_root->data[4], box_root->data[5]};

        corners[0] = mMat * v0;
        corners[1] = mMat * v1;
        corners[2] = mMat * v2;
        corners[3] = mMat * v3;
        corners[4] = mMat * v4;
        corners[5] = mMat * v5;
        corners[6] = mMat * v6;
        corners[7] = mMat * v7;

        // Contain the transformed mesh template root AABB.
        for(unsigned int i=0; i<8; i++) {
            LOGI("corners[%d]: %f, %f, %f.", i, corners[i].x, corners[i].y, corners[i].z);
            box_world.contain(corners[i]);
        }

        // Calculate the center of this object.
        Vertex center = box_world.getCenter();
        c[0] = center.x;
        c[1] = center.y;
        c[2] = center.z;
	}
public:
    // The box that wraps the transformed bvh of this object.
	AABB box_world;
	float c[3];

	void transformAndRefitAABB(const Matrix4 &new_mMat) {
        memcpy(&mMat, &new_mMat, sizeof(Matrix4));
        // Calculate the inverse of the model matrix, for ray traversal.
        inverseMatrix(mMat, &inv_mMat);

        printMatrix(inv_mMat);

        // Refit the world AABB using the new transform.
        refitWorldAABB();
        LOGI("RigidOBject world AABB:");
        box_world.print();
    }

    void supplyMeshArray(RigidMeshTemplate* _meshArray) {
        meshArray = _meshArray;
        LOGI("Mesh array supplied to object.");
    }

    // Constructor
    RigidObject() {}
    RigidObject(const int _meshTmpl) : meshTmpl(_meshTmpl) {
        LOGI("RigidObject constructor, meshID: %d", meshTmpl);
    }
    // Destructor
    ~RigidObject() {}

    float intersect(const Ray &ray, int* nearestTri) {
        float t = FLT_MAX;
        if(box_world.intersect(ray)) { // TODO this might be redundant.
            // Transform the ray with the inverse of the model transform.
            Ray trRay = transformRay(inv_mMat, &ray);
//            Ray trRay = ray;
//            LOGI("Ray: dir: %f, %f, %f, orig: %f, %f, %f", ray.dir.x, ray.dir.y, ray.dir.z, ray.orig.x, ray.orig.y, ray.orig.z);
//            LOGI("trRay: dir: %f, %f, %f, orig: %f, %f, %f", trRay.dir.x, trRay.dir.y, trRay.dir.z, trRay.orig.x, trRay.orig.y, trRay.orig.z);

            int _nearestTri = -1;
            t = meshArray[meshTmpl].traverse(trRay, &_nearestTri);

            if(_nearestTri >= 0) {
                *nearestTri = _nearestTri;
            }

            // TODO: Recalculate t using the model transform somehow.

        }

        return t;
    }

    void containSelfInAABB(AABB &box) {
        box.contain(box_world);
    }
};

#endif // _RIGIDOBJECT_H
