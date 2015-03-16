#ifndef _TRIANGLE_H
#define _TRIANGLE_H

#include <Structures.h>
#include <AABB.h>

struct Triangle {
//	int v1_idx;
//	int v2_idx;
//	int v3_idx;
    Vertex v1;
    Vertex v2;
    Vertex v3;
	int mat;
	float c[3];
	Vertex n;
	float n_inv_norm;
	int n_calculated; // 0 or 1

	void containSelfInAABB(AABB &box) {
        box.contain(v1);
        box.contain(v2);
        box.contain(v3);
    }

    float intersect(const Ray &ray) {
        Vertex vec1 = {v2.x - v1.x, v2.y - v1.y, v2.z - v1.z};
        Vertex vec2 = {v3.x - v1.x, v3.y - v1.y, v3.z - v1.z};

        Vertex D = ray.dir;
        Vertex P;
        cross(D, vec2, P);

        float det = dot(vec1, P);
        if(det > -0.000001f && det < 0.000001f) return -1.0f;
        float inv_det = 1.0f / det;

        Vertex O = ray.orig;
        Vertex T = O - v1;

        float u = dot(T, P) * inv_det;
        if(u < 0.0f || u > 1.0f) return -1.0f;

        Vertex Q;
        cross(T, vec1, Q);

        float v = dot(D, Q) * inv_det;
        if(v < 0.0f || u + v > 1.0f) return -1.0f;

        return dot(vec2, Q) * inv_det; // t
    }
};

#endif // _TRIANGLE_H
