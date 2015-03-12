#ifndef _AABB
#define _AABB

#include <Structures.h>

#define MIN_X   0
#define MIN_Y   1
#define MIN_Z   2
#define MAX_X   3
#define MAX_Y   4
#define MAX_Z   5

class AABB {
public:
    float data[6];
public:
    AABB() {
        data[MIN_X] = FLT_MAX;
        data[MIN_Y] = FLT_MAX;
        data[MIN_Z] = FLT_MAX;
        data[MAX_X] = -FLT_MAX;
        data[MAX_Y] = -FLT_MAX;
        data[MAX_Z] = -FLT_MAX;
    }

    AABB(const float &minX, const float &maxX,
         const float &minY, const float &maxY,
         const float &minZ, const float &maxZ) {
    	data[MIN_X] = minX;
		data[MIN_Y] = minY;
		data[MIN_Z] = minZ;
		data[MAX_X] = maxX;
		data[MAX_Y] = maxY;
		data[MAX_Z] = maxZ;
    }

    void set(const float &minX, const float &maxX,
             const float &minY, const float &maxY,
             const float &minZ, const float &maxZ) {
       	data[MIN_X] = minX;
   		data[MIN_Y] = minY;
   		data[MIN_Z] = minZ;
   		data[MAX_X] = maxX;
   		data[MAX_Y] = maxY;
   		data[MAX_Z] = maxZ;
    }

    void reset() {
    	data[MIN_X] = FLT_MAX;
		data[MIN_Y] = FLT_MAX;
		data[MIN_Z] = FLT_MAX;
		data[MAX_X] = -FLT_MAX;
		data[MAX_Y] = -FLT_MAX;
		data[MAX_Z] = -FLT_MAX;
    }

    void contain(const AABB &box) {
        if (box.data[MIN_X] < data[MIN_X]) data[MIN_X] = box.data[MIN_X];
        if (box.data[MIN_Y] < data[MIN_Y]) data[MIN_Y] = box.data[MIN_Y];
        if (box.data[MIN_Z] < data[MIN_Z]) data[MIN_Z] = box.data[MIN_Z];
        if (box.data[MAX_X] > data[MAX_X]) data[MAX_X] = box.data[MAX_X];
        if (box.data[MAX_Y] > data[MAX_Y]) data[MAX_Y] = box.data[MAX_Y];
        if (box.data[MAX_Z] > data[MAX_Z]) data[MAX_Z] = box.data[MAX_Z];
    }

    void contain(const Vertex &vertex) {
    	if (vertex.x < data[MIN_X]) data[MIN_X] = vertex.x;
		if (vertex.y < data[MIN_Y]) data[MIN_Y] = vertex.y;
		if (vertex.z < data[MIN_Z]) data[MIN_Z] = vertex.z;
		if (vertex.x > data[MAX_X]) data[MAX_X] = vertex.x;
		if (vertex.y > data[MAX_Y]) data[MAX_Y] = vertex.y;
		if (vertex.z > data[MAX_Z]) data[MAX_Z] = vertex.z;
    }
/*
    void AABB::containPrimitivesFromTo(const vector<Primitive*>& primitiveList, int inclFrom, int notInclTo) {
        Primitive* pr;
        AABB box;
        for (int i = inclFrom; i < notInclTo; i++) {
            pr = primitiveList[i];
            box.set(pr->minx, pr->maxx, pr->miny, pr->maxy, pr->minz, pr->maxz);
            contain(box);
        }
    }

    void AABB::containPrimitives(const std::vector<Primitive*>& primitiveList) {
        Primitive* pr;
        AABB box;
        for (int i = 0; i < primitiveList.size(); i++) {
            pr = primitiveList[i];
            box.set(pr->minx, pr->maxx, pr->miny, pr->maxy, pr->minz, pr->maxz);
            contain(box);
        }
    }

    void AABB::containPrimitive(Primitive* pr) {
        AABB box;
        box.set(pr->minx, pr->maxx, pr->miny, pr->maxy, pr->minz, pr->maxz);
        contain(box);
    }
*/
    int getWidestAxis() {
    	float dx = data[MAX_X] - data[MIN_X];
    	float dy = data[MAX_Y] - data[MIN_Y];
    	float dz = data[MAX_Z] - data[MIN_Z];
    	return dx > dy ? (dx > dz ? 0 : 2) : (dy > dz ? 1 : 2);
    }

    float getXsize() {
        return data[MAX_X] - data[MIN_X];
    }
    float getYsize() {
        return data[MAX_Y] - data[MIN_Y];
    }
    float getZsize() {
        return data[MAX_Z] - data[MIN_Z];
    }

    Vertex getCenter() {
        Vertex center = {(data[MAX_X] + data[MIN_X]) / 2.0f,
                         (data[MAX_Y] + data[MIN_Y]) / 2.0f,
                         (data[MAX_Z] + data[MIN_Z]) / 2.0f};
        return center;
    }

    float getArea() {
    	float dx = data[MAX_X] - data[MIN_X];
		float dy = data[MAX_Y] - data[MIN_Y];
		float dz = data[MAX_Z] - data[MIN_Z];
        return 2.0f * (dx*dy + dx*dz + dy*dz);
    }

    bool checkRayIntersect(const float &fromX, const float &fromY, const float &fromZ,
    		 	 	 	 	 	 const float &x, const float &y, const float &z) {
        const float dirEps = 1e-8f;
        float minx, maxx, tmin, tmax;
        float x_recip = 1.0f / x;
        float y_recip = 1.0f / y;
        float z_recip = 1.0f / z;

        if (fabs(x) < dirEps) {
            if (data[MIN_X] < fromX && data[MAX_X] > fromX) {
                minx = -FLT_MAX;
                maxx = FLT_MAX;
            }
            else
                return false;
        }
        else {
            float t1 = (data[MIN_X] - fromX) * x_recip;
            float t2 = (data[MAX_X] - fromX) * x_recip;
            if (t1 < t2) {
                minx = t1;
                maxx = t2;
            }
            else {
                minx = t2;
                maxx = t1;
            }
            if (maxx < 0.0f)
                return false;
        }

        tmin = minx;
        tmax = maxx;

        if (fabs(y) < dirEps) {
            if (data[MIN_Y] < fromY && data[MAX_Y] > fromY) {
                minx = -FLT_MAX;
                maxx = FLT_MAX;
            }
            else
                return false;
        }
        else {
            float t1 = (data[MIN_Y] - fromY) * y_recip;
            float t2 = (data[MAX_Y] - fromY) * y_recip;
            if (t1 < t2) {
                minx = t1;
                maxx = t2;
            }
            else {
                minx = t2;
                maxx = t1;
            }
            if (maxx < 0.0f)
                return false;
        }

        if (minx > tmin)
            tmin = minx;
        if (maxx < tmax)
            tmax = maxx;

        if (fabs(z) < dirEps) {
            if (data[MIN_Z] < fromZ && data[MAX_Z] > fromZ) {
                minx = -FLT_MAX;
                maxx = FLT_MAX;
            }
            else
                return false;
        }
        else {
            float t1 = (data[MIN_Z] - fromZ) * z_recip;
            float t2 = (data[MAX_Z] - fromZ) * z_recip;
            if (t1 < t2) {
                minx = t1;
                maxx = t2;
            }
            else {
                minx = t2;
                maxx = t1;
            }
            if (maxx < 0.0f)
                return false;
        }

        if (minx > tmin)
            tmin = minx;

        if (maxx < tmax)
            tmax = maxx;

        return tmin <= tmax; // yes, intersection was found
    }

    bool checkRayIntersect_new(const float &fromX, const float &fromY, const float &fromZ,
							   const float &inv_x, const float &inv_y, const float &inv_z) {
    	// lb is the corner of AABB with minimal coordinates - left bottom, rt is maximal corner
    	// r.org is origin of ray
    	float t1 = (data[MIN_X] - fromX)*inv_x;
    	float t2 = (data[MAX_X] - fromX)*inv_x;
    	float t3 = (data[MIN_Y] - fromY)*inv_y;
    	float t4 = (data[MAX_Y] - fromY)*inv_y;
    	float t5 = (data[MIN_Z] - fromZ)*inv_z;
    	float t6 = (data[MAX_Z] - fromZ)*inv_z;

    	float tmin = fmax(fmax(fmin(t1, t2), fmin(t3, t4)), fmin(t5, t6));
    	float tmax = fmin(fmin(fmax(t1, t2), fmax(t3, t4)), fmax(t5, t6));

    	// if tmax < 0, ray (line) is intersecting AABB, but whole AABB is behind us
    	if (tmax < 0)
    	{
    	    //t = tmax;
    	    //return -1.0f;
    		return false;
    	}

    	// if tmin > tmax, ray doesn't intersect AABB
    	if (tmin > tmax)
    	{
    	    //t = tmax;
    	    //return -1.0f;
    		return false;
    	}

    	//t = tmin;
    	//return tmin;
    	return true;
    }

    bool checkRayIntersect_new_t(const float &fromX, const float &fromY, const float &fromZ,
							   const float &inv_x, const float &inv_y, const float &inv_z, float* t) {
    	// lb is the corner of AABB with minimal coordinates - left bottom, rt is maximal corner
    	// r.org is origin of ray
    	float t1 = (data[MIN_X] - fromX)*inv_x;
    	float t2 = (data[MAX_X] - fromX)*inv_x;
    	float t3 = (data[MIN_Y] - fromY)*inv_y;
    	float t4 = (data[MAX_Y] - fromY)*inv_y;
    	float t5 = (data[MIN_Z] - fromZ)*inv_z;
    	float t6 = (data[MAX_Z] - fromZ)*inv_z;

    	float tmin = fmax(fmax(fmin(t1, t2), fmin(t3, t4)), fmin(t5, t6));
    	float tmax = fmin(fmin(fmax(t1, t2), fmax(t3, t4)), fmax(t5, t6));

    	// if tmax < 0, ray (line) is intersecting AABB, but whole AABB is behind us
    	if (tmax < 0)
    	{
    	    *t = tmax;
    	    //return -1.0f;
    		return false;
    	}

    	// if tmin > tmax, ray doesn't intersect AABB
    	if (tmin > tmax)
    	{
            *t = tmax;
    	    //return -1.0f;
    		return false;
    	}

    	*t = tmin;
    	//return tmin;
    	return true;
    }
};

#endif
