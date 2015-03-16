#ifndef MATRIX_H
#define MATRIX_H

#include <string>
#include <Structures.h>

//struct Matrix4 {
//    float row0[4];
//    float row1[4];
//    float row2[4];
//    float row3[4];
//};

// Column-major.
struct Matrix4 {
    float cols[16];
};

inline Vertex operator * (const Matrix4 &mat, const Vertex &v3) {
    Vertex res;
    res.x = mat.cols[0] * v3.x + mat.cols[4] * v3.y + mat.cols[8] * v3.z + mat.cols[12];
    res.y = mat.cols[1] * v3.x + mat.cols[5] * v3.y + mat.cols[9] * v3.z + mat.cols[13];
    res.z = mat.cols[2] * v3.x + mat.cols[6] * v3.y + mat.cols[10] * v3.z + mat.cols[14];
    return res;
}

void inverseMatrix(const Matrix4 &m, Matrix4 *resultMat);
void printMatrix(const Matrix4 &m);

const Matrix4 identityM = {{1.0f, 0.0f, 0.0f, 0.0f,
                          0.0f, 1.0f, 0.0f, 0.0f,
                          0.0f, 0.0f, 1.0f, 0.0f,
                          0.0f, 0.0f, 0.0f, 1.0f}};

const Matrix4 identityM2 = {{2.0f, 0.0f, 0.0f, 0.0f,
                          0.0f, 2.0f, 0.0f, 0.0f,
                          0.0f, 0.0f, 2.0f, 0.0f,
                          3.0f, 0.0f, 0.0f, 1.0f}};

inline Ray transformRay(const Matrix4 &mat, const Ray* ray) {
    Ray trRay;
    Vertex dir = ray->dir;
    Vertex orig = ray->orig;

    // Ray direction only gets rotated and scaled.
    dir.x = mat.cols[0] * ray->dir.x + mat.cols[4] * ray->dir.y + mat.cols[8] * ray->dir.z;
    dir.y = mat.cols[1] * ray->dir.x + mat.cols[5] * ray->dir.y + mat.cols[9] * ray->dir.z;
    dir.z = mat.cols[2] * ray->dir.x + mat.cols[6] * ray->dir.y + mat.cols[10] * ray->dir.z;

    // Ray origin gets rotated, scaled and translated.
    orig.x = mat.cols[0] * ray->orig.x + mat.cols[4] * ray->orig.y + mat.cols[8] * ray->orig.z + mat.cols[12];
    orig.y = mat.cols[1] * ray->orig.x + mat.cols[5] * ray->orig.y + mat.cols[9] * ray->orig.z + mat.cols[13];
    orig.z = mat.cols[2] * ray->orig.x + mat.cols[6] * ray->orig.y + mat.cols[10] * ray->orig.z + mat.cols[14];

//    normalize(dir);

    // Recalculate the inverse direction of the ray.
    memcpy(&(trRay.dir), &dir, sizeof(Vertex));
    memcpy(&(trRay.orig), &orig, sizeof(Vertex));
    trRay.inv_dir.x = 1.0f / dir.x;
    trRay.inv_dir.y = 1.0f / dir.y;
    trRay.inv_dir.z = 1.0f / dir.z;

    return trRay;
}

#endif // MATRIX_H
