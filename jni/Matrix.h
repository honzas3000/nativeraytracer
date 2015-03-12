#ifndef MATRIX_H
#define MATRIX_H

#include <Structures.h>

struct Matrix4 {
    float row0[4];
    float row1[4];
    float row2[4];
    float row3[4];
};

inline Vertex operator * (const Matrix4 &m4, const Vertex &v3) {
    Vertex res;
    res.x = m4.row0[0]*v3.x + m4.row0[1]*v3.y + m4.row0[2]*v3.z + m4.row0[3];
    res.y = m4.row1[0]*v3.x + m4.row1[1]*v3.y + m4.row1[2]*v3.z + m4.row1[3];
    res.z = m4.row2[0]*v3.x + m4.row2[1]*v3.y + m4.row2[2]*v3.z + m4.row2[3];
    return res;
}

#endif // MATRIX_H
