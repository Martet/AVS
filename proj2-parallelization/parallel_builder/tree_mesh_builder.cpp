/**
 * @file    tree_mesh_builder.cpp
 *
 * @author  MARTIN ZMITKO <xlogin00@stud.fit.vutbr.cz>
 *
 * @brief   Parallel Marching Cubes implementation using OpenMP tasks + octree early elimination
 *
 * @date    15 December 2023
 **/

#include <iostream>
#include <math.h>
#include <limits>

#include "tree_mesh_builder.h"

TreeMeshBuilder::TreeMeshBuilder(unsigned gridEdgeSize)
    : BaseMeshBuilder(gridEdgeSize, "Octree")
{

}

unsigned TreeMeshBuilder::evaluateCube(const Vec3_t<float> &cubeOffset, const unsigned a, const ParametricScalarField &field)
{
    const unsigned newA = a / 2;
    const Vec3_t<float> center(
        (cubeOffset.x + newA) * mGridResolution,
        (cubeOffset.y + newA) * mGridResolution,
        (cubeOffset.z + newA) * mGridResolution
    );
    if(evaluateFieldAt(center, field) > mIsoLevel + a * mGridResolution * sqrt(3) / 2.0)
        return 0;

    if(a > 1)
    {
        unsigned totalTriangles = 0;
        for(unsigned i = 0; i < 8; i++)
        {
            #pragma omp task shared(totalTriangles)
            {
                const Vec3_t<float> childOffset(
                    cubeOffset.x + sc_vertexNormPos[i].x * newA,
                    cubeOffset.y + sc_vertexNormPos[i].y * newA,
                    cubeOffset.z + sc_vertexNormPos[i].z * newA
                );
                unsigned childTriangles = evaluateCube(childOffset, newA, field);
                #pragma omp atomic update
                totalTriangles += childTriangles;
            }
        }
        #pragma omp taskwait
        return totalTriangles;
    } else {
        return buildCube(cubeOffset, field);
    }
}

unsigned TreeMeshBuilder::marchCubes(const ParametricScalarField &field)
{
    unsigned totalTriangles;
    #pragma omp parallel
    #pragma omp single nowait
    totalTriangles = evaluateCube(Vec3_t<float>(), mGridSize, field);
    return totalTriangles;
}

float TreeMeshBuilder::evaluateFieldAt(const Vec3_t<float> &pos, const ParametricScalarField &field)
{
    // NOTE: This method is called from "buildCube(...)"!

    // 1. Store pointer to and number of 3D points in the field
    //    (to avoid "data()" and "size()" call in the loop).
    const Vec3_t<float> *pPoints = field.getPoints().data();
    const unsigned count = unsigned(field.getPoints().size());

    float value = std::numeric_limits<float>::max();

    // 2. Find minimum square distance from points "pos" to any point in the
    //    field.
    for(unsigned i = 0; i < count; ++i)
    {
        float distanceSquared  = (pos.x - pPoints[i].x) * (pos.x - pPoints[i].x);
        distanceSquared       += (pos.y - pPoints[i].y) * (pos.y - pPoints[i].y);
        distanceSquared       += (pos.z - pPoints[i].z) * (pos.z - pPoints[i].z);

        // Comparing squares instead of real distance to avoid unnecessary
        // "sqrt"s in the loop.
        value = std::min(value, distanceSquared);
    }

    // 3. Finally take square root of the minimal square distance to get the real distance
    return sqrt(value);
}

void TreeMeshBuilder::emitTriangle(const BaseMeshBuilder::Triangle_t &triangle)
{
    // NOTE: This method is called from "buildCube(...)"!

    // Store generated triangle into vector (array) of generated triangles.
    // The pointer to data in this array is return by "getTrianglesArray(...)" call
    // after "marchCubes(...)" call ends.
    #pragma omp critical(emitTriangle)
    mTriangles.push_back(triangle);
}
