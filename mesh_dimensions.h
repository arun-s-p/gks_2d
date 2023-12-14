// mesh_dimensions.h

#ifndef MESH_DIMENSIONS_H
#define MESH_DIMENSIONS_H

struct MeshDimensions {
    int NX, NY;    // Mesh size
    int IL, JL;    // Lower bounds in X and Y directions
    int IE, JE;    // Upper bounds in X and Y directions
    int IB, JB;    // Bounds for ghost cells in X and Y directions
};

#endif // MESH_DIMENSIONS_H
