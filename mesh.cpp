#include <iostream>
#include <cmath>

// Include necessary headers (DIMS, MESH_VAR, FLO_PARAM) if defined in separate files

// Assuming these variables are global or defined in another header file
extern int NX, NY, IL, JL, IE, JE, IB, JB;
extern float VTHICK0;
extern float DX, DY;
extern float *XX, *YY, **VOL;

void ALLOC_MESH() {
    // Implementation of ALLOC_MESH function if needed
}

void MESH() {
    // Include DIMS, MESH_VAR, FLO_PARAM if not defined globally
    // Using the provided global variables, adjust as needed
    
    float ETAXL, ETAXR;
    float ETAY, pi;
    int I, J;
    
    ALLOC_MESH();

    pi = 4.0 * atan(1.0);
    DX = 20.0 * VTHICK0 / (float(NX) - 1.0);
    XX[2] = 0.0; // Assuming XX is 1-indexed in Fortran
    for (I = 4; I <= IL; ++I)
        XX[I] = XX[I - 1] + DX;

    DY = 20.0 * VTHICK0 / (float(NY) - 1.0);
    YY[2] = 0.0; // Assuming YY is 1-indexed in Fortran
    for (J = 4; J <= JL; ++J)
        YY[J] = YY[J - 1] + DY;

    XX[1] = 2.0 * XX[2] - XX[3];
    XX[0] = 2.0 * XX[2] - XX[4];
    XX[IE - 1] = 2.0 * XX[IL - 1] - XX[IL - 2];
    XX[IB - 1] = 2.0 * XX[IL - 1] - XX[IL - 3];

    YY[1] = 2.0 * YY[2] - YY[3];
    YY[0] = 2.0 * YY[2] - YY[4];
    YY[JE - 1] = 2.0 * YY[JL - 1] - YY[JL - 2];
    YY[JB - 1] = 2.0 * YY[JL - 1] - YY[JL - 3];

    for (J = 0; J < JE; ++J) {
        for (I = 0; I < IE; ++I) {
            VOL[I][J] = (XX[I + 1] - XX[I]) * (YY[J + 1] - YY[J]);
        }
    }

    // Output results
    for (I = 0; I < IB; ++I)
        std::cout << I + 1 << " " << XX[I] << std::endl;

    for (J = 0; J < JB; ++J)
        std::cout << J + 1 << " " << YY[J] << std::endl;
}

