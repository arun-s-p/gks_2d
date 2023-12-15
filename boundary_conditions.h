#include <iostream>

// Include necessary headers (DIMS, FLO_VAR) if defined in separate files

// Assuming these variables are global or defined in another header file
extern float W[4][IL][JL];
extern int IL, JL, IE, JE;

void BCONDS() {
    // Include DIMS, FLO_VAR if not defined globally
    // Using the provided global variables, adjust as needed

    // LOCAL VARIABLES
    int I, J;
    int JHALO, JINT;

    // SET SLIP CONDITIONS AT THE BOTTOM AND AT THE TOP
    for (I = 2; I <= IL - 1; ++I) {
        JHALO = 2; JINT = 3;
        for (int N = 0; N < 4; ++N) {
            W[N][I][JHALO] = W[N][I][JINT];
        }
        W[2][I][JHALO] = -W[2][I][JINT];

        JHALO = 1; JINT = 4;
        for (int N = 0; N < 4; ++N) {
            W[N][I][JHALO] = W[N][I][JINT];
        }
        W[2][I][JHALO] = -W[2][I][JINT];

        JHALO = JL; JINT = JL - 1;
        for (int N = 0; N < 4; ++N) {
            W[N][I][JHALO] = W[N][I][JINT];
        }
        W[2][I][JHALO] = -W[2][I][JINT];

        JHALO = JE; JINT = JL - 2;
        for (int N = 0; N < 4; ++N) {
            W[N][I][JHALO] = W[N][I][JINT];
        }
        W[2][I][JHALO] = -W[2][I][JINT];
    }

    // SET PERIODIC BOUNDARY CONDITIONS ON LEFT AND RIGHT
    for (J = 3; J <= JL - 1; ++J) {
        for (int N = 0; N < 4; ++N) {
            W[N][2][J] = W[N][IL - 1][J];
            W[N][1][J] = W[N][IL - 2][J];
            W[N][IL][J] = W[N][3][J];
            W[N][IE][J] = W[N][4][J];
        }
    }
}

