#include <cmath>

// Assuming these variables are global or defined in another header file
extern float W[4][IL][JL];
extern float W0[4][IL][JL];
extern float DW[4][IL][JL];
extern float DTL[IL][JL];
extern float VOL[IL][JL];
extern float GAMMA;
extern int IL, JL;

void DERIVS_BGK();  // Declaration for the DERIVS_BGK function
void BCONDS();      // Declaration for the BCONDS function

void UPDATE() {
    // LOCAL VARIABLES
    float DU, DV, DP;
    float DT;
    float EL, RTMAX, RTRMS;
    float WC[4];
    float GM1;
    float DX, DY, DXG;
    float MEANOMEGA, MAXOMEGA, OMEGA, RELVOR;
    int I, J, N;
    int IMAX, JMAX;

    GM1 = GAMMA - 1.0;

    // Save the current solution
    for (J = 0; J < JL; ++J) {
        for (I = 0; I < IL; ++I) {
            for (N = 0; N < 4; ++N) {
                W0[N][I][J] = W[N][I][J];
            }
        }
    }

    // Call the DERIVS_BGK function
    DERIVS_BGK();

    RTMAX = 0.0;
    RTRMS = 0.0;
    IMAX = 3;
    JMAX = 3;

    for (J = 2; J < JL - 1; ++J) {
        for (I = 2; I < IL - 1; ++I) {
            EL = W0[1][I][J];
            DT = DTL[I][J] / VOL[I][J];

            // COMPUTE CONSERVATIVE VARIABLES FROM PRIMITIVE VARIABLES
            WC[0] = W0[0][I][J];
            WC[1] = W0[0][I][J] * W0[1][I][J];
            WC[2] = W0[0][I][J] * W0[2][I][J];
            WC[3] = W0[3][I][J] / GM1 + 0.5 * W0[0][I][J] * (pow(W0[1][I][J], 2) + pow(W0[2][I][J], 2));

            // UPDATE CONSERVATIVE VARIABLES
            for (N = 0; N < 4; ++N) {
                WC[N] = WC[N] - DT * DW[N][I][J];
            }

            // RECOMPUTE PRIMITIVE VARIABLES FROM CONSERVATIVE VARIABLES
            W[0][I][J] = WC[0];
            W[1][I][J] = WC[1] / WC[0];
            W[2][I][J] = WC[2] / WC[0];
            W[3][I][J] = GM1 * (WC[3] - 0.5 * (pow(WC[1], 2) + pow(WC[2], 2)) / WC[0]);

            // Track maximum and RMS values
            RTMAX = std::max(RTMAX, std::abs(W[3][I][J]));
            RTRMS += pow(W[3][I][J], 2);
        }
    }

    RTRMS = std::sqrt(RTRMS / ((JL - 3) * (IL - 3)));

    // Call the BCONDS function
    BCONDS();
}

