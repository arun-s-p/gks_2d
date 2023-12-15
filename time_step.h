#include <cmath>

// Assuming these variables are global or defined in another header file
extern float W[4][IL][JL];
extern float XX[IL], YY[JL];
extern float VOL[IL][JL];
extern float RLV[IL][JL];
extern float DTL[IL][JL];
extern float DTMIN, CFL, GAMMA, PRN;
extern int IL, JL, IE, JE, KVIS;

void TSTEP() {
    // LOCAL VARIABLES
    float DX, DY, DS;
    float A, PRG, PA, RA, CS, QS;
    int I, J, ramp, order;

    for (J = 2; J <= JL - 1; ++J) {
        for (I = 2; I <= IL; ++I) {
            PA = W[3][I - 1][J] + W[3][I][J];
            RA = W[0][I - 1][J] + W[0][I][J];
            DS = YY[J + 1] - YY[J];

            CS = std::sqrt(GAMMA * PA / RA) * DS;
            QS = 0.5 * (W[1][I - 1][J] + W[1][I][J]) * DS;
            A = std::abs(QS) + CS;

            DTL[I - 1][J] += A;
            DTL[I][J] += A;
        }
    }

    // SPECTRAL RADIUS IN J DIRECTION
    for (J = 2; J <= JL; ++J) {
        for (I = 2; I <= IL - 1; ++I) {
            PA = W[3][I][J - 1] + W[3][I][J];
            RA = W[0][I][J - 1] + W[0][I][J];
            DS = XX[I + 1] - XX[I];

            CS = std::sqrt(GAMMA * PA / RA) * DS;
            QS = 0.5 * (W[2][I][J - 1] + W[2][I][J]) * DS;
            A = std::abs(QS) + CS;

            DTL[I][J - 1] += A;
            DTL[I][J] += A;
        }
    }

    // IN CASE OF VISCOUS FLOW TAKE INTO ACCOUNT VISCOSITY
    if (KVIS > 0) {
        PRG = GAMMA / PRN;

        // SPECTRAL RADIUS IN I DIRECTION
        for (J = 2; J <= JL - 1; ++J) {
            for (I = 2; I <= IL; ++I) {
                DTL[I - 1][J] += 2.0 * PRG * RLV[I - 1][J] / (W[0][I - 1][J] * (XX[I] - XX[I - 1]));
                DTL[I][J] += 2.0 * PRG * RLV[I][J] / (W[0][I][J] * (XX[I + 1] - XX[I]));
            }
        }

        // SPECTRAL RADIUS IN J DIRECTION
        for (J = 2; J <= JL; ++J) {
            for (I = 2; I <= IL - 1; ++I) {
                DTL[I][J - 1] += 2.0 * PRG * RLV[I][J - 1] / (W[0][I][J - 1] * (YY[J] - YY[J - 1]));
                DTL[I][J] += 2.0 * PRG * RLV[I][J] / (W[0][I][J] * (YY[J + 1] - YY[J]));
            }
        }
    }

    // DIVIDE BY RESPECTIVE VOLUMES
    for (J = 2; J <= JL - 1; ++J) {
        for (I = 2; I <= IL - 1; ++I) {
            DTL[I][J] = 4.0 * VOL[I][J] / DTL[I][J];
        }
    }

    // SET TIME STEP IN HALO CELLS
    for (J = 2; J <= JL; ++J) {
        DTL[1][J] = DTL[4][J];
        DTL[2][J] = DTL[3][J];
        DTL[IL][J] = DTL[IL - 1][J];
        DTL[IE][J] = DTL[IL - 2][J];
    }

    for (I = 2; I <= IL; ++I) {
        DTL[I][1] = DTL[I][4];
        DTL[I][2] = DTL[I][3];
        DTL[I][JL] = DTL[I][JL - 1];
        DTL[I][JE] = DTL[I][JL - 2];
    }

    DTL = CFL * DTL;

    // FIND MINIMUM TIME STEP
    DTMIN = DTL[2][2];
    for (J = 2; J <= JL - 1; ++J) {
        for (I = 2; I <= IL - 1; ++I) {
            if (DTL[I][J] < DTMIN) DTMIN = DTL[I][J];
        }
    }

    // CONSTANT TIME STEP
    for (J = 0; J <= JE; ++J) {
        for (I = 0; I <= IE; ++I) {
            DTL[I][J] = DTMIN;
        }
    }
}

