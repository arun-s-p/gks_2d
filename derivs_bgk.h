#include <cmath>
#include <iostream>

// Assuming these variables are global or defined in another header file
extern int IL, JL;
extern double **X, **VOL;
extern double **W, **WX;
extern double **DW, **FW, **VW;
extern double **DTMIN;
extern double GAMMA, PRN, RMU0;

void BGKFLUX(const double *WL, const double *WR, const double *W1P, const double *W2P,
             double DSL, double DSR, double DTMIN, double RMU0, double *FC, double GAMMA, double PRN);

void DERIVS_BGK() {
    const double EPS = 1.0E-7;
    double DWL, DWR, DSL, DSR, DX, DY, XY;
    double FC[4], FD[4], FV[4], FCP[4], FDP[4], FVP[4];
    double WL[4], WR[4], WLP[4], WRP[4], W1P[4], W2P[4];
    double DWT0[4];
    double MU;
    int I, J, N, IP1, JP1;
    double PLW[2], PRW[2], betaL[2], betaR[2], alphaL[2], alphaR[2];
    double omegaL[2], omegaR[2];
    double SMV = 1E-6;

    DW = 0.0;
    FW = 0.0;
    VW = 0.0;
    WX = 0.0;
    DWT0 = 0.0;

    // Find flux contributions in I direction
    for (J = 1; J <= JL - 1; ++J) {
        for (I = 1; I <= IL - 1; ++I) {
            DY = X[2][I][J + 1] - X[2][I][J];
            DSL = 0.5 * (X[1][I + 1][J] - X[1][I][J]);
            DSR = 0.5 * (X[1][I + 2][J] - X[1][I + 1][J]);

            // WENO interpolation
            for (N = 0; N < 4; ++N) {
                // Reconstruct WL and WR for N = 1 to 4 from W - X direction I index
                PLW[0] = 0.5 * (W[N][I][J] + W[N][I + 1][J]);
                PLW[1] = -0.5 * W[N][I - 1][J] + 1.5 * W[N][I][J];

                PRW[0] = 1.5 * W[N][I + 1][J] - 0.5 * W[N][I + 2][J];
                PRW[1] = 0.5 * (W[N][I][J] + W[N][I + 1][J]);

                betaL[0] = pow(W[N][I + 1][J] - W[N][I][J], 2);
                betaL[1] = pow(W[N][I][J] - W[N][I - 1][J], 2);

                betaR[0] = pow(W[N][I + 2][J] - W[N][I + 1][J], 2);
                betaR[1] = pow(W[N][I + 1][J] - W[N][I][J], 2);

                alphaL[0] = 2.0 / (3.0 * pow(SMV + betaL[0], 2));
                alphaL[1] = 1.0 / (3.0 * pow(SMV + betaL[1], 2));

                alphaR[0] = 1.0 / (3.0 * pow(SMV + betaR[0], 2));
                alphaR[1] = 2.0 / (3.0 * pow(SMV + betaR[1], 2));

                omegaL[0] = alphaL[0] / (alphaL[0] + alphaL[1]);
                omegaL[1] = alphaL[1] / (alphaL[0] + alphaL[1]);

                omegaR[0] = alphaR[0] / (alphaR[0] + alphaR[1]);
                omegaR[1] = alphaR[1] / (alphaR[0] + alphaR[1]);

                WL[N] = omegaL[0] * PLW[0] + omegaL[1] * PLW[1];
                WR[N] = omegaR[0] * PRW[0] + omegaR[1] * PRW[1];
            }

            W1P[0] = W[0][I][J];
            W1P[1] = W[2][I][J];
            W1P[2] = -W[1][I][J];
            W1P[3] = W[3][I][J];

            W2P[0] = W[0][I + 1][J];
            W2P[1] = W[2][I + 1][J];
            W2P[2] = -W[1][I + 1][J];
            W2P[3] = W[3][I + 1][J];

            BGKFLUX(WL, WR, W1P, W2P, DSL, DSR, DTMIN, RMU0, FC, GAMMA, PRN);

            FC[0] *= DY;

            for (N = 0; N < 4; ++N) {
                // Accumulate complete convective flux
                DW[N][I][J] += FC[N];
                DW[N][I + 1][J] -= FC[N];
            }
        }
    }

    // Find flux contributions in J direction
    for (J = 1; J <= JL - 1; ++J) {
        for (I = 1; I <= IL - 1; ++I) {
            DX = X[1][I + 1][J] - X[1][I][J];
            DSL = 0.5 * (X[2][I][J + 1] - X[2][I][J]);
            DSR = 0.5 * (X[2][I][J + 2] - X[2][I][J + 1]);

            // WENO interpolation
            for (N = 0; N < 4; ++N) {
                // Reconstruct WL and WR for N = 1 to 4 from W - Y direction J index
                PLW[0] = 0.5 * (W[N][I][J] + W[N][I][J + 1]);
                PLW[1] = -0.5 * W[N][I][J - 1] + 1.5 * W[N][I][J];

                PRW[0] = 1.5 * W[N][I][J + 1] - 0.5 * W[N][I][J + 2];
                PRW[1] = 0.5 * (W[N][I][J] + W[N][I][J + 1]);

                betaL[0] = pow(W[N][I][J + 1] - W[N][I][J], 2);
                betaL[1] = pow(W[N][I][J] - W[N][I][J - 1], 2);

                betaR[0] = pow(W[N][I][J + 2] - W[N][I][J + 1], 2);
                betaR[1] = pow(W[N][I][J + 1] - W[N][I][J], 2);

                alphaL[0] = 2.0 / (3.0 * pow(SMV + betaL[0], 2));
                alphaL[1] = 1.0 / (3.0 * pow(SMV + betaL[1], 2));

                alphaR[0] = 1.0 / (3.0 * pow(SMV + betaR[0], 2));
                alphaR[1] = 2.0 / (3.0 * pow(SMV + betaR[1], 2));

                omegaL[0] = alphaL[0] / (alphaL[0] + alphaL[1]);
                omegaL[1] = alphaL[1] / (alphaL[0] + alphaL[1]);

                omegaR[0] = alphaR[0] / (alphaR[0] + alphaR[1]);
                omegaR[1] = alphaR[1] / (alphaR[0] + alphaR[1]);

                WL[N] = omegaL[0] * PLW[0] + omegaL[1] * PLW[1];
                WR[N] = omegaR[0] * PRW[0] + omegaR[1] * PRW[1];
            }

            if (J == 1 || J == JL - 1) {
                for (N = 0; N < 4; ++N)
                    WL[N] = 0.5 * (W[N][I][J] + W[N][I][J + 1]);

                for (N = 0; N < 4; ++N)
                    WR[N] = WL[N];
            }

            // Rotate the variables
            WLP[0] = WL[0];
            WLP[1] = WL[2];
            WLP[2] = -WL[1];
            WLP[3] = WL[3];

            WRP[0] = WR[0];
            WRP[1] = WR[2];
            WRP[2] = -WRP[1];
            WRP[3] = WR[3];

            W1P[0] = W[0][I][J];
            W1P[1] = W[2][I][J];
            W1P[2] = -W[1][I][J];
            W1P[3] = W[3][I][J];

            W2P[0] = WR[0][I][J + 1];
            W2P[1] = WR[2][I][J + 1];
            W2P[2] = -WRP[1][I][J + 1];
            W2P[3] = WR[3][I][J + 1];

            BGKFLUX(WLP, WRP, W1P, W2P, DSL, DSR, DTMIN, RMU0, FCP, GAMMA, PRN);

            FC[0] = FCP[0]*DX;
            FC[1] = -FCP[2]*DX;
            FC[2] = FCP[1]*DX;
            FC[3] = FCP[3]*DX;

            for (N = 0; N < 4; ++N) {
                // Accumulate complete convective flux
                DW[N][I][J] += FC[N];
                DW[N][I + 1][J] -= FC[N];
            }
        }
    }

}
