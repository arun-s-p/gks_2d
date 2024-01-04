#include "derivs_bgk.h"

void UpdateField() {
    // local variables
    float du, dv, dp;
    float dt;
    float wc[4];
    float gm1;
    int i, j, n;
    int imax, jmax;

    time_t tstart, tend;
    time(&tstart);
    gm1 = *gam - 1.0;

    // save the current solution
    for (n = 0; n < 4; ++n) {
        for (j = 0; j < *jl; ++j) {
            for (i = 0; i < *il; ++i) {
                w0[i][j][n] = w[i][j][n];
            }
        }
    }
   

    // call the derivs_bgk function
    derivs_bgk();

    for (j = 2; j < *jl - 1; ++j) {
        for (i = 2; i < *il - 1; ++i) {
            dt = *dtmin / volume[i][j];

            // compute conservative variables from primitive variables
            wc[0] = w0[i][j][0];
            wc[1] = w0[i][j][0] * w0[i][j][1];
            wc[2] = w0[i][j][0] * w0[i][j][2];
            wc[3] = w0[i][j][3] / gm1 + 0.5 * w0[i][j][0] * (pow(w0[i][j][1], 2) + pow(w0[i][j][2], 2));

            // update conservative variables
            for (n = 0; n < 4; ++n) {
                wc[n] = wc[n] - dt * dw[i][j][n];
            }

            // recompute primitive variables from conservative variables
            w[i][j][0] = wc[0];
            w[i][j][1] = wc[1] / wc[0];
            w[i][j][2] = wc[2] / wc[0];
            w[i][j][3] = gm1 * (wc[3] - 0.5 * (pow(wc[1], 2) + pow(wc[2], 2)) / wc[0]);

        }
    }

    time(&tend);
    t_update += (tend - tstart);
    // call the bconds function
    ApplyBoundaryConditions();
}

