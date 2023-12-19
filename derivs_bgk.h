#include "bgkflux.h"

void derivs_bgk() {
    const float eps = 1.0e-7;
    float dsl, dsr, dx, dy;
    std::vector<float> fc(4), fcp(4);
    float wl[4], wr[4], wlp[4], wrp[4], w1p[4], w2p[4];
    // std::vector<float>  dwt0;
    float mu;
    int i, j, n, ip1, jp1;
    float plw[2], prw[2], betal[2], betar[2], alphal[2], alphar[2];
    float omegal[2], omegar[2];
    float smv = 1e-6;

    ResetVariables(ie, je, 4, 2);
    // dwt0.resize(4, 0.0);

    // find flux contributions in i direction
    for (j = 1; j < jl - 1; ++j) {
        for (i = 1; i < il - 1; ++i) {
            dy = coords[i][j + 1][1] - coords[i][j][1];
            dsl = 0.5 * (coords[i + 1][j][0] - coords[i][j][0]);
            dsr = 0.5 * (coords[i + 2][j][0] - coords[i + 1][j][0]);

            // weno interpolation
            for (n = 0; n < 4; ++n) {
                // reconstruct wl and wr for n = 1 to 4 from w - x direction i index
                plw[0] = 0.5 * (w[i][j][n] + w[i + 1][j][n]);
                plw[1] = -0.5 * w[i - 1][j][n] + 1.5 * w[i][j][n];

                prw[0] = 1.5 * w[i + 1][j][n] - 0.5 * w[i + 2][j][n];
                prw[1] = 0.5 * (w[i][j][n] + w[i + 1][j][n]);

                betal[0] = pow(w[i + 1][j][n] - w[i][j][n], 2);
                betal[1] = pow(w[i][j][n] - w[i - 1][j][n], 2);

                betar[0] = pow(w[i + 2][j][n] - w[i + 1][j][n], 2);
                betar[1] = pow(w[i + 1][j][n] - w[i][j][n], 2);

                alphal[0] = 2.0 / (3.0 * pow(smv + betal[0], 2));
                alphal[1] = 1.0 / (3.0 * pow(smv + betal[1], 2));

                alphar[0] = 1.0 / (3.0 * pow(smv + betar[0], 2));
                alphar[1] = 2.0 / (3.0 * pow(smv + betar[1], 2));

                omegal[0] = alphal[0] / (alphal[0] + alphal[1]);
                omegal[1] = alphal[1] / (alphal[0] + alphal[1]);

                omegar[0] = alphar[0] / (alphar[0] + alphar[1]);
                omegar[1] = alphar[1] / (alphar[0] + alphar[1]);

                wl[n] = omegal[0] * plw[0] + omegal[1] * plw[1];
                wr[n] = omegar[0] * prw[0] + omegar[1] * prw[1];
            }

            w1p[0] = w[i][j][0];
            w1p[1] = w[i][j][1];
            w1p[2] = w[i][j][2];
            w1p[3] = w[i][j][3];

            w2p[0] = w[i + 1][j][0];
            w2p[1] = w[i + 1][j][1];
            w2p[2] = w[i + 1][j][2];
            w2p[3] = w[i + 1][j][3];

            bgkflux(wl, wr, w1p, w2p, dsl, dsr, dtmin, rmu0, fc, gam, prandtl);

            fc[0] *= dy;

            for (n = 0; n < 4; ++n) {
                // accumulate complete convective flux
                dw[i][j][n] += fc[n];
                dw[i + 1][j][n] -= fc[n];
            }
        }
    }

    // find flux contributions in j direction
    /*for (j = 1; j < jl - 1; ++j) {
        for (i = 1; i < il - 1; ++i) {
            dx = coords[i + 1][j][0] - coords[i][j][0];
            dsl = 0.5 * (coords[i][j + 1][1] - coords[i][j][1]);
            dsr = 0.5 * (coords[i][j + 2][1] - coords[i][j + 1][1]);

            // weno interpolation
            for (n = 0; n < 4; ++n) {
                // reconstruct wl and wr for n = 1 to 4 from w - y direction j index
                plw[0] = 0.5 * (w[i][j][n] + w[i][j + 1][n]);
                plw[1] = -0.5 * w[i][j - 1][n] + 1.5 * w[i][j][n];

                prw[0] = 1.5 * w[i][j + 1][n] - 0.5 * w[i][j + 2][n];
                prw[1] = 0.5 * (w[i][j][n] + w[i][j + 1][n]);

                betal[0] = pow(w[i][j + 1][n] - w[i][j][n], 2);
                betal[1] = pow(w[i][j][n] - w[i][j - 1][n], 2);

                betar[0] = pow(w[i][j + 2][n] - w[i][j + 1][n], 2);
                betar[1] = pow(w[i][j + 1][n] - w[i][j][n], 2);

                alphal[0] = 2.0 / (3.0 * pow(smv + betal[0], 2));
                alphal[1] = 1.0 / (3.0 * pow(smv + betal[1], 2));

                alphar[0] = 1.0 / (3.0 * pow(smv + betar[0], 2));
                alphar[1] = 2.0 / (3.0 * pow(smv + betar[1], 2));

                omegal[0] = alphal[0] / (alphal[0] + alphal[1]);
                omegal[1] = alphal[1] / (alphal[0] + alphal[1]);

                omegar[0] = alphar[0] / (alphar[0] + alphar[1]);
                omegar[1] = alphar[1] / (alphar[0] + alphar[1]);

                wl[n] = omegal[0] * plw[0] + omegal[1] * plw[1];
                wr[n] = omegar[0] * prw[0] + omegar[1] * prw[1];
            }

            if (j == 1 || j == jl - 2) {
                for (n = 0; n < 4; ++n)
                    wl[n] = 0.5 * (w[i][j][n] + w[i][j + 1][n]);
                    wr[n] = wl[n];
            }

            // rotate the variables
            wlp[0] = wl[0];
            wlp[1] = wl[2];
            wlp[2] = -wl[1];
            wlp[3] = wl[3];

            wrp[0] = wr[0];
            wrp[1] = wr[2];
            wrp[2] = -wr[1];
            wrp[3] = wr[3];

            w1p[0] = w[i][j][0];
            w1p[1] = w[i][j][2];
            w1p[2] = -w[i][j][1];
            w1p[3] = w[i][j][3];

            w2p[0] = w[i][j + 1][0];
            w2p[1] = w[i][j + 1][2];
            w2p[2] = -w[i][j + 1][1];
            w2p[3] = w[i][j + 1][3];

            // bgkflux(wlp, wrp, w1p, w2p, dsl, dsr, dtmin, rmu0, fcp, gam, prandtl);

            fc[0] = fcp[0]*dx;
            fc[1] = -fcp[2]*dx;
            fc[2] = fcp[1]*dx;
            fc[3] = fcp[3]*dx;

            for (n = 0; n < 4; ++n) {
                // accumulate complete convective flux
                dw[i][j][n] += fc[n];
                dw[i + 1][j][n] -= fc[n];
            }
        }
    }*/

}
