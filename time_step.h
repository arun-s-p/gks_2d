__global__ void CalculateDT() {
    // local variables
    float ds;
    float a, prg, pa, ra, cs, qs;
    int i, j;

    // time_t tstart, tend;
    // time(&tstart);

    // spectral radius in i direction
    for (j = 2; j < *jl - 1; ++j) {
        for (i = 2; i < *il; ++i) {
            pa = w[i - 1][j][3] + w[i][j][3];
            ra = w[i - 1][j][0] + w[i][j][0];
            ds = *dy; //y[j + 1] - y[j];

            cs = std::sqrt(*gam * pa / ra) * ds;
            qs = 0.5 * (w[i - 1][j][1] + w[i][j][1]) * ds;
            a = std::abs(qs) + cs;

            dtl[i - 1][j] += a;
            dtl[i][j] += a;
        }
    }

/*    // spectral radius in j direction
    for (j = 2; j < jl; ++j) {
        for (i = 2; i < il - 1; ++i) {
            pa = w[i][j - 1][3] + w[i][j][3];
            ra = w[i][j - 1][0] + w[i][j][0];
            ds = dx; //x[i + 1] - x[i];

            cs = std::sqrt(gam * pa / ra) * ds;
            qs = 0.5 * (w[2][i][j - 1] + w[2][i][j]) * ds;
            a = std::abs(qs) + cs;

            dtl[i][j - 1] += a;
            dtl[i][j] += a;
        }
    }
*/
/*    // in case of viscous flow take into account viscosity
    if (kvis > 0) {
        prg = gam / prandtl;

        // spectral radius in i direction
        for (j = 2; j < jl - 1; ++j) {
            for (i = 2; i < il; ++i) {
                dtl[i - 1][j] += 2.0 * prg * rlv[i - 1][j] / (w[i - 1][j][0] * (x[i] - x[i - 1]));
                dtl[i][j] += 2.0 * prg * rlv[i][j] / (w[i][j][0] * (x[i + 1] - x[i]));
            }
        }

        // spectral radius in j direction
        for (j = 2; j < jl; ++j) {
            for (i = 2; i < il - 1; ++i) {
                dtl[i][j - 1] += 2.0 * prg * rlv[i][j - 1] / (w[i][j - 1][0] * (y[j] - y[j - 1]));
                dtl[i][j] += 2.0 * prg * rlv[i][j] / (w[i][j][0] * (y[j + 1] - y[j]));
            }
        }
    }
*/
/*    // divide by respective volumes
    for (j = 2; j < jl - 1; ++j) {
        for (i = 2; i < il - 1; ++i) {
            dtl[i][j] = 4.0 * volume[i][j] / dtl[i][j];
        }
    }

    // find minimum time step
    dtmin = dtl[2][2];
    for (j = 2; j < jl - 1; ++j) {
        for (i = 2; i < il - 1; ++i) {
            if (dtl[i][j] < dtmin) dtmin = dtl[i][j];
        }
    }

    // constant time step
    for (j = 0; j < je; ++j) {
        for (i = 0; i < ie; ++i) {
            dtl[i][j] = dtmin*cfl;
        }
    }
    dtmin *= cfl;
*/
    // time(&tend);
    // t_tstep += (tend - tstart);
}
