
void ApplyBoundaryConditions() {

    // local variables
    int i, j;
    int halo, interior;

    // set slip conditions at the bottom and at the top
    for (i = 2; i <= *il - 1; ++i) {
        halo = 1; interior = 2;
        for (int n = 0; n < 4; ++n) {
            w[i][halo][n] = w[i][interior][n];
        }
        w[i][halo][2] = -w[i][interior][2];

        halo = 0; interior = 3;
        for (int n = 0; n < 4; ++n) {
            w[i][halo][n] = w[i][interior][n];
        }
        w[i][halo][2] = -w[i][interior][2];

        halo = *jl-1; interior = *jl - 2;
        for (int n = 0; n < 4; ++n) {
            w[i][halo][n] = w[i][interior][n];
        }
        w[i][halo][2] = -w[i][interior][2];

        halo = *jl; interior = *jl - 3;
        for (int n = 0; n < 4; ++n) {
            w[i][halo][n] = w[i][interior][n];
        }
        w[i][halo][2] = -w[i][interior][2];
    }

    // set periodic boundary conditions on left and right
    for (j = 2; j < *jl - 1; ++j) {
        for (int n = 0; n < 4; ++n) {
            w[1][j][n] = w[*il - 2][j][n];
            w[0][j][n] = w[*il - 3][j][n];
            w[*il-1][j][n] = w[2][j][n];
            w[*il][j][n] = w[3][j][n];
        }
    }
}

