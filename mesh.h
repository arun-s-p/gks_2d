void CartesianMesh2D() {

    AllocateMeshVariables(ib,jb,2);

    float pi;
    int i, j;

    pi = 4.0 * atan(1.0);
    dx = 20.0 * vthick0 / (float(nx) - 1.0);
    for (i = 0; i < ib; ++i)
        x[i] = (i-2)*dx;

    dy = 20.0 * vthick0 / (float(ny) - 1.0);
    for (j = 0; j < jb; ++j)
        y[j] = (j-2)*dy;

    int ncells =0;
    for (j = 0; j < je; ++j) {
        float deltaY = y[j + 1] - y[j];
        for (i = 0; i < ie; ++i) {
            float deltaX = x[i + 1] - x[i];
            volume[i][j] = deltaX * deltaY;
            ncells++;
        }
    }

    for ( j = 0; j < jb; ++j)
    {
        for ( i = 0; i < ib; ++i)
        {
            coords[i][j][0] = x[i];
            coords[i][j][1] = y[j];
        }
        
    }

    std::cout <<"Cartesian mesh genarated with "<<ncells<<" cells (including ghost cells) \n";
}

