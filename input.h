void InputParameters() {
    // read grid limits
    nx = 128;
    ny = 151;
    il = nx + 2;
    jl = ny + 2;
    ie = nx + 3;
    ib = nx + 4;
    je = ny + 3;
    jb = ny + 4;

    // read flow parameters
    kvis = 1;
    gam = 1.667;
    prandtl = 1.0;
    rho0 = 1;
    rmu0 = 2.2360612501607223e-004;
    p0 = 10.0;
    c0 = std::sqrt(gam * p0 / rho0);

    rm = 0.2;
    u0 = 2.0 * rm * c0;
    v0 = 0.0;

    // reynolds number
    re = 200.0;
    rlen = 1.0 / 14.0;
    vthick0 = 2.0 * re * rmu0 / (u0 * rho0);
    mthick0 = 4.0 * vthick0;

    // read run parameters
    ncyc = 1000;
    cfl = 0.8;

    // read file names
    iread = 0;  // iread = 1 for restart
    fname0 = "mlayer.rst";  // restart file
}

