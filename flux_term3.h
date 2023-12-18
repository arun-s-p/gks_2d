for (int i = 0; i < 5; ++i) {
    dw0l[i] = (ww0[i] - w1[i]) / dsl;
    dw0r[i] = (w2[i] - ww0[i]) / dsr;
}

    std::vector<float> ac0l(5), ac0r(5); 

    dxe3(ux0, vx0, wx0, lambda0, dw0l, ac0l, gam);
    dxe3(ux0, vx0, wx0, lambda0, dw0r, ac0r, gam);

    float laminv0 = 1.0 / lambda0;
    float up0[7], um0[7], uf0[7], vf0[7],wf0[7];

    up0[0] = 0.5 * erfc(-ux0 * sqrt(lambda0));
    up0[1] = ux0 * up0[0] + 0.50 * exp(-lambda0 * ux0 * ux0) / sqrt(lambda0 * flo_pi);

    um0[0] = 0.5 * erfc(ux0 * sqrt(lambda0));
    um0[1] = ux0 * um0[0] - 0.50 * exp(-lambda0 * ux0 * ux0) / sqrt(lambda0 * flo_pi);

    uf0[0] = 1.0;      uf0[1] = ux0 * uf0[0];
    vf0[0] = 1.0;      vf0[1] = vx0 * vf0[0];
    wf0[0] = 1.0;      wf0[1] = wx0 * wf0[0];

    for (int i = 2; i < 7; ++i) {
        up0[i] = ux0 * up0[i-1] + (i - 1.0) * 0.5 * up0[i-2] * laminv0;
        um0[i] = ux0 * um0[i-1] + (i - 1.0) * 0.5 * um0[i-2] * laminv0;
        uf0[i] = ux0 * uf0[i-1] + (i - 1.0) * 0.5 * uf0[i-2] * laminv0;
        vf0[i] = vx0 * vf0[i-1] + (i - 1.0) * 0.5 * vf0[i-2] * laminv0;
        wf0[i] = wx0 * wf0[i-1] + (i - 1.0) * 0.5 * wf0[i-2] * laminv0;
    }

    float xf02  = 0.5 * ck * laminv0;
    float xf04  = 0.25 * ck * (ck + 2.0) * pow(laminv0, 2);
    float xvw02 = xf02 + vf0[2] + wf0[2];
    float term3l[5], term3r[5], term3[5];

    term3l[0] = ac0l[0] * up0[0] + ac0l[1] * up0[1] + ac0l[2] * up0[0] * vf0[0] + ac0l[3] * up0[0] * wf0[0] + ac0l[4] * (up0[2] + up0[0] * xvw02);
    term3l[1] = ac0l[0] * up0[1] + ac0l[1] * up0[2] + ac0l[2] * up0[1] * vf0[0] + ac0l[3] * up0[1] * wf0[0] + ac0l[4] * (up0[3] + up0[1] * xvw02);
    term3l[2] = ac0l[0] * up0[0] * vf0[0] + ac0l[1] * up0[1] * vf0[0] + ac0l[2] * up0[0] * vf0[1] + ac0l[3] * up0[0] * vf0[0] * wf0[0] + ac0l[4] * (up0[2] * vf0[0] + up0[0] * vf0[2] + up0[0] * vf0[0] * wf0[1] + up0[0] * vf0[0] * xf02);
    term3l[3] = ac0l[0] * up0[0] * wf0[0] + ac0l[1] * up0[1] * wf0[0] + ac0l[2] * up0[0] * vf0[0] * wf0[0] + ac0l[3] * up0[0] * wf0[1] + ac0l[4] * (up0[2] * wf0[0] + up0[0] * wf0[2] + up0[0] * wf0[0] * vf0[1] + up0[0] * wf0[0] * xf02);
    term3l[4] = ac0l[0] * (up0[2] + up0[0] * xvw02) + ac0l[1] * (up0[3] + up0[1] * xvw02) + ac0l[2] * (up0[2] * vf0[0] + up0[0] * vf0[2] + up0[0] * vf0[0] * wf0[1] + up0[0] * vf0[0] * xf02) + ac0l[3] * (up0[2] * wf0[0] + up0[0] * wf0[2] + up0[0] * wf0[0] * vf0[1] + up0[0] * wf0[0] * xf02) + ac0l[4] * (up0[4] + up0[0] * vf0[3] + up0[0] * wf0[3] + up0[0] * xf04 + 2.0 * up0[2] * xvw02 + 2.0 * up0[0] * (vf0[1] * xf02 + wf0[1] * xf02 + vf0[1] * wf0[1]));
    term3l[4] *= 0.5;

    term3r[0] = ac0r[0] * um0[0] + ac0r[1] * um0[1] + ac0r[2] * um0[0] * vf0[0] + ac0r[3] * um0[0] * wf0[0] + ac0r[4] * (um0[2] + um0[0] * xvw02);
    term3r[1] = ac0r[0] * um0[1] + ac0r[1] * um0[2] + ac0r[2] * um0[1] * vf0[0] + ac0r[3] * um0[1] * wf0[0] + ac0r[4] * (um0[3] + um0[1] * xvw02);
    term3r[2] = ac0r[0] * um0[0] * vf0[0] + ac0r[1] * um0[1] * vf0[0] + ac0r[2] * um0[0] * vf0[1] + ac0r[3] * um0[0] * vf0[0] * wf0[0] + ac0r[4] * (um0[2] * vf0[0] + um0[0] * vf0[2] + um0[0] * vf0[0] * wf0[1] + um0[0] * vf0[0] * xf02);
    term3r[3] = ac0r[0] * um0[0] * wf0[0] + ac0r[1] * um0[1] * wf0[0] + ac0r[2] * um0[0] * vf0[0] * wf0[0] + ac0r[3] * um0[0] * wf0[1] + ac0r[4] * (um0[2] * wf0[0] + um0[0] * wf0[2] + um0[0] * wf0[0] * vf0[1] + um0[0] * wf0[0] * xf02);
    term3r[4] = ac0r[0] * (um0[2] + um0[0] * xvw02) + ac0r[1] * (um0[3] + um0[1] * xvw02) + ac0r[2] * (um0[2] * vf0[0] + um0[0] * vf0[2] + um0[0] * vf0[0] * wf0[1] + um0[0] * vf0[0] * xf02) + ac0r[3] * (um0[2] * wf0[0] + um0[0] * wf0[2] + um0[0] * wf0[0] * vf0[1] + um0[0] * wf0[0] * xf02) + ac0r[4] * (um0[4] + um0[0] * vf0[3] + um0[0] * wf0[3] + um0[0] * xf04 + 2.0 * um0[2] * xvw02 + 2.0 * um0[0] * (vf0[1] * xf02 + wf0[1] * xf02 + vf0[1] * wf0[1]));
    term3r[4] *= 0.5;

    for (int i = 0; i < 5; ++i) {
        term3[i] = term3l[i] + term3r[i];
    }

