    float slopa[5], aterml[5], atermr[5], term4[5];

    aterml[0] = atl[0] * upl[0] + atl[1] * upl[1] + atl[2] * upl[0] * vfl[1] +
                atl[3] * upl[0] * wfl[1] + atl[4] * (upl[2] + upl[0] * xvwl2);
    aterml[1] = atl[0] * upl[1] + atl[1] * upl[2] + atl[2] * upl[1] * vfl[1] +
                atl[3] * upl[1] * wfl[1] + atl[4] * (upl[3] + upl[1] * xvwl2);
    aterml[2] = atl[0] * upl[0] * vfl[1] + atl[1] * upl[1] * vfl[1] +
                atl[2] * upl[0] * vfl[2] + atl[3] * upl[0] * vfl[1] * wfl[1] +
                atl[4] * (upl[2] * vfl[1] + upl[0] * vfl[3] +
                          upl[0] * vfl[1] * wfl[2] + upl[0] * vfl[1] * xvwl2);
    aterml[3] = atl[0] * upl[0] * wfl[1] + atl[1] * upl[1] * wfl[1] +
                atl[2] * upl[0] * vfl[1] * wfl[1] + atl[3] * upl[0] * wfl[2] +
                atl[4] * (upl[2] * wfl[1] + upl[0] * wfl[3] +
                          upl[0] * wfl[1] * vfl[2] + upl[0] * wfl[1] * xvwl2);
    aterml[4] = atl[0] * (upl[2] + upl[0] * xvwl2) +
                atl[1] * (upl[3] + upl[1] * xvwl2) +
                atl[2] * (upl[2] * vfl[1] + upl[0] * vfl[3] +
                          upl[0] * vfl[1] * wfl[2] + upl[0] * vfl[1] * xvwl2) +
                atl[3] * (upl[2] * wfl[1] + upl[0] * wfl[3] +
                          upl[0] * wfl[1] * vfl[2] + upl[0] * wfl[1] * xvwl2) +
                atl[4] * (upl[4] + upl[0] * vfl[4] + upl[0] * wfl[4] +
                          upl[0] * xvwl2 * xvwl2 +
                          2. * upl[2] * xvwl2 +
                          2. * upl[0] * (vfl[2] * xvwl2 + wfl[2] * xvwl2 +
                                       vfl[2] * wfl[2]));
    aterml[4] *= 0.5;

    atermr[0] = atr[0] * umr[0] + atr[1] * umr[1] + atr[2] * umr[0] * vfr[1] +
                atr[3] * umr[0] * wfr[1] + atr[4] * (umr[2] + umr[0] * xvwr2);
    atermr[1] = atr[0] * umr[1] + atr[1] * umr[2] + atr[2] * umr[1] * vfr[1] +
                atr[3] * umr[1] * wfr[1] + atr[4] * (umr[3] + umr[1] * xvwr2);
    atermr[2] = atr[0] * umr[0] * vfr[1] + atr[1] * umr[1] * vfr[1] +
                atr[2] * umr[0] * vfr[2] + atr[3] * umr[0] * vfr[1] * wfr[1] +
                atr[4] * (umr[2] * vfr[1] + umr[0] * vfr[3] +
                          umr[0] * vfr[1] * wfr[2] + umr[0] * vfr[1] * xvwr2);
    atermr[3] = atr[0] * umr[0] * wfr[1] + atr[1] * umr[1] * wfr[1] +
                atr[2] * umr[0] * vfr[1] * wfr[1] + atr[3] * umr[0] * wfr[2] +
                atr[4] * (umr[2] * wfr[1] + umr[0] * wfr[3] +
                          umr[0] * wfr[1] * vfr[2] + umr[0] * wfr[1] * xvwr2);
    atermr[4] = atr[0] * (umr[2] + umr[0] * xvwr2) +
                atr[1] * (umr[3] + umr[1] * xvwr2) +
                atr[2] * (umr[2] * vfr[1] + umr[0] * vfr[3] +
                          umr[0] * vfr[1] * wfr[2] + umr[0] * vfr[1] * xvwr2) +
                atr[3] * (umr[2] * wfr[1] +umr[0] * wfr[3] +
                          umr[0] * wfr[1] * vfr[2] + umr[0] * wfr[1] * xvwr2) +
                atr[4] * (umr[4] + umr[0] * vfr[4] + umr[0] * wfr[4] +
                          umr[0] * xvwr2 * xvwr2 +
                          2. * umr[2] * xvwr2 +
                          2. * umr[0] * (vfr[2] * xvwr2 + wfr[2] * xvwr2 +
                                       vfr[2] * wfr[2]));
    atermr[4] *= 0.5;

    float gam0 = tau * (dt - alpha4);
    float gam2 = -gam0 - alpha5;
    float gam5 = tau * alpha4;

    for (int i = 0; i < 5; ++i) {
        slopa[i] = (term3[i] * gam2 + term1[i] - gam5 * (aterml[i] + atermr[i])) / gam0;
    }

    std::vector<float> aa(5);
    dxe3(ux0, vx0, wx0, lambda0, slopa, aa, gam);

    term4[0] = aa[0] * uf0[0] + aa[1] * uf0[1] + aa[2] * uf0[0] * vf0[1] +
               aa[3] * uf0[0] * wf0[1] + aa[4] * (uf0[2] + uf0[0] * xvw02);
    term4[1] = aa[0] * uf0[1] + aa[1] * uf0[2] + aa[2] * uf0[1] * vf0[1] +
               aa[3] * uf0[1] * wf0[1] + aa[4] * (uf0[3] + uf0[1] * xvw02);
    term4[2] = aa[0] * uf0[0] * vf0[1] + aa[1] * uf0[1] * vf0[1] +
               aa[2] * uf0[0] * vf0[2] + aa[3] * uf0[0] * vf0[1] * wf0[1] +
               aa[4] * (uf0[2] * vf0[1] + uf0[0] * vf0[3] +
                        uf0[0] * vf0[1] * wf0[2] + uf0[0] * vf0[1] * xvw02);
    term4[3] = aa[0] * uf0[0] * wf0[1] + aa[1] * uf0[1] * wf0[1] +
               aa[2] * uf0[0] * vf0[1] * wf0[1] + aa[3] * uf0[0] * wf0[2] +
               aa[4] * (uf0[2] * wf0[1] + uf0[0] * wf0[3] +
                        uf0[0] * wf0[1] * vf0[2] + uf0[0] * wf0[1] * xvw02);
    term4[4] = aa[0] * (uf0[2] + uf0[0] * xvw02) +
               aa[1] * (uf0[3] + uf0[1] * xvw02) +
               aa[2] * (uf0[2] * vf0[1] + uf0[0] * vf0[3] +
                        uf0[0] * vf0[1] * wf0[2] + uf0[0] * vf0[1] * xvw02) +
               aa[3] * (uf0[2] * wf0[1] + uf0[0] * wf0[3] +
                        uf0[0] * wf0[1] * vf0[2] + uf0[0] * wf0[1] * xvw02) +
               aa[4] * (uf0[4] + uf0[0] * vf0[4] + uf0[0] * wf0[4] +
                        uf0[0] * xvw02 * xvw02 +
                        2. * uf0[2] * xvw02 +
                        2. * uf0[0] * (vf0[2] * xvw02 + wf0[2] * xvw02 +
                                       vf0[2] * wf0[2]));
    term4[4] *= 0.5;

