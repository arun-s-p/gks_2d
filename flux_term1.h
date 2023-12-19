// --- FIND NON-EQUILIBRIUM TERM COEFFS (SIMPLIFICATION POSSIBLE..LATER)

std::vector<float> atl(5), atr(5);
float rhsl[5], rhsr[5];
float extau = exp(-dt / tau);
float alpha4 = tau * (1.0 - extau);
float alpha5 = tau * (dt * extau - alpha4);
float alpha1 = dt - alpha4;
float alpha6 = tau * (dt * extau - 2.0 * alpha4);
float alpha7 = tau * alpha4;
float term1l[5], term1r[5], term1[5];

rhsl[0] = acl[0] * ufl[1] + acl[1] * ufl[2] + acl[2] * ufl[1] * vfl[1] +
          acl[3] * ufl[1] * wfl[1] + acl[4] * (ufl[3] + ufl[1] * xvwl2);
rhsl[1] = acl[0] * ufl[2] + acl[1] * ufl[3] + acl[2] * ufl[2] * vfl[1] +
          acl[3] * ufl[2] * wfl[1] + acl[4] * (ufl[4] + ufl[2] * xvwl2);
rhsl[2] = acl[0] * ufl[1] * vfl[1] + acl[1] * ufl[2] * vfl[1] +
          acl[2] * ufl[1] * vfl[2] + acl[3] * ufl[1] * vfl[1] * wfl[1] +
          acl[4] * (ufl[3] * vfl[1] + ufl[1] * vfl[3] +
                    ufl[1] * vfl[1] * wfl[2] + ufl[1] * vfl[1] * xfl2);
rhsl[3] = acl[0] * ufl[1] * wfl[1] + acl[1] * ufl[2] * wfl[1] +
          acl[2] * ufl[1] * vfl[1] * wfl[1] + acl[3] * ufl[1] * wfl[2] +
          acl[4] * (ufl[3] * wfl[1] + ufl[1] * wfl[3] +
                    ufl[1] * wfl[1] * vfl[2] + ufl[1] * wfl[1] * xfl2);
rhsl[4] =
    acl[0] * (ufl[3] + ufl[1] * xvwl2) + acl[1] * (ufl[4] + ufl[2] * xvwl2) +
    acl[2] * (ufl[3] * vfl[1] + ufl[1] * vfl[3] + ufl[1] * vfl[1] * wfl[2] +
              ufl[1] * vfl[1] * xfl2) +
    acl[3] * (ufl[3] * wfl[1] + ufl[1] * wfl[3] + ufl[1] * wfl[1] * vfl[2] +
              ufl[1] * wfl[1] * xfl2) +
    acl[4] * (ufl[5] + ufl[1] * vfl[4] + ufl[1] * wfl[4] + ufl[1] * xfl4 +
              2.0 * ufl[3] * xvwl2 +
              2.0 * ufl[1] * (vfl[2] * xfl2 + wfl[2] * xfl2 + vfl[2] * wfl[2]));
rhsl[4] *= 0.5;

rhsr[0] = acr[0] * ufr[1] + acr[1] * ufr[2] + acr[2] * ufr[1] * vfr[1] +
          acr[3] * ufr[1] * wfr[1] + acr[4] * (ufr[3] + ufr[1] * xvwr2);
rhsr[1] = acr[0] * ufr[2] + acr[1] * ufr[3] + acr[2] * ufr[2] * vfr[1] +
          acr[3] * ufr[2] * wfr[1] + acr[4] * (ufr[4] + ufr[2] * xvwr2);
rhsr[2] = acr[0] * ufr[1] * vfr[1] + acr[1] * ufr[2] * vfr[1] +
          acr[2] * ufr[1] * vfr[2] + acr[3] * ufr[1] * vfr[1] * wfr[1] +
          acr[4] * (ufr[3] * vfr[1] + ufr[1] * vfr[3] +
                    ufr[1] * vfr[1] * wfr[2] + ufr[1] * vfr[1] * xfr2);
rhsr[3] = acr[0] * ufr[1] * wfr[1] + acr[1] * ufr[2] * wfr[1] +
          acr[2] * ufr[1] * vfr[1] * wfr[1] + acr[3] * ufr[1] * wfr[2] +
          acr[4] * (ufr[3] * wfr[1] + ufr[1] * wfr[3] +
                    ufr[1] * wfr[1] * vfr[2] + ufr[1] * wfr[1] * xfr2);
rhsr[4] =
    acr[0] * (ufr[3] + ufr[1] * xvwr2) + acr[1] * (ufr[4] + ufr[2] * xvwr2) +
    acr[2] * (ufr[3] * vfr[1] + ufr[1] * vfr[3] + ufr[1] * vfr[1] * wfr[2] +
              ufr[1] * vfr[1] * xfr2) +
    acr[3] * (ufr[3] * wfr[1] + ufr[1] * wfr[3] + ufr[1] * wfr[1] * vfr[2] +
              ufr[1] * wfr[1] * xfr2) +
    acr[4] * (ufr[5] + ufr[1] * vfr[4] + ufr[1] * wfr[4] + ufr[1] * xfr4 +
              2.0 * ufr[3] * xvwr2 +
              2.0 * ufr[1] * (vfr[2] * xfr2 + wfr[2] * xfr2 + vfr[2] * wfr[2]));
rhsr[5] *= 0.5;

for (int i = 0; i < 5; ++i) {
  rhsl[i] = -rhsl[i];
  rhsr[i] = -rhsr[i];
}

dxe3(uxl, vxl, wxl, lambdal, rhsl, atl, gam);
dxe3(uxr, vxr, wxr, lambdar, rhsr, atr, gam);

term1l[0] = acl[0] * upl[1] + acl[1] * upl[2] + acl[2] * upl[1] * vfl[1] +
            acl[3] * upl[1] * wfl[1] + acl[4] * (upl[3] + upl[1] * xvwl2);
term1l[1] = acl[0] * upl[2] + acl[1] * upl[3] + acl[2] * upl[2] * vfl[1] +
            acl[3] * upl[2] * wfl[1] + acl[4] * (upl[4] + upl[2] * xvwl2);
term1l[2] = acl[0] * upl[1] * vfl[1] + acl[1] * upl[2] * vfl[1] +
            acl[2] * upl[1] * vfl[2] + acl[3] * upl[1] * vfl[1] * wfl[1] +
            acl[4] * (upl[3] * vfl[1] + upl[1] * vfl[3] +
                      upl[1] * vfl[1] * wfl[2] + upl[1] * vfl[1] * xfl2);
term1l[3] = acl[0] * upl[1] * wfl[1] + acl[1] * upl[2] * wfl[1] +
            acl[2] * upl[1] * vfl[1] * wfl[1] + acl[3] * upl[1] * wfl[2] +
            acl[4] * (upl[3] * wfl[1] + upl[1] * wfl[3] +
                      upl[1] * wfl[1] * vfl[2] + upl[1] * wfl[1] * xfl2);
term1l[4] =
    acl[0] * (upl[3] + upl[1] * xvwl2) + acl[1] * (upl[4] + upl[2] * xvwl2) +
    acl[2] * (upl[3] * vfl[1] + upl[1] * vfl[3] + upl[1] * vfl[1] * wfl[2] +
              upl[1] * vfl[1] * xfl2) +
    acl[3] * (upl[3] * wfl[1] + upl[1] * wfl[3] + upl[1] * wfl[1] * vfl[2] +
              upl[1] * wfl[1] * xfl2) +
    acl[4] * (upl[5] + upl[1] * vfl[4] + upl[1] * wfl[4] + upl[1] * xfl4 +
              2.0 * upl[3] * xvwl2 +
              2.0 * upl[1] * (vfl[2] * xfl2 + wfl[2] * xfl2 + vfl[2] * wfl[2]));
term1l[4] *= 0.5;

term1r[0] = acr[0] * umr[1] + acr[1] * umr[2] + acr[2] * umr[1] * vfr[1] +
            acr[3] * umr[1] * wfr[1] + acr[4] * (umr[3] + umr[1] * xvwr2);
term1r[1] = acr[0] * umr[2] + acr[1] * umr[3] + acr[2] * umr[2] * vfr[1] +
            acr[3] * umr[2] * wfr[1] + acr[4] * (umr[4] + umr[2] * xvwr2);
term1r[2] = acr[0] * umr[1] * vfr[1] + acr[1] * umr[2] * vfr[1] +
            acr[2] * umr[1] * vfr[2] + acr[3] * umr[1] * vfr[1] * wfr[1] +
            acr[4] * (umr[3] * vfr[1] + umr[1] * vfr[3] +
                      umr[1] * vfr[1] * wfr[2] + umr[1] * vfr[1] * xfr2);
term1r[3] = acr[0] * umr[1] * wfr[1] + acr[1] * umr[2] * wfr[1] +
            acr[2] * umr[1] * vfr[1] * wfr[1] + acr[3] * umr[1] * wfr[2] +
            acr[4] * (umr[3] * wfr[1] + umr[1] * wfr[3] +
                      umr[1] * wfr[1] * vfr[2] + umr[1] * wfr[1] * xfr2);
term1r[4] =
    acr[0] * (umr[3] + umr[1] * xvwr2) + acr[1] * (umr[4] + umr[2] * xvwr2) +
    acr[2] * (umr[3] * vfr[1] + umr[1] * vfr[3] + umr[1] * vfr[1] * wfr[2] +
              umr[1] * vfr[1] * xfr2) +
    acr[3] * (umr[3] * wfr[1] + umr[1] * wfr[3] + umr[1] * wfr[1] * vfr[2] +
              umr[1] * wfr[1] * xfr2) +
    acr[4] * (umr[5] + umr[1] * vfr[4] + umr[1] * wfr[4] + umr[1] * xfr4 +
              2.0 * umr[3] * xvwr2 +
              2.0 * umr[1] * (vfr[2] * xfr2 + wfr[2] * xfr2 + vfr[2] * wfr[2]));
term1r[4] *= 0.5;

for (int i = 0; i < 5; ++i) {
  term1[i] = alpha6 * (term1l[i] + term1r[i]);
}
