float term2l[5], term2r[5], term2[5], dw0l[5], dw0r[5];

term2l[0] =
    rhol * alpha4 * upl[0] +
    alpha6 * (acl[0] * upl[1] + acl[1] * upl[2] + acl[2] * upl[1] * vfl[0] +
              acl[3] * upl[1] * wfl[0] + acl[4] * (upl[3] + upl[1] * xvwl2));
term2l[1] =
    rhol * alpha4 * upl[1] +
    alpha6 * (acl[0] * upl[2] + acl[1] * upl[3] + acl[2] * upl[2] * vfl[0] +
              acl[3] * upl[2] * wfl[0] + acl[4] * (upl[4] + upl[2] * xvwl2));
term2l[2] =
    rhol * alpha4 * upl[0] * vfl[0] +
    alpha6 * (acl[0] * upl[1] * vfl[0] + acl[1] * upl[2] * vfl[0] +
              acl[2] * upl[1] * vfl[1] + acl[3] * upl[1] * vfl[0] * wfl[0] +
              acl[4] * (upl[3] * vfl[0] + upl[1] * (vfl[2] + wfl[2] * xfl2)));
term2l[3] =
    rhol * alpha4 * upl[0] * wfl[0] +
    alpha6 * (acl[0] * upl[1] * wfl[0] + acl[1] * upl[2] * wfl[0] +
              acl[2] * upl[1] * vfl[0] * wfl[0] + acl[3] * upl[1] * wfl[1] +
              acl[4] * (upl[3] * wfl[0] + upl[1] * (wfl[2] + vfl[2] * xfl2)));
term2l[4] =
    rhol * alpha4 * 0.5 * (upl[2] + upl[0] * xvwl2) +
    0.5 * alpha6 *
        (acl[0] * (upl[3] + upl[1] * xvwl2) +
         acl[1] * (upl[4] + upl[2] * xvwl2) +
         acl[2] * (upl[3] * vfl[0] + upl[1] * vfl[2] +
                   upl[1] * vfl[0] * wfl[1] + upl[1] * vfl[0] * xfl2) +
         acl[3] * (upl[3] * wfl[0] + upl[1] * wfl[2] +
                   upl[1] * wfl[0] * vfl[1] + upl[1] * wfl[0] * xfl2) +
         acl[4] * (upl[5] + upl[1] * vfl[3] + upl[1] * wfl[3] + upl[1] * xfl4 +
                   2.0 * upl[3] * xvwl2 +
                   2.0 * upl[1] *
                       (vfl[1] * xfl2 + wfl[1] * xfl2 + vfl[1] * wfl[1])));
term2l[4] *= 0.5;

term2r[0] =
    rhor * alpha4 * umr[0] +
    alpha6 * (acr[0] * umr[1] + acr[1] * umr[2] + acr[2] * umr[1] * vfr[0] +
              acr[3] * umr[1] * wfr[0] + acr[4] * (umr[3] + umr[1] * xvwr2));
term2r[1] =
    rhor * alpha4 * umr[1] +
    alpha6 * (acr[0] * umr[2] + acr[1] * umr[3] + acr[2] * umr[2] * vfr[0] +
              acr[3] * umr[2] * wfr[0] + acr[4] * (umr[4] + umr[2] * xvwr2));
term2r[2] =
    rhor * alpha4 * umr[0] * vfr[0] +
    alpha6 * (acr[0] * umr[1] * vfr[0] + acr[1] * umr[2] * vfr[0] +
              acr[2] * umr[1] * vfr[1] + acr[3] * umr[1] * vfr[0] * wfr[0] +
              acr[4] * (umr[3] * vfr[0] + umr[1] * (vfr[2] + wfr[2] * xfr2)));
term2r[3] =
    rhor * alpha4 * umr[0] * wfr[0] +
    alpha6 * (acr[0] * umr[1] * wfr[0] + acr[1] * umr[2] * wfr[0] +
              acr[2] * umr[1] * vfr[0] * wfr[0] + acr[3] * umr[1] * wfr[1] +
              acr[4] * (umr[3] * wfr[0] + umr[1] * (wfr[2] + vfr[2] * xfr2)));
term2r[4] =
    rhor * alpha4 * 0.5 * (umr[2] + umr[0] * xvwr2) +
    0.5 * alpha6 *
        (acr[0] * (umr[3] + umr[1] * xvwr2) +
         acr[1] * (umr[4] + umr[2] * xvwr2) +
         acr[2] * (umr[3] * vfr[0] + umr[1] * vfr[2] +
                   umr[1] * vfr[0] * wfr[1] + umr[1] * vfr[0] * xfr2) +
         acr[3] * (umr[3] * wfr[0] + umr[1] * wfr[2] +
                   umr[1] * wfr[0] * vfr[1] + umr[1] * wfr[0] * xfr2) +
         acr[4] * (umr[5] + umr[1] * vfr[3] + umr[1] * wfr[3] + umr[1] * xfr4 +
                   2.0 * umr[3] * xvwr2 +
                   2.0 * umr[1] *
                       (vfr[1] * xfr2 + wfr[1] * xfr2 + vfr[1] * wfr[1])));
term2r[4] *= 0.5;

for (int i = 0; i < 5; ++i) {
  term2[i] = term2l[i] + term2r[i];
}
