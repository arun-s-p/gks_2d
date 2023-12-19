float term6l[5], term6r[5], term6[5];

term6l[0] = atl[0] * upl[0] + atl[1] * upl[1] + atl[2] * upl[0] * vfl[1] +
            atl[3] * upl[0] * wfl[1] + atl[4] * (upl[2] + upl[0] * xvwl2);
term6l[1] = atl[0] * upl[1] + atl[1] * upl[2] + atl[2] * upl[1] * vfl[1] +
            atl[3] * upl[1] * wfl[1] + atl[4] * (upl[3] + upl[1] * xvwl2);
term6l[2] = atl[0] * upl[0] * vfl[1] + atl[1] * upl[1] * vfl[1] +
            atl[2] * upl[0] * vfl[2] + atl[3] * upl[0] * vfl[1] * wfl[1] +
            atl[4] * (upl[2] * vfl[1] + upl[0] * vfl[3] +
                      upl[0] * vfl[1] * wfl[2] + upl[0] * vfl[1] * xvwl2);
term6l[3] = atl[0] * upl[0] * wfl[1] + atl[1] * upl[1] * wfl[1] +
            atl[2] * upl[0] * vfl[1] * wfl[1] + atl[3] * upl[0] * wfl[2] +
            atl[4] * (upl[2] * wfl[1] + upl[0] * wfl[3] +
                      upl[0] * wfl[1] * vfl[2] + upl[0] * wfl[1] * xvwl2);
term6l[4] =
    atl[0] * (upl[2] + upl[0] * xvwl2) + atl[1] * (upl[3] + upl[1] * xvwl2) +
    atl[2] * (upl[2] * vfl[1] + upl[0] * vfl[3] + upl[0] * vfl[1] * wfl[2] +
              upl[0] * vfl[1] * xvwl2) +
    atl[3] * (upl[2] * wfl[1] + upl[0] * wfl[3] + upl[0] * wfl[1] * vfl[2] +
              upl[0] * wfl[1] * xvwl2) +
    atl[4] *
        (upl[4] + upl[0] * vfl[4] + upl[0] * wfl[4] + upl[0] * xvwl2 * xvwl2 +
         2. * upl[2] * xvwl2 +
         2. * upl[0] * (vfl[2] * xvwl2 + wfl[2] * xvwl2 + vfl[2] * wfl[2]));
term6l[5] = 0.5 * term6l[5];

term6r[0] = atr[0] * umr[0] + atr[1] * umr[1] + atr[2] * umr[0] * vfr[1] +
            atr[3] * umr[0] * wfr[1] + atr[4] * (umr[2] + umr[0] * xvwr2);
term6r[1] = atr[0] * umr[1] + atr[1] * umr[2] + atr[2] * umr[1] * vfr[1] +
            atr[3] * umr[1] * wfr[1] + atr[4] * (umr[3] + umr[1] * xvwr2);
term6r[2] = atr[0] * umr[0] * vfr[1] + atr[1] * umr[1] * vfr[1] +
            atr[2] * umr[0] * vfr[2] + atr[3] * umr[0] * vfr[1] * wfr[1] +
            atr[4] * (umr[2] * vfr[1] + umr[0] * vfr[3] +
                      umr[0] * vfr[1] * wfr[2] + umr[0] * vfr[1] * xvwr2);
term6r[3] = atr[0] * umr[0] * wfr[1] + atr[1] * umr[1] * wfr[1] +
            atr[2] * umr[0] * vfr[1] * wfr[1] + atr[3] * umr[0] * wfr[2] +
            atr[4] * (umr[2] * wfr[1] + umr[0] * wfr[3] +
                      umr[0] * wfr[1] * vfr[2] + umr[0] * wfr[1] * xvwr2);
term6r[4] =
    atr[0] * (umr[2] + umr[0] * xvwr2) + atr[1] * (umr[3] + umr[1] * xvwr2) +
    atr[2] * (umr[2] * vfr[1] + umr[0] * vfr[3] + umr[0] * vfr[1] * wfr[2] +
              umr[0] * vfr[1] * xvwr2) +
    atr[3] * (umr[2] * wfr[1] + umr[0] * wfr[3] + umr[0] * wfr[1] * vfr[2] +
              umr[0] * wfr[1] * xvwr2) +
    atr[4] *
        (umr[4] + umr[0] * vfr[4] + umr[0] * wfr[4] + umr[0] * xvwr2 * xvwr2 +
         2. * umr[2] * xvwr2 +
         2. * umr[0] * (vfr[2] * xvwr2 + wfr[2] * xvwr2 + vfr[2] * wfr[2]));
term6r[5] = 0.5 * term6r[5];

// summing term6l and term6r
for (int i = 0; i < 5; ++i) {
  term6[i] = term6l[i] + term6r[i];
}
