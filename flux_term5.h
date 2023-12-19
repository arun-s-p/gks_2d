float term5l[5], term5r[5], term5[5];

// calculation for term5l
term5l[0] = ac0l[0] * up0[1] + ac0l[1] * up0[2] + ac0l[2] * up0[1] * vf0[1] +
            ac0l[3] * up0[1] * wf0[1] + ac0l[4] * (up0[3] + up0[1] * xvw02);
term5l[1] = ac0l[0] * up0[2] + ac0l[1] * up0[3] + ac0l[2] * up0[2] * vf0[1] +
            ac0l[3] * up0[2] * wf0[1] + ac0l[4] * (up0[4] + up0[2] * xvw02);
term5l[2] = ac0l[0] * up0[1] * vf0[1] + ac0l[1] * up0[2] * vf0[1] +
            ac0l[2] * up0[1] * vf0[2] + ac0l[3] * up0[1] * vf0[1] * wf0[1] +
            ac0l[4] * (up0[3] * vf0[1] + up0[1] * vf0[3] +
                       up0[1] * vf0[1] * wf0[2] + up0[1] * vf0[1] * xvw02);
term5l[3] = ac0l[0] * up0[1] * wf0[1] + ac0l[1] * up0[2] * wf0[1] +
            ac0l[2] * up0[1] * vf0[1] * wf0[1] + ac0l[3] * up0[1] * wf0[2] +
            ac0l[4] * (up0[3] * wf0[1] + up0[1] * wf0[3] +
                       up0[1] * wf0[1] * vf0[2] + up0[1] * wf0[1] * xvw02);
term5l[4] =
    ac0l[0] * (up0[3] + up0[1] * xvw02) + ac0l[1] * (up0[4] + up0[2] * xvw02) +
    ac0l[2] * (up0[3] * vf0[1] + up0[1] * vf0[3] + up0[1] * vf0[1] * wf0[2] +
               up0[1] * vf0[1] * xvw02) +
    ac0l[3] * (up0[3] * wf0[1] + up0[1] * wf0[3] + up0[1] * wf0[1] * vf0[2] +
               up0[1] * wf0[1] * xvw02) +
    ac0l[4] *
        (up0[5] + up0[1] * vf0[4] + up0[1] * wf0[4] + up0[1] * xvw02 * xvw02 +
         2. * up0[3] * xvw02 +
         2. * up0[1] * (vf0[2] * xvw02 + wf0[2] * xvw02 + vf0[2] * wf0[2]));
term5l[4] *= 0.5;

// calculation for term5r
term5r[0] = ac0r[0] * up0[1] + ac0r[1] * up0[2] + ac0r[2] * up0[1] * vf0[1] +
            ac0r[3] * up0[1] * wf0[1] + ac0r[4] * (up0[3] + up0[1] * xvw02);
term5r[1] = ac0r[0] * up0[2] + ac0r[1] * up0[3] + ac0r[2] * up0[2] * vf0[1] +
            ac0r[3] * up0[2] * wf0[1] + ac0r[4] * (up0[4] + up0[2] * xvw02);
term5r[2] = ac0r[0] * up0[1] * vf0[1] + ac0r[1] * up0[2] * vf0[1] +
            ac0r[2] * up0[1] * vf0[2] + ac0r[3] * up0[1] * vf0[1] * wf0[1] +
            ac0r[4] * (up0[3] * vf0[1] + up0[1] * vf0[3] +
                       up0[1] * vf0[1] * wf0[2] + up0[1] * vf0[1] * xvw02);
term5r[3] = ac0r[0] * up0[1] * wf0[1] + ac0r[1] * up0[2] * wf0[1] +
            ac0r[2] * up0[1] * vf0[1] * wf0[1] + ac0r[3] * up0[1] * wf0[2] +
            ac0r[4] * (up0[3] * wf0[1] + up0[1] * wf0[3] +
                       up0[1] * wf0[1] * vf0[2] + up0[1] * wf0[1] * xvw02);
term5r[4] =
    ac0r[0] * (up0[3] + up0[1] * xvw02) + ac0r[1] * (up0[4] + up0[2] * xvw02) +
    ac0r[2] * (up0[3] * vf0[1] + up0[1] * vf0[3] + up0[1] * vf0[1] * wf0[2] +
               up0[1] * vf0[1] * xvw02) +
    ac0r[3] * (up0[3] * wf0[1] + up0[1] * wf0[3] + up0[1] * wf0[1] * vf0[2] +
               up0[1] * wf0[1] * xvw02) +
    ac0r[4] *
        (up0[5] + up0[1] * vf0[4] + up0[1] * wf0[4] + up0[1] * xvw02 * xvw02 +
         2. * up0[3] * xvw02 +
         2. * up0[1] * (vf0[2] * xvw02 + wf0[2] * xvw02 + vf0[2] * wf0[2]));
term5r[4] *= 0.5;

// summing term5l and term5r
for (int i = 0; i < 5; ++i) {
  term5[i] = term5l[i] + term5r[i];
}
