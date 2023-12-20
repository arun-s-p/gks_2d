#include <math.h>

void dxe3(float u, float v, float w, float lambda, float dw[5],
          std::vector<float> &a, float gam) {
  float dd, cc, bb, aa;

  dd = 2.0 * dw[4] -
       (u * u + v * v + w * w + 1.0 / (gam - 1.0) / lambda) * dw[0];
  cc = dw[3] - w * dw[0];
  bb = dw[2] - v * dw[0];
  aa = dw[1] - u * dw[0];

  a[4] = (gam - 1.0) * lambda * lambda *
         (dd - 2.0 * u * aa - 2.0 * v * bb - 2.0 * w * cc);
  a[1] = 2.0 * lambda * (aa - u * a[4] / lambda);
  a[2] = 2.0 * lambda * (bb - v * a[4] / lambda);
  a[3] = 2.0 * lambda * (cc - w * a[4] / lambda);
  a[0] = dw[0] - a[1] * u - a[2] * v - a[3] * w -
         a[4] * (u * u + v * v + w * w + 1.0 / (gam - 1.0) / lambda);
}

void bgkflux(const float wlp[4], const float wrp[4], const float w1p[4],
             const float w2p[4], const float dsl, const float dsr,
             const float dt, const float mu, std::vector<float> &fs,
             const float gam, const float prn) {
  // local variables
  float wl[6], wr[6], w1[6], w2[6];
  float rhol, rhor, uxl, uxr, vxl, vxr, wxl, wxr, el, er, pl, pr;
  float flo_pi = 4.0 * atan(1.0);

  float upl[7], umr[7], ufl[7], ufr[7], vfl[7], vfr[7], wfl[7], wfr[7];
  float ww0[5];

  std::vector<float> acl(5), acr(5);
  float dwl[5], dwr[5];

  time_t tstart, tend;
  time(&tstart);
  // extract variables from input arrays
  wl[0] = wlp[0];
  wr[0] = wrp[0];
  wl[1] = wlp[0] * wlp[1];
  wr[1] = wrp[0] * wrp[1];
  wl[2] = wlp[0] * wlp[2];
  wr[2] = wrp[0] * wrp[2];
  wl[3] = 0.0;
  wr[3] = 0.0;
  wl[4] = (wlp[3] / (gam - 1.0)) +
          0.5 * wlp[0] * (wlp[1] * wlp[1] + wlp[2] * wlp[2]);
  wr[4] = (wrp[3] / (gam - 1.0)) +
          0.5 * wrp[0] * (wrp[1] * wrp[1] + wrp[2] * wrp[2]);

  w1[0] = w1p[0];
  w2[0] = w2p[0];
  w1[1] = w1p[0] * w1p[1];
  w2[1] = w2p[0] * w2p[1];
  w1[2] = w1p[0] * w1p[2];
  w2[2] = w2p[0] * w2p[2];
  w1[3] = 0.0;
  w2[3] = 0.0;
  w1[4] = (w1p[3] / (gam - 1.0)) +
          0.5 * w1p[0] * (w1p[1] * w1p[1] + w1p[2] * w1p[2]);
  w2[4] = (w2p[3] / (gam - 1.0)) +
          0.5 * w2p[0] * (w2p[1] * w2p[1] + w2p[2] * w2p[2]);

  rhol = wl[0];
  rhor = wr[0];
  uxl = wl[1] / wl[0];
  uxr = wr[1] / wr[0];
  vxl = wl[2] / wl[0];
  vxr = wr[2] / wr[0];
  wxl = wl[3] / wl[0];
  wxr = wr[3] / wr[0];
  el = wl[4];
  er = wr[4];

  pl = (gam - 1.0) * (el - 0.5 * wl[0] * (uxl * uxl + vxl * vxl + wxl * wxl));
  pr = (gam - 1.0) * (er - 0.5 * wr[0] * (uxr * uxr + vxr * vxr + wxr * wxr));

  // --------------------------------------------------------------
  // compute lambda = m / (2 k t) = rho / (2 p)
  // this is the damping factor in the maxwellian distribution
  // --------------------------------------------------------------
  float lambdal = 0.5 * wl[0] / pl;
  float laminvl = 1.0 / lambdal;
  float lambdar = 0.5 * wr[0] / pr;
  float laminvr = 1.0 / lambdar;

  // ----------------------------------------------------------------
  // moments
  // f- (-inf->inf) m- (-inf->0) p- (0->inf)
  // ----------------------------------------------------------------

  upl[0] = 0.5 * erfc(-uxl * std::sqrt(lambdal));
  upl[1] =
      uxl * upl[0] + 0.50 * exp(-lambdal * uxl * uxl) / std::sqrt(lambdal * flo_pi);
  umr[0] = 0.5 * erfc(uxr * std::sqrt(lambdar));
  umr[1] =
      uxr * umr[0] - 0.50 * exp(-lambdar * uxr * uxr) / std::sqrt(lambdar * flo_pi);

  ufl[0] = 1.0;
  ufl[1] = uxl * ufl[0];
  ufr[0] = 1.0;
  ufr[1] = uxr * ufr[0];

  vfl[0] = 1.0;
  vfl[1] = vxl * vfl[0];
  vfr[0] = 1.0;
  vfr[1] = vxr * vfr[0];

  wfl[0] = 1.0;
  wfl[1] = wxl * wfl[0];
  wfr[0] = 1.0;
  wfr[1] = wxr * wfr[0];

  for (int i = 2; i <= 6; ++i) {
    upl[i] = uxl * upl[i - 1] + (i - 1.0) * 0.5 * upl[i - 2] * laminvl;
    ufl[i] = uxl * ufl[i - 1] + (i - 1.0) * 0.5 * ufl[i - 2] * laminvl;
    vfl[i] = vxl * vfl[i - 1] + (i - 1.0) * 0.5 * vfl[i - 2] * laminvl;
    wfl[i] = wxl * wfl[i - 1] + (i - 1.0) * 0.5 * wfl[i - 2] * laminvl;

    umr[i] = uxr * umr[i - 1] + (i - 1.0) * 0.5 * umr[i - 2] * laminvr;
    ufr[i] = uxr * ufr[i - 1] + (i - 1.0) * 0.5 * ufr[i - 2] * laminvr;
    vfr[i] = vxr * vfr[i - 1] + (i - 1.0) * 0.5 * vfr[i - 2] * laminvr;
    wfr[i] = wxr * wfr[i - 1] + (i - 1.0) * 0.5 * wfr[i - 2] * laminvr;
  }

  // --- contribution from internal energy (here just rotation)
  float ck = 1.0;
  float xfl2 = 0.5 * ck * laminvl;
  float xfl4 = 0.25 * ck * (ck + 2.0) * pow(laminvl, 2);
  float xvwl2 = xfl2 + vfl[2] + wfl[2];
  float xfr2 = 0.5 * ck * laminvr;
  float xfr4 = 0.25 * ck * (ck + 2.0) * pow(laminvr, 2);
  float xvwr2 = xfr2 + vfr[2] + wfr[2];

  // --- eqn 4.18, collapsing left and right at the cell interface
  ww0[0] = rhol * upl[0] + rhor * umr[0];
  ww0[1] = rhol * upl[1] + rhor * umr[1];
  ww0[2] = rhol * upl[0] * vfl[1] + rhor * umr[0] * vfr[1];
  ww0[3] = rhol * upl[0] * wfl[1] + rhor * umr[0] * wfr[1];
  ww0[4] =
      (rhol * (upl[2] + upl[0] * xvwl2) + rhor * (umr[2] + umr[0] * xvwr2)) *
      0.5;

  float ux0 = ww0[1] / ww0[0];
  float vx0 = ww0[2] / ww0[0];
  float wx0 = ww0[3] / ww0[0];

  float pp0 = (gam - 1.0) *
              (ww0[4] - 0.5 * ww0[0] * (ux0 * ux0 + vx0 * vx0 + wx0 * wx0));
  float lambda0 = 0.5 * std::abs(ww0[0] / pp0);

  // viscous
  float vis = mu / pp0;
  float eps = 2.0;
  float tau = vis + dt * std::abs(pl - pr) / (pl + pr);

  // uncomment these for inviscid
  // double taui = 5.0 * (fabs((rhol / lambdal) - (rhor / lambdar))) /
  // (fabs((rhol / lambdal) + (rhor / lambdar))); tau = 0.05 * dt + dt *
  // fmin(1.0, taui);

  for (int i = 0; i < 5; ++i) {
    dwl[i] = (wl[i] - w1[i]) / dsl;
    dwr[i] = (w2[i] - wr[i]) / dsr;
  }

  dxe3(uxl, vxl, wxl, lambdal, dwl, acl, gam);
  dxe3(uxr, vxr, wxr, lambdar, dwr, acr, gam);

#include "flux_term1.h"

#include "flux_term2.h"

#include "flux_term3.h"

#include "flux_term4.h"

#include "flux_term5.h"

#include "flux_term6.h"

  // Perform the calculations
  float alpha3 = 0.5 * dt * dt - tau * dt + tau * tau * (1.0 - extau);
  float alpha2 = tau * (-dt + alpha4) - alpha5;

  fs[0] = alpha1 * ww0[0] * uf0[0] 
          + alpha2 * term5[0] 
          + alpha3 * term4[0] 
          + term2[0] 
          - alpha7 * term6[0];
  fs[1] = alpha1 * ww0[0] * uf0[1]
          + alpha2 * term5[1] 
          + alpha3 * term4[1] 
          + term2[1] 
          - alpha7 * term6[1];
  fs[2] = alpha1 * ww0[0] * uf0[0] * vf0[0]
          + alpha2 * term5[2] 
          + alpha3 * term4[2] 
          + term2[2] 
          - alpha7 * term6[2];
  fs[3] = 0.5 * alpha1 * ww0[0] * (uf0[2] + uf0[1] * xvw02)
          + alpha2 * term5[4] 
          + alpha3 * term4[4] 
          + term2[4] 
          - alpha7 * term6[4];

  fs[0] /= dt;
  fs[1] /= dt;
  fs[2] /= dt;
  fs[3] /= dt;

  time(&tend);
  t_flux += (tend - tstart);
}
