#include <cmath>

void DXE3(float U, float V, float W, float LAMBDA, float DW[5], float A[5], float GAM) {
    float DD, CC, BB, AA;

    DD = 2.0 * DW[4] - (U * U + V * V + W * W + 1.0 / (GAM - 1.0) / LAMBDA) * DW[0];
    CC = DW[3] - W * DW[0];
    BB = DW[2] - V * DW[0];
    AA = DW[1] - U * DW[0];

    A[4] = (GAM - 1.0) * LAMBDA * LAMBDA * (DD - 2.0 * U * AA - 2.0 * V * BB - 2.0 * W * CC);
    A[1] = 2.0 * LAMBDA * (AA - U * A[4] / LAMBDA);
    A[2] = 2.0 * LAMBDA * (BB - V * A[4] / LAMBDA);
    A[3] = 2.0 * LAMBDA * (CC - W * A[4] / LAMBDA);
    A[0] = DW[0] - A[1] * U - A[2] * V - A[3] * W - A[4] * (U * U + V * V + W * W + 1.0 / (GAM - 1.0) / LAMBDA);
}

void BGKFLUX(const double WLP[4], const double WRP[4], const double W1P[4], const double W2P[4],
             const double DSL, const double DSR, const double DT, const double MU,
             const double GAMMA, const double PRN, double FS[4]) {
    // Local variables
    double WL[6], WR[6], W1[6], W2[6];
    double RHOL, RHOR, UXL, UXR, VXL, VXR, WXL, WXR, EL, ER, PL, PR;
    double FLO_PI = 4.0 * atan(1.0);

    // Extract variables from input arrays
    WL[0] = WLP[0];    WR[0] = WRP[0];
    WL[1] = WLP[0] * WLP[1];    WR[1] = WRP[0] * WRP[1];
    WL[2] = WLP[0] * WLP[2];    WR[2] = WRP[0] * WRP[2];
    WL[3] = 0.0;        WR[3] = 0.0;
    WL[4] = (WLP[3] / (GAMMA - 1.0)) + 0.5 * WLP[0] * (WLP[1] * WLP[1] + WLP[2] * WLP[2]);
    WR[4] = (WRP[3] / (GAMMA - 1.0)) + 0.5 * WRP[0] * (WRP[1] * WRP[1] + WRP[2] * WRP[2]);

    W1[0] = W1P[0];    W2[0] = W2P[0];
    W1[1] = W1P[0] * W1P[1];    W2[1] = W2P[0] * W2P[1];
    W1[2] = W1P[0] * W1P[2];    W2[2] = W2P[0] * W2P[2];
    W1[3] = 0.0;        W2[3] = 0.0;
    W1[4] = (W1P[3] / (GAMMA - 1.0)) + 0.5 * W1P[0] * (W1P[1] * W1P[1] + W1P[2] * W1P[2]);
    W2[4] = (W2P[3] / (GAMMA - 1.0)) + 0.5 * W2P[0] * (W2P[1] * W2P[1] + W2P[2] * W2P[2]);

    RHOL = WL[0];    RHOR = WR[0];
    UXL = WL[1] / WL[0];    UXR = WR[1] / WR[0];
    VXL = WL[2] / WL[0];    VXR = WR[2] / WR[0];
    WXL = WL[3] / WL[0];    WXR = WR[3] / WR[0];
    EL = WL[4];    ER = WR[4];

    PL = (GAMMA - 1.0) * (EL - 0.5 * WL[0] * (UXL * UXL + VXL * VXL + WXL * WXL));
    PR = (GAMMA - 1.0) * (ER - 0.5 * WR[0] * (UXR * UXR + VXR * VXR + WXR * WXR));

    // Continue from the previous C++ code

// --------------------------------------------------------------
// COMPUTE LAMBDA = m / (2 k T) = RHO / (2 P)
// THIS IS THE DAMPING FACTOR IN THE MAXWELLIAN DISTRIBUTION
// --------------------------------------------------------------
double LAMBDAL = 0.5 * WL[0] / PL;
double LAMINVL = 1.0 / LAMBDAL;
double LAMBDAR = 0.5 * WR[0] / PR;
double LAMINVR = 1.0 / LAMBDAR;

// ----------------------------------------------------------------
// MOMENTS
// F- (-inf->inf) M- (-inf->0) P- (0->inf)
// ----------------------------------------------------------------
double UPL[7], UMR[7], UFL[7], UFR[7], VFL[7], VFR[7], WFL[7], WFR[7];

UPL[0] = 0.5 * erfc(-UXL * sqrt(LAMBDAL));
UPL[1] = UXL * UPL[0] + 0.50 * exp(-LAMBDAL * UXL * UXL) / sqrt(LAMBDAL * FLO_PI);
UMR[0] = 0.5 * erfc(UXR * sqrt(LAMBDAR));
UMR[1] = UXR * UMR[0] - 0.50 * exp(-LAMBDAR * UXR * UXR) / sqrt(LAMBDAR * FLO_PI);

UFL[0] = 1.0;                   UFL[1] = UXL * UFL[0];
UFR[0] = 1.0;                   UFR[1] = UXR * UFR[0];

VFL[0] = 1.0;                   VFL[1] = VXL * VFL[0];
VFR[0] = 1.0;                   VFR[1] = VXR * VFR[0];

WFL[0] = 1.0;                   WFL[1] = WXL * WFL[0];
WFR[0] = 1.0;                   WFR[1] = WXR * WFR[0];

for (int I = 2; I <= 6; ++I) {
    UPL[I] = UXL * UPL[I - 1] + (I - 1.0) * 0.5 * UPL[I - 2] * LAMINVL;
    UFL[I] = UXL * UFL[I - 1] + (I - 1.0) * 0.5 * UFL[I - 2] * LAMINVL;
    VFL[I] = VXL * VFL[I - 1] + (I - 1.0) * 0.5 * VFL[I - 2] * LAMINVL;
    WFL[I] = WXL * WFL[I - 1] + (I - 1.0) * 0.5 * WFL[I - 2] * LAMINVL;

    UMR[I] = UXR * UMR[I - 1] + (I - 1.0) * 0.5 * UMR[I - 2] * LAMINVR;
    UFR[I] = UXR * UFR[I - 1] + (I - 1.0) * 0.5 * UFR[I - 2] * LAMINVR;
    VFR[I] = VXR * VFR[I - 1] + (I - 1.0) * 0.5 * VFR[I - 2] * LAMINVR;
    WFR[I] = WXR * WFR[I - 1] + (I - 1.0) * 0.5 * WFR[I - 2] * LAMINVR;
}

// --- CONTRIBUTION FROM INTERNAL ENERGY (HERE JUST ROTATION)
double XFL2 = 0.5 * CK * LAMINVL;
double XFL4 = 0.25 * CK * (CK + 2.0) * pow(LAMINVL, 2);
double XVWL2 = XFL2 + VFL[2] + WFL[2];

double XFR2 = 0.5 * CK * LAMINVR;
double XFR4 = 0.25 * CK * (CK + 2.0) * pow(LAMINVR, 2);
double XVWR2 = XFR2 + VFR[2] + WFR[2];

// --- Eqn 4.18, collapsing left and right at the cell interface
double WW0[6];
WW0[1] = RHOL * UPL[0] + RHOR * UMR[0];
WW0[2] = RHOL * UPL[1] + RHOR * UMR[1];
WW0[3] = RHOL * UPL[0] * VFL[1] + RHOR * UMR[0] * VFR[1];
WW0[4] = RHOL * UPL[0] * WFL[1] + RHOR * UMR[0] * WFR[1];
WW0[5] = (RHOL * (UPL[2] + UPL[0] * XVWL2) + RHOR * (UMR[2] + UMR[0] * XVWR2)) * 0.5;

double UX0 = WW0[2] / WW0[1];
double VX0 = WW0[3] / WW0[1];
double WX0 = WW0[4] / WW0[1];

double PP0 = GM1 * (WW0[5] - 0.5 * WW0[1] * (UX0 * UX0 + VX0 * VX0 + WX0 * WX0));
double LAMBDA0 = 0.5 * WW0[1] / PP0;

// viscous
double VIS = MU / PP0;
double eps = 2.0;
double TAU = VIS + DT * fabs(PL - PR) / (PL + PR);

// Uncomment these for inviscid
// double TAUI = 5.0 * (fabs((RHOL / LAMBDAL) - (RHOR / LAMBDAR))) / (fabs((RHOL / LAMBDAL) + (RHOR / LAMBDAR)));
// TAU = 0.05 * DT + DT * fmin(1.0, TAUI);

double DWL[6], DWR[6], ACL[6], ACR[6];

for (int I = 1; I <= 5; ++I) {
    DWL[I] = (WL[I] - W1[I]) / DSL;
    DWR[I] = (W2[I] - WR[I]) / DSR;
}

DXE3(UXL, VXL, WXL, LAMBDAL, DWL, ACL, GAMMA);
DXE3(UXR, VXR, WXR, LAMBDAR, DWR, ACR, GAMMA);

// --- FIND NON-EQUILIBRIUM TERM COEFFS (SIMPLIFICATION POSSIBLE..LATER)

double RHSL[6], RHSR[6], ATL[6], ATR[6];

RHSL[1] = ACL[0] * UFL[1] + ACL[1] * UFL[2] + ACL[2] * UFL[1] * VFL[1] +
          ACL[3] * UFL[1] * WFL[1] + ACL[4] * (UFL[3] + UFL[1] * XVWL2);
RHSL[2] = ACL[0] * UFL[2] + ACL[1] * UFL[3] + ACL[2] * UFL[2] * VFL[1] +
          ACL[3] * UFL[2] * WFL[1] + ACL[4] * (UFL[4] + UFL[2] * XVWL2);
RHSL[3] = ACL[0] * UFL[1] * VFL[1] + ACL[1] * UFL[2] * VFL[1] +
          ACL[2] * UFL[1] * VFL[2] + ACL[3] * UFL[1] * VFL[1] * WFL[1] +
          ACL[4] * (UFL[3] * VFL[1] + UFL[1] * VFL[3] + UFL[1] * VFL[1] * WFL[2] +
                    UFL[1] * VFL[1] * XFL2);
RHSL[4] = ACL[0] * UFL[1] * WFL[1] + ACL[1] * UFL[2] * WFL[1] +
          ACL[2] * UFL[1] * VFL[1] * WFL[1] + ACL[3] * UFL[1] * WFL[2] +
          ACL[4] * (UFL[3] * WFL[1] + UFL[1] * WFL[3] + UFL[1] * WFL[1] * VFL[2] +
                    UFL[1] * WFL[1] * XFL2);
RHSL[5] = (ACL[0] * (UFL[3] + UFL[1] * XVWL2) +
           ACL[1] * (UFL[4] + UFL[2] * XVWL2) +
           ACL[2] * (UFL[3] * VFL[1] + UFL[1] * VFL[3] + UFL[1] * VFL[1] * WFL[2] +
                    UFL[1] * VFL[1] * XFL2) +
           ACL[3] * (UFL[3] * WFL[1] + UFL[1] * WFL[3] + UFL[1] * WFL[1] * VFL[2] +
                    UFL[1] * WFL[1] * XFL2) +
           ACL[4] * (UFL[5] + UFL[1] * VFL[4] + UFL[1] * WFL[4] + UFL[1] * XFL4 +
                    2.0 * UFL[3] * XVWL2 +
                    2.0 * UFL[1] * (VFL[2] * XFL2 + WFL[2] * XFL2 + VFL[2] * WFL[2])))
          * 0.5;

RHSL[5] *= 0.5;

RHSR[1] = ACR[0] * UFR[1] + ACR[1] * UFR[2] + ACR[2] * UFR[1] * VFR[1] +
          ACR[3] * UFR[1] * WFR[1] + ACR[4] * (UFR[3] + UFR[1] * XVWR2);
RHSR[2] = ACR[0] * UFR[2] + ACR[1] * UFR[3] + ACR[2] * UFR[2] * VFR[1] +
          ACR[3] * UFR[2] * WFR[1] + ACR[4] * (UFR[4] + UFR[2] * XVWR2);
RHSR[3] = ACR[0] * UFR[1] * VFR[1] + ACR[1] * UFR[2] * VFR[1] +
          ACR[2] * UFR[1] * VFR[2] + ACR[3] * UFR[1] * VFR[1] * WFR[1] +
          ACR[4] * (UFR[3] * VFR[1] + UFR[1] * VFR[3] + UFR[1] * VFR[1] * WFR[2] +
                    UFR[1] * VFR[1] * XFR2);
RHSR[4] = ACR[0] * UFR[1] * WFR[1] + ACR[1] * UFR[2] * WFR[1] +
          ACR[2] * UFR[1] * VFR[1] * WFR[1] + ACR[3] * UFR[1] * WFR[2] +
          ACR[4] * (UFR[3] * WFR[1] + UFR[1] * WFR[3] + UFR[1] * WFR[1] * VFR[2] +
                    UFR[1] * WFR[1] * XFR2);
RHSR[5] = (ACR[0] * (UFR[3] + UFR[1] * XVWR2) +
           ACR[1] * (UFR[4] + UFR[2] * XVWR2) +
           ACR[2] * (UFR[3] * VFR[1] + UFR[1] * VFR[3] + UFR[1] * VFR[1] * WFR[2] +
                    UFR[1] * VFR[1] * XFR2) +
           ACR[3] * (UFR[3] * WFR[1] + UFR[1] * WFR[3] + UFR[1] * WFR[1] * VFR[2] +
                    UFR[1] * WFR[1] * XFR2) +
           ACR[4] * (UFR[5] + UFR[1] * VFR[4] + UFR[1] * WFR[4] + UFR[1] * XFR4 +
                    2.0 * UFR[3] * XVWR2 +
                    2.0 * UFR[1] * (VFR[2] * XFR2 + WFR[2] * XFR2 + VFR[2] * WFR[2])))
          * 0.5;

RHSR[5] *= 0.5;

for (int I = 0; I < 5; ++I) {
    RHSL[I] = -RHSL[I] / RHOL;
    RHSR[I] = -RHSR[I] / RHOR;
}

DXE3(UXL, VXL, WXL, LAMBDAL, RHSL, ATL, GAMMA);
DXE3(UXR, VXR, WXR, LAMBDAR, RHSR, ATR, GAMMA);

double EXTAU = exp(-DT / TAU);

double ALPHA4 = TAU * (1.0 - EXTAU);
double ALPHA5 = TAU * (DT * EXTAU - ALPHA4);
double ALPHA1 = DT - ALPHA4;

double ALPHA6 = TAU * (DT * EXTAU - 2.0 * ALPHA4);
double ALPHA7 = TAU * ALPHA4;

double TERM1L[6], TERM1R[6], TERM1[6];

TERM1L[0] = ACL[0] * UPL[1] + ACL[1] * UPL[2] + ACL[2] * UPL[1] * VFL[1] +
            ACL[3] * UPL[1] * WFL[1] + ACL[4] * (UPL[3] + UPL[1] * XVWL2);
TERM1L[1] = ACL[0] * UPL[2] + ACL[1] * UPL[3] + ACL[2] * UPL[2] * VFL[1] +
            ACL[3] * UPL[2] * WFL[1] + ACL[4] * (UPL[4] + UPL[2] * XVWL2);
TERM1L[2] = ACL[0] * UPL[1] * VFL[1] + ACL[1] * UPL[2] * VFL[1] +
            ACL[2] * UPL[1] * VFL[2] + ACL[3] * UPL[1] * VFL[1] * WFL[1] +
            ACL[4] * (UPL[3] * VFL[1] + UPL[1] * VFL[3] + UPL[1] * VFL[1] * WFL[2] +
                      UPL[1] * VFL[1] * XFL2);
TERM1L[3] = ACL[0] * UPL[1] * WFL[1] + ACL[1] * UPL[2] * WFL[1] +
            ACL[2] * UPL[1] * VFL[1] * WFL[1] + ACL[3] * UPL[1] * WFL[2] +
            ACL[4] * (UPL[3] * WFL[1] + UPL[1] * WFL[3] + UPL[1] * WFL[1] * VFL[2] +
                      UPL[1] * WFL[1] * XFL2);
TERM1L[4] = (ACL[0] * (UPL[3] + UPL[1] * XVWL2) +
             ACL[1] * (UPL[4] + UPL[2] * XVWL2) +
             ACL[2] * (UPL[3] * VFL[1] + UPL[1] * VFL[3] + UPL[1] * VFL[1] * WFL[2] +
                      UPL[1] * VFL[1] * XFL2) +
             ACL[3] * (UPL[3] * WFL[1] + UPL[1] * WFL[3] + UPL[1] * WFL[1] * VFL[2] +
                      UPL[1] * WFL[1] * XFL2) +
             ACL[4] * (UPL[5] + UPL[1] * VFL[4] + UPL[1] * WFL[4] + UPL[1] * XFL4 +
                      2.0 * UPL[3] * XVWL2 +
                      2.0 * UPL[1] * (VFL[2] * XFL2 + WFL[2] * XFL2 + VFL[2] * WFL[2])))
            * 0.5;

TERM1L[4] *= 0.5;

TERM1R[0] = ACR[0] * UMR[1] + ACR[1] * UMR[2] + ACR[2] * UMR[1] * VFR[1] +
            ACR[3] * UMR[1] * WFR[1] + ACR[4] * (UMR[3] + UMR[1] * XVWR2);
TERM1R[1] = ACR[0] * UMR[2] + ACR[1] * UMR[3] + ACR[2] * UMR[2] * VFR[1] +
            ACR[3] * UMR[2] * WFR[1] + ACR[4] * (UMR[4] + UMR[2] * XVWR2);
TERM1R[2] = ACR[0] * UMR[1] * VFR[1] + ACR[1] * UMR[2] * VFR[1] +
            ACR[2] * UMR[1] * VFR[2] + ACR[3] * UMR[1] * VFR[1] * WFR[1] +
            ACR[4] * (UMR[3] * VFR[1] + UMR[1] * VFR[3] + UMR[1] * VFR[1] * WFR[2] +
                      UMR[1] * VFR[1] * XFR2);
TERM1R[3] = ACR[0] * UMR[1] * WFR[1] + ACR[1] * UMR[2] * WFR[1] +
            ACR[2] * UMR[1] * VFR[1] * WFR[1] + ACR[3] * UMR[1] * WFR[2] +
            ACR[4] * (UMR[3] * WFR[1] + UMR[1] * WFR[3] + UMR[1] * WFR[1] * VFR[2] +
                      UMR[1] * WFR[1] * XFR2);
TERM1R[4] = (ACR[0] * (UMR[3] + UMR[1] * XVWR2) +
             ACR[1] * (UMR[4] + UMR[2] * XVWR2) +
             ACR[2] * (UMR[3] * VFR[1] + UMR[1] * VFR[3] + UMR[1] * VFR[1] * WFR[2] +
                      UMR[1] * VFR[1] * XFR2) +
             ACR[3] * (UMR[3] * WFR[1] + UMR[1] * WFR[3] + UMR[1] * WFR[1] * VFR[2] +
                      UMR[1] * WFR[1] * XFR2) +
             ACR[4] * (UMR[5] + UMR[1] * VFR[4] + UMR[1] * WFR[4] + UMR[1] * XFR4 +
                      2.0 * UMR[3] * XVWR2 +
                      2.0 * UMR[1] * (VFR[2] * XFR2 + WFR[2] * XFR2 + VFR[2] * WFR[2])))
            * 0.5;

TERM1R[4] *= 0.5;

for (int I = 0; I < 5; ++I) {
    TERM1[I] = ALPHA6 * (TERM1L[I] + TERM1R[I]);
}

double TERM2L[6], TERM2R[6], TERM2[6], DW0L[6], DW0R[6];

TERM2L[0] = RHOL * ALPHA4 * UPL[0] + ALPHA6 * (ACL[0] * UPL[1] + ACL[1] * UPL[2] +
            ACL[2] * UPL[1] * VFL[0] + ACL[3] * UPL[1] * WFL[0] + ACL[4] * (UPL[3] + UPL[1] * XVWL2));
TERM2L[1] = RHOL * ALPHA4 * UPL[1] + ALPHA6 * (ACL[0] * UPL[2] + ACL[1] * UPL[3] +
            ACL[2] * UPL[2] * VFL[0] + ACL[3] * UPL[2] * WFL[0] + ACL[4] * (UPL[4] + UPL[2] * XVWL2));
TERM2L[2] = RHOL * ALPHA4 * UPL[0] * VFL[0] + ALPHA6 * (ACL[0] * UPL[1] * VFL[0] + ACL[1] * UPL[2] * VFL[0] +
            ACL[2] * UPL[1] * VFL[1] + ACL[3] * UPL[1] * VFL[0] * WFL[0] +
            ACL[4] * (UPL[3] * VFL[0] + UPL[1] * (VFL[2] + WFL[2] * XFL2)));
TERM2L[3] = RHOL * ALPHA4 * UPL[0] * WFL[0] + ALPHA6 * (ACL[0] * UPL[1] * WFL[0] + ACL[1] * UPL[2] * WFL[0] +
            ACL[2] * UPL[1] * VFL[0] * WFL[0] + ACL[3] * UPL[1] * WFL[1] +
            ACL[4] * (UPL[3] * WFL[0] + UPL[1] * (WFL[2] + VFL[2] * XFL2)));
TERM2L[4] = RHOL * ALPHA4 * 0.5 * (UPL[2] + UPL[0] * XVWL2) + 0.5 * ALPHA6 * (
            ACL[0] * (UPL[3] + UPL[1] * XVWL2) + ACL[1] * (UPL[4] + UPL[2] * XVWL2) +
            ACL[2] * (UPL[3] * VFL[0] + UPL[1] * VFL[2] + UPL[1] * VFL[0] * WFL[1] +
                      UPL[1] * VFL[0] * XFL2) +
            ACL[3] * (UPL[3] * WFL[0] + UPL[1] * WFL[2] + UPL[1] * WFL[0] * VFL[1] +
                      UPL[1] * WFL[0] * XFL2) +
            ACL[4] * (UPL[5] + UPL[1] * VFL[3] + UPL[1] * WFL[3] + UPL[1] * XFL3 +
                      2.0 * UPL[3] * XVWL2 +
                      2.0 * UPL[1] * (VFL[1] * XFL2 + WFL[1] * XFL2 + VFL[1] * WFL[1])));
TERM2L[4] *= 0.5;

TERM2R[0] = RHOR * ALPHA4 * UMR[0] + ALPHA6 * (ACR[0] * UMR[1] + ACR[1] * UMR[2] +
            ACR[2] * UMR[1] * VFR[0] + ACR[3] * UMR[1] * WFR[0] + ACR[4] * (UMR[3] + UMR[1] * XVWR2));
TERM2R[1] = RHOR * ALPHA4 * UMR[1] + ALPHA6 * (ACR[0] * UMR[2] + ACR[1] * UMR[3] +
            ACR[2] * UMR[2] * VFR[0] + ACR[3] * UMR[2] * WFR[0] + ACR[4] * (UMR[4] + UMR[2] * XVWR2));
TERM2R[2] = RHOR * ALPHA4 * UMR[0] * VFR[0] + ALPHA6 * (ACR[0] * UMR[1] * VFR[0] + ACR[1] * UMR[2] * VFR[0] +
            ACR[2] * UMR[1] * VFR[1] + ACR[3] * UMR[1] * VFR[0] * WFR[0] +
            ACR[4] * (UMR[3] * VFR[0] + UMR[1] * (VFR[2] + WFR[2] * XFR2)));
TERM2R[3] = RHOR * ALPHA4 * UMR[0] * WFR[0] + ALPHA6 * (ACR[0] * UMR[1] * WFR[0] + ACR[1] * UMR[2] * WFR[0] +
            ACR[2] * UMR[1] * VFR[0] * WFR[0] + ACR[3] * UMR[1] * WFR[1] +
            ACR[4] * (UMR[3] * WFR[0] + UMR[1] * (WFR[2] + VFR[2] * XFR2)));
TERM2R[4] = RHOR * ALPHA4 * 0.5 * (UMR[2] + UMR[0] * XVWR2) + 0.5 * ALPHA6 * (
            ACR[0] * (UMR[3] + UMR[1] * XVWR2) + ACR[1] * (UMR[4] + UMR[2] * XVWR2) +
            ACR[2] * (UMR[3] * VFR[0] + UMR[1] * VFR[2] + UMR[1] * VFR[0] * WFR[1] +
                      UMR[1] * VFR[0] * XFR2) +
            ACR[3] * (UMR[3] * WFR[0] + UMR[1] * WFR[2] + UMR[1] * WFR[0] * VFR[1] +
                      UMR[1] * WFR[0] * XFR2) +
            ACR[4] * (UMR[5] + UMR[1] * VFR[3] + UMR[1] * WFR[3] + UMR[1] * XFR3 +
                      2.0 * UMR[3] * XVWR2 +
                      2.0 * UMR[1] * (VFR[1] * XFR2 + WFR[1] * XFR2 + VFR[1] * WFR[1])));
TERM2R[4] *= 0.5;

for (int I = 0; I < 5; ++I) {
    TERM2[I] = TERM2L[I] + TERM2R[I];
}

for (int I = 0; I < 5; ++I) {
    DW0L[I] = (WW0[I] - W1[I]) / DSL;
    DW0R[I] = (W2[I] - WW0[I]) / DSR;
}
    double LAMINV0 = 1.0 / LAMBDA0;

    UP0[0] = 0.5 * erfc(-UX0 * sqrt(LAMBDA0));
    UP0[1] = UX0 * UP0[0] + 0.50 * exp(-LAMBDA0 * UX0 * UX0) / sqrt(LAMBDA0 * M_PI);

    UM0[0] = 0.5 * erfc(UX0 * sqrt(LAMBDA0));
    UM0[1] = UX0 * UM0[0] - 0.50 * exp(-LAMBDA0 * UX0 * UX0) / sqrt(LAMBDA0 * M_PI);

    UF0[0] = 1.0;      UF0[1] = UX0 * UF0[0];
    VF0[0] = 1.0;      VF0[1] = VX0 * VF0[0];
    WF0[0] = 1.0;      WF0[1] = WX0 * WF0[0];

    for (int I = 2; I < 7; ++I) {
        UP0[I] = UX0 * UP0[I-1] + (I - 1.0) * 0.5 * UP0[I-2] * LAMINV0;
        UM0[I] = UX0 * UM0[I-1] + (I - 1.0) * 0.5 * UM0[I-2] * LAMINV0;
        UF0[I] = UX0 * UF0[I-1] + (I - 1.0) * 0.5 * UF0[I-2] * LAMINV0;
        VF0[I] = VX0 * VF0[I-1] + (I - 1.0) * 0.5 * VF0[I-2] * LAMINV0;
        WF0[I] = WX0 * WF0[I-1] + (I - 1.0) * 0.5 * WF0[I-2] * LAMINV0;
    }

    double XF02  = 0.5 * CK * LAMINV0;
    double XF04  = 0.25 * CK * (CK + 2.0) * pow(LAMINV0, 2);
    double XVW02 = XF02 + VF0[2] + WF0[2];

    TERM3L[0] = AC0L[0] * UP0[0] + AC0L[1] * UP0[1] + AC0L[2] * UP0[0] * VF0[0] + AC0L[3] * UP0[0] * WF0[0] + AC0L[4] * (UP0[2] + UP0[0] * XVW02);
    TERM3L[1] = AC0L[0] * UP0[1] + AC0L[1] * UP0[2] + AC0L[2] * UP0[1] * VF0[0] + AC0L[3] * UP0[1] * WF0[0] + AC0L[4] * (UP0[3] + UP0[1] * XVW02);
    TERM3L[2] = AC0L[0] * UP0[0] * VF0[0] + AC0L[1] * UP0[1] * VF0[0] + AC0L[2] * UP0[0] * VF0[1] + AC0L[3] * UP0[0] * VF0[0] * WF0[0] + AC0L[4] * (UP0[2] * VF0[0] + UP0[0] * VF0[2] + UP0[0] * VF0[0] * WF0[1] + UP0[0] * VF0[0] * XF02);
    TERM3L[3] = AC0L[0] * UP0[0] * WF0[0] + AC0L[1] * UP0[1] * WF0[0] + AC0L[2] * UP0[0] * VF0[0] * WF0[0] + AC0L[3] * UP0[0] * WF0[1] + AC0L[4] * (UP0[2] * WF0[0] + UP0[0] * WF0[2] + UP0[0] * WF0[0] * VF0[1] + UP0[0] * WF0[0] * XF02);
    TERM3L[4] = AC0L[0] * (UP0[2] + UP0[0] * XVW02) + AC0L[1] * (UP0[3] + UP0[1] * XVW02) + AC0L[2] * (UP0[2] * VF0[0] + UP0[0] * VF0[2] + UP0[0] * VF0[0] * WF0[1] + UP0[0] * VF0[0] * XF02) + AC0L[3] * (UP0[2] * WF0[0] + UP0[0] * WF0[2] + UP0[0] * WF0[0] * VF0[1] + UP0[0] * WF0[0] * XF02) + AC0L[4] * (UP0[4] + UP0[0] * VF0[3] + UP0[0] * WF0[3] + UP0[0] * XF04 + 2.0 * UP0[2] * XVW02 + 2.0 * UP0[0] * (VF0[1] * XF02 + WF0[1] * XF02 + VF0[1] * WF0[1]));
    TERM3L[4] *= 0.5;

    TERM3R[0] = AC0R[0] * UM0[0] + AC0R[1] * UM0[1] + AC0R[2] * UM0[0] * VF0[0] + AC0R[3] * UM0[0] * WF0[0] + AC0R[4] * (UM0[2] + UM0[0] * XVW02);
    TERM3R[1] = AC0R[0] * UM0[1] + AC0R[1] * UM0[2] + AC0R[2] * UM0[1] * VF0[0] + AC0R[3] * UM0[1] * WF0[0] + AC0R[4] * (UM0[3] + UM0[1] * XVW02);
    TERM3R[2] = AC0R[0] * UM0[0] * VF0[0] + AC0R[1] * UM0[1] * VF0[0] + AC0R[2] * UM0[0] * VF0[1] + AC0R[3] * UM0[0] * VF0[0] * WF0[0] + AC0R[4] * (UM0[2] * VF0[0] + UM0[0] * VF0[2] + UM0[0] * VF0[0] * WF0[1] + UM0[0] * VF0[0] * XF02);
    TERM3R[3] = AC0R[0] * UM0[0] * WF0[0] + AC0R[1] * UM0[1] * WF0[0] + AC0R[2] * UM0[0] * VF0[0] * WF0[0] + AC0R[3] * UM0[0] * WF0[1] + AC0R[4] * (UM0[2] * WF0[0] + UM0[0] * WF0[2] + UM0[0] * WF0[0] * VF0[1] + UM0[0] * WF0[0] * XF02);
    TERM3R[4] = AC0R[0] * (UM0[2] + UM0[0] * XVW02) + AC0R[1] * (UM0[3] + UM0[1] * XVW02) + AC0R[2] * (UM0[2] * VF0[0] + UM0[0] * VF0[2] + UM0[0] * VF0[0] * WF0[1] + UM0[0] * VF0[0] * XF02) + AC0R[3] * (UM0[2] * WF0[0] + UM0[0] * WF0[2] + UM0[0] * WF0[0] * VF0[1] + UM0[0] * WF0[0] * XF02) + AC0R[4] * (UM0[4] + UM0[0] * VF0[3] + UM0[0] * WF0[3] + UM0[0] * XF04 + 2.0 * UM0[2] * XVW02 + 2.0 * UM0[0] * (VF0[1] * XF02 + WF0[1] * XF02 + VF0[1] * WF0[1]));
    TERM3R[4] *= 0.5;

    for (int I = 0; I < 5; ++I) {
        TERM3[I] = TERM3L[I] + TERM3R[I];
    }

    ATERML[0] = ATL[0] * UPL[0] + ATL[1] * UPL[1] + ATL[2] * UPL[0] * VFL[1] +
                ATL[3] * UPL[0] * WFL[1] + ATL[4] * (UPL[2] + UPL[0] * XVWL2);
    ATERML[1] = ATL[0] * UPL[1] + ATL[1] * UPL[2] + ATL[2] * UPL[1] * VFL[1] +
                ATL[3] * UPL[1] * WFL[1] + ATL[4] * (UPL[3] + UPL[1] * XVWL2);
    ATERML[2] = ATL[0] * UPL[0] * VFL[1] + ATL[1] * UPL[1] * VFL[1] +
                ATL[2] * UPL[0] * VFL[2] + ATL[3] * UPL[0] * VFL[1] * WFL[1] +
                ATL[4] * (UPL[2] * VFL[1] + UPL[0] * VFL[3] +
                          UPL[0] * VFL[1] * WFL[2] + UPL[0] * VFL[1] * XVWL2);
    ATERML[3] = ATL[0] * UPL[0] * WFL[1] + ATL[1] * UPL[1] * WFL[1] +
                ATL[2] * UPL[0] * VFL[1] * WFL[1] + ATL[3] * UPL[0] * WFL[2] +
                ATL[4] * (UPL[2] * WFL[1] + UPL[0] * WFL[3] +
                          UPL[0] * WFL[1] * VFL[2] + UPL[0] * WFL[1] * XVWL2);
    ATERML[4] = ATL[0] * (UPL[2] + UPL[0] * XVWL2) +
                ATL[1] * (UPL[3] + UPL[1] * XVWL2) +
                ATL[2] * (UPL[2] * VFL[1] + UPL[0] * VFL[3] +
                          UPL[0] * VFL[1] * WFL[2] + UPL[0] * VFL[1] * XVWL2) +
                ATL[3] * (UPL[2] * WFL[1] + UPL[0] * WFL[3] +
                          UPL[0] * WFL[1] * VFL[2] + UPL[0] * WFL[1] * XVWL2) +
                ATL[4] * (UPL[4] + UPL[0] * VFL[4] + UPL[0] * WFL[4] +
                          UPL[0] * XVWL2 * XVWL2 +
                          2. * UPL[2] * XVWL2 +
                          2. * UPL[0] * (VFL[2] * XVWL2 + WFL[2] * XVWL2 +
                                       VFL[2] * WFL[2]));
    ATERML[4] *= 0.5;

    ATERMR[0] = ATR[0] * UMR[0] + ATR[1] * UMR[1] + ATR[2] * UMR[0] * VFR[1] +
                ATR[3] * UMR[0] * WFR[1] + ATR[4] * (UMR[2] + UMR[0] * XVWR2);
    ATERMR[1] = ATR[0] * UMR[1] + ATR[1] * UMR[2] + ATR[2] * UMR[1] * VFR[1] +
                ATR[3] * UMR[1] * WFR[1] + ATR[4] * (UMR[3] + UMR[1] * XVWR2);
    ATERMR[2] = ATR[0] * UMR[0] * VFR[1] + ATR[1] * UMR[1] * VFR[1] +
                ATR[2] * UMR[0] * VFR[2] + ATR[3] * UMR[0] * VFR[1] * WFR[1] +
                ATR[4] * (UMR[2] * VFR[1] + UMR[0] * VFR[3] +
                          UMR[0] * VFR[1] * WFR[2] + UMR[0] * VFR[1] * XVWR2);
    ATERMR[3] = ATR[0] * UMR[0] * WFR[1] + ATR[1] * UMR[1] * WFR[1] +
                ATR[2] * UMR[0] * VFR[1] * WFR[1] + ATR[3] * UMR[0] * WFR[2] +
                ATR[4] * (UMR[2] * WFR[1] + UMR[0] * WFR[3] +
                          UMR[0] * WFR[1] * VFR[2] + UMR[0] * WFR[1] * XVWR2);
    ATERMR[4] = ATR[0] * (UMR[2] + UMR[0] * XVWR2) +
                ATR[1] * (UMR[3] + UMR[1] * XVWR2) +
                ATR[2] * (UMR[2] * VFR[1] + UMR[0] * VFR[3] +
                          UMR[0] * VFR[1] * WFR[2] + UMR[0] * VFR[1] * XVWR2) +
                ATR[3] * (UMR[2] * WFR[1] +UMR[0] * WFR[3] +
                          UMR[0] * WFR[1] * VFR[2] + UMR[0] * WFR[1] * XVWR2) +
                ATR[4] * (UMR[4] + UMR[0] * VFR[4] + UMR[0] * WFR[4] +
                          UMR[0] * XVWR2 * XVWR2 +
                          2. * UMR[2] * XVWR2 +
                          2. * UMR[0] * (VFR[2] * XVWR2 + WFR[2] * XVWR2 +
                                       VFR[2] * WFR[2]));
    ATERMR[4] *= 0.5;

    double GAM0 = TAU * (DT - ALPHA4);
    double GAM2 = -GAM0 - ALPHA5;
    double GAM5 = TAU * ALPHA4;

    for (int I = 0; I < 5; ++I) {
        SLOPA[I] = (TERM3[I] * GAM2 + TERM1[I] - GAM5 * (ATERML[I] + ATERMR[I])) / GAM0;
    }

    // Assuming DXE3 function is available
    DXE3(UX0, VX0, WX0, LAMBDA0, SLOPA, AA, GAMMA);

    TERM4[0] = AA[0] * UF0[0] + AA[1] * UF0[1] + AA[2] * UF0[0] * VF0[1] +
               AA[3] * UF0[0] * WF0[1] + AA[4] * (UF0[2] + UF0[0] * XVW02);
    TERM4[1] = AA[0] * UF0[1] + AA[1] * UF0[2] + AA[2] * UF0[1] * VF0[1] +
               AA[3] * UF0[1] * WF0[1] + AA[4] * (UF0[3] + UF0[1] * XVW02);
    TERM4[2] = AA[0] * UF0[0] * VF0[1] + AA[1] * UF0[1] * VF0[1] +
               AA[2] * UF0[0] * VF0[2] + AA[3] * UF0[0] * VF0[1] * WF0[1] +
               AA[4] * (UF0[2] * VF0[1] + UF0[0] * VF0[3] +
                        UF0[0] * VF0[1] * WF0[2] + UF0[0] * VF0[1] * XVW02);
    TERM4[3] = AA[0] * UF0[0] * WF0[1] + AA[1] * UF0[1] * WF0[1] +
               AA[2] * UF0[0] * VF0[1] * WF0[1] + AA[3] * UF0[0] * WF0[2] +
               AA[4] * (UF0[2] * WF0[1] + UF0[0] * WF0[3] +
                        UF0[0] * WF0[1] * VF0[2] + UF0[0] * WF0[1] * XVW02);
    TERM4[4] = AA[0] * (UF0[2] + UF0[0] * XVW02) +
               AA[1] * (UF0[3] + UF0[1] * XVW02) +
               AA[2] * (UF0[2] * VF0[1] + UF0[0] * VF0[3] +
                        UF0[0] * VF0[1] * WF0[2] + UF0[0] * VF0[1] * XVW02) +
               AA[3] * (UF0[2] * WF0[1] + UF0[0] * WF0[3] +
                        UF0[0] * WF0[1] * VF0[2] + UF0[0] * WF0[1] * XVW02) +
               AA[4] * (UF0[4] + UF0[0] * VF0[4] + UF0[0] * WF0[4] +
                        UF0[0] * XVW02 * XVW02 +
                        2. * UF0[2] * XVW02 +
                        2. * UF0[0] * (VF0[2] * XVW02 + WF0[2] * XVW02 +
                                       VF0[2] * WF0[2]));
    TERM4[4] *= 0.5;

    // Calculation for TERM5L
    TERM5L[0] = AC0L[0] * UP0[1] + AC0L[1] * UP0[2] + AC0L[2] * UP0[1] * VF0[1] +
                AC0L[3] * UP0[1] * WF0[1] + AC0L[4] * (UP0[3] + UP0[1] * XVW02);
    TERM5L[1] = AC0L[0] * UP0[2] + AC0L[1] * UP0[3] + AC0L[2] * UP0[2] * VF0[1] +
                AC0L[3] * UP0[2] * WF0[1] + AC0L[4] * (UP0[4] + UP0[2] * XVW02);
    TERM5L[2] = AC0L[0] * UP0[1] * VF0[1] + AC0L[1] * UP0[2] * VF0[1] +
                AC0L[2] * UP0[1] * VF0[2] + AC0L[3] * UP0[1] * VF0[1] * WF0[1] +
                AC0L[4] * (UP0[3] * VF0[1] + UP0[1] * VF0[3] +
                           UP0[1] * VF0[1] * WF0[2] + UP0[1] * VF0[1] * XVW02);
    TERM5L[3] = AC0L[0] * UP0[1] * WF0[1] + AC0L[1] * UP0[2] * WF0[1] +
                AC0L[2] * UP0[1] * VF0[1] * WF0[1] + AC0L[3] * UP0[1] * WF0[2] +
                AC0L[4] * (UP0[3] * WF0[1] + UP0[1] * WF0[3] +
                           UP0[1] * WF0[1] * VF0[2] + UP0[1] * WF0[1] * XVW02);
    TERM5L[4] = AC0L[0] * (UP0[3] + UP0[1] * XVW02) +
                AC0L[1] * (UP0[4] + UP0[2] * XVW02) +
                AC0L[2] * (UP0[3] * VF0[1] + UP0[1] * VF0[3] +
                           UP0[1] * VF0[1] * WF0[2] + UP0[1] * VF0[1] * XVW02) +
                AC0L[3] * (UP0[3] * WF0[1] + UP0[1] * WF0[3] +
                           UP0[1] * WF0[1] * VF0[2] + UP0[1] * WF0[1] * XVW02) +
                AC0L[4] * (UP0[5] + UP0[1] * VF0[4] + UP0[1] * WF0[4] +
                           UP0[1] * XVW02 * XVW02 +
                           2. * UP0[3] * XVW02 +
                           2. * UP0[1] * (VF0[2] * XVW02 + WF0[2] * XVW02 +
                                          VF0[2] * WF0[2]));
    TERM5L[4] *= 0.5;

    // Calculation for TERM5R
    TERM5R[0] = AC0R[0] * UP0[1] + AC0R[1] * UP0[2] + AC0R[2] * UP0[1] * VF0[1] +
                AC0R[3] * UP0[1] * WF0[1] + AC0R[4] * (UP0[3] + UP0[1] * XVW02);
    TERM5R[1] = AC0R[0] * UP0[2] + AC0R[1] * UP0[3] + AC0R[2] * UP0[2] * VF0[1] +
                AC0R[3] * UP0[2] * WF0[1] + AC0R[4] * (UP0[4] + UP0[2] * XVW02);
    TERM5R[2] = AC0R[0] * UP0[1] * VF0[1] + AC0R[1] * UP0[2] * VF0[1] +
                AC0R[2] * UP0[1] * VF0[2] + AC0R[3] * UP0[1] * VF0[1] * WF0[1] +
                AC0R[4] * (UP0[3] * VF0[1] + UP0[1] * VF0[3] +
                           UP0[1] * VF0[1] * WF0[2] + UP0[1] * VF0[1] * XVW02);
    TERM5R[3] = AC0R[0] * UP0[1] * WF0[1] + AC0R[1] * UP0[2] * WF0[1] +
                AC0R[2] * UP0[1] * VF0[1] * WF0[1] + AC0R[3] * UP0[1] * WF0[2] +
                AC0R[4] * (UP0[3] * WF0[1] + UP0[1] * WF0[3] +
                           UP0[1] * WF0[1] * VF0[2] + UP0[1] * WF0[1] * XVW02);
    TERM5R[4] = AC0R[0] * (UP0[3] + UP0[1] * XVW02) +
                AC0R[1] * (UP0[4] + UP0[2] * XVW02) +
                AC0R[2] * (UP0[3] * VF0[1] + UP0[1] * VF0[3] +
                           UP0[1] * VF0[1] * WF0[2] + UP0[1] * VF0[1] * XVW02) +
                AC0R[3] * (UP0[3] * WF0[1] + UP0[1] * WF0[3] +
                           UP0[1] * WF0[1] * VF0[2] + UP0[1] * WF0[1] * XVW02) +
                AC0R[4] * (UP0[5] + UP0[1] * VF0[4] + UP0[1] * WF0[4] +
                           UP0[1] * XVW02 * XVW02 +
                           2. * UP0[3] * XVW02 +
                           2. * UP0[1] * (VF0[2] * XVW02 + WF0[2] * XVW02 +
                                          VF0[2] * WF0[2]));
    TERM5R[4] *= 0.5;

    // Summing TERM5L and TERM5R
    for (int i = 0; i < 5; ++i) {
        TERM5[i] = TERM5L[i] + TERM5R[i];
    }

    TERM6L[0] = ATL[0] * UPL[0] + ATL[1] * UPL[1] + ATL[2] * UPL[0] * VFL[1] +
                ATL[3] * UPL[0] * WFL[1] + ATL[4] * (UPL[2] + UPL[0] * XVWL2);
    TERM6L[1] = ATL[0] * UPL[1] + ATL[1] * UPL[2] + ATL[2] * UPL[1] * VFL[1] +
                ATL[3] * UPL[1] * WFL[1] + ATL[4] * (UPL[3] + UPL[1] * XVWL2);
    TERM6L[2] = ATL[0] * UPL[0] * VFL[1] + ATL[1] * UPL[1] * VFL[1] +
                ATL[2] * UPL[0] * VFL[2] + ATL[3] * UPL[0] * VFL[1] * WFL[1] +
                ATL[4] * (UPL[2] * VFL[1] + UPL[0] * VFL[3] +
                          UPL[0] * VFL[1] * WFL[2] +
                          UPL[0] * VFL[1] * XVWL2);
    TERM6L[3] = ATL[0] * UPL[0] * WFL[1] + ATL[1] * UPL[1] * WFL[1] +
                ATL[2] * UPL[0] * VFL[1] * WFL[1] + ATL[3] * UPL[0] * WFL[2] +
                ATL[4] * (UPL[2] * WFL[1] + UPL[0] * WFL[3] +
                          UPL[0] * WFL[1] * VFL[2] +
                          UPL[0] * WFL[1] * XVWL2);
    TERM6L[4] = ATL[0] * (UPL[2] + UPL[0] * XVWL2) +
                ATL[1] * (UPL[3] + UPL[1] * XVWL2) +
                ATL[2] * (UPL[2] * VFL[1] + UPL[0] * VFL[3] +
                          UPL[0] * VFL[1] * WFL[2] +
                          UPL[0] * VFL[1] * XVWL2) +
                ATL[3] * (UPL[2] * WFL[1] + UPL[0] * WFL[3] +
                          UPL[0] * WFL[1] * VFL[2] +
                          UPL[0] * WFL[1] * XVWL2) +
                ATL[4] * (UPL[4] + UPL[0] * VFL[4] + UPL[0] * WFL[4] +
                          UPL[0] * XVWL2 * XVWL2 +
                          2. * UPL[2] * XVWL2 +
                          2. * UPL[0] * (VFL[2] * XVWL2 + WFL[2] * XVWL2 +
                                        VFL[2] * WFL[2]));
    TERM6L[5] = 0.5 * TERM6L[5];

    TERM6R[0] = ATR[0] * UMR[0] + ATR[1] * UMR[1] + ATR[2] * UMR[0] * VFR[1] +
                ATR[3] * UMR[0] * WFR[1] + ATR[4] * (UMR[2] + UMR[0] * XVWR2);
    TERM6R[1] = ATR[0] * UMR[1] + ATR[1] * UMR[2] + ATR[2] * UMR[1] * VFR[1] +
                ATR[3] * UMR[1] * WFR[1] + ATR[4] * (UMR[3] + UMR[1] * XVWR2);
    TERM6R[2] = ATR[0] * UMR[0] * VFR[1] + ATR[1] * UMR[1] * VFR[1] +
                ATR[2] * UMR[0] * VFR[2] + ATR[3] * UMR[0] * VFR[1] * WFR[1] +
                ATR[4] * (UMR[2] * VFR[1] + UMR[0] * VFR[3] +
                          UMR[0] * VFR[1] * WFR[2] +
                          UMR[0] * VFR[1] * XVWR2);
    TERM6R[3] = ATR[0] * UMR[0] * WFR[1] + ATR[1] * UMR[1] * WFR[1] +
                ATR[2] * UMR[0] * VFR[1] * WFR[1] + ATR[3] * UMR[0] * WFR[2] +
                ATR[4] * (UMR[2] * WFR[1] + UMR[0] * WFR[3] +
                          UMR[0] * WFR[1] * VFR[2] +
                          UMR[0] * WFR[1] * XVWR2);
    TERM6R[4] = ATR[0] * (UMR[2] + UMR[0] * XVWR2) +
                ATR[1] * (UMR[3] + UMR[1] * XVWR2) +
                ATR[2] * (UMR[2] * VFR[1] + UMR[0] * VFR[3] +
                          UMR[0] * VFR[1] * WFR[2] +
                          UMR[0] * VFR[1] * XVWR2) +
                ATR[3] * (UMR[2] * WFR[1] + UMR[0] * WFR[3] +
                          UMR[0] * WFR[1] * VFR[2] +
                          UMR[0] * WFR[1] * XVWR2) +
                ATR[4] * (UMR[4] + UMR[0] * VFR[4] + UMR[0] * WFR[4] +
                          UMR[0] * XVWR2 * XVWR2 +
                          2. * UMR[2] * XVWR2 +
                          2. * UMR[0] * (VFR[2] * XVWR2 + WFR[2] * XVWR2 +
                                        VFR[2] * WFR[2]));
    TERM6R[5] = 0.5 * TERM6R[5];

    // Summing TERM6L and TERM6R
    for (int i = 0; i < 5; ++i) {
        TERM6[i] = TERM6L[i] + TERM6R[i];
    }

    // Perform the calculations
    ALPHA3 = 0.5 * DT * DT - TAU * DT + TAU * TAU * (1.0 - EXTAU);
    ALPHA2 = TAU * (-DT + ALPHA4) - ALPHA5;

    FS[0] = ALPHA1 * WW0[0] * UF0[0] + ALPHA2 * TERM5[0] + ALPHA3 * TERM4[0] +
            TERM2[0] - ALPHA7 * TERM6[0];
    FS[1] = ALPHA1 * WW0[0] * UF0[1] + ALPHA2 * TERM5[1] + ALPHA3 * TERM4[1] +
            TERM2[1] - ALPHA7 * TERM6[1];
    FS[2] = ALPHA1 * WW0[0] * UF0[0] * VF0[0] + ALPHA2 * TERM5[2] + ALPHA3 * TERM4[2] +
            TERM2[2] - ALPHA7 * TERM6[2];
    FS[3] = 0.5 * ALPHA1 * WW0[0] * (UF0[2] + UF0[0] * XVW02) + ALPHA2 * TERM5[4] +
            ALPHA3 * TERM4[4] + TERM2[4] - ALPHA7 * TERM6[4];

    FS[0] /= DT
    FS[1] /= DT
    FS[2] /= DT
    FS[3] /= DT

}

