// Include necessary headers if defined in separate files

// Assuming these variables are global or defined in another header file
extern float GAMMA, RM, RHO0, P0, C0, U0, V0, HH0, EI0, PTOT, CVTTOT;
extern float RE, PRN, PRT, T0, RMU0;
extern float CFL, PWENO;
extern int KVIS;
extern int SELECTED_SCHEME, WENO_FLAG;
extern const int UPS_SCHEME, BGK_SCHEME, CUSP_SCHEME;
extern float gammaw[4];
extern float MTHICK0, VTHICK0;

// Equivalent structure to encapsulate flow parameters
struct FlowParameters {
    float GAMMA, RM, RHO0, P0, C0, U0, V0, HH0, EI0, PTOT, CVTTOT;
    float RE, PRN, PRT, T0, RMU0;
    float CFL, PWENO;
    int KVIS;
    int SELECTED_SCHEME, WENO_FLAG;
    const int UPS_SCHEME, BGK_SCHEME, CUSP_SCHEME;
    float gammaw[4];
    float MTHICK0, VTHICK0;
};

// Function to initialize flow parameters
FlowParameters InitializeFlowParameters() {
    FlowParameters flowParams;
    flowParams.GAMMA = 0.0;
    flowParams.RM = 0.0;
    flowParams.RHO0 = 0.0;
    flowParams.P0 = 0.0;
    // ... (initialize other parameters similarly)

    return flowParams;
}

