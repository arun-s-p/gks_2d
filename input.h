#include <iostream>
#include <cmath>
#include <string>

// Include necessary headers (DIMS, FLO_VAR, MESH_VAR, FLO_PARAM, TIME_VAR, IODEFS) if defined in separate files

// Assuming these variables are global or defined in another header file
extern int NX, NY, IL, JL, IE, IB, JE, JB;
extern int SELECTED_SCHEME;
extern int KVIS;
extern float GAMMA, PRN, RHO0, RMU0, P0, C0, RM, U0, V0, RE, RLEN, VTHICK0, MTHICK0;
extern int NCYC;
extern float CFL;
extern int IREAD;
extern std::string FNAME0, FNAME1, FNAME2;

void INPUT() {
    // Include DIMS, FLO_VAR, MESH_VAR, FLO_PARAM, TIME_VAR, IODEFS if not defined globally
    // Using the provided global variables, adjust as needed

    // READ GRID LIMITS
    NX = 128;
    NY = 151;
    IL = NX + 2;
    JL = NY + 2;
    IE = NX + 3;
    IB = NX + 4;
    JE = NY + 3;
    JB = NY + 4;

    SELECTED_SCHEME = BGK_SCHEME;

    // READ FLOW PARAMETERS
    KVIS = 1;    // CHANGE THIS TO READ STATEMENT
    GAMMA = 1.667;    // CHANGE THIS TO READ STATEMENT
    PRN = 1.0;    // CHANGE THIS TO READ STATEMENT
    RHO0 = 1;    // CHANGE THIS TO READ STATEMENT
    RMU0 = 2.2360612501607223E-004;    // CHANGE THIS TO READ STATEMENT
    P0 = 10.0;
    C0 = std::sqrt(GAMMA * P0 / RHO0);

    RM = 0.2;    // CHANGE THIS TO READ STATEMENT
    U0 = 2.0 * RM * C0;
    V0 = 0.0;

    // REYNOLDS NUMBER
    RE = 200.0;
    RLEN = 1.0 / 14.0;    // CHANGE THIS TO READ STATEMENT
    VTHICK0 = 2.0 * RE * RMU0 / (U0 * RHO0);
    MTHICK0 = 4.0 * VTHICK0;

    // READ RUN PARAMETERS
    NCYC = 1000;
    CFL = 0.6;    // CHANGE THIS TO READ STATEMENT

    // READ FILE NAMES
    IREAD = 0;  // IREAD = 1 FOR RESTART
    FNAME0 = "MLAYER.rst";  // RESTART FILE
}

