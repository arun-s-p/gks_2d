#include <iostream>
#include <fstream>
#include <iomanip>

// Include necessary headers (DIMS, FLO_VAR, MESH_VAR, FLO_PARAM, IODEFS, TIME_VAR) if defined in separate files

// Assuming these variables are global or defined in another header file
extern float W[4][IL][JL];
extern float X[3][IL][JL];
extern float XX[IL], YY[JL];
extern float VTHICK0, TIME;
extern int NX, NY, JL, IL, IE, JE;
extern int IOUT, IFLO;
extern std::string FNAME0, FNAME1, FNAME2;

void OUTPUT() {
    // Include DIMS, FLO_VAR, MESH_VAR, FLO_PARAM, IODEFS, TIME_VAR if not defined globally
    // Using the provided global variables, adjust as needed

    // LOCAL VARIABLES
    float XC, YC, DU, DV, OMEGA, EDDY_T, UTHICK, MT_T;
    float URA, RRA, UFAVRE, MTHICK, VTHICK, DUDY;
    float UAVG[JL];
    int I, J, N, IR;

    std::string jobname = "kh";
    EDDY_T = TIME / (1.0 / 14.0);

    IFLO = 18;
    ISIM = 28;

    std::ofstream outfile(FNAME1.c_str());

    // DUMP SOLUTION FOR RESTARTING
    for (J = 0; J < JE; ++J) {
        for (I = 0; I < IE; ++I) {
            N = static_cast<int>(EDDY_T) + RE;
            outfile << std::scientific << std::setw(18) << W[0][I][J]
                    << std::scientific << std::setw(18) << W[1][I][J]
                    << std::scientific << std::setw(18) << W[2][I][J]
                    << std::scientific << std::setw(18) << W[3][I][J] << "\n";
        }
    }
    outfile << TIME << "\n";
    outfile.close();

    // WRITE SOLUTION IN TECPLOT FORMAT
    std::ofstream tecplotFile(FNAME1.c_str());
    tecplotFile << "TITLE = BLAYER NEW SCHEME\n";
    tecplotFile << "VARIABLES = \"X\" \"Y\" \"RHO\" \"U\" \"V\" \"PRESSURE\" \"VORTICITY\"\n";
    tecplotFile << "Zone I=" << NX - 1 << ", J=" << NY - 1 << ", F=POINT\n";

    for (J = 2; J <= JL - 1; ++J) {
        for (I = 2; I <= IL - 1; ++I) {
            XC = 0.5 * (XX[I] + XX[I + 1]);
            YC = 0.5 * (YY[J] + YY[J + 1]);

            // COMPUTE VORTICITY
            float DX = X[0][I + 1][J] - X[0][I - 1][J];
            float DY = X[1][I][J + 1] - X[1][I][J - 1];

            DU = W[1][I][J + 1] - W[1][I][J - 1];
            DV = W[2][I + 1][J] - W[2][I - 1][J];
            OMEGA = DV / DX - DU / DY;

            tecplotFile << std::scientific << std::setw(18) << XC << std::scientific << std::setw(18) << YC
                        << std::scientific << std::setw(18) << W[0][I][J] << std::scientific << std::setw(18) << W[1][I][J]
                        << std::scientific << std::setw(18) << W[2][I][J] << std::scientific << std::setw(18) << W[3][I][J]
                        << std::scientific << std::setw(18) << OMEGA << "\n";
        }
    }
    tecplotFile.close();
}
