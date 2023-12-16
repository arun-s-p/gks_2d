#include <iostream>
#include <cmath>

#include "io_variables.h"
#include "flow_parameters.h"
#include "time_variables.h"
#include "mesh_variables.h"

#include "input.h"
#include "mesh.h"

int main() {
    int nstart;
    double mom, tstop, tscale, viztime, vtime;
    char fname[20];

    std::cout << "Reading inputs...\n";
    input();
    // INIT();
    std::cout << "Reference velocity = " << u0 << 
              "\nInitial vorticity thickness = " << vthick0 << 
              "\nReynolds number = " << re << 
              "\nMach number = " << rm << std::endl;

    std::cout << "Generating mesh...\n";
    mesh();

    // counter = 0;
    // viztime = 0.0;
    // vtime = 1.0;
    // tstop = 20.0;
    // counter++;

    // if (IREAD == 0) {
    //     TIME = 0.0;
    //     NSTART = 100;

    //     // Do NSTART cycles to get a good estimate of initial pressure
    //     for (int CYC = 1; CYC <= NSTART; ++CYC) {
    //         TSTEP();
    //         DTL[0][0] = DTMIN / static_cast<double>(NSTART);
    //         std::cout << DTMIN << std::endl;
    //         DTMIN = DTMIN / static_cast<double>(NSTART);

    //         UPDATE();
    //         TIME += DTMIN;
    //     }
    // }

    // viztime = viztime + vtime;

    // // Actual time stepping starts here
    // for (int CYC = 1; CYC <= NCYC; ++CYC) {
    //     TSTEP();
    //     UPDATE();
    //     TIME += DTMIN;

    //     if (CYC % 50 == 0)
    //         std::cout << CYC << " " << DTMIN << " " << TIME << std::endl;

    //     tscale = TIME * 0.5 * U0 / VTHICK0;

    //     if (tscale > viztime) {
    //         OUTPUT();
    //         std::cout << "----------writing---------- t = " << TIME << std::endl;
    //         counter++;
    //         viztime = viztime + vtime;
    //     }

    //     if (tscale > tstop)
    //         break;
    // }

    // OUTPUT();
    DeallocateMeshVariables();

    return 0;
}

