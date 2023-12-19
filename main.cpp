#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "io_variables.h"
#include "flow_parameters.h"
#include "time_variables.h"
#include "mesh_variables.h"

#include "input.h"
#include "mesh.h"
#include "initialize.h"
#include "time_step.h"
#include "update.h"

int main() {
    int nstart;
    double mom, tstop, tscale, viztime, vtime;
    char fname[20];

    std::cout << "Reading inputs...\n";
    InputParameters();
    // INIT();
    std::cout << "Reference velocity = " << u0 << 
              "\nInitial vorticity thickness = " << vthick0 << 
              "\nReynolds number = " << re << 
              "\nMach number = " << rm << std::endl;

    std::cout << "Generating mesh...\n";
    CartesianMesh2D();

    InitializeField();

    counter = 1;
    viztime = 0.0;
    vtime = 1.0;
    tstop = 20.0;

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

    viztime += vtime;

    // Actual time stepping starts here
    for ( cyc = 1; cyc <= ncyc; ++cyc) {
        CalculateDT();
        UpdateField();
        simtime += dtmin;

        // if (CYC % 50 == 0)
            // std::cout << CYC << " " << DTMIN << " " << TIME << std::endl;

        tscale = simtime * 0.5 * u0 / vthick0;

        if (tscale > viztime) {
            WriteVTK(counter);
            std::cout << "----------writing---------- t = " << simtime << std::endl;
            std::cout << dtmin << std::endl;
            counter++;
            viztime += vtime;
        }

        // if (tscale > tstop)
            // break;
    }

    // OUTPUT();
    Finalize();

    return 0;
}

