#include <fstream>

#include "flow_variables.h"
#include "solution_variables.h"

#include "boundary_conditions.h"
#include "output.h"

void InitializeField() {

    AllocateFlowVariables(ie, je, 4, 2);
    AllocateSolutionVariables(ie, je, 4);

    // local variables
    float eps, xc, yc, ey, ee;
    int i, j;

    float pi = 4.0 * atan(1.0);

    /*if (iread == 1) {  // read from previous solution for restarting
        std::ifstream infile(fname0.c_str());
        for (j = 0; j < je; ++j) {
            for (i = 0; i < ie; ++i) {
                infile >> w[0][i][j] >> w[1][i][j] >> w[2][i][j] >> w[3][i][j];
                rlv[i][j] = rmu0;
            }
        }
        infile >> time;
        infile.close();
    }*/

    if (iread == 0) {  // initialize flow field
        for (j = 0; j < je; ++j) {
            for (i = 0; i < ie; ++i) {
                xc = 0.5 * (coords[i][j][0] + coords[i + 1][j][0]);
                yc = 0.5 * (coords[i][j][1] + coords[i][j + 1][1]) - 10.0 * vthick0;
                ey = 2.0 * yc / rlen;
                ee = exp(-ey * ey);

                w[i][j][0] = rho0;
                w[i][j][1] = 0.5*u0*std::tanh(2.*yc/vthick0) -
                             ee*.05*yc*20.*vthick0/(2.*pi*10.)*cos(4.*pi*xc/(20.*vthick0))*exp(-pow(yc,2)/10.) -
                             ee*.025*yc*20.*vthick0/(pi*10.)*cos(2.*pi*xc/(20.*vthick0))*exp(-pow(yc,2)/10.);
                w[i][j][2] = 0.0;
                w[i][j][3] = p0;
                rlv[i][j] = rmu0;
            }
        }
    }

    ApplyBoundaryConditions();
    WriteVTK(0);
    std::cout << "Mixing layer initialised and inital field written to VTK \n";
}

void Finalize() {
    std::cout << "Clearing allocated memory for ";
    DeallocateMeshVariables();
    std::cout << "mesh variables ";
    DeallocateFlowVariables(ie, je, 4);
    std::cout << ", flow variables ";
    DeallocateSolutionVariables();
    std::cout << "& solution variables";
}
