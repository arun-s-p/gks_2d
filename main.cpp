#include <iostream>
#include <cmath>

// Assuming these variables are global or defined in another header file
extern int IREAD, NCYC, counter;
extern double U0, VTHICK0, RE, RM;
extern double TIME, DTMIN;
extern double **DTL;
extern double **W, **W0, **DW;
extern double **XX, **YY;
extern double RLEN;

void INPUT();
void MESH();
void INIT();
void TSTEP();
void UPDATE();
void OUTPUT();

int main() {
    int NSTART;
    double mom, tstop, tscale, viztime, vtime;
    char FNAME[20];

    INPUT();
    MESH();
    INIT();

    std::cout << "U0 = " << U0 << " Vt = " << VTHICK0 << " RE = " << RE << " M = " << RM << std::endl;

    counter = 0;
    viztime = 0.0;
    vtime = 1.0;
    tstop = 20.0;
    counter++;

    if (IREAD == 0) {
        TIME = 0.0;
        NSTART = 100;

        // Do NSTART cycles to get a good estimate of initial pressure
        for (int CYC = 1; CYC <= NSTART; ++CYC) {
            TSTEP();
            DTL[0][0] = DTMIN / static_cast<double>(NSTART);
            std::cout << DTMIN << std::endl;
            DTMIN = DTMIN / static_cast<double>(NSTART);

            UPDATE();
            TIME += DTMIN;
        }
    }

    viztime = viztime + vtime;

    // Actual time stepping starts here
    for (int CYC = 1; CYC <= NCYC; ++CYC) {
        TSTEP();
        UPDATE();
        TIME += DTMIN;

        if (CYC % 50 == 0)
            std::cout << CYC << " " << DTMIN << " " << TIME << std::endl;

        tscale = TIME * 0.5 * U0 / VTHICK0;

        if (tscale > viztime) {
            OUTPUT();
            std::cout << "----------writing---------- t = " << TIME << std::endl;
            counter++;
            viztime = viztime + vtime;
        }

        if (tscale > tstop)
            break;
    }

    OUTPUT();

    return 0;
}

