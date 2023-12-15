#include <vector>

// Assuming these variables are global or defined in another header file
extern float TIME, DTMIN, viztime, vtime;
extern std::vector<std::vector<float>> DTL;
extern int CYC, NCYC, counter;

// Equivalent structure to encapsulate time variables
struct TimeVariables {
    float TIME, DTMIN, viztime, vtime;
    std::vector<std::vector<float>> DTL;
    int CYC, NCYC, counter;
};

// Function to allocate memory for time variables
void AllocateTimeVariables(int nx, int ny) {
    DTL.resize(nx, std::vector<float>(ny, 0.0));
}

// Function to deallocate memory for time variables
void DeallocateTimeVariables() {
    DTL.clear();
}

// Function to initialize time variables
TimeVariables InitializeTimeVariables() {
    TimeVariables timeVars;
    timeVars.TIME = 0.0;
    timeVars.DTMIN = 0.0;
    timeVars.viztime = 0.0;
    timeVars.vtime = 0.0;
    // ... (initialize other variables similarly)

    return timeVars;
}

