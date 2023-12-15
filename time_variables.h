#include <vector>

float simtime, dtmin, viztime, vtime;
std::vector<std::vector<float>> dtl;
int cyc, ncyc, counter;

// function to allocate memory for time variables
void AllocateTimeVariables(int nx, int ny) {
    dtl.resize(nx, std::vector<float>(ny, 0.0));
}

// function to deallocate memory for time variables
void DeallocateTimeVariables() {
    dtl.clear();
}
