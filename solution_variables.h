#include <vector>

std::vector<std::vector<std::vector<float>>> dw, fw, vw;
std::vector<std::vector<std::vector<float>>> w0;

// Function to allocate memory for solution variables
void AllocateSolutionVariables(int nx, int ny, int nvar) {
    dw.resize(nx, std::vector<std::vector<float>>(ny, std::vector<float>(nvar, 0.0)));
    fw.resize(nx, std::vector<std::vector<float>>(ny, std::vector<float>(nvar, 0.0)));
    vw.resize(nx, std::vector<std::vector<float>>(ny, std::vector<float>(nvar, 0.0)));
    w0.resize(nx, std::vector<std::vector<float>>(ny, std::vector<float>(nvar, 0.0)));
    std::cout <<" solution variables\n";
}

// Function to deallocate memory for solution variables
void DeallocateSolutionVariables() {
    dw.clear();
    fw.clear();
    vw.clear();
    w0.clear();
}
