#include <vector>

std::vector<std::vector<std::vector<float>>> dw; 
std::vector<std::vector<std::vector<float>>> w0;

// Function to allocate memory for solution variables
void AllocateSolutionVariables(int nx, int ny, int nvar) {
    dw.resize(nx, std::vector<std::vector<float>>(ny, std::vector<float>(nvar, 0.0)));
    w0.resize(nx, std::vector<std::vector<float>>(ny, std::vector<float>(nvar, 0.0)));
    std::cout <<" solution variables\n";
}

void ResetVariables(int nx, int ny, int nvar, int ndim) {
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nvar; ++k) {
                dw[i][j][k] = 0.0;
                for (int l = 0; l < ndim; l++) {
                    wx[i][j][k][l] = 0.0;
                }
            }
        }
    }
}

// Function to deallocate memory for solution variables
void DeallocateSolutionVariables() {
    dw.clear();
    w0.clear();
}
