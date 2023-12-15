#include <vector>

// Assuming these variables are global or defined in another header file
extern std::vector<float> XX, YY;
extern std::vector<std::vector<std::vector<float>>> X;
extern std::vector<std::vector<float>> VOL;
extern float RLEN, DX, DY;

// Equivalent structure to encapsulate mesh variables
struct MeshVariables {
  std::vector<float> XX, YY;
  std::vector<std::vector<std::vector<float>>> X;
  std::vector<std::vector<float>> VOL;
  float RLEN, DX, DY;
};

// Function to allocate memory for mesh variables
void AllocateMeshVariables(int nx, int ny, int ndim) {
  XX.resize(nx);
  YY.resize(ny);

  X.resize(nx);
  for (int i = 0; i < nx; ++i) {
    X[i].resize(ny);
    for (int j = 0; j < ny; ++j) {
      X[i][j].resize(ndim, 0.0);
    }
  }

  VOL.resize(nx, std::vector<float>(ny, 0.0));
}

// Function to deallocate memory for mesh variables
void DeallocateMeshVariables() {
  XX.clear();
  YY.clear();
  X.clear();
  VOL.clear();
}

// Function to initialize mesh variables
MeshVariables InitializeMeshVariables() {
  MeshVariables meshVars;
  meshVars.RLEN = 0.0;
  meshVars.DX = 0.0;
  meshVars.DY = 0.0;
  // ... (initialize other variables similarly)

  return meshVars;
}
