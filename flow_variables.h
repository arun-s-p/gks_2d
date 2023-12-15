#include <vector>

// Assuming these variables are global or defined in another header file
extern float ***W;
extern float ****WX;
extern float **RLV;
extern float RHOM, UM;

// Equivalent structure to encapsulate dynamically allocated variables
struct FlowVariables {
  float ***W;
  float ****WX;
  float **RLV;
  float RHOM;
  float UM;
};

// Function to allocate memory for flow variables
void AllocateFlowVariables(int nx, int ny, int nvar, int ndim) {
  W = new float **[nx];
  WX = new float ***[nx];
  RLV = new float *[nx];

  for (int i = 0; i < nx; ++i) {
    W[i] = new float *[ny];
    WX[i] = new float **[ny];
    RLV[i] = new float[nvar];

    for (int j = 0; j < ny; ++j) {
      W[i][j] = new float[nvar];
      WX[i][j] = new float *[nvar];

      for (int k = 0; k < nz; ++k) {
        W[i][j][k] = 0.0;
        WX[i][j][k] = new float[ndim];
      }
    }
  }
}

// Function to deallocate memory for flow variables
void DeallocateFlowVariables(int nx, int ny, int nvar) {
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      delete[] W[i][j];
      delete[] WX[i][j];

      for (int k = 0; k < nvar; ++k) {
        delete[] WX[i][j][k];
      }
    }

    delete[] W[i];
    delete[] WX[i];
    delete[] RLV[i];
  }

  delete[] W;
  delete[] WX;
  delete[] RLV;
}

// Equivalent structure to encapsulate statically allocated variables
struct FlowConstants {
  float RHOM;
  float UM;
};

// Function to initialize flow constants
FlowConstants InitializeFlowConstants() {
  FlowConstants flowConstants;
  flowConstants.RHOM = 0.0;
  flowConstants.UM = 0.0;
  return flowConstants;
}
