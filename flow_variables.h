float ***w;
float ****wx;
float **rlv;
float rhom, um;

// Function to allocate memory for flow variables
void AllocateFlowVariables(int nx, int ny, int nvar, int ndim) {
  cudaMallocManaged(&w, nx * sizeof(float**));
  // w = new float **[nx];
  // wx = new float ***[nx];
  rlv = new float *[nx];

  for (int i = 0; i < nx; ++i) {
    cudaMallocManaged(&w[i], ny * sizeof(float*));
    // w[i] = new float *[ny];
    // wx[i] = new float **[ny];
    rlv[i] = new float[ny];

    for (int j = 0; j < ny; ++j) {
      // w[i][j] = new float[nvar];
      // wx[i][j] = new float *[nvar];
      cudaMallocManaged(&w[i][j], nvar * sizeof(float));
      for (int k = 0; k < nvar; ++k) {
        w[i][j][k] = 0.0;
        // wx[i][j][k] = new float[ndim];
      }
    }
  }

  std::cout << "Memory allocated for flow variables, ";
}

// Function to deallocate memory for flow variables
void DeallocateFlowVariables(int nx, int ny, int nvar) {
  for (int i = 0; i < nx; ++i) {
    // for (int j = 0; j < ny; ++j) {
    //    delete[] w[i][j];

      // for (int k = 0; k < nvar; ++k) {
      //   delete[] wx[i][j][k];
      // }
      // delete[] wx[i][j];
    // }

    // delete[] w[i];
    // delete[] wx[i];
    delete[] rlv[i];
  }

  // delete[] w;
  // delete[] wx;
  delete[] rlv;

  cudaFree(w);
}

