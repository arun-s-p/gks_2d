int nx, ny;    // mesh size
int *il, *jl;    // lower bounds in x and y directions
int ie, je;    // upper bounds in x and y directions
int ib, jb;    // bounds for ghost cells in x and y directions

std::vector< float > x, y;
std::vector< std::vector< std::vector< float > > > coords;
std::vector< std::vector< float > > volume;
float *dx, *dy;
float rlen;

// function to allocate memory for mesh variables
void AllocateMeshVariables(int nx, int ny, int ndim) {
  cudaMallocManaged(&dx, sizeof(float));
  cudaMallocManaged(&dy, sizeof(float));
  cudaMallocManaged(&il, sizeof(int));
  cudaMallocManaged(&jl, sizeof(int));
  x.resize(nx);
  y.resize(ny);

  coords.resize(nx);
  for (int i = 0; i < nx; ++i) {
    coords[i].resize(ny);
    for (int j = 0; j < ny; ++j) {
      coords[i][j].resize(ndim, 0.0);
    }
  }

  volume.resize(nx-1, std::vector<float>(ny-1, 0.0));
}

// function to deallocate memory for mesh variables
void DeallocateMeshVariables() {
  x.clear();
  y.clear();
  coords.clear();
  volume.clear();

  cudaFree(dx);
  cudaFree(dy);
  cudaFree(il);
  cudaFree(jl);
}
