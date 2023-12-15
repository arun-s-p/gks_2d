int nx, ny;    // mesh size
int il, jl;    // lower bounds in x and y directions
int ie, je;    // upper bounds in x and y directions
int ib, jb;    // bounds for ghost cells in x and y directions

std::vector<float> xx, yy;
std::vector<std::vector<std::vector<float>>> x;
std::vector<std::vector<float>> volume;
float rlen, dx, dy;

// function to allocate memory for mesh variables
void AllocateMeshVariables(int nx, int ny, int ndim) {
  xx.resize(nx);
  yy.resize(ny);

  x.resize(nx);
  for (int i = 0; i < nx; ++i) {
    x[i].resize(ny);
    for (int j = 0; j < ny; ++j) {
      x[i][j].resize(ndim, 0.0);
    }
  }

  volume.resize(nx, std::vector<float>(ny, 0.0));
}

// function to deallocate memory for mesh variables
void DeallocateMeshVariables() {
  xx.clear();
  yy.clear();
  x.clear();
  volume.clear();
}
