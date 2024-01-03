float **dtl;
float *dtmin;
float simtime, viztime, vtime;
// std::vector< std::vector< float > > dtl;
int cyc, ncyc, counter;

// function to allocate memory for time variables
void AllocateTimeVariables(int nx, int ny) {
    // dtl.resize(nx, std::vector<float>(ny, 0.0));
    cudaMallocManaged(&dtmin, sizeof(float));
    cudaMalloc((void***)&dtl, nx * sizeof(float*));
    for (int i = 0; i < nx; ++i) {
        cudaMalloc((void**)&dtl , ny *  sizeof(float));
    }
}

// function to deallocate memory for time variables
void DeallocateTimeVariables() {
    // dtl.clear();
    cudaFree(dtl);
    cudaFree(dtmin);
}
