void WriteVTK(int counter) {

    // local variables
    float xc, yc, du, dv, omega;
    int i, j, np;

    std::string jobname = "kh";

    // write solution in tecplot format
    std::ostringstream fname1_stream;
    fname1_stream << jobname << std::setfill('0') << std::setw(4) << counter << ".vtk";
    fname1 = fname1_stream.str();
    std::ofstream vtkFile(fname1.c_str());

    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "VTK from C++\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET STRUCTURED_GRID\n";
    vtkFile << "DIMENSIONS " << nx-1 << " " << ny-1 << " 1\n";
    np = (nx-1)*(ny-1);
    vtkFile << "POINTS " << np << " float\n";

    // Write coordinates
    for (int j = 2; j < jl-1; ++j) {
        for (int i = 2; i < il-1; ++i) {
            xc = 0.5 * (x[i] + x[i + 1]);
            yc = 0.5 * (y[j] + y[j + 1]);
            vtkFile << xc << " " << yc << " " << 0.0 << "\n";
        }
    }

    vtkFile << "POINT_DATA " << np << "\n";
    vtkFile << "SCALARS RHO float\n";
    vtkFile << "LOOKUP_TABLE default\n";

    for (int j = 2; j < jl - 1; ++j) {
        for (int i = 2; i < il - 1; ++i) {
            vtkFile << w[i][j][0] << "\n";
        }
    }

    // vtkFile << "POINT_DATA " << np << "\n";
    vtkFile << "SCALARS U float\n";
    vtkFile << "LOOKUP_TABLE default\n";

    for (int j = 2; j < jl - 1; ++j) {
        for (int i = 2; i < il - 1; ++i) {
            vtkFile << w[i][j][1] << "\n";
        }
    }

    // vtkFile << "POINT_DATA " << np << "\n";
    vtkFile << "SCALARS V float\n";
    vtkFile << "LOOKUP_TABLE default\n";

    for (int j = 2; j < jl - 1; ++j) {
        for (int i = 2; i < il - 1; ++i) {
            vtkFile << w[i][j][2] << "\n";
        }
    }

    // vtkFile << "POINT_DATA " << np << "\n";
    vtkFile << "SCALARS Pressure float\n";
    vtkFile << "LOOKUP_TABLE default\n";

    for (int j = 2; j < jl - 1; ++j) {
        for (int i = 2; i < il - 1; ++i) {
            vtkFile << w[i][j][3] << "\n";
        }
    }

    vtkFile.close();
}

void RestartDump(){
    
    int i, j;
    std::ofstream outfile(fname0.c_str());

    // dump solution for restarting
    for (j = 0; j < je; ++j) {
        for (i = 0; i < ie; ++i) {
            outfile << std::scientific << std::setw(18) 
                    << w[i][j][0] << ' ' << w[i][j][1] << ' ' 
                    << w[i][j][2] << ' ' << w[i][j][3] << '\n';
        }
    }
    outfile << time << "\n";
    outfile.close();
}
