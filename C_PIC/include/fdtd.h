#ifndef FDTD
#define FDTD

    #include "simulation.hpp"


    extern void applyBoundaryTop(double* B, double* E, int numcells);
    extern void applyBoundaryBottom(double* B, double* E, int numcells);
    extern void applyBoundaryLeft(double* B, double* E, int numcells);
    extern void applyBoundaryRight(double* B, double* E, int numcells);

    extern void start_loop(Grid* grid, double dt, int nt);
    extern void start_loop1DFDTD(Grid* grid, double dt, int nt);
    extern void start_loop2DFDTD(Grid* grid, double dt, int nt);
    extern void start_loop3DFDTD(Grid* grid, double dt, int nt);
    extern void updateH2D_FirstTMz(Grid* grid);
    extern void updateH2DTMz(Grid* grid);
    extern void updateE2DTMz(Grid* grid);
    extern void updateE3D(Grid* grid);
    extern void updateH3D(Grid* grid);
    extern void write_data(double* arr, int* dims, char* target, int dim);

#endif