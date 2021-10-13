#include "push.hu"
#include "grid.h"
#include "rng.hpp"
#include <hdf5/serial/hdf5.h>
RandGen* initializeRNG(rng_t j)
{
    RandGen* rng = new RandGen;
    rng->u = j ^ rng->v;
    rng->v = rng->u;
    rng->w = rng->v;
    return rng;
};
RandGen* initializeRNG(rng_t j, double mu, double sigma)
{
    RandGen* rng = new RandGen;
    rng->u = j ^ rng->v;
    rng->v = rng->u;
    rng->w = rng->v;
    rng->mu = mu;
    rng->sigma = sigma;
    return rng;
};
int main(int argc, char** argv)
{
    double m = 1e-28;
    Grid* grid = new_grid(10, 10, 10, 10, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, m, 1);
    int nx = grid->nx;
    int ny = grid->ny;
    int nz = grid->nz;
    field_t b[] = {0, 0, 1e3};
    setup_grid_constb(grid, unif, unif, 1e3*m,500*m, b);
    Grid* grid_cu;
    Particle* particle_cu;
    field_t *e_cu, *b_cu;
    cudaMalloc(&grid_cu, sizeof(Grid));
    cudaMalloc(&particle_cu, sizeof(Particle)*grid->numparticles);
    cudaMalloc(&e_cu, sizeof(field_t)*nx*ny*nz*3);
    cudaMalloc(&b_cu, sizeof(field_t)*nx*ny*nz*3);
	printf("%d\n",sizeof(Particle));
    cudaMemcpy(grid_cu, grid, sizeof(Grid), cudaMemcpyHostToDevice);
    cudaMemcpy(e_cu, grid->e_field, sizeof(field_t)*nx*ny*nz*3, cudaMemcpyHostToDevice);
    cudaMemcpy(particle_cu, grid->particles, sizeof(Particle)*grid->numparticles, cudaMemcpyHostToDevice);
    cudaMemcpy(b_cu, grid->b_field, sizeof(field_t)*nx*ny*nz*3, cudaMemcpyHostToDevice);
    double time = 1;
    double step = 1e-5;
    
    printf("\nREE %e %e %e\n", (grid->b_field)[0], (grid->b_field)[1], (grid->b_field)[2]);
    int nt = time/step;
    part_t* logger_cu, *logger;
    logger = (part_t*)malloc(sizeof(part_t)*grid->numparticles*3*(nt));
    cudaMalloc(&logger_cu, sizeof(part_t)*grid->numparticles*3*(nt));
    
    for(int i = 1; i <= nt; i++)
    {
        push_stage(particle_cu, e_cu, b_cu,grid_cu, 5e-6, 1, logger_cu, i);
    }
    cudaMemcpy(grid, grid_cu, sizeof(Grid), cudaMemcpyDeviceToHost);
    cudaMemcpy(grid->e_field, e_cu, sizeof(field_t)*nx*ny*nz*3, cudaMemcpyDeviceToHost);
    cudaMemcpy(grid->particles, particle_cu, sizeof(Particle)*grid->numparticles, cudaMemcpyDeviceToHost);
    cudaMemcpy(grid->b_field, b_cu, sizeof(field_t)*nx*ny*nz*3, cudaMemcpyDeviceToHost);
    cudaMemcpy(logger, logger_cu, sizeof(part_t)*grid->numparticles*3*(nt), cudaMemcpyDeviceToHost);
    cudaFree(grid_cu);
    cudaFree(logger_cu);
    cudaFree(e_cu);
    cudaFree(b_cu);
    cudaFree(particle_cu);
    hid_t file, set, space;
    herr_t err;
    hsize_t dims[] = {(nt), grid->numparticles, 3};
    file = H5Fcreate("positions.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    space = H5Screate_simple(3, dims, NULL);
    set = H5Dcreate2(file, "positions", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    err = H5Dwrite(set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, logger);
    err = H5Dclose(set);
    err = H5Sclose(space);
    err = H5Fclose(file);
    free_grid(grid);
    free(logger);
    return 0;
}
