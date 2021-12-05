#include "push.hu"
#include "grid.h"
#include "rng.hpp"
#include <chrono>
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
//#define PRINTPROG
void test_push(int n, int p)
{
    double m = 1e-28;
    Grid* grid = new_grid(n, n, n, p, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, m, 1);
    int nx = grid->nx;
    int ny = grid->ny;
    int nz = grid->nz;
    field_t b[] = {0, 0, 1e3};
    setup_grid_constb(grid, unif, unif, 1e3*m,500*m, b);
    Grid* grid_cu;
    Particle* particle_cu;
    field_t *e_cu, *b_cu;
    cudaError_t err = cudaMalloc(&grid_cu, sizeof(Grid));
    CHECKERR(err, __LINE__)
    err = cudaMalloc(&particle_cu, sizeof(Particle)*grid->numparticles);
    CHECKERR(err, __LINE__)
    err = cudaMalloc(&e_cu, sizeof(field_t)*nx*ny*nz*3);
    CHECKERR(err, __LINE__)
    err = cudaMalloc(&b_cu, sizeof(field_t)*nx*ny*nz*3);
    CHECKERR(err, __LINE__)
    err = cudaMemcpy(grid_cu, grid, sizeof(Grid), cudaMemcpyHostToDevice);
    CHECKERR(err, __LINE__)
    err = cudaMemcpy(e_cu, grid->e_field, sizeof(field_t)*nx*ny*nz*3, cudaMemcpyHostToDevice);
    CHECKERR(err, __LINE__)
    err = cudaMemcpy(particle_cu, grid->particles, sizeof(Particle)*grid->numparticles, cudaMemcpyHostToDevice);
    CHECKERR(err, __LINE__)
    err = cudaMemcpy(b_cu, grid->b_field, sizeof(field_t)*nx*ny*nz*3, cudaMemcpyHostToDevice);
    CHECKERR(err, __LINE__)
    double time = 1;
    double step = 1e-5;
    
    int nt = 1000;
    part_t* logger_cu, *logger;
    #ifdef LOG
        logger = (part_t*)malloc(sizeof(part_t)*grid->numparticles*3*(nt));
        cudaMalloc(&logger_cu, sizeof(part_t)*grid->numparticles*3*(nt));
    #endif
    int mthreads = max_threads;
    int ppt = grid->numparticles/((int)max_threads) + 1;
    int blocks = fmin(grid->numparticles/threads_per_block + 1, SMs*blocks_per_sm);
    int interval = nt/10;
    auto begin = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= nt; i++)
    {
        push_stage(particle_cu, e_cu, b_cu,grid_cu, 5e-6, ppt, logger_cu, blocks, i, 1);
        #ifdef PRINTPROG
            if(i % interval == 0)
            {
                printf("%d%c done\n", (int)(ceil((double)i/nt * 100)), '%');
            }
        #endif
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf("%d, %d, %e\n", n, grid->numparticles, elapsed.count()*1e-9);
    cudaDeviceSynchronize();
    begin = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= nt; i++)
    {
        push_stage(grid->particles, grid->e_field, grid->b_field,grid, 5e-6, ppt, logger, blocks, i, 0);
        #ifdef PRINTPROG
            if(i % interval == 0)
            {
                printf("%d%c done\n", (int)(ceil((double)i/nt * 100)), '%');
            }
        #endif
    }
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf("%d, %d, %e\n", n, grid->numparticles, elapsed.count()*1e-9);
    err = cudaMemcpy(grid, grid_cu, sizeof(Grid), cudaMemcpyDeviceToHost);
    CHECKERR(err, __LINE__)
    err = cudaMemcpy(grid->e_field, e_cu, sizeof(field_t)*nx*ny*nz*3, cudaMemcpyDeviceToHost);
    CHECKERR(err, __LINE__)
    err = cudaMemcpy(grid->particles, particle_cu, sizeof(Particle)*grid->numparticles, cudaMemcpyDeviceToHost);
    CHECKERR(err, __LINE__)
    err = cudaMemcpy(grid->b_field, b_cu, sizeof(field_t)*nx*ny*nz*3, cudaMemcpyDeviceToHost);
    CHECKERR(err, __LINE__)
    #ifdef LOG
        err = cudaMemcpy(logger, logger_cu, sizeof(part_t)*grid->numparticles*3*(nt), cudaMemcpyDeviceToHost);
        CHECKERR(err, __LINE__)
        cudaFree(logger_cu);
    #endif
    err = cudaFree(grid_cu);
    cudaFree(e_cu);
    cudaFree(b_cu);
    cudaFree(particle_cu);

    #ifdef LOG
        hid_t file, set, space;
        herr_t err_h5;
        hsize_t dims[] = {(nt), grid->numparticles, 3};
        file = H5Fcreate("positions.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        space = H5Screate_simple(3, dims, NULL);
        set = H5Dcreate2(file, "positions", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        err_h5 = H5Dwrite(set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, logger);
        err_h5 = H5Dclose(set);
        err_h5 = H5Sclose(space);
        err_h5 = H5Fclose(file);
        free(logger);
    #endif
    free_grid(grid);
}
int main(int argc, char** argv)
{
    //#define BENCHMARK
    #ifdef BENCHMARK
    int trials = 4;
    int n, p;
    p = 1000000;
    for(n = 50; n <= 800; n+=50)
    {
        test_push(n, p);
    }
    n = 800;
    for(p = 50000; p <= 10000000; p += 50000)
    {
        test_push(n, p);
    }
    #endif
    #ifndef BENCHMARK
        test_push(10, 100000000);
    #endif
    return 0;
}
