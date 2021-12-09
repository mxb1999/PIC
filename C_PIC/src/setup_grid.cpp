#include "fdtd.h"
#include "simulation.hpp"
#include "push.hu"
#include "rng.hpp"
#define mu0 1.256637062e-6
void dist_particles(Grid* grid, distribution position_distribution, distribution momentum_distribution, const double temperature, double variance);

void initialize(Grid* grid, double sigma_m, double sigma_e, double* dt, distribution position_distribution, distribution momentum_distribution, const double temperature, double variance)
{
    int nx = grid->nx;
    int ny = grid->ny;
    int nz = grid->nz;
    double dx = (grid->xlims[1]-grid->xlims[0])/nx;
    double dy = (grid->ylims[1]-grid->ylims[0])/ny;
    double dz = (grid->zlims[1]-grid->zlims[0])/nz;
    double ratio = 1/sqrt(3);
    *dt = dx*ratio/c;
    initializeE(grid, sigma_e, *dt);
    initializeH(grid, sigma_m, *dt);
    dist_particles(grid, position_distribution, momentum_distribution, temperature, variance);
};
/*
void free_grid(Grid* grid)
{
    delete [] grid->particles;
    delete [] grid->hx;
    delete [] grid->e;
    delete [] grid->eself;
    delete [] grid->hself;
    delete [] grid->old_e;
    delete [] grid->old_h;
    delete [] grid->j;
    delete grid;
}*/
void dist_particles(Grid* grid, distribution position_distribution, distribution momentum_distribution, const double temperature, double variance) {
    double posvariancex = fmin(fabs(grid->xlims[0]), fabs(grid->xlims[1]));
    double posvariancey = fmin(fabs(grid->ylims[0]), fabs(grid->ylims[1]));
    double posvariancez = fmin(fabs(grid->zlims[0]), fabs(grid->zlims[1]));
    double xspan = grid->xlims[1] - grid->xlims[0];
    double yspan = grid->ylims[1] - grid->ylims[0];
    double zspan = grid->zlims[1] - grid->zlims[0];
    double xmin = grid->xlims[0];
    double ymin = grid->ylims[0];
    double zmin = grid->zlims[0];
    int nump = grid->num_particles;
    printf("Mass %e\n", M_I);
    switch(position_distribution)
    {
        case gaussian:
            {
            RandGen* gen_positionx = initializeRNG(412<<3, 0, posvariancex);
            RandGen* gen_positiony = initializeRNG(341<<2, 0, posvariancey);
            RandGen* gen_positionz = initializeRNG(234<<1, 0, posvariancez);
            for(int i = 0; i < nump; i++)
            {
                Particle* p = grid->particles + i;
                p->x = normal_deviation(gen_positionx);
                p->y = normal_deviation(gen_positiony);
                p->z = normal_deviation(gen_positionz);
                printf("%e %e %e\n", p->x, p->y, p->z);
            };
            free(gen_positionx);
            free(gen_positiony);
            free(gen_positionz);
            break;}
        case unif:{
            RandGen* generator = initializeRNG(412<<3);
            for(int i = 0; i < nump; i++)
            {
                Particle* p = grid->particles + i;
                rng_t randx = int64RNG(generator);
                rng_t randy = int64RNG(generator);
                rng_t randz = int64RNG(generator);
                p->x = xspan*((double)randx/ULLONG_MAX)+xmin;
                p->y = yspan*((double)randy/ULLONG_MAX)+ymin;
                p->z = zspan*((double)randz/ULLONG_MAX)+zmin;
            };}
        break;
    }
    switch(momentum_distribution)
    {
        case gaussian:{
            RandGen* gen_x = initializeRNG(412<<3, temperature, variance);
            RandGen* gen_y = initializeRNG(341<<2, temperature, variance);
            RandGen* gen_z = initializeRNG(234<<1, temperature, variance);
            for(int i = 0; i < nump; i++)
            {
                Particle* p = grid->particles + i;
                p->px = normal_deviation(gen_x);
                p->py = normal_deviation(gen_y);
                p->pz = normal_deviation(gen_z);
            };
            free(gen_x);
            free(gen_y);
            free(gen_z);
            break;}
        case unif:
            {RandGen* generator = initializeRNG(13123);
            for(int i = 0; i < nump; i++)
            {
                Particle* p = grid->particles + i;
                rng_t randx = int64RNG(generator);
                rng_t randy = int64RNG(generator);
                rng_t randz = int64RNG(generator);

                double px = ((double)randx/ULLONG_MAX);
                double py = ((double)randy/ULLONG_MAX);
                double pz = ((double)randz/ULLONG_MAX);
                double norm = temperature/sqrt(px*px + py*py + pz*pz);
                px*=norm;
                py*=norm;
                pz*=norm;
                p->px = px;
                p->py = py;
                p->pz = pz;
            };}
        break;
    }
}
void setup_grid_constb(Grid* grid, distribution position_distribution, distribution momentum_distribution, const double temperature, double variance, field_t b[3])
{


    field_t* hx = grid->hx;
    field_t* hy = grid->hy;
    field_t* hz = grid->hz;
    int nx = grid->nx;
    int ny = grid->ny;
    int nz = grid->nz;
    for(int i = 0; i < nx; i++) {
        for(int j = 0; j < ny; j++) {
            for(int k = 0; k < nz; k++) {
                MU(i, j, k) = mu0;
            }
        }
    }
    for(int i = 0; i < nx-1; i++) {
        for(int j = 0; j < ny-1; j++) {
            for(int k = 0; k < nz; k++) {
                Hz(i, j, k) = b[2]/MU(i, j, k);
            }
        }
    }
    for(int i = 0; i < nx-1; i++) {
        for(int j = 0; j < ny; j++) {
            for(int k = 0; k < nz-1; k++) {
                Hy(i, j, k) = b[1]/MU(i, j, k);
            }
        }
    }
    for(int i = 0; i < nx; i++) {
        for(int j = 0; j < ny-1; j++) {
            for(int k = 0; k < nz-1; k++) {
                Hx(i, j, k) = b[0]/MU(i, j, k);
            }
        }
    }
};
Grid* define_grid(double* dimsx, double* dimsy, double* dimsz,
                  int nx, int ny, int nz,
                  int num_particles,
                  int permeable_bc){
    Grid* grid = (Grid*)malloc(sizeof(Grid));
    IFNOMEMRET(grid);
    int buff = 4*permeable_bc;
    grid->num_particles=num_particles;
    grid->xlims[0] = dimsx[0];
    grid->xlims[1] = dimsx[1];
    grid->ylims[0] = dimsy[0];
    grid->ylims[1] = dimsy[1];
    grid->zlims[0] = dimsz[0];
    grid->zlims[1] = dimsz[1];
    grid->nx = nx;
    grid->ny = ny;
    grid->nz = nz;
    grid->dx = (XMAX-XMIN)/(double)nx;
    grid->dy = (YMAX-YMIN)/(double)ny;
    grid->dz = (ZMAX-ZMIN)/(double)nz;
    /*grid->ex = (double*)calloc(sizeof(double), (nx-1)*ny*nz);
    grid->ey = (double*)calloc(sizeof(double), nx*(ny-1)*nz);
    grid->ez = (double*)calloc(sizeof(double), nx*ny*(nz-1));
    grid->hx = (double*)calloc(sizeof(double), nx*(ny-1)*(nz-1));
    grid->hy = (double*)calloc(sizeof(double), (nx-1)*ny*(nz-1));
    grid->hz = (double*)calloc(sizeof(double), (nx-1)*(ny-1)*nz);*/
    grid->exself = (double*)calloc(sizeof(double), (nx-1)*ny*nz);
    grid->exhyc = (double*)calloc(sizeof(double), (nx-1)*ny*nz);
    grid->exhzc = (double*)calloc(sizeof(double), (nx-1)*ny*nz);
    grid->eyself = (double*)calloc(sizeof(double), nx*(ny-1)*nz);
    grid->eyhxc = (double*)calloc(sizeof(double), nx*(ny-1)*nz);
    grid->eyhzc = (double*)calloc(sizeof(double), nx*(ny-1)*nz);
    grid->ezself = (double*)calloc(sizeof(double), nx*ny*(nz-1));
    grid->ezhyc = (double*)calloc(sizeof(double), nx*ny*(nz-1));
    grid->ezhxc = (double*)calloc(sizeof(double), nx*ny*(nz-1));
    grid->hxself = (double*)calloc(sizeof(double), nx*(ny-1)*(nz-1));
    grid->hxeyc = (double*)calloc(sizeof(double), nx*(ny-1)*(nz-1));
    grid->hxezc = (double*)calloc(sizeof(double), nx*(ny-1)*(nz-1));
    grid->hyself = (double*)calloc(sizeof(double), (nx-1)*ny*(nz-1));
    grid->hyexc = (double*)calloc(sizeof(double), (nx-1)*ny*(nz-1));
    grid->hyezc = (double*)calloc(sizeof(double), (nx-1)*ny*(nz-1));
    grid->hzself = (double*)calloc(sizeof(double), (nx-1)*(ny-1)*nz);
    grid->hzeyc = (double*)calloc(sizeof(double), (nx-1)*(ny-1)*nz);
    grid->hzexc = (double*)calloc(sizeof(double), (nx-1)*(ny-1)*nz);

    grid->oldEx = (double*)calloc(sizeof(double), 12*ny*nz);
    grid->oldEy = (double*)calloc(sizeof(double), 12*nx*nz);
    grid->oldEz = (double*)calloc(sizeof(double), 12*ny*nx);
    grid->oldHx = (double*)calloc(sizeof(double), 12*ny*nz);
    grid->oldHy = (double*)calloc(sizeof(double), 12*nx*nz);
    grid->oldHz = (double*)calloc(sizeof(double), 12*ny*nx);
    /*grid->jx = (double*)calloc(sizeof(double), nx*ny*nz);
    grid->jy = (double*)calloc(sizeof(double), nx*ny*nz);
    grid->jz = (double*)calloc(sizeof(double), nx*ny*nz);
    grid->rho = (double*)calloc(sizeof(double), nx*ny*nz);*/
    grid->mu = (double*)calloc(sizeof(double), nx*ny*nz);
    grid->eps = (double*)calloc(sizeof(double), nx*ny*nz);
    cudaMallocManaged(&grid->particles, sizeof(Particle)*num_particles);
    cudaMallocManaged(&grid->ex, sizeof(field_t)*(nx-1)*ny*nz);
    cudaMallocManaged(&grid->jx, sizeof(field_t)*nx*ny*nz);
    cudaMallocManaged(&grid->jy, sizeof(field_t)*nx*ny*nz);
    cudaMallocManaged(&grid->jz, sizeof(field_t)*nx*ny*nz);
    cudaMallocManaged(&grid->rho, sizeof(field_t)*nx*ny*nz);

    cudaMallocManaged(&grid->ey, sizeof(field_t)*nx*(ny-1)*nz);
    cudaMallocManaged(&grid->ez, sizeof(field_t)*nx*ny*(nz-1));
    cudaMallocManaged(&grid->hx, sizeof(field_t)*nx*(ny-1)*(nz-1));
    cudaMallocManaged(&grid->hy, sizeof(field_t)*(nx-1)*ny*(nz-1));
    cudaMallocManaged(&grid->hz, sizeof(field_t)*(nx-1)*(ny-1)*nz);
    return grid;
};
