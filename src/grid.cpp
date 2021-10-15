#include "grid.h"
#include "push.hu"
#include "rng.hpp"

Grid* new_grid(const int nx, const int ny, const int nz, const int numparticles, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double mass, double charge)
{
    Grid* result = new Grid;
    space_t dx, dy, dz;
    dx = (xmax-xmin)/nx;
    dy = (ymax-ymin)/ny;
    dz = (zmax-zmin)/nz;
    *result = {NULL, nx, ny, nz, numparticles, dx, dy, dz, xmin, xmax, ymin, ymax, zmin, zmax, NULL, NULL, mass, charge};
    result->particles = new Particle[numparticles];
    result->e_field = new field_t[nx*ny*nz*3]{0.0};
    result->b_field = new field_t[nx*ny*nz*3]{0.0};
    return result;
};
void free_grid(Grid* grid)
{
    delete [] grid->particles;
    delete [] grid->b_field;
    delete [] grid->e_field;
    delete grid;
}
void setup_grid_constb(Grid* grid, distribution position_distribution, distribution momentum_distribution, const double temperature, double variance, field_t b[3])
{
    double posvariancex = fmin(fabs(grid->xlims[0]), fabs(grid->xlims[1]));
    double posvariancey = fmin(fabs(grid->ylims[0]), fabs(grid->ylims[1]));
    double posvariancez = fmin(fabs(grid->zlims[0]), fabs(grid->zlims[1]));
    double xspan = grid->xlims[1] - grid->xlims[0];
    double yspan = grid->ylims[1] - grid->ylims[0];
    double zspan = grid->zlims[1] - grid->zlims[0];
    double xmin = grid->xlims[0];
    double ymin = grid->ylims[0];
    double zmin = grid->zlims[0];
    int nump = grid->numparticles;
    switch(position_distribution)
    {
        case gaussian:
            {RandGen* gen_positionx = initializeRNG(412<<3, 0, posvariancex);
            RandGen* gen_positiony = initializeRNG(341<<2, 0, posvariancey);
            RandGen* gen_positionz = initializeRNG(234<<1, 0, posvariancez);
            for(int i = 0; i < nump; i++)
            {
                Particle* p = grid->particles + i;
                p->x = normal_deviation(gen_positionx);
                p->y = normal_deviation(gen_positiony);
                p->z = normal_deviation(gen_positionz);
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
    field_t* bfield = grid->b_field;
    int nx = grid->nx;
    int ny = grid->ny;
    int nz = grid->nz;
    for(int i = 0; i < nx*ny*nz*3; i += 3)
    {
        bfield[i] = b[0];
        bfield[i+1] = b[1];
        bfield[i+2] = b[2];
    }
};