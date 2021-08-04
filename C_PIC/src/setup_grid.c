#include "fdtd.h"


void initialize(Grid* grid, double sigma_m, double sigma_e, double dt)
{
    int nx = grid->nx;
    int ny = grid->ny;
    int nz = grid->nz;
    double dx = (grid->xlims[1]-grid->xlims[0])/nx;
    double dy = (grid->ylims[1]-grid->ylims[0])/ny;
    double dz = (grid->zlims[1]-grid->zlims[0])/nz;
    double ratio = 1/sqrt(3);
    dt = dx*ratio/c;
    initializeE(grid, sigma_e, dt);
    initializeH(grid, sigma_m, dt);

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
    grid->ex = (double*)calloc(sizeof(double), (nx-1)*ny*nz);
    grid->ey = (double*)calloc(sizeof(double), nx*(ny-1)*nz);
    grid->ez = (double*)calloc(sizeof(double), nx*ny*(nz-1));
    grid->hx = (double*)calloc(sizeof(double), nx*(ny-1)*(nz-1));
    grid->hy = (double*)calloc(sizeof(double), (nx-1)*ny*(nz-1));
    grid->hz = (double*)calloc(sizeof(double), (nx-1)*(ny-1)*nz);
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
    grid->jx = (double*)calloc(sizeof(double), nx*ny*nz);
    grid->jy = (double*)calloc(sizeof(double), nx*ny*nz);
    grid->jz = (double*)calloc(sizeof(double), nx*ny*nz);
    grid->rho = (double*)calloc(sizeof(double), nx*ny*nz);
    grid->mu = (double*)calloc(sizeof(double), nx*ny*nz);
    grid->eps = (double*)calloc(sizeof(double), nx*ny*nz);
    grid->particles = (double*)calloc(sizeof(Particle), num_particles);
    return grid;
};
            