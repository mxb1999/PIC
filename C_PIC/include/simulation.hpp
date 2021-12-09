/*
    Define the overarching structure of the simulation in C
    - Define the initial geometry of the fields and the sources
        - In initial sim 'antennae' used to generate EM waves, external B field applied w/ helmholtz coils
        - Need to differentiate between E sources, EM sources, and B sources.
    -




*/
#include <hdf5/serial/hdf5.h>
#include <stdio.h>
#include <cuda_runtime.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#ifndef SIMULATION
    #define SIMULATION
    typedef double part_t; // type to experiment with mixed precision down the line

    enum distribution {gaussian, unif};
    typedef struct
    {
        part_t x, y, z;
        part_t px, py, pz;
    } Particle;
    typedef struct
    {
        double xlims[2], ylims[2], zlims[2];
        int nx, ny, nz, num_particles;
        double dx, dy, dz;
        //Electric field components
        double *ex, *ey, *ez;
        //Auxiliary magnetic field components
        double *hx, *hy, *hz;
        //Precalculated update coefficients for E updates
        double *exself, *eyself, *ezself;
        double *exhyc, *exhzc, *eyhxc, *eyhzc, *ezhxc, *ezhyc;
        //Precalculated update coefficients for H updates
        double *hxself, *hyself, *hzself;
        double *hxeyc, *hxezc, *hyexc, *hyezc, *hzexc, *hzeyc;
        double *mu, *eps;
        //Arrays used to store old values for ABC calculation
        double *oldHx, *oldHy, *oldHz;
        double *oldEx, *oldEy, *oldEz;
        Particle* particles;
        double *jx, *jy, *jz, *rho; //Current and charge density arrays for the mesh
        long long flags;
        double mass_p;
        double q_p;
    }Grid;
    typedef double space_t; // same as above, but for spatial variables
    typedef double field_t;
    #ifndef M_PI
        #define M_PI 3.14159265358979323846
    #endif
    #define GRID grid //the name used to refer to the grid in all functions
    #define EXSIZE (NX-1)*(NY)*(NZ)
    #define EYSIZE (NX)*(NY-1)*(NZ)
    #define EZSIZE (NX)*(NY)*(NZ-1)
    #define BXSIZE (NX)*(NY-1)*(NZ-1)
    #define BYSIZE (NX-1)*(NY)*(NZ-1)
    #define BZSIZE (NX-1)*(NY-1)*(NZ)
    #define M_I GRID->mass_p
    #define Q_I GRID->q_p
    #define PARTICLES GRID->particles
    #define Jx(i, j, k) ACCESS3D(GRID->jx, i, j, k, ny, nz)
    #define Jy(i, j, k) ACCESS3D(GRID->jy, i, j, k, ny, nz)
    #define Jz(i, j, k) ACCESS3D(GRID->jz, i, j, k, ny, nz)
    #define RHO(i, j, k) ACCESS3D(GRID->rho, i, j, k, ny, nz)
    #define Jx_cu(i, j, k) ACCESS3D(GRID.jx, i, j, k, ny, nz)
    #define Jy_cu(i, j, k) ACCESS3D(GRID.jy, i, j, k, ny, nz)
    #define Jz_cu(i, j, k) ACCESS3D(GRID.jz, i, j, k, ny, nz)
    #define RHO_cu(i, j, k) ACCESS3D(GRID.rho, i, j, k, ny, nz)
    #define ACCESS3D(array, i, j, k, size2, size3) (array)[((i)*(size2)+(j))*(size3)+(k)]
    #define Ex(i, j, k) ACCESS3D(GRID->ex, i, j, k, ny, nz)
    #define Ey(i, j, k) ACCESS3D(GRID->ey, i, j, k, ny-1, nz)
    #define Ez(i, j, k) ACCESS3D(GRID->ez, i, j, k, ny, nz-1)
    #define Hx(i, j, k) ACCESS3D(GRID->hx, i, j, k, ny-1, nz-1)
    #define Hy(i, j, k) ACCESS3D(GRID->hy, i, j, k, ny, nz-1)
    #define Hz(i, j, k) ACCESS3D(GRID->hz, i, j, k, ny-1, nz)
    #define Ex_cu(i, j, k) ACCESS3D(GRID.ex, i, j, k, ny, nz)
    #define Ey_cu(i, j, k) ACCESS3D(GRID.ey, i, j, k, ny-1, nz)
    #define Ez_cu(i, j, k) ACCESS3D(GRID.ez, i, j, k, ny, nz-1)
    #define Hx_cu(i, j, k) ACCESS3D(GRID.hx, i, j, k, ny-1, nz-1)
    #define Hy_cu(i, j, k) ACCESS3D(GRID.hy, i, j, k, ny, nz-1)
    #define Hz_cu(i, j, k) ACCESS3D(GRID.hz, i, j, k, ny-1, nz)
    #define CExSelf(i, j, k) ACCESS3D(GRID->exself, i, j, k, ny, nz)
    #define CEySelf(i, j, k) ACCESS3D(GRID->eyself, i, j, k, ny-1, nz)
    #define CEzSelf(i, j, k) ACCESS3D(GRID->ezself, i, j, k, ny, nz-1)
    #define CHxSelf(i, j, k) ACCESS3D(GRID->hxself, i, j, k, ny-1, nz-1)
    #define CHySelf(i, j, k) ACCESS3D(GRID->hyself, i, j, k, ny, nz-1)
    #define CHzSelf(i, j, k) ACCESS3D(GRID->hzself, i, j, k, ny-1, nz)
    #define CExhy(i, j, k) ACCESS3D(GRID->exhyc, i, j, k, ny, nz)
    #define CEyhz(i, j, k) ACCESS3D(GRID->eyhzc, i, j, k, ny-1, nz)
    #define CEzhx(i, j, k) ACCESS3D(GRID->ezhxc, i, j, k, ny, nz-1)
    #define CEzhy(i, j, k) ACCESS3D(GRID->ezhyc, i, j, k, ny, nz-1)
    #define CExhz(i, j, k) ACCESS3D(GRID->exhzc, i, j, k, ny, nz)
    #define CEyhx(i, j, k) ACCESS3D(GRID->eyhxc, i, j, k, ny-1, nz)

    #define CHxey(i, j, k) ACCESS3D(GRID->hxeyc, i, j, k, ny-1, nz-1)
    #define CHyex(i, j, k) ACCESS3D(GRID->hyexc, i, j, k, ny, nz-1)
    #define CHxez(i, j, k) ACCESS3D(GRID->hxezc, i, j, k, ny-1, nz-1)
    #define CHyez(i, j, k) ACCESS3D(GRID->hyezc, i, j, k, ny, nz-1)
    #define CHzey(i, j, k) ACCESS3D(GRID->hzeyc, i, j, k, ny-1, nz)
    #define CHzex(i, j, k) ACCESS3D(GRID->hzexc, i, j, k, ny-1, nz)

    #define EX_MIN_OLD(i, j) (GRID->oldEx + ((i)*nz + j)*12)//accessing at nx-1
    #define EX_MAX_OLD(i, j) (GRID->oldEx + ((i)*nz + j)*12 + 6)
    #define EY_MIN_OLD(i, j) (GRID->oldEy + ((i)*nz + j)*12)//accessing at nx-1
    #define EY_MAX_OLD(i, j) (GRID->oldEy + ((i)*nz + j)*12 + 6)
    #define EZ_MIN_OLD(i, j) (GRID->oldEz + ((i)*nz + j)*12)//accessing at nx-1
    #define EZ_MAX_OLD(i, j) (GRID->oldEz + ((i)*nz + j)*12 + 6)

    #define HX_MIN_OLD(i, j) (GRID->oldHx + ((i)*nz + j)*12)//accessing at nx-1
    #define HX_MAX_OLD(i, j) (GRID->oldHx + ((i)*nz + j)*12 + 6)
    #define HY_MIN_OLD(i, j) (GRID->oldHy + ((i)*nz + j)*12)//accessing at nx-1
    #define HY_MAX_OLD(i, j) (GRID->oldHy + ((i)*nz + j)*12 + 6)
    #define HZ_MIN_OLD(i, j) (GRID->oldHz + ((i)*nz + j)*12)//accessing at nx-1
    #define HZ_MAX_OLD(i, j) (GRID->oldHz + ((i)*nz + j)*12 + 6)

    #define MU(i,j,k) ACCESS3D(GRID->mu, i, j, k, ny, nz)
    #define EPSILON(i,j,k) ACCESS3D(GRID->eps, i, j, k, ny, nz)

    #define NX GRID->nx
    #define NY GRID->ny
    #define NZ GRID->nz
    #define DX GRID->dx
    #define DY GRID->dy
    #define DZ GRID->dz
    #define NP GRID->num_particles

    #define XMIN GRID->xlims[0]
    #define XMAX GRID->xlims[1]
    #define YMIN GRID->ylims[0]
    #define YMAX GRID->ylims[1]
    #define ZMIN GRID->zlims[0]
    #define ZMAX GRID->zlims[1]

    #define DEFINE_GRID_CONSTANTS {\
        nx = NX, ny = NY, nz = NZ;\
        dx = DX, dy = (YMAX-YMIN)/(double)ny, dz = (ZMAX-ZMIN)/(double)nz;\
    };

    #define X_PARTITION 1
    #define Y_PARTITION 2
    #define Z_PARTITION 3
    #define CELL(i,j,k,d,ny,nz) ((((i)*ny+(j))*nz+(k))*3 + d)

    #define ACCESS4D(array, i, j, k, d, size2, size3, size4) array[(((i)*size2+(j))*size3+(k))*size4+(d)]

    #define CURLXNEXT(vector, i, j, k, nx, ny, nz, dy, dz, curl){\
        double dvz_dy = (vector[CELL(i,j+1,k,2,ny,nz)] - vector[CELL(i,j,k,2,ny,nz)])/dy;\
        double dvy_dz = (vector[CELL(i,j,k+1,1,ny,nz)] - vector[CELL(i,j,k,1,ny,nz)])/dz;\
        curl[0] = dvz_dy-dvy_dz;\
    };
        //printf("%e %f %f\n", dvz_dx, vector[CELL(i+1,j,k,2,ny,nz)], vector[CELL(i,j,k,2,ny,nz)]);\
        //printf("X %f %f\n", dvz_dy, dvy_dz);
    #define CURLYNEXT(vector, i, j, k, nx, ny, nz, dx, dz, curl){\
        double dvz_dx = (vector[CELL(i+1,j,k,2,ny,nz)] - vector[CELL(i,j,k,2,ny,nz)])/dx;\
        double dvx_dz = (vector[CELL(i,j,k+1,0,ny,nz)] - vector[CELL(i,j,k,0,ny,nz)])/dz;\
        curl[1] = dvx_dz-dvz_dx;\
    };
        //printf("Y %f %f\n", dvx_dz, dvz_dx);
    #define CURLZNEXT(vector, i, j, k, nx, ny, nz, dx, dy, curl){\
        double dvy_dx = (vector[CELL(i+1,j,k,1,ny,nz)] - vector[CELL(i,j,k,1,ny,nz)])/dx;\
        double dvx_dy = (vector[CELL(i,j+1,k,0,ny,nz)] - vector[CELL(i,j,k,0,ny,nz)])/dy;\
        curl[2] = dvy_dx-dvx_dy;\
    };
        //printf("Z %f %f\n", dvy_dx, dvx_dy);

    #define CURLNEXT(vector, i, j, k, nx, ny, nz, dx, dy, dz, curl){\
        CURLXNEXT(vector, i, j, k, nx, ny, nz, dy, dz, curl);\
        CURLYNEXT(vector, i, j, k, nx, ny, nz, dx, dz, curl);\
        CURLZNEXT(vector, i, j, k, nx, ny, nz, dx, dy, curl);\
    };
    #define CURLXPREV(vector, i, j, k, nx, ny, nz, dy, dz, curl){\
        double dvz_dy = (vector[CELL(i,j,k,2,ny,nz)] - vector[CELL(i,j-1,k,2,ny,nz)])/dy;\
        double dvy_dz = (vector[CELL(i,j,k,1,ny,nz)] - vector[CELL(i,j,k-1,1,ny,nz)])/dz;\
        curl[0] = dvz_dy-dvy_dz;\
    };
        //printf("%e %f %f\n", dvz_dx, vector[CELL(i+1,j,k,2,ny,nz)], vector[CELL(i,j,k,2,ny,nz)]);\
        //printf("X %f %f\n", dvz_dy, dvy_dz);
    #define CURLYPREV(vector, i, j, k, nx, ny, nz, dx, dz, curl){\
        double dvz_dx = (vector[CELL(i,j,k,2,ny,nz)] - vector[CELL(i-1,j,k,2,ny,nz)])/dx;\
        double dvx_dz = (vector[CELL(i,j,k,0,ny,nz)] - vector[CELL(i,j,k-1,0,ny,nz)])/dz;\
        curl[1] = dvx_dz-dvz_dx;\
    };
        //printf("Y %f %f\n", dvx_dz, dvz_dx);
    #define CURLZPREV(vector, i, j, k, nx, ny, nz, dx, dy, curl){\
        double dvy_dx = (vector[CELL(i,j,k,1,ny,nz)] - vector[CELL(i-1,j,k,1,ny,nz)])/dx;\
        double dvx_dy = (vector[CELL(i,j,k,0,ny,nz)] - vector[CELL(i,j-1,k,0,ny,nz)])/dy;\
        curl[2] = dvy_dx-dvx_dy;\
    };
        //printf("Z %f %f\n", dvy_dx, dvx_dy);

    void free_grid(Grid* grid);
    void setup_grid_constb(Grid* grid, distribution position_distribution, distribution momentum_distribution, const double temperature, double variance, field_t b[3]);
    #define CURLPREV(vector, i, j, k, nx, ny, nz, dx, dy, dz, curl){\
        CURLXPREV(vector, i, j, k, nx, ny, nz, dy, dz, curl);\
        CURLYPREV(vector, i, j, k, nx, ny, nz, dx, dz, curl);\
        CURLZPREV(vector, i, j, k, nx, ny, nz, dx, dy, curl);\
    };
    extern Grid* define_grid(double* dimsx, double* dimsy, double* dimsz,
                            int nx, int ny, int nz,
                            int num_particles,
                            int permeable_bc);

    extern Grid* split_grid(Grid* original, int partition_dim, int partition_point);
    #define IFNOMEMRET(pointer) {\
        if(pointer == NULL)\
        {\
            printf("No memory available for allocation\n");\
            return NULL;\
        }\
    };
    #define CROSSX(a, b) a[2]*b[3] - a[3]*b[2]
    #define CROSSY(a, b) a[3]*b[1] - a[1]*b[3]
    #define CROSSZ(a, b) a[1]*b[2] - a[2]*b[1]
    #define CROSS(a, b) {CROSSX(a, b), CROSSY(a, b), CROSSZ(a, b)}
    #define IFNOMEM(pointer) {\
        if(pointer == NULL)\
        {\
            return;\
        }\
    };
    #define DYNAMIC_ARRAY(pointer, size, number){\
        pointer = (double*)malloc(size*number);\
        IFNOMEM(pointer);\
    };
    #define VOIDIFNULL(pointer){\
        if(pointer == NULL)\
        {\
            return;\
        }\
    };
    #define DOT(a, b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])
    #define MAG(a) sqrt(DOT(a, a))
    extern void initializeE(Grid* grid, double sigma_e, double dt);
    extern void initializeH(Grid* grid, double sigma_e, double dt);
    void writeArr(void* arr, int type, char* filename, char* name, int dimnum, int* dimlist);
void initialize(Grid* grid, double sigma_m, double sigma_e, double* dt, distribution position_distribution, distribution momentum_distribution, const double temperature, double variance);
    extern void interpolate_E(Grid* grid, int ix, int iy, int iz, double* target);
    extern void interpolate_B(Grid* grid, int ix, int iy, int iz, double* target);
    Grid* new_grid(const int nx, const int ny, const int nz, const int num_particles, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double mass, double charge);
    extern void update_from_particles(Grid* grid, double dt);
    void current_deposition(Grid* grid, double dt);
    void gpu_particle_push(Grid* grid, double step, part_t* logger, int stepnum);
void solve_fields_gpu(Grid* grid, double dt, double sigma_e, double sigma_m, double epsilon, double mu);

#define c  299792458.0
#define eps0  8.85418782e-12
#define mu0  1.25663706e-6
#define freq  14.0e6
#define pi  3.14159265
#endif