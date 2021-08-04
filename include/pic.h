
/*
    Header file describing basic PIC data structures and constants needed for up to three-dimensional code

*/
#ifndef PIC_HEAD


    #define PIC_HEAD

    #include <math.h>
    #include <stdio.h>
    #include <stdlib.h>
    #define IFEMPTY(object) {\
        if(object == NULL)\
        {\
            return NULL;\
        }\
    };
    #define NORM(vec) sqrt(pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2))
    #define DOUBLEARR(pointer, size) pointer=(double*)malloc(sizeof(double)*size)
    #define X(particle) particle->x[0]
    #define Y(particle) particle->x[1]
    #define Z(particle) particle->x[2]
    #define CROSS(a, b) {\
                            (a[1]*b[2] - a[2]*b[1])\
                            (a[2]*b[0] - a[0]*b[2])\
                            (a[0]*b[1] - a[1]*b[0])\
                            };
    #define PX(particle) particle->p[0]
    #define PY(particle) particle->p[1]
    #define PZ(particle) particle->p[2]
    #define XZone(particle, xmin, xmax, nx) (int)((xmax/(particle->x[0]-xmin)*nx))

    typedef struct
    {
        double p[3]; //store three dimensions of momentum
        double x[3]; //store three dimensions of position, but store as float
                    // express position within current zone
    }Particle;


    #define NX(grid) grid->dims[0]
    #define NY(grid) grid->dims[1]
    #define NZ(grid) grid->dims[2]
    #define n(grid) grid->n
    #define j(grid) grid->j
    #define rho(grid) grid->rho
    #define E(grid) grid->E
    #define B(grid) grid->B
    typedef struct
    {
        // Store the physical variables for mesh (or some subset of it)
        // Need to track: particle number, charge density, current density, Electric field, Magnetic field
        int dims[3];
        int nump;
        double minmax[6];
        double* n;
        double* rho;
        double* j;
        double* E;
        double* B;
    }Grid;

    extern void push(Particle** particles, Grid* grid, double dt);
    extern void interpolate(Particle** particles, Grid* grid, double dt);
    extern void solveFields(Particle** particles, Grid* grid, double dt);
    extern void updateParticles(Particle** particles, Grid* grid, double dt);
    
    extern void initialize(Particle** particles, Grid* grid);
    extern void pic_loop();

    //Constants
    const double sigma = 1.7e-4;
    const double e0 =8.85418782e-12;
    const double me =9.10938356e-31;
    const double pi =3.14159265359;
    const double kb= 1.3806485279e-16;   //Boltzmann constant in erg/K
    const double kb2= 1.3806485279e-23;   //Boltzmann constant in J/K
    const double ec= 1.60217662e-19;
    const double c= 29979245800.0;              // Speed of light in cm/s
    const double estat=4.80320427e-10; 	       // electron charge in statC
    const double Z = 3.1;                        // ionization state
    const double Te = 2.0e3*11604.5052;          // Temperature of electron in K
    const double Te_eV = 2.0e3;
    const double Ti = 1.0e3*11604.5052;          // Temperature of ion in K
    
#endif