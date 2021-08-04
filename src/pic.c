#include "pic.h"

double mi_kg;	   // Mass of ion in kg
double mi;          // Mass of ion in g
Particle* new_empty_particle()
{
    Particle* result = (Particle*)malloc(sizeof(Particle));
    IFEMPTY(result);
    return result;
}

Particle* new_particle(double* position, double* momentum)
{
    Particle* result = (Particle*)malloc(sizeof(Particle));
    IFEMPTY(result);
    X(result) = position[0];
    Y(result) = position[1];
    Z(result) = position[2];
    PX(result) = momentum[0];
    PY(result) = momentum[1];
    PZ(result) = momentum[2];
    return result;
}

Grid* new_grid(int* dims, double* minmax, int numParticles)
{
    Grid* result = (Grid*)malloc(sizeof(Grid));
    IFEMPTY(result);
    int nx = dims[0];
    int ny = dims[1];
    int nz = dims[2];
    NX(result) = nx;
    NY(result) = ny;
    NZ(result) = nz;
    result->nump = numParticles;
    result->minmax[0] = minmax[0]; // min x
    result->minmax[1] = minmax[1]; // max x
    result->minmax[2] = minmax[2]; // min y
    result->minmax[3] = minmax[3]; // max y
    result->minmax[4] = minmax[4]; // min z
    result->minmax[5] = minmax[5]; // max z
    DOUBLEARR(B(result), nx*3);
    DOUBLEARR(E(result), nx*3);
    DOUBLEARR(n(result), nx);
    DOUBLEARR(rho(result), nx);
    DOUBLEARR(j(result), nx*3);
    return result;
}

void push(Particle** particles, Grid* grid, double dt)
{
    //use Boris algorithm to solve for the velocities at the split step
    int nx = grid->dims[0];
    double xmin = grid->minmax[0];
    double xmax = grid->minmax[1];
    int particle_count = grid->nump;
    double q = Z*ec;
    for(int i = 0; i < particle_count; i++)
    {
        Particle* particle = particles[i];
        double speed = NORM(particle->p)/mi_kg;
        double gamma = 1/sqrt(1-pow(speed/c,2));
        double* p = particle->p;
        printf("%d : %f %f %f\n", i, p[0], p[1], p[2]);
        X(particle) += p[0]/(gamma*mi_kg) * dt; 
        int ix = XZone(particle, xmin, xmax, nx);
        double epsilon = q/(2)*E(grid)[ix];
        double p_minus[3] = {p[0]+epsilon*dt, p[1]+epsilon*dt,p[2]+epsilon*dt};
        double speed_minus = NORM(p_minus);
        double gamma_minus = 1/sqrt(1-pow(speed_minus/c,2));
        double BX = B(grid)[ix*3];
        double BY = B(grid)[ix*3+1];
        double BZ = B(grid)[ix*3+2];
        double bMag = sqrt(pow(BX,2) + pow(BY,2) + pow(BZ,2));
        double theta = q*dt*bMag/gamma_minus;
        double tantheta = tan(theta/2);
        double t_vec[3] = {tantheta/bMag*BX, tantheta/bMag*BY, tantheta/bMag*BZ};
        double crossprod[3] = CROSS(p_minus, t_vec);
        double pprime[3] =  {p_minus[0]+crossprod[0], p_minus[1]+crossprod[1], p_minus[2]+crossprod[2]};
        double tmag = NORM(t_vec);
        double t_const = 2/(1+pow(tmag,2));
        double pprime_cross_t[3] = CROSS(pprime, t_vec);
        double pplus[3] = {p_minus[0]+pprime_cross_t[0]*t_const, p_minus[1]+pprime_cross_t[1]*t_const, p_minus[2]+pprime_cross_t[2]*t_const};
        p[0] = pplus[0]+epsilon*dt;
        p[1] = pplus[1]+epsilon*dt;
        p[2] = pplus[2]+epsilon*dt;
        printf("\t%f %f %f\n", i, p[0], p[1], p[2]);
    }


};
void interpolate(Particle** particles, Grid* grid, double dt)
{

};
void solveFields(Particle** particles, Grid* grid, double dt)
{

};
void updateParticles(Particle** particles, Grid* grid, double dt)
{

};

void initialize(Particle** particles, Grid* grid)
{

};
#define NPC 5
#define NUM_X 20
void pic_loop()
{
    double mi_kg = 10230.0*me;	   // Mass of ion in kg
    double mi = 10230*(1.0e3*me);          // Mass of ion in g
    int num_particles = NPC*NUM_X;
    Particle** p = (Particle**)malloc(sizeof(Particle*)*num_particles);
    int dims[3] = {NUM_X,NUM_X,NUM_X};
    double minmax[6] = {-5e-4,5e-4,-5e-4,5e-4,-5e-4,5e-4};
    int* dimaddr = &dims[0];
    Grid* grid = new_grid(dimaddr, &minmax[0], num_particles);
    printf("Success\n");
};



int main(int argc, char** argv)
{
    
    pic_loop();
    return 0;
}
