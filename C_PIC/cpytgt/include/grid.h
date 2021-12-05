/** Matthew Burns
 *
 * Store the representation of a 3D spatial grid
 * Store particle positions and velocities
 * Store field information
 *  */
#ifndef GRID
#define GRID
typedef double part_t; // type to experiment with mixed precision down the line
typedef struct
{
    part_t x, y, z;
    part_t px, py, pz;
} Particle;

typedef double space_t; // same as above, but for spatial variables
typedef double field_t;
enum distribution {gaussian, unif};
typedef struct
{
    //needed for particle push
    Particle* particles; //array of particles in the grid
    int nx, ny, nz, num_particles; //number of spatial zones in each dimension, number of particles in grid
    space_t dx, dy, dz, xlims[2], ylims[2], zlims[2]; //width of zone in each dimension, spatial limits (max and min) in each dimension
    field_t* b_field; //Magnetic field, flattened 3D array
    field_t* e_field; //Electric field, flattened 3D array
    part_t mass_p;
    part_t q_p;
} Grid;

Grid* new_grid(const int nx, const int ny, const int nz, const int num_particles, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double mass, double charge);
void free_grid(Grid* grid);

void setup_grid_constb(Grid* grid, distribution position_distribution, distribution momentum_distribution, const double temperature, double variance, field_t b[3]);
#endif