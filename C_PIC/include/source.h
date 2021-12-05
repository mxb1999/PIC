/*Define data structures necessary for EM sources
Define different functions for adding E and B to the grid
*/
#include "simulation.hpp"
typedef struct {
    int num_segments; //store the number of wire segments
    int isfunc; //store if the current is a function pointer or a static value
    double* x; //store the x extrema of the wire segments: 2*num_segments
    double* y; //store the y extrema of the wire segments: 2*num_segments
    double* z; //store the z extrema of the wire segments: 2*num_segments
    double* current; //either a pointer to J or a pointer to a function which returns a time dependent J
} Wire;
typedef struct
{

} Coil;
#define END(wire, segment) {wire->x[segment*2+1],\
                            wire->y[segment*2+1],\
                            wire->z[segment*2+1],};
#define START(wire, segment) {wire->x[segment*2],\
                              wire->y[segment*2],\
                              wire->z[segment*2],};

#define SETEND(wire, segment, xmax, ymax, zmax) {\
    wire->x[segment*2+1] = xmax;\
    wire->y[segment*2+1] = ymax;\
    wire->z[segment*2+1] = zmax;\
};

#define SETBEGINNING(wire, segment, xmin, ymin, zmin) {\
    wire->x[segment*2+1] = xmin;\
    wire->y[segment*2+1] = ymin;\
    wire->z[segment*2+1] = zmin;\
};


extern double* apply_wire(Wire* source, Grid* GRID);
extern double* apply_particle(Particle* source, Grid* GRID);