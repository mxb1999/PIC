#ifndef FIELDGATHER
    #define FIELDGATHER
    #include "simulation.hpp"
    //Implement a momentum conserving PIC scheme->interpolate fields from staggered grid to nodes
    void gather_e(Grid* grid, field_t* target_e);
    void gather_b(Grid* grid, field_t* target_b);
#endif