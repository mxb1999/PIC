/*
    Definitions for faraday cage surrounding the PFRC apparatus, enforcing BCs


*/

typedef struct
{
    double xlim[2];
    double ylim[2];
    double zlim[2];
    int xlimcell[2];
    int ylimcell[2];
    int zlimcell[2];
} Wire;

typedef struct
{
    Wire* wires;

} FaradayCage;


