/*Necessary Aspects of Problem Geometry:
- Define min and max Z for each coil
- Define internal boundary conditions for E and B
*/

typedef struct 
{
    double z;
    double dz;
    double radius_inner;
    double radius_outer;
} FRC;

extern double* internal_BC_E; //Set the internal boundary conditions for E within the rings: 4D: 3, nx, ny, nz (though the last cell of each is not used along its own dimension)-> E[0] == nx-1, ny, nz
extern double* internal_BC_B; //Set the internal boundary conditions for B within the rings: 4D: 3, nx, ny, nz (though the last cells of other 2 dimensions are not allocated)-> B[0] == nx, ny-1, nz-1
extern double* radial_epsilon;



extern FRC* get_rings(int numrings,
                      double* zlocations,
                      double* dz,
                      double* inner_r,
                      double* outer_r);