#include "simulation.hpp"

//assume linear interpolation -> Energy Conserving algorithm
void interpolate_E(Grid* grid, int ix, int iy, int iz, double* position, double* target)
{
    int nx, ny, nz;
    double dx, dy, dz;
    DEFINE_GRID_CONSTANTS;
    //assume energy conserving
    double xdiff = (target[0]- ix*dx)/dx;
    double ydiff = (target[1]-iy*dy)/dy;
    double zdiff = (target[2] - iz*dz)/dz;
    double a[] = {1-xdiff, xdiff, 1-ydiff, ydiff, 1-zdiff, zdiff};
    double SE[] = {a[2*1] * a[2*2],     a[2*0] * a[2*2]    ,   a[2*0] * a[2*1],
                   0,                   a[2*0+1] * a[2*2]  , a[2*0+1] * a[2*1],
                   a[2*1+1] * a[2*2],   0                  , a[2*0] * a[2*1+1],
                   0,                   0                  ,   a[2*0+1] * a[2*1+1],
                   a[2*1] * a[2*2+1],   a[2*0] * a[2*2+1]  , 0,
                   0,                   a[2*0+1] * a[2*2+1], 0,
                   a[2*1+1] * a[2*2+1], 0,                   0,
                   0,                   0,                   0};
    //iterate for E in X direction
    target[0] += SE[0*3+0]*Ex(ix, iy, iz);
    target[0] += SE[2*3+0]*Ex(ix, iy+1, iz);
    target[0] += SE[4*3+0]*Ex(ix, iy, iz+1);
    target[0] += SE[6*3+0]*Ex(ix, iy+1, iz+1);
    target[1] += SE[0*3+1]*Ey(ix, iy, iz);
    target[1] += SE[1*3+1]*Ey(ix+1, iy, iz);
    target[1] += SE[4*3+1]*Ey(ix, iy, iz+1);
    target[1] += SE[5*3+1]*Ey(ix+1, iy, iz+1);
    target[2] += SE[0*3+2]*Ez(ix, iy, iz);
    target[2] += SE[1*3+2]*Ez(ix+1, iy, iz);
    target[2] += SE[2*3+2]*Ez(ix, iy+1, iz);
    target[2] += SE[3*3+2]*Ez(ix+1, iy+1, iz);
};
void interpolate_B(Grid* grid, int ix, int iy, int iz, double* position, double* target)
{
    int nx, ny, nz;
    double dx, dy, dz;
    DEFINE_GRID_CONSTANTS;
    //assume energy conserving
    double xdiff = (target[0]- ix*dx)/dx;
    double ydiff = (target[1]-iy*dy)/dy;
    double zdiff = (target[2] - iz*dz)/dz;
    double a[] = {1-xdiff, xdiff, 1-ydiff, ydiff, 1-zdiff, zdiff};
    double SB[] = {a[2*0]  , a[2*1]  , a[2*2],
                   a[2*0+1], 0       , 0,
                   0       , a[2*1+1], 0,
                   0       , 0       , 0,
                   0       , 0       , 0,
                   0       , 0       , a[2*2+1],
                   0       , 0       , 0,
                   0       , 0       , 0};
    //add all B components
    target[0] += SB[0*3+0]*Bx(ix, iy, iz);
    target[0] += SB[1*3+0]*Bx(ix+1, iy, iz);
    target[1] += SB[0*3+1]*By(ix, iy, iz);
    target[1] += SB[2*3+1]*By(ix, iy+1, iz);
    target[2] += SB[0*3+2]*Bz(ix, iy, iz);
    target[2] += SB[5*3+2]*Bz(ix, iy, iz+1);
};