#include "source.h"

static void get_discrete_position(Grid* grid, Particle* p, int* position)
{
    VOIDIFNULL(p);
    VOIDIFNULL(position);
    int nx, ny, nz;
    double dx, dy, dz;
    DEFINE_GRID_CONSTANTS;
    position[0] = (p->x-XMIN)/dx;
    position[1] = (p->y-YMIN)/dy;
    position[2] = (p->z-ZMIN)/dz;
}


void gatherCurrents(Grid* grid, double dt)
{
    int nx, ny, nz, np = NP;
    double dx, dy, dz;
    DEFINE_GRID_CONSTANTS;
    //interpolate particles to aligned grid
    int i;
    for(i = 0; i < np; i++)
    {
        Particle* p = grid->particles + i;
        double position[] = {p->x, p->y, p->z};
        double velocity[] = {p->px/M_I, p->py/M_I, p->pz/M_I};
        double ion_charge = Q_I;
        double vx = velocity[0], vy = velocity[1], vz = velocity[2];
        int indices[3];
        get_discrete_position(grid, p, indices);
        int ix = indices[0], iy = indices[1], iz = indices[2];
        double xdiff = (position[0]- ix*dx)/dx;
        double ydiff = (position[1]-iy*dy)/dy;
        double zdiff = (position[2] - iz*dz)/dz;
        double a[] = {1-xdiff, xdiff, 1-ydiff, ydiff, 1-zdiff, zdiff};
        RHO(ix+0,iy+0,iz+0) += ion_charge*a[3*0+0]*a[3*1+0]*a[3*2+0];
        RHO(ix+1,iy+0,iz+0) += ion_charge*a[3*0+1]*a[3*1+0]*a[3*2+0];
        RHO(ix+0,iy+1,iz+0) += ion_charge*a[3*0+0]*a[3*1+1]*a[3*2+0];
        RHO(ix+0,iy+0,iz+1) += ion_charge*a[3*0+0]*a[3*1+0]*a[3*2+1];
        RHO(ix+1,iy+1,iz+0) += ion_charge*a[3*0+1]*a[3*1+1]*a[3*2+0];
        RHO(ix+1,iy+0,iz+1) += ion_charge*a[3*0+1]*a[3*1+0]*a[3*2+1];
        RHO(ix+0,iy+1,iz+1) += ion_charge*a[3*0+0]*a[3*1+1]*a[3*2+1];
        RHO(ix+1,iy+1,iz+1) += ion_charge*a[3*0+1]*a[3*1+1]*a[3*2+1];

        Jx(ix+0,iy+0,iz+0) = vx*RHO(ix+0, iy+0, iz+0);
        Jx(ix+1,iy+0,iz+0) = vx*RHO(ix+1, iy+0, iz+0);
        Jx(ix+0,iy+1,iz+0) = vx*RHO(ix+0, iy+1, iz+0);
        Jx(ix+0,iy+0,iz+1) = vx*RHO(ix+0, iy+0, iz+1);
        Jx(ix+1,iy+1,iz+0) = vx*RHO(ix+1, iy+1, iz+0);
        Jx(ix+1,iy+0,iz+1) = vx*RHO(ix+1, iy+0, iz+1);
        Jx(ix+0,iy+1,iz+1) = vx*RHO(ix+0, iy+1, iz+1);
        Jx(ix+1,iy+1,iz+1) = vx*RHO(ix+1, iy+1, iz+1);

        Jy(ix+0,iy+0,iz+0) = vy*RHO(ix+0, iy+0, iz+0);
        Jy(ix+1,iy+0,iz+0) = vy*RHO(ix+1, iy+0, iz+0);
        Jy(ix+0,iy+1,iz+0) = vy*RHO(ix+0, iy+1, iz+0);
        Jy(ix+0,iy+0,iz+1) = vy*RHO(ix+0, iy+0, iz+1);
        Jy(ix+1,iy+1,iz+0) = vy*RHO(ix+1, iy+1, iz+0);
        Jy(ix+1,iy+0,iz+1) = vy*RHO(ix+1, iy+0, iz+1);
        Jy(ix+0,iy+1,iz+1) = vy*RHO(ix+0, iy+1, iz+1);
        Jy(ix+1,iy+1,iz+1) = vy*RHO(ix+1, iy+1, iz+1);

        Jz(ix+0,iy+0,iz+0) = vz*RHO(ix+0, iy+0, iz+0);
        Jz(ix+1,iy+0,iz+0) = vz*RHO(ix+1, iy+0, iz+0);
        Jz(ix+0,iy+1,iz+0) = vz*RHO(ix+0, iy+1, iz+0);
        Jz(ix+0,iy+0,iz+1) = vz*RHO(ix+0, iy+0, iz+1);
        Jz(ix+1,iy+1,iz+0) = vz*RHO(ix+1, iy+1, iz+0);
        Jz(ix+1,iy+0,iz+1) = vz*RHO(ix+1, iy+0, iz+1);
        Jz(ix+0,iy+1,iz+1) = vz*RHO(ix+0, iy+1, iz+1);
        Jz(ix+1,iy+1,iz+1) = vz*RHO(ix+1, iy+1, iz+1);
    }
    //interpolate to E nodes
    int i, j, k;
    for(i = 0; i < nx; i++)
    {
        for(j = 0; j < ny; j++)
        {
            for(k = 0; k < nz; k++)
            {
                if(i < nx-1)
                {
                    double jxavg = (Jx(i,j,k)+Jx(i+1,j,k))/2;
                    double epsavg = (EPSILON(i,j,k) + EPSILON(i+1,j,k))/2;
                    Ex(i,j,k) -= jxavg*dt/epsavg;
                }
                if(j < ny-1)
                {
                    double jyavg = (Jy(i,j,k)+Jy(i,j+1,k))/2;
                    double epsavg = (EPSILON(i,j,k) + EPSILON(i,j+1,k))/2;
                    Ey(i,j,k) -= jyavg*dt/epsavg;
                }
                if(k < nz-1)
                {
                    double jzavg = (Jz(i,j,k)+Jz(i,j,k+1))/2;
                    double epsavg = (EPSILON(i,j,k) + EPSILON(i,j,k+1))/2;
                    Ez(i,j,k) -= jzavg*dt/epsavg;
                }
            }
        }
    }
}