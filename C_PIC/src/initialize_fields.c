#include "fdtd.h"

void initializeE(Grid* grid, double sigma_e, double dt)
{
    int nx, ny, nz;
    double dx, dy, dz;
    DEFINE_GRID_CONSTANTS;
    double epsilon = eps0;
    for(int i = 0; i < nx-1; i++)
    {
        for(int j = 0; j < ny; j++)
        {
            for(int k = 0; k < nz; k++)
            {
                CExSelf(i,j,k) = (1 - sigma_e*dt/(2*epsilon))/(1 + sigma_e*dt/(2*epsilon));
                CExhy(i,j,k) = (dt/(epsilon*dz))*(1/(1+sigma_e*dt/(2*epsilon)));
                CExhz(i,j,k) = (dt/(epsilon*dy))*(1/(1+sigma_e*dt/(2*epsilon)));
            }
        }
    }
    printf("%p\n", grid->eyself);
    for(int i = 0; i < nx; i++)
    {
        for(int j = 0; j < ny-1; j++)
        {
            for(int k = 0; k < nz; k++)
            {
                //printf("%p\n", &(ACCESS3D(grid->eyself, i,j,k,ny-1, nz)));
                CEySelf(i,j,k) = (1 - sigma_e*dt/(2*epsilon))/(1 + sigma_e*dt/(2*epsilon));
                //getchar();
                CEyhx(i,j,k) = (dt/(epsilon*dz))*(1/(1+sigma_e*dt/(2*epsilon)));
                CEyhz(i,j,k) = (dt/(epsilon*dx))*(1/(1+sigma_e*dt/(2*epsilon)));
            }
        }
    }
    for(int i = 0; i < nx; i++)
    {
        for(int j = 0; j < ny; j++)
        {
            for(int k = 0; k < nz-1; k++)
            {
                CEzSelf(i,j,k) = (1 - sigma_e*dt/(2*epsilon))/(1 + sigma_e*dt/(2*epsilon));
                CEzhy(i,j,k) = (dt/(epsilon*dx))*(1/(1+sigma_e*dt/(2*epsilon)));
                CEzhx(i,j,k) = (dt/(epsilon*dy))*(1/(1+sigma_e*dt/(2*epsilon)));
                
            }
        }
    }
}
void initializeH(Grid* grid, double sigma_m, double dt)
{
        int nx, ny, nz;
    double dx, dy, dz;
    DEFINE_GRID_CONSTANTS;
    double mu = mu0;
    for(int i = 0; i < nx; i++)
    {
        for(int j = 0; j < ny-1; j++)
        {
            for(int k = 0; k < nz-1; k++)
            {
                CHxSelf(i,j,k) = (1 - sigma_m*dt/(2*mu))/(1 + sigma_m*dt/(2*mu));
                CHxey(i,j,k) = (dt/(mu*dz))*(1/(1+sigma_m*dt/(2*mu)));
                CHxez(i,j,k) = (dt/(mu*dy))*(1/(1+sigma_m*dt/(2*mu)));
            }
        }
    }
    for(int i = 0; i < nx-1; i++)
    {
        for(int j = 0; j < ny; j++)
        {
            for(int k = 0; k < nz-1; k++)
            {
                CHySelf(i,j,k) = (1 - sigma_m*dt/(2*mu))/(1 + sigma_m*dt/(2*mu));
                CHyex(i,j,k) = (dt/(mu*dz))*(1/(1+sigma_m*dt/(2*mu)));
                CHyez(i,j,k) = (dt/(mu*dx))*(1/(1+sigma_m*dt/(2*mu)));
            }
        }
    }
    for(int i = 0; i < nx-1; i++)
    {
        for(int j = 0; j < ny-1; j++)
        {
            for(int k = 0; k < nz; k++)
            {
                CHzSelf(i,j,k) = (1 - sigma_m*dt/(2*mu))/(1 + sigma_m*dt/(2*mu));
                CHzey(i,j,k) = (dt/(mu*dx))*(1/(1+sigma_m*dt/(2*mu)));
                CHzex(i,j,k) = (dt/(mu*dy))*(1/(1+sigma_m*dt/(2*mu)));
            }
        }
    }
}