#include "fdtd.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "push.h"


#define SIZE 200

void test_push(int n, int p)
{
    double m = 1e-28;
    double dims[] = {-5e-5, 5e-5};
    double dimsy[] = {-5e-5, 5e-5};
    double dimsz[] = {-5e-5, 5e-5};
    int gridres = 100;
    Grid* grid = define_grid(dims, dimsy, dimsz, gridres,gridres,gridres,100, 0);
    grid->mass_p = m;
    grid->q_p = 1e-19;
    int nx = grid->nx;
    int ny = grid->ny;
    int nz = grid->nz;
    field_t b[] = {0, 0, 1e3};
    setup_grid_constb(grid, unif, unif, 1e3*m,500*m, b);
    for(int i = 0; i < 1e8; i++) {
        particle_push(grid, 1e-5, NULL, i);

    }
    //free_grid(grid);
}
void updateH2D_FirstTMz(Grid* grid)
{
    int nx = NX, ny = NY, nz = NZ;
    int ly = ny-1;
    int lx = nx-1;
    for(int i = 0; i < nx - 1; i++)
    {
        Hy(i, ly, 0) = Hy(i, ly, 0)*CHySelf(i, ly, 0) + CHyez(i,ly,0)*(Ez(i+1,ly,0) - Ez(i,ly,0));
    }
    for(int j = 0; j < ny - 1; j++)
    {
        Hx(lx, j, 0) = Hx(lx, j, 0)*CHxSelf(lx, j, 0) - CHyez(lx,j,0)*(Ez(lx,j+1,0) - Ez(lx,j,0));
    }
}
void updateH2DTMz(Grid* grid)
{
    int nx = NX, ny = NY, nz = NZ;
    for(int i = 0; i < nx - 1; i++)
    {
        for(int j = 0; j < ny - 1; j++)
        {
            Hy(i, j, 0) = Hy(i, j, 0)*CHySelf(i, j, 0) + CHyez(i,j,0)*(Ez(i+1,j,0) - Ez(i,j,0));
            Hx(i, j, 0) = Hx(i, j, 0)*CHxSelf(i, j, 0) - CHyez(i,j,0)*(Ez(i,j+1,0) - Ez(i,j,0));
        }
    }
    updateH2D_FirstTMz(grid);
}
void updateE2DTMz(Grid* grid)
{
    int nx = NX, ny = NY, nz = NZ;
    for(int i = 1; i < nx; i++)
    {
        for(int j = 1; j < ny; j++)
        {
            Ez(i, j, 0) = CEzSelf(i,j,0)*Ez(i, j, 0) + CEzhy(i, j, 0)*(Hy(i, j, 0) - Hy(i - 1, j, 0)) - CEzhx(i, j, 0)*(Hx(i, j, 0) - Hx(i, j - 1, 0));
        }
    }

}
void updateE3D(Grid* grid)
{
    int nx = NX, ny = NY, nz = NZ;
    for(int i = 0; i < nx-1; i++)
    {
        for(int j = 0; j < ny-1; j++)
        {
            for(int k = 0; k < nz-1; k++)
            {
                if(j > 0 && k > 0)
                {
                    Ex(i, j, k) = Ex(i, j, k) + CExhy(i,j,k)*((Hz(i,j,k) - Hz(i,j-1,k)) - (Hy(i, j, k) - Hy(i, j, k-1)));
                }
                if(i > 0 && k > 0)
                {
                    Ey(i, j, k) = Ey(i, j, k) + CEyhx(i,j,k)*((Hx(i,j,k) - Hx(i,j,k-1)) - (Hz(i,j,k) - Hz(i-1,j,k)));
                }
                if(j > 0 && i > 0)
                {
                    Ez(i, j, k) = Ez(i, j, k) + CEzhx(i,j,k)*((Hy(i,j,k)-Hy(i-1,j,k)) - (Hx(i,j,k) - Hx(i,j-1,k)));
                }
            }
        }
    }
}
void updateH3D(Grid* grid)
{
    int nx = NX, ny = NY, nz = NZ;
    for(int i = 0; i < nx-1; i++)
    {
        int lastx = i == nx-2;
        for(int j = 0; j < ny-1; j++)
        {
            int lasty = i == ny-2;
            for(int k = 0; k < nz-1; k++)
            {
                int lastz = i == nz-2;
                Hx(i, j, k) = CHxSelf(i, j, k)*Hx(i, j, k) + CHxey(i, j, k)*((Ey(i, j, k+1) - Ey(i, j, k)) - (Ez(i, j+1, k) - Ez(i, j, k)));
                Hy(i, j, k) = CHySelf(i, j, k)*Hy(i, j, k) + CHyex(i, j, k)*((Ez(i+1, j, k) - Ez(i, j, k)) - (Ex(i, j, k+1) - Ex(i, j, k)));
                Hz(i, j, k) = CHzSelf(i, j, k)*Hz(i, j, k) + CHzex(i, j, k)*((Ex(i, j+1, k) - Ex(i, j, k)) - (Ey(i+1, j, k) - Ey(i, j ,k)));
                if(lastx) {
                    Hx(nx-1, j, k) = CHxSelf(nx-1, j, k)*Hx(nx-1, j, k) + CHxey(nx-1, j, k)*((Ey(nx-1, j, k+1) - Ey(nx-1, j, k)) - (Ez(nx-1, j+1, k) - Ez(nx-1, j, k)));
                }
                if(lasty) {
                    Hy(i, ny-1, k) = CHySelf(i, ny-1, k)*Hy(i, ny-1, k) + CHyex(i, ny-1, k)*((Ez(i+1, ny-1, k) - Ez(i, ny-1, k)) - (Ex(i, ny-1, k+1) - Ex(i, ny-1, k)));
                }
                if(lastz) {
                    Hz(i, j, nz-1) = CHzSelf(i, j, nz-1)*Hz(i, j, nz-1) + CHzex(i, j, nz-1)*((Ex(i, j+1, nz-1) - Ex(i, j, nz-1)) - (Ey(i+1, j, nz-1) - Ey(i, j ,nz-1)));
                }

            }
        }
    }
    /*
    for(int j = 0; j < ny-1; j++)
    {
        for(int k = 0; k < nz-1; k++)
        {
            Hx(nx-1, j, k) = CHxSelf(nx-1, j, k)*Hx(nx-1, j, k) + CHxey(nx-1, j, k)*((Ey(nx-1, j, k+1) - Ey(nx-1, j, k)) - (Ez(nx-1, j+1, k) - Ez(nx-1, j, k)));
        }
    }
    for(int i = 0; i < nx-1; i++)
    {
        for(int j = 0; j < ny-1; j++)
        {
            Hz(i, j, nz-1) = CHzSelf(i, j, nz-1)*Hz(i, j, nz-1) + CHzex(i, j, nz-1)*((Ex(i, j+1, nz-1) - Ex(i, j, nz-1)) - (Ey(i+1, j, nz-1) - Ey(i, j ,nz-1)));
        }
    }
    for(int i = 0; i < nx-1; i++)
    {
        for(int k = 0; k < nz-1; k++)
        {
            Hy(i, ny-1, k) = CHySelf(i, ny-1, k)*Hy(i, ny-1, k) + CHyex(i, ny-1, k)*((Ez(i+1, ny-1, k) - Ez(i, ny-1, k)) - (Ex(i, ny-1, k+1) - Ex(i, ny-1, k)));
        }
    }*/
}



void write_data(double* arr, int* dims, char* target, int dim)
{
    //Assume a 2D array for now, keep it simple and do CSV output
    FILE* csvfile = fopen(target, "w");
    for(int i = 0; i < dims[0]; i++)
    {
        char stringrep[50];
        for(int j = 0; j < dims[1]; j++)
        {
            snprintf(stringrep, 50, "%0.10e", arr[(i*dims[1]+j)*3+dim]);
            if(j == dims[1]-1)
            {
                break;
            }
            fprintf(csvfile, "%s, ", stringrep);
        }
        fprintf(csvfile, "%s", stringrep);
        fprintf(csvfile, "\n");
    }
    fclose(csvfile);
};



int main()
{
    /*double dims[] = {-1, 1};
    double dimsy[] = {-1, 1};
    double dimsz[] = {-1, 1};
    int gridres = 10;
    double m = 1e-28;
    Grid* grid = define_grid(dims, dimsy, dimsz, gridres,gridres,gridres,100, 0);
    grid->mass_p = m;
    grid->q_p = 1e-19;
    int nx = grid->nx;
    int ny = grid->ny;
    int nz = grid->nz;
    field_t b[] = {0, 0, 1e-10};
    setup_grid_constb(grid, unif, unif, 1e3*m,500*m, b);

    double constE = 1e-3;
    double constB = 1;
    double dt = 4e-12;
    int nt = 1000;
    //initialize(grid, 0.0, 0.0, dt);
    start_loop3DFDTD(grid,  dt, nt);
    return 0;*/
    double dims[] = {-1, 1};
    double m = 1e-28;

    double dimsy[] = {-1, 1};
    double dimsz[] = {-1, 1};
    int gridres = 100;
    Grid* grid = define_grid(dims, dimsy, dimsz, gridres,gridres,gridres,1000000, 0);

    double constE = 1e-3;
    double constB = 1;
    double dt = 1e-13;
    M_I = m;
    Q_I = 1e-19;
    initialize(grid, 0.0, 0.0, &dt, unif, unif, 1e3*m,500*m);
    double time = dt*1000;
    printf("%e %e\n", time, dt);
    start_loop3DFDTD(grid, dt, time);
    return 0;
};
