#include "fdtd.h"
#include "push.h"
static double add_ricker(double time, double location, double dt, Grid* grid)
{
    double dx = (XMAX-XMIN)/NX;
    double cdtds = c*dt/dx;
    double ppw = 5.0;
    double arg = M_PI*((cdtds*time-location)/ppw - 1.0);
    arg = arg * arg;
    return (1.0 - 2.0*arg)*exp(-1*arg);

}
void start_loop1DFDTD(Grid* grid, double dt, int nt)
{
    FILE* snapshot =fopen("dataout.csv", "w");
    int nx, ny, nz;
    double dx, dy, dz;
    DEFINE_GRID_CONSTANTS;
    double multE = pow(c, 2)*dt/dx;
    double multB = dt/dx;
    int dims[2] = {1,nx};
    int midx = nx/2;
    double jmult = pow(dx,2)*mu0*pow(c, 2)*dt;
    double mult = 377.0;
    double sourceterm = 1e3;//current in the source wire
    double time = 0;
    double omega = 1/(2*M_PI);
    double T = 1/freq;
    int frame = 1;
    for(int t = 0; t < nt; t++)
    {
        Hy(nx-1, 0, 0) = Hy(nx-2, 0, 0);
        for(int i = 0; i < nx-1; i++)
        {
            Hy(i, 0, 0) = Hy(i, 0, 0) + (Ez(i+1, 0, 0)-Ez(i, 0, 0))/mult;
        }
        Ez(0,0,0) = Ez(1,0,0);// + mult*(hy[0]-hy[nx-1]);
        for(int i = 1; i < nx; i++)
        {
            Ez(i,0,0) = Ez(i,0,0) + mult*(Hy(i, 0, 0)-Hy(i-1, 0, 0));
        }
        Ez(50,0,0) += exp(-1*(t-30.)*(t-30.)/100.);
        time += dt;
        if (t % 10 == 0) {
            for (int q = 0; q < nx-1; q++)
                fprintf(snapshot, "%e, ", Ez(q,0,0));
            fprintf(snapshot, "%e\n", Ez(nx-1,0,0));
        }
    }
    fclose(snapshot);
}

void start_loop2DFDTD(Grid* grid, double dt, int nt)
{
    int nx, ny, nz;
    double dx, dy, dz;
    DEFINE_GRID_CONSTANTS;
    printf("Hello %e\n", c*dt/dx);
    double ratio = 1/sqrt(3);
    dt = dx*ratio/c;
    for(int time = 0; time < nt; time++)
    {
        updateH2DTMz(grid);
        updateE2DTMz(grid);
        Ez(nx/2, ny/2, 0) += add_ricker(time*dt, 0.0, dt, grid);
        //if (time % 10 == 0) {
            char filename[100];
            sprintf(filename, "output/dataout%d.csv", time/10);
            FILE* snapshot =fopen(filename, "w");
            for(int i = 0; i < ny; i++)
            {
                for(int j = 0; j < nx-1; j++)
                {
                    fprintf(snapshot, "%e, ", Ez(j,i,0));
                }
                fprintf(snapshot, "%e\n", Ez(nx-1,i,0));
            }
            fclose(snapshot);
        //}
    }


}

void start_loop3DFDTD(Grid* grid, double dt, int nt)
{
    int nx, ny, nz;
    double dx, dy, dz;
    hid_t file, set, space;
    herr_t err_h5;
    hsize_t dims[] = {(nt), grid->num_particles, 3};
    DEFINE_GRID_CONSTANTS;
    double ratio = 1/sqrt(3);
    dt = sqrt(pow(dx,2) + pow(dy,2)+pow(dz,2))/(c);
    printf("dt %e\n", dt);
    part_t* logger = (part_t*)malloc(sizeof(part_t)*grid->num_particles*3*(nt));
    for(int time = 0; time < nt; time++)
    {
        printf("%e\n", Hz(nx/2, ny/2, nz/2));
        updateH3D(grid);
        updateE3D(grid);
        //particle_push(grid, dt, logger, time);
        Hz(nx/2, ny/2, nz/2) += add_ricker(time*dt, 0.0, dt, grid);
        if (time % 10 == 0) {
            char filename[100];
            sprintf(filename, "output/dataout%d.csv", time/10);
            FILE* snapshot =fopen(filename, "w");
            for(int i = 0; i < ny; i++)
            {
                for(int j = 0; j < nz-1; j++)
                {
                    fprintf(snapshot, "%e, ", Hz((nx)/2,i,j));
                }
                fprintf(snapshot, "%e\n", Hz((nx)/2,i,nz-1));
            }
            fclose(snapshot);
        }
    }
    file = H5Fcreate("positions.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    space = H5Screate_simple(3, dims, NULL);
    set = H5Dcreate2(file, "positions", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    err_h5 = H5Dwrite(set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, logger);
    err_h5 = H5Dclose(set);
    err_h5 = H5Sclose(space);
    err_h5 = H5Fclose(file);
    free(logger);
}