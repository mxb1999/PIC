#include "simulation.hpp"
__global__
void updateE3D(Grid grid, double dt, double sigma_e, double epsilon)
{
    int ix = blockDim.x*blockIdx.x + threadIdx.x;
    int iy = blockDim.y*blockIdx.y + threadIdx.y;
    int iz = blockDim.z*blockIdx.z + threadIdx.z;
    int nx = grid.nx, ny = grid.ny, nz = grid.nz;
    if(ix >= nx - 1 || iy >= ny - 1 || iz >= nz - 1) {
        return;
    }

    int lastx = ix == nx-2;
    int lasty = iy == ny-2;
    int lastz = iz == nz-2;
    int xind = ix + lastx;
    int yind = iy + lasty;
    int zind = iz + lastz;
    double mult = dt/eps0;
    double mult2 = 4/(grid.dx*grid.dx);
    //assume regular spaceing
    field_t eSelf = (1 - sigma_e*dt/(2*epsilon))/(1 + sigma_e*dt/(2*epsilon));
    field_t eother = (dt/(epsilon*grid.dx))*(1/(1+sigma_e*dt/(2*epsilon)));
    Ex_cu(ix, iy, iz) += (Jx_cu(ix + 1, iy, iz) + Jx_cu(ix, iy, iz))/2*mult + (RHO_cu(ix + 1, iy, iz) + RHO_cu(ix, iy, iz))/2*mult2;
    Ey_cu(ix, iy, iz) += (Jy_cu(ix, iy + 1, iz) + Jy_cu(ix, iy, iz))/2*mult + (RHO_cu(ix, iy + 1, iz) + RHO_cu(ix, iy, iz))/2*mult2;
    Ez_cu(ix, iy, iz) += (Jz_cu(ix, iy, iz + 1) + Jz_cu(ix, iy, iz))/2*mult + (RHO_cu(ix, iy, iz + 1) + RHO_cu(ix, iy, iz))/2*mult2;
    if(lastx) {
        Ey_cu(xind, iy, iz) += (Jy_cu(xind, iy + 1, iz) + Jy_cu(xind, iy, iz))/2*mult + (RHO_cu(xind, iy + 1, iz) + RHO_cu(xind, iy, iz))/2*mult2;
        Ez_cu(xind, iy, iz) += (Jz_cu(xind, iy, iz + 1) + Jz_cu(xind, iy, iz))/2*mult + (RHO_cu(xind, iy, iz + 1) + RHO_cu(xind, iy, iz))/2*mult2;
    }
    if(lasty) {
        Ex_cu(ix, yind, iz) += (Jy_cu(ix + 1, yind, iz) + Jy_cu(ix, yind, iz))/2*mult + (RHO_cu(ix + 1, yind, iz) + RHO_cu(ix, yind, iz))/2*mult2;
        Ez_cu(ix, yind, iz) += (Jz_cu(ix, yind, iz + 1) + Jz_cu(ix, yind, iz))/2*mult + (RHO_cu(ix, yind, iz + 1) + RHO_cu(ix, yind, iz))/2*mult2;
    }
    if(lastz) {
        Ex_cu(ix, iy, zind) += (Jx_cu(ix + 1, iy, zind) + Jx_cu(ix, iy, zind))/2*mult + (RHO_cu(ix + 1, iy, zind) + RHO_cu(ix, iy, zind))/2*mult2;
        Ey_cu(ix, iy, zind) += (Jy_cu(ix, iy + 1, zind) + Jy_cu(ix, iy, zind))/2*mult + (RHO_cu(ix, iy + 1, zind) + RHO_cu(ix, iy, zind))/2*mult2;
    }
    if(iy > 0 && iz > 0)
    {
        Ex_cu(ix, iy, iz) = Ex_cu(ix, iy, iz) + eother*((Hz_cu(ix,iy,iz) - Hz_cu(ix,iy-1,iz)) - (Hy_cu(ix, iy, iz) - Hy_cu(ix, iy, iz-1)));
    }
    if(ix > 0 && iz > 0)
    {
        Ey_cu(ix, iy, iz) = Ey_cu(ix, iy, iz) + eother*((Hx_cu(ix,iy,iz) - Hx_cu(ix,iy,iz-1)) - (Hz_cu(ix,iy,iz) - Hz_cu(ix-1,iy,iz)));
    }
    if(iy > 0 && ix > 0)
    {
        Ez_cu(ix, iy, iz) = Ez_cu(ix, iy, iz) + eother*((Hy_cu(ix,iy,iz)-Hy_cu(ix-1,iy,iz)) - (Hx_cu(ix,iy,iz) - Hx_cu(ix,iy-1,iz)));
    }
    /*if(lastx) {
        Ey_cu(nx-1, iy, iz) = Ey_cu(nx-1, iy, iz) + eother*((Hx_cu(nx-1,iy,iz) - Hx_cu(nx-1,iy,iz-1)) - (Hz_cu(nx-1,iy,iz) - Hz_cu(nx-2,iy,iz)));
        Ez_cu(nx-1, iy, iz) = Ez_cu(nx-1, iy, iz) + eother*((Hy_cu(nx-1,iy,iz)-Hy_cu(nx-1-1,iy,iz)) - (Hx_cu(nx-1,iy,iz) - Hx_cu(nx-1,iy-1,iz)));
        if(lasty) {
            Ez_cu(nx-1, ny-1, iz) = Ez_cu(nx-1, ny-1, iz) + eother*((Hy_cu(nx-1,ny-1,iz)-Hy_cu(nx-1-1,ny-1,iz)) - (Hx_cu(nx-1,ny-1,iz) - Hx_cu(nx-1,ny-1-1,iz)));
        }
        if(lastz) {
            Ey_cu(nx-1, iy, nz-1) = Ey_cu(nx-1, iy, nz-1) + eother*((Hx_cu(nx-1,iy,nz-1) - Hx_cu(nx-1,iy,nz-2)) - (Hz_cu(nx-1,iy,nz-1) - Hz_cu(nx-1-1,iy,nz-1)));
        }
    }else if(lasty) {
        Ez_cu(ix, ny-1, iz) = Ez_cu(ix, ny-1, iz) + eother*((Hy_cu(ix,ny-1,iz)-Hy_cu(ix-1,ny-1,iz)) - (Hx_cu(ix,ny-1,iz) - Hx_cu(ix,ny-1-1,iz)));
        Ex_cu(ix, ny-1, iz) = Ex_cu(ix, ny-1, iz) + eother*((Hz_cu(ix,ny-1,iz) - Hz_cu(ix,ny-1-1,iz)) - (Hy_cu(ix, ny-1, iz) - Hy_cu(ix, ny-1, iz-1)));
        if(lastz) {
            Ex_cu(ix, ny-1, nz-1) = Ex_cu(ix, ny-1, nz-1) + eother*((Hz_cu(ix,ny-1,nz-1) - Hz_cu(ix,ny-1-1,nz-1)) - (Hy_cu(ix, ny-1, nz-1) - Hy_cu(ix, ny-1, nz-1-1)));
        }
    }else if(lastz) {
        Ey_cu(ix, iy, nz-1) = Ey_cu(ix, iy, nz-1) + eother*((Hx_cu(ix,iy,nz-1) - Hx_cu(ix,iy,nz-1-1)) - (Hz_cu(ix,iy,nz-1) - Hz_cu(ix-1,iy,nz-1)));
        Ex_cu(ix, iy, nz-1) = Ex_cu(ix, iy, nz-1) + eother*((Hz_cu(ix,iy,nz-1) - Hz_cu(ix,iy-1,nz-1)) - (Hy_cu(ix, iy, nz-1) - Hy_cu(ix, iy, nz-1-1)));

    }*/
}
__global__
void updateH3D(Grid grid, double dt, double sigma_m, double mu)
{
    int ix = blockDim.x*blockIdx.x + threadIdx.x;
    int iy = blockDim.y*blockIdx.y + threadIdx.y;
    int iz = blockDim.z*blockIdx.z + threadIdx.z;
    int nx = grid.nx, ny = grid.ny, nz = grid.nz;
    if(ix >= nx - 1 || iy >= ny - 1 || iz >= nz - 1) {
        return;
    }
    field_t hSelf = (1 - sigma_m*dt/(2*mu))/(1 + sigma_m*dt/(2*mu));
    field_t hother = (dt/(mu*grid.dx))*(1/(1+sigma_m*dt/(2*mu)));
    int lastx = ix == nx-2;
    int lasty = iy == ny-2;
    int lastz = iz == nz-2;
    Hx_cu(ix, iy, iz) = hSelf*Hx_cu(ix, iy, iz) + hother*((Ey_cu(ix, iy, iz+1) - Ey_cu(ix, iy, iz)) - (Ez_cu(ix, iy+1, iz) - Ez_cu(ix, iy, iz)));
    Hy_cu(ix, iy, iz) = hSelf*Hy_cu(ix, iy, iz) + hother*((Ez_cu(ix+1, iy, iz) - Ez_cu(ix, iy, iz)) - (Ex_cu(ix, iy, iz+1) - Ex_cu(ix, iy, iz)));
    Hz_cu(ix, iy, iz) = hSelf*Hz_cu(ix, iy, iz) + hother*((Ex_cu(ix, iy+1, iz) - Ex_cu(ix, iy, iz)) - (Ey_cu(ix+1, iy, iz) - Ey_cu(ix, iy ,iz)));

    //conditionally update
    if(lastx) {
        Hx_cu(nx-1, iy, iz) = hSelf*Hx_cu(nx-1, iy, iz) + hother*((Ey_cu(nx-1, iy, iz+1) - Ey_cu(nx-1, iy, iz)) -\
            (Ez_cu(nx-1, iy+1, iz) - Ez_cu(nx-1, iy, iz)));
    }
    if(lasty) {
        Hy_cu(ix, ny-1, iz) = hSelf*Hy_cu(ix, ny-1, iz) + hother*((Ez_cu(ix+1, ny-1, iz) - Ez_cu(ix, ny-1, iz)) - \
            (Ex_cu(ix, ny-1, iz+1) - Ex_cu(ix, ny-1, iz)));
    }
    if(lastz) {
        Hz_cu(ix, iy, nz-1) = hSelf*Hz_cu(ix, iy, nz-1) + hother*((Ex_cu(ix, iy+1, nz-1) - Ex_cu(ix, iy, nz-1)) - \
            (Ey_cu(ix+1, iy, nz-1) - Ey_cu(ix, iy ,nz-1)));
    }
    /*
    Hx_cu(nx-1, iy, iz) = lastx*(hSelf*Hx_cu(nx-1, iy, iz) + hother*((Ey_cu(nx-1, iy, iz+1) - Ey_cu(nx-1, iy, iz)) -\
            (Ez_cu(nx-1, iy+1, iz) - Ez_cu(nx-1, iy, iz))));
    Hy_cu(ix, ny-1, iz) = lasty*(hSelf*Hy_cu(ix, ny-1, iz) + hother*((Ez_cu(ix+1, ny-1, iz) - Ez_cu(ix, ny-1, iz)) - \
            (Ex_cu(ix, ny-1, iz+1) - Ex_cu(ix, ny-1, iz))));
    Hz_cu(ix, iy, nz-1) = lastz*(hSelf*Hz_cu(ix, iy, nz-1) + hother*((Ex_cu(ix, iy+1, nz-1) - Ex_cu(ix, iy, nz-1)) - \
            (Ey_cu(ix+1, iy, nz-1) - Ey_cu(ix, iy ,nz-1))));
    */
    /*
    for(int iy = 0; iy < ny-1; iy++)
    {
        for(int iz = 0; iz < nz-1; iz++)
        {
            Hx(nx-1, iy, iz) = hSelf*Hx(nx-1, iy, iz) + hother*((Ey_cu(nx-1, iy, iz+1) - Ey_cu(nx-1, iy, iz)) - (Ez_cu(nx-1, iy+1, iz) - Ez_cu(nx-1, iy, iz)));
        }
    }
    for(int ix = 0; ix < nx-1; ix++)
    {
        for(int iy = 0; iy < ny-1; iy++)
        {
            Hz(ix, iy, nz-1) = hSelf*Hz(ix, iy, nz-1) + hother*((Ex_cu(ix, iy+1, nz-1) - Ex_cu(ix, iy, nz-1)) - (Ey_cu(ix+1, iy, nz-1) - Ey_cu(ix, iy ,nz-1)));
        }
    }
    for(int ix = 0; ix < nx-1; ix++)
    {
        for(int iz = 0; iz < nz-1; iz++)
        {
            Hy(ix, ny-1, iz) = hSelf*Hy(ix, ny-1, iz) + hother*((Ez_cu(ix+1, ny-1, iz) - Ez_cu(ix, ny-1, iz)) - (Ex_cu(ix, ny-1, iz+1) - Ex_cu(ix, ny-1, iz)));
        }
    }*/
}
#define THREADS 8
void solve_fields_gpu(Grid* grid, double dt, double sigma_e, double sigma_m, double epsilon, double mu){
    int numblocks = NX/THREADS+1;
    dim3 blocks(numblocks, numblocks, numblocks);
    dim3 threads(THREADS, THREADS, THREADS);
    updateH3D<<<blocks, threads>>>(*grid, dt, sigma_m, mu);
    cudaError_t stat = cudaGetLastError();
    if(stat != cudaSuccess) {
        printf("%d %s\n", __LINE__, cudaGetErrorString(stat));
        getchar();
    }
    updateE3D<<<blocks, threads>>>(*grid, dt, sigma_e, epsilon);
    stat = cudaGetLastError();

    cudaDeviceSynchronize();
    memset(grid->rho, 0, sizeof(double)*NX*NY*NZ);
    memset(grid->jx, 0, sizeof(double)*NX*NY*NZ);
    memset(grid->jy, 0, sizeof(double)*NX*NY*NZ);
    memset(grid->jz, 0, sizeof(double)*NX*NY*NZ);
    if(stat != cudaSuccess) {
        printf("%d %s\n",__LINE__, cudaGetErrorString(stat));
        getchar();
    }


}