#include "simulation.hpp"
#define TRILIN_INTERP(x, y, z, xa, ya, za) {\
    xd = 1-fabs((x)-(xa))/dx;\
    yd = 1-fabs((y)-(ya))/dy;\
    zd = 1-fabs((z)-(za))/dz;\
    weight = xd*yd*zd;\
};

//using trilinear interpolation to map currents
void current_deposition(Grid* grid, double dt) {
    Particle* allparticles = PARTICLES;
    double dx = DX, dy = DY, dz = DZ;
    int nx = NX, ny = NY, nz = NZ;
    double xmin = XMIN, ymin = YMIN, zmin = ZMIN;

    //DEFINE_GRID_CONSTANTS;
    int num_particles = NP;
    double m = M_I;
    double q = Q_I;
    /*for(i = 0; i < num_particles; i++) {
        Particle p = allparticles[i];
        double x = p.x, y = p.y, z = p.z;
        double vx = p.px/m, vy = p.py/m, vz = p.pz/m;
        int ix = (int)((x-xmin)/dx), iy = (int)((y-ymin)/dy), iz = (int)((z-zmin)/dz);
        int not_lastx = (ix != nx-1);
        int not_lasty = (iy != ny-1);
        int not_lastz = (iz != nz-1);
        double xg[] = {dx*ix + xmin, dx*(ix+1) + xmin};
        double yg[] = {dy*iy + ymin, dy*(iy+1) + ymin};
        double zg[] = {dz*iz + zmin, dz*(iz+1) + zmin};
        int ox = ix, oy = iy, oz = iz;
        double xd, yd, zd, weight;
        double jx_p = vx*q;
        double jy_p = vy*q;
        double jz_p = vz*q;
        //ix iy iz
        TRILIN_INTERP(x, y, z, xg[ox-ix], yg[oy-iy], zg[oz-iz]);

        Jx(ox, oy, oz) += jx_p*weight;
        Jy(ox, oy, oz) += jy_p*weight;
        Jz(ox, oy, oz) += jz_p*weight;
        RHO(ox, oy, oz) += weight*q;
        //ix+1 iy iz
        ox = ix + not_lastx;
        TRILIN_INTERP(x, y, z, xg[ox-ix], yg[oy-iy], zg[oz-iz]);
        Jx(ox, oy, oz) += jx_p*weight*(not_lastx);
        Jy(ox, oy, oz) += jy_p*weight*(not_lastx);
        Jz(ox, oy, oz) += jz_p*weight*(not_lastx);
        RHO(ox, oy, oz) += weight*q;
        //ix+1 iy+1 iz
        oy = iy + 1;
        TRILIN_INTERP(x, y, z, xg[ox-ix], yg[oy-iy], zg[oz-iz]);
        Jx(ox, oy, oz) += jx_p*weight*(not_lastx)*(not_lasty);
        Jy(ox, oy, oz) += jy_p*weight*(not_lastx)*(not_lasty);
        Jz(ox, oy, oz) += jz_p*weight*(not_lastx)*(not_lasty);
        RHO(ox, oy, oz) += weight*q;
        //ix+1 iy iz+1
        oy = iy;
        oz = iz + not_lastz;
        TRILIN_INTERP(x, y, z, xg[ox-ix], yg[oy-iy], zg[oz-iz]);
        Jx(ox, oy, oz) += jx_p*weight*(not_lastx)*(not_lastz);
        Jy(ox, oy, oz) += jy_p*weight*(not_lastx)*(not_lastz);
        Jz(ox, oy, oz) += jz_p*weight*(not_lastx)*(not_lastz);
        RHO(ox, oy, oz) += weight*q;
        //ix iy+1 iz+1
        ox = ix;
        oy = iy + not_lasty;
        TRILIN_INTERP(x, y, z, xg[ox-ix], yg[oy-iy], zg[oz-iz]);
        Jx(ox, oy, oz) += jx_p*weight*(not_lastz)*(not_lasty);
        Jy(ox, oy, oz) += jy_p*weight*(not_lastz)*(not_lasty);
        Jz(ox, oy, oz) += jz_p*weight*(not_lastz)*(not_lasty);
        RHO(ox, oy, oz) += weight*q;
        //ix+1 iy+1 iz+1
        ox = ix + not_lastx;
        TRILIN_INTERP(x, y, z, xg[ox-ix], yg[oy-iy], zg[oz-iz]);
        Jx(ox, oy, oz) += jx_p*weight*(not_lastx)*(not_lasty)*(not_lastz);
        Jy(ox, oy, oz) += jy_p*weight*(not_lastx)*(not_lasty)*(not_lastz);
        Jz(ox, oy, oz) += jz_p*weight*(not_lastx)*(not_lasty)*(not_lastz);
        RHO(ox, oy, oz) += weight*q;
    }*/
    double mult2 = 4/(dx*dx);
    #pragma omp parallel for num_threads(10)
    for(int i = 0; i < nx-1; i++) {
        int lastx = (i == nx-2);
        for(int j = 0; j < ny-1; j++) {
            int lasty = (j == ny-2);
            for(int k = 0; k < nz-1; k++) {
                int lastz = (k == nz-2);
                double eps = eps0;//EPSILON(i, j, k);
                double mult = dt/eps;
                int xind = i + lastx;
                int yind = j + lasty;
                int zind = k + lastz;

                Ex(i, j, k) += (Jx(i + 1, j, k) + Jx(i, j, k))/2*mult + (RHO(i + 1, j, k) + RHO(i, j, k))/2*mult2;
                Ey(i, j, k) += (Jy(i, j + 1, k) + Jy(i, j, k))/2*mult + (RHO(i, j + 1, k) + RHO(i, j, k))/2*mult2;
                Ez(i, j, k) += (Jz(i, j, k + 1) + Jz(i, j, k))/2*mult + (RHO(i, j, k + 1) + RHO(i, j, k))/2*mult2;
                if(lastx) {
                    Ey(xind, j, k) += (Jy(xind, j + 1, k) + Jy(xind, j, k))/2*mult + (RHO(xind, j + 1, k) + RHO(xind, j, k))/2*mult2;
                    Ez(xind, j, k) += (Jz(xind, j, k + 1) + Jz(xind, j, k))/2*mult + (RHO(xind, j, k + 1) + RHO(xind, j, k))/2*mult2;
                }
                if(lasty) {
                    Ex(i, yind, k) += (Jy(i + 1, yind, k) + Jy(i, yind, k))/2*mult + (RHO(i + 1, yind, k) + RHO(i, yind, k))/2*mult2;
                    Ez(i, yind, k) += (Jz(i, yind, k + 1) + Jz(i, yind, k))/2*mult + (RHO(i, yind, k + 1) + RHO(i, yind, k))/2*mult2;
                }
                if(lastz) {
                    Ex(i, j, zind) += (Jx(i + 1, j, zind) + Jx(i, j, zind))/2*mult + (RHO(i + 1, j, zind) + RHO(i, j, zind))/2*mult2;
                    Ey(i, j, zind) += (Jy(i, j + 1, zind) + Jy(i, j, zind))/2*mult + (RHO(i, j + 1, zind) + RHO(i, j, zind))/2*mult2;
                }

            }
        }
    }
    //getchar();


}