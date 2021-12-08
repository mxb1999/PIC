#include "push.h"
//#define LOG
void interpolate_E(Grid* grid, int ix, int iy, int iz, part_t* position, double* target)
{
    int nx, ny, nz;
    double dx, dy, dz;
    DEFINE_GRID_CONSTANTS;
    //assume energy conserving
    part_t xdiff = (position[0]- ix*dx)/dx;
    part_t ydiff = (position[1]-iy*dy)/dy;
    part_t zdiff = (position[2] - iz*dz)/dz;
    int not_lastx = ix != nx-1;
    int not_lasty = iy != ny-1;
    int not_lastz = iz != nz-1;
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
    target[0] += SE[0*3+0]*Ex(ix- !not_lastx, iy, iz);
    target[0] += SE[2*3+0]*Ex(ix- !not_lastx, iy+not_lasty, iz)*not_lasty;
    target[0] += SE[4*3+0]*Ex(ix- !not_lastx, iy, iz+not_lastz)*not_lastz;
    target[0] += SE[6*3+0]*Ex(ix- !not_lastx, iy+not_lasty, iz+not_lastz)*not_lasty*not_lastz;
    target[1] += SE[0*3+1]*Ey(ix, iy- !not_lasty, iz);
    target[1] += SE[1*3+1]*Ey(ix+not_lastx, iy- !not_lasty, iz)*not_lastx;
    target[1] += SE[4*3+1]*Ey(ix, iy- !not_lasty, iz+not_lastz)*not_lastz;
    target[1] += SE[5*3+1]*Ey(ix+not_lastx, iy- !not_lasty, iz+not_lastz)*not_lastx*not_lastz;
    target[2] += SE[0*3+2]*Ez(ix, iy, iz-!not_lastz);
    target[2] += SE[1*3+2]*Ez(ix+not_lastx, iy, iz-!not_lastz)*not_lastx;
    target[2] += SE[2*3+2]*Ez(ix, iy+not_lasty, iz-!not_lastz)*not_lasty;
    target[2] += SE[3*3+2]*Ez(ix+not_lastx, iy+not_lasty, iz-!not_lastz)*not_lastx*not_lasty;
};
void interpolate_B(Grid* grid, int ix, int iy, int iz, double* position, double* target)
{
    int nx, ny, nz;
    double dx, dy, dz;
    DEFINE_GRID_CONSTANTS;
    int not_lastx = ix != nx-1;
    int not_lasty = iy != ny-1;
    int not_lastz = iz != nz-1;
    //assume energy conserving
    part_t xdiff = (position[0]- ix*dx)/dx;
    part_t ydiff = (position[1]-iy*dy)/dy;
    part_t zdiff = (position[2] - iz*dz)/dz;
    double a[] = {1-xdiff, xdiff, 1-ydiff, ydiff, 1-zdiff, zdiff};
    double SB[] = {a[2*0]  , a[2*1]  , a[2*2],
                   a[2*0+1], 0       , 0,
                   0       , a[2*1+1], 0,
                   0       , 0       , 0,
                   0       , 0       , 0,
                   0       , 0       , a[2*2+1],
                   0       , 0       , 0,
                   0       , 0       , 0};
    double mu_loc = mu0;//MU(ix, iy, iz);
    //printf("Target %e %e %e\n", target[0], target[1], target[2]);
    //add all B components
    target[0] += SB[0*3+0]*Hx(ix, iy-!not_lasty, iz-!not_lastz)*mu_loc;
    target[0] += SB[1*3+0]*Hx(ix+not_lastx, iy-!not_lasty, iz-!not_lastz)*mu_loc*not_lastx;
    target[1] += SB[0*3+1]*Hy(ix-!not_lastx, iy, iz-!not_lastz)*mu_loc;
    target[1] += SB[2*3+1]*Hy(ix-!not_lastx, iy+not_lasty, iz-!not_lastz)*mu_loc*not_lasty;
    target[2] += SB[0*3+2]*Hz(ix-!not_lastx, iy-!not_lasty, iz)*mu_loc;
    target[2] += SB[5*3+2]*Hz(ix-!not_lastx, iy-!not_lasty, iz+not_lastz)*mu_loc*not_lastz;
};
void particle_push(Grid* grid, double step, part_t* logger, int stepnum) {
    //update each particle according to the fields at its current position
    //will want to seperate the field gather and the particle push steps later, for now focus on basic procedure
    int p_index;//base particle
    int numparticles = NP;
    part_t m = M_I;
    space_t dx, dy, dz;
    dx = DX;
    dy = DY;
    dz = DZ;
    int nx = NX, ny = NY, nz = NZ;
    space_t xmin = XMIN, ymin = YMIN, zmin = ZMIN;
    part_t q = Q_I;
    Particle* p_arr = PARTICLES;
    for(p_index = 0; p_index < numparticles; p_index++) {

        Particle* p = &p_arr[p_index];
        part_t px, py, pz, x, y, z;
        px = p->px;
        py = p->py;
        pz = p->pz;
        x = p->x;
        y = p->y;
        z = p->z;
        double next_x = x + (px/m)*step;
        double next_y = y + (py/m)*step;
        double next_z = z + (pz/m)*step;
        int condx_max, condx_min, condy_max, condy_min, condz_max, condz_min;
        condx_max = (next_x >= grid->xlims[1]);
        condx_min = (next_x <= grid->xlims[0]);
        condy_max = (next_y >= grid->ylims[1]);
        condy_min = (next_y <= grid->ylims[0]);
        condz_max = (next_z >= grid->zlims[1]);
        condz_min = (next_z <= grid->zlims[0]);
        x = next_x*(!condx_max && !condx_min) + condx_max*(grid->xlims[0] + 1e-10) + condx_min*(grid->xlims[1] - 1e-10);
        y = next_y*(!condy_max && !condy_min) + condy_max*(grid->ylims[0] + 1e-10) + condy_min*(grid->ylims[1] - 1e-10);
        z = next_z*(!condz_max && !condz_min) + condz_max*(grid->zlims[0] + 1e-10) + condz_min*(grid->zlims[1] - 1e-10);
        p->x = x;
        p->y = y;
        p->z = z;
        field_t elocal[] = {0.0, 0.0, 0.0};
        field_t blocal[] = {0.0, 0.0, 0.0};
        part_t pos[3] = {x, y, z};

        int ix, iy, iz;
        ix = (int)((x - xmin)/dx);
        iy = (int)((y - ymin)/dy);
        iz = (int)((z - zmin)/dz);
        interpolate_E(grid, ix, iy, iz, pos, elocal);
        interpolate_B(grid, ix, iy, iz, pos, blocal);
        part_t qconst = q*step/2;
        part_t p_temp[] = {
            px+elocal[0]*qconst,
            py+elocal[1]*qconst,
            pz+elocal[2]*qconst
        };
        part_t t_vec[] = {
            blocal[0]*qconst,
            blocal[1]*qconst,
            blocal[2]*qconst
        };

        part_t p_prime[] = {
            p_temp[0] + (p_temp[1]*t_vec[2] - p_temp[2]*t_vec[1]),
            p_temp[1] + (p_temp[2]*t_vec[0] - p_temp[0]*t_vec[2]),
            p_temp[2] + (p_temp[0]*t_vec[1] - p_temp[1]*t_vec[0])
        };
        part_t tmag = sqrt(t_vec[0]*t_vec[0] + t_vec[1]*t_vec[1] + t_vec[2]*t_vec[2]);
        part_t tconst = 2/(tmag*tmag + 1);
        p_temp[0] += tconst*(p_prime[1]*t_vec[2] - p_prime[2]*t_vec[1]);
        p_temp[1] += tconst*(p_prime[2]*t_vec[0] - p_prime[0]*t_vec[2]);
        p_temp[2] += tconst*(p_prime[0]*t_vec[1] - p_prime[1]*t_vec[0]);
        //particle push with Boris step
        p->px = p_temp[0] + elocal[0]*qconst;
        p->py = p_temp[1] + elocal[1]*qconst;
        p->pz = p_temp[2] + elocal[2]*qconst;
        #ifdef LOG
            logger[(stepnum*numparticles + p_index)*3] = x;
            logger[(stepnum*numparticles + p_index)*3 + 1] = y;
            logger[(stepnum*numparticles + p_index)*3 + 2] = z;
        #endif

    }
}