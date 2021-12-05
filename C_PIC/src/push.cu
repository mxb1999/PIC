#include "push.hu"
#include "simulation.hpp"
inline __device__
void interpolate_E(Grid* grid, int ix, int iy, int iz, part_t position[3], field_t* target)
{
    int nx, ny, nz;
    space_t dx, dy, dz;
    DEFINE_GRID_CONSTANTS;
    //assume energy conserving
    space_t xdiff = (target[0]- ix*dx)/dx;
    space_t ydiff = (target[1]-iy*dy)/dy;
    space_t zdiff = (target[2] - iz*dz)/dz;
    space_t a[] = {1-xdiff, xdiff, 1-ydiff, ydiff, 1-zdiff, zdiff};
    space_t SE[] = {a[2*1] * a[2*2],     a[2*0] * a[2*2]    ,   a[2*0] * a[2*1],
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
inline __device__
void interpolate_B(Grid* grid, int ix, int iy, int iz, part_t position[3], field_t* target)
{
    int nx, ny, nz;
    space_t dx, dy, dz;
    DEFINE_GRID_CONSTANTS;
    //assume energy conserving
    space_t xdiff = (target[0]- ix*dx)/dx;
    space_t ydiff = (target[1]-iy*dy)/dy;
    space_t zdiff = (target[2] - iz*dz)/dz;
    space_t a[] = {1-xdiff, xdiff, 1-ydiff, ydiff, 1-zdiff, zdiff};
    space_t SB[] = {a[2*0]  , a[2*1]  , a[2*2],
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
inline __device__
void get_field_cu(Particle* p, Grid* grid, field_t* efield[3], field_t* bfield[3], part_t x, part_t y, part_t z, space_t dx, space_t dy, space_t dz, int nx, int ny, int nz, space_t xmin, space_t ymin, space_t zmin)
{
    //keep it simple for now, just populate with the field recorded in nearest node
    int ix, iy, iz;
    int offset;
    ix = (int)(x - xmin)/dx;
    iy = (int)(y - ymin)/dy;
    iz = (int)(z - zmin)/dz;
    offset = ((ix*ny + iy)*(nz) + iz)*3;
    //want to interpolate with 2nd order spline
    //W(x) = 1/dx^3(-1*x^2 + 3/4*dx^2) if 0<=x<=dx/2
    //W(x) = 1/dx^3*1/8(2x - 3*dx)^2 if dx/2<=x<=3/2dx

    field_t* blocalx = bfield + offset;
    field_t* elocaly = efield + offset;
    bfield[0] = blocal[0];
    bfield[1] = blocal[1];
    bfield[2] = blocal[2];

    efield[0] = elocal[0];
    efield[1] = elocal[1];
    efield[2] = elocal[2];
};
void get_field(Particle* p, Grid* grid, field_t* efield, field_t* bfield, field_t* e, field_t* b)
{
    //keep it simple for now, just populate with the field recorded in nearest node
    int ix, iy, iz;
    int offset;
    space_t dx, dy, dz;
    dx = grid->dx;
    dy = grid->dy;
    dz = grid->dz;
    ix = (int)(p->x - grid->xlims[0])/dx;
    iy = (int)(p->y - grid->ylims[0])/dy;
    iz = (int)(p->z - grid->zlims[0])/dz;
    offset = ((ix*grid->ny + iy)*(grid->nz) + iz)*3;
    field_t* blocal = bfield + offset;
    field_t* elocal = efield + offset;
    b[0] = blocal[0];
    b[1] = blocal[1];
    b[2] = blocal[2];

    e[0] = elocal[0];
    e[1] = elocal[1];
    e[2] = elocal[2];
};

__global__
void particle_push_cu(Particle* p_arr, field_t* efield, field_t* bfield, Grid* grid, double step, int ppt, part_t* logger, int stepnum)
{
    //update each particle according to the fields at its current position
    //will want to seperate the field gather and the particle push steps later, for now focus on basic procedure
    int p_index = (blockIdx.x*blockDim.x + threadIdx.x)*ppt;//base particle
    int maxnum = NP;
    if(p_index >= maxnum)
    {
        return;
    }
    maxnum = (p_index + ppt >= maxnum) ? maxnum - 1 : p_index + ppt;
    part_t m = grid->mass_p;
    space_t dx, dy, dz;
    dx = DX;
    dy = DY;
    dz = DZ;
    int nx = NX, ny = NY, nz = NZ;
    space_t xmin = XMIN, ymin = YMIN, zmin = ZMIN;
    part_t q = Q_I;
    for(int i = p_index; i < maxnum; i++)//for each particle assigned to this thread
    {
        Particle* p = p_arr + i;
        part_t px, py, pz, x, y, z;
        px = p->px;
        py = p->py;
        pz = p->pz;
        x = p->x;
        y = p->y;
        z = p->z;
        double next_x = x + px*step/m;
        double next_y = y + py*step/m;
        double next_z = z + pz*step/m;
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
        field_t elocal[3];
        field_t blocal[3];
        part_t pos[3] = {x, y, z};
        int ix, iy, iz;
        ix = (int)(x - xmin)/dx;
        iy = (int)(y - ymin)/dy;
        iz = (int)(z - zmin)/dz;
        interpolate_E(grid, ix, iy, iz, pos, elocal);
        interpolate_B(grid, ix, iy, iz, pos, blocal);
        part_t qconst = q*step/2;
        part_t p_temp[3] = {
            px+elocal[0]*qconst,
            py+elocal[1]*qconst,
            pz+elocal[2]*qconst
        };
        part_t t_vec[3] = {
            blocal[0]*qconst,
            blocal[1]*qconst,
            blocal[2]*qconst
        };
        part_t p_prime[3] = {
            p_temp[0] + (p_temp[1]*t_vec[2] - p_temp[2]*t_vec[1]),
            p_temp[1] + (p_temp[2]*t_vec[0] - p_temp[0]*t_vec[2]),
            p_temp[2] + (p_temp[0]*t_vec[1] - p_temp[1]*t_vec[0])
        };
        part_t tmag = sqrt(t_vec[0]*t_vec[0] + t_vec[1]*t_vec[1] + t_vec[2]*t_vec[2]);
        part_t tconst = 2/(tmag*tmag + 1);
      /*  t_vec[0] *= tconst;
        t_vec[1] *= tconst;
        t_vec[2] *= tconst;*/
        p_temp[0] += tconst*(p_prime[1]*t_vec[2] - p_prime[2]*t_vec[1]);
        p_temp[1] += tconst*(p_prime[2]*t_vec[0] - p_prime[0]*t_vec[2]);
        p_temp[2] += tconst*(p_prime[0]*t_vec[1] - p_prime[1]*t_vec[0]);
        /*
        part_t p_plus[3] = {
            p_minus[0] + (p_prime[1]*t_vec[2] - p_prime[2]*t_vec[1]),
            p_minus[1] + (p_prime[2]*t_vec[0] - p_prime[0]*t_vec[2]),
            p_minus[2] + (p_prime[0]*t_vec[1] - p_prime[1]*t_vec[0])
        };*/
        //particle push with Boris step
        p->px = p_temp[0] + elocal[0]*qconst;
        p->py = p_temp[1] + elocal[1]*qconst;
        p->pz = p_temp[2] + elocal[2]*qconst;
        #ifdef LOG
            logger[((stepnum - 1)*grid->num_particles + i)*3] = p->x;
            logger[((stepnum - 1)*grid->num_particles + i)*3 + 1] = p->y;
            logger[((stepnum - 1)*grid->num_particles + i)*3 + 2] = p->z;
        #endif

    }
};//main pushing kernel
void particle_push(Grid* grid, double step, part_t* logger, int stepnum)
{
    //update each particle according to the fields at its current position
    //will want to seperate the field gather and the particle push steps later, for now focus on basic procedure
    Particle* p_arr = grid->particles;
    field_t* efield = grid->e_field;
    field_t* bfield = grid->b_field;
    part_t m = grid->mass_p;
    part_t q = grid->q_p;
    int num_p = grid->num_particles;
    for(int i = 0; i < num_p; i++)//for each particle assigned to this thread
    {
        Particle* p = p_arr + i;
        part_t px, py, pz, x, y, z;
        px = p->px;
        py = p->py;
        pz = p->pz;
        x = p->x;
        y = p->y;
        z = p->z;
        double next_x = x + px*step/m;
        double next_y = y + py*step/m;
        double next_z = z + pz*step/m;
        int condx_max, condx_min, condy_max, condy_min, condz_max, condz_min;
        condx_max = (next_x >= grid->xlims[1]);
        condx_min = (next_x <= grid->xlims[0]);
        condy_max = (next_y >= grid->ylims[1]);
        condy_min = (next_y <= grid->ylims[0]);
        condz_max = (next_z >= grid->zlims[1]);
        condz_min = (next_z <= grid->zlims[0]);
        p->x = next_x*(!condx_max && !condx_min) + condx_max*(grid->xlims[0] + 1e-10) + condx_min*(grid->xlims[1] - 1e-10);
        p->y = next_y*(!condy_max && !condy_min) + condy_max*(grid->ylims[0] + 1e-10) + condy_min*(grid->ylims[1] - 1e-10);
        p->z = next_z*(!condz_max && !condz_min) + condz_max*(grid->zlims[0] + 1e-10) + condz_min*(grid->zlims[1] - 1e-10);
        field_t elocal[3];
        field_t blocal[3];
        get_field(p, grid, efield, bfield, elocal, blocal);
        part_t qconst = q*step/2;
        part_t p_minus[3] = {
            px+elocal[0]*qconst,
            py+elocal[1]*qconst,
            pz+elocal[2]*qconst
        };
        part_t t_vec[3] = {
            blocal[0]*qconst,
            blocal[1]*qconst,
            blocal[2]*qconst
        };
        part_t p_prime[3] = {
            p_minus[0] + (p_minus[1]*t_vec[2] - p_minus[2]*t_vec[1]),
            p_minus[1] + (p_minus[2]*t_vec[0] - p_minus[0]*t_vec[2]),
            p_minus[2] + (p_minus[0]*t_vec[1] - p_minus[1]*t_vec[0])
        };
        part_t tmag = sqrt(t_vec[0]*t_vec[0] + t_vec[1]*t_vec[1] + t_vec[2]*t_vec[2]);
        part_t tconst = 2/(tmag*tmag + 1);
        t_vec[0] *= tconst;
        t_vec[1] *= tconst;
        t_vec[2] *= tconst;
        part_t p_plus[3] = {
            p_minus[0] + (p_prime[1]*t_vec[2] - p_prime[2]*t_vec[1]),
            p_minus[1] + (p_prime[2]*t_vec[0] - p_prime[0]*t_vec[2]),
            p_minus[2] + (p_prime[0]*t_vec[1] - p_prime[1]*t_vec[0])
        };
        //particle push with Boris step
        p->px = p_plus[0] + elocal[0]*qconst;
        p->py = p_plus[1] + elocal[1]*qconst;
        p->pz = p_plus[2] + elocal[2]*qconst;
        double pnorm = sqrt(px*px + py*py + pz*pz);
        px = p->px;
        py = p->py;
        pz = p->pz;
        double pnorm2 = sqrt(px*px + py*py + pz*pz);
        #ifdef LOG
            logger[((stepnum - 1)*grid->num_particles + i)*3] = p->x;
            logger[((stepnum - 1)*grid->num_particles + i)*3 + 1] = p->y;
            logger[((stepnum - 1)*grid->num_particles + i)*3 + 2] = p->z;
        #endif

    }
};//main pushing kernel
void push_stage(Particle* p, Grid* grid, double step, int ppt, part_t* logger, int blocks, int stepnum, int opt)
{
    field_t* efield, *bfield;
    cudaMalloc(efield, sizeof(field_t)*NX*NY*NZ*3);
    cudaMalloc(bfield, sizeof(field_t)*NX*NY*NZ*3);
    cudaMemset(efield, 0, )
    if(opt == 1)
    {

        particle_push_cu<<<blocks, threads_per_block>>>(p, efield, bfield, grid, step, ppt, logger, stepnum);
        cudaDeviceSynchronize();
        cudaError_t err = cudaGetLastError();
        if(err != cudaSuccess)
        {
            printf("ERROR %s\n", cudaGetErrorString(err));
        }
    }else
    {
        particle_push(grid, step, logger, stepnum);
    }
};