#include "push.hu"
#include "grid.h"

inline __device__
void get_field_cu(Particle* p, Grid* grid, field_t* efield, field_t* bfield, field_t* e, field_t* b, part_t x, part_t y, part_t z, space_t dx, space_t dy, space_t dz, int nx, int ny, int nz, space_t xmin, space_t ymin, space_t zmin)
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
    int maxnum = grid->num_particles;
    if(p_index >= maxnum)
    {
        return;
    }
    maxnum = (p_index + ppt >= maxnum) ? maxnum - 1 : p_index + ppt;
    part_t m = grid->mass_p;
    part_t q = grid->q_p;
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
        p->x = next_x*(!condx_max && !condx_min) + condx_max*(grid->xlims[0] + 1e-10) + condx_min*(grid->xlims[1] - 1e-10);
        p->y = next_y*(!condy_max && !condy_min) + condy_max*(grid->ylims[0] + 1e-10) + condy_min*(grid->ylims[1] - 1e-10);
        p->z = next_z*(!condz_max && !condz_min) + condz_max*(grid->zlims[0] + 1e-10) + condz_min*(grid->zlims[1] - 1e-10);
        field_t elocal[3];
        field_t blocal[3];
        get_field_cu(p, grid, efield, bfield, elocal, blocal);
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
void push_stage(Particle* p, field_t* efield, field_t* bfield, Grid* grid, double step, int ppt, part_t* logger, int blocks, int stepnum, int opt)
{
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