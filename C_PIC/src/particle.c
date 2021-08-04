#include "source.h"

static void get_discrete_position(Grid* grid, Particle* p, int* position)
{
    VOIDIFNULL(p);
    VOIDIFNULL(position);
    int nx, ny, nz;
    double dx, dy, dz;
    DEFINE_GRID_CONSTANTS;
    position[0] = (p->x[0]-XMIN)/dx;
    position[1] = (p->x[1]-YMIN)/dy;
    position[2] = (p->x[2]-ZMIN)/dz;
}
void push(Grid* grid, double dt)
{
    int np = NP;
    int i;
    for(i = 0; i < np; i++)
    {
        Particle* p = grid->particles + i;
        double* v = p->v;
        int position[3];
        get_discrete_position(grid, p, position);
        double applied_e[3], applied_b[3];
        interpolate_E(grid, position[0], position[1], position[2], applied_e);
        interpolate_B(grid, position[0], position[1], position[2], applied_b);
        double vnorm = sqrt(pow(v[0],2) + pow(v[1],2) + pow(v[2],2))/M_I;

        double gamma_mi = 1/sqrt(1-pow(vnorm/c,2));
        double u[3] = {gamma_mi*v[0], gamma_mi*v[1], gamma_mi*v[2]};
        double u_prime[3];
        double vcrossB[3] = CROSS(v, applied_b);
        double const1 = Q_I*dt/M_I;
        u_prime[0] = u[0] + const1*(applied_e[0] + vcrossB[0]/2);
        u_prime[1] = u[1] + const1*(applied_e[1] + vcrossB[1]/2);
        u_prime[2] = u[2] + const1*(applied_e[2] + vcrossB[2]/2);
        double gamma_prime = sqrt(1 + DOT(u_prime,u_prime)/pow(c, 2));
        double tau[3];
        tau[0] = const1/2*applied_b[0];
        tau[1] = const1/2*applied_b[1];
        tau[2] = const1/2*applied_b[2];
        double w[3] = DOT(u_prime, tau);
        double taunormsq = DOT(tau, tau);
        double sigma = (pow(gamma_prime, 2)-taunormsq)/2;
        double next_gamma = sqrt(sigma + sqrt(pow(sigma,2) + (taunormsq + DOT(w, w))));
        double t_vec[3] = {tau[0]/next_gamma, tau[1]/next_gamma, tau[2]/next_gamma};
        double denom = (1+DOT(t_vec,t_vec))*next_gamma;
        double u_dot_t = DOT(u_prime, t_vec);
        double u_cross_t[3] = CROSS(u_prime, t_vec);
        double v_next[3] = {(u[0] + u_dot_t*t_vec[0]+u_cross_t[0])/denom,\
                            (u[1] + u_dot_t*t_vec[1]+u_cross_t[1])/denom,\
                            (u[2] + u_dot_t*t_vec[2]+u_cross_t[2])/denom};
        p->v[0] = v_next[0];
        p->v[1] = v_next[1];
        p->v[2] = v_next[2];
    }
    
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
        double* position = p->x;
        double* velocity = p->v;
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