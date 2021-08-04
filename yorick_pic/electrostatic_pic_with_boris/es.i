// 2D Electrostatic PIC
// Flow of solar wind around a charged plate


extern EPSO, QE, den, A, n0, phi0, phi_p, Te, box

extern x_var,y_var,phi,efx,efy,den,solvertol;

#include "potsolv.i"
#include "myplots.i"

#include "palettes/pal.i"
#include "palettes/pal2.i"
#include "palettes/colorbar.i"



/* Define constants */

EPS0 = 8.85418782e-12	// permittivity of free space, Farad m^-1 = A^2 s^4 kg^-1 m^-3
QE = 1.602176662e-19	// electron charge, C, also converts eV to J (1.60218e-19)
K = 1.38064853e-23	// boltzmann constant, J K^-1
AMU = 1.66053904e-27	// atomic mass unit, kg
ME = 9.10938356e-31	// electron mass, kg
MP = 1.672621898e-27	// proton mass, kg
MN = 1.674927471e-27	// neutron mass, kg
M = 32*AMU		// ion mass (32*AMU is molecular oxygen, O2)

/* Simulation variables */

n0 = 1e12 		// density, # m^-3
phi0 = 0		// reference potential, V
Te = 1			// electron temperature, eV
Ti = 0.1		// ion temperature, eV
v_drift = 7.0e3		// ion injection velocity, m/s
phi_p = -5		// wall potential, V

solvertol = 1.0e-2	// solver tolerance

setB=1
		// A good scan is Bx= 0, 0.02, 0.1, 1.0
bx = 1.0			// set initial magnetic field strength


/* Plasma parameters */
lD = sqrt(EPS0*Te/(n0*QE))	// Debye length, m
vth = sqrt(2*QE*Ti/M)	// thermal velocity with Ti in eV

/* Simulation domain */
myn = 2;		// multiplier on resolution
nx = 17*myn		// number of nodes in x
ny = 11*myn		// number of nodes in y
ts = 250		// number of time steps
dh = lD/myn		// cell size

np_insert = (ny-1)*15	// insert 15 particles per cell 

nn = nx*ny		// total number of nodes
dt = myn*0.1*dh/v_drift	// time step, moves 0.1 of a cell at vdrift
Lx = (nx-1)*dh		// domain length in x
Ly = (ny-1)*dh		// domain length in y

/* Construct x and y arrays for plotting */
x_var = array(0.0,nx,ny);y_var = array(0.0,nx,ny);
for (xx=1; xx<=nx; ++xx){; y_var(xx,) = span(0.0,Ly,ny);};
for (yy=1; yy<=ny; ++yy){; x_var(,yy) = span(0.0,Lx,nx);};


/* Specify dimensions of charged plate */
box = array(0,2,2);
box(1,:) = [floor(nx/3),floor(nx/3)+(2*myn)]
box(2,:) = [1,floor(ny/3)]

/* Making an object domain of 0s and 1s for visualization of the plate */
object = array(0,nx,ny);
for (j=box(2,1);j<=box(2,2);++j){
	object(box(1,1):box(1,2),j) = array(1.0,box(1,2)-box(1,1)+1,1)
}

/* Calculate the specific weight */
flux = n0*v_drift*Ly	// flux of entering particles
npt = flux*dt		// number of real particles created per timestep
spwt = npt/np_insert	// specific weight, real particles per macroparticle
mp_q = 1		// macroparticle charge
max_part = int(3.0e6)	// buffer size


/* Allocate particle array */
part_x = array(0.0,max_part,2);		// particle positions
part_v = array(0.0,max_part,2);		// particle velocities

/* Set up multiplication matrix for potential solver
	Here we setup the finite difference stencil	*/

A = array(0.0,nn,nn);	// allocates empty nn*nn matrix

/* Set regular stencil on internal nodes */
for (j=2; j<=ny-1; ++j){
	for(i=2; i<=nx-1; ++i){
		u = (j-1)*nx+i
		A(u,u) = -4/(dh*dh)	// phi(i,j)
		A(u,u-1) = 1/(dh*dh)	// phi(i-1,j)
		A(u,u+1) = 1/(dh*dh)	// phi(i+1,j)
		A(u,u-nx) = 1/(dh*dh)	// phi(i,j-1)
		A(u,u+nx) = 1/(dh*dh)	// phi(i,j+1)
	}
}

/* Neumann boundary on y=0 */
for (i=1;i<=nx;++i){
	u=i
	A(u,u) = -1/dh;			// phi(i,j)
	A(u,u+nx) = 1/dh;		// phi(i,j+1)
}

/* Neumann boundary on y=Ly */
for (i=1;i<=nx;++i){
	u=(ny-1)*nx+i
	A(u,u-nx) = 1/dh		// phi(i,j-1)
	A(u,u) = -1/dh			// phi(i,j)
}	

/* Neumann boundary on x=Lx */
for (j=1; j<=ny; ++j){
	u=(j-1)*nx+nx
	A(u,:) = array(0.0,nn)		// clear row
	A(u,u-1) = 1/dh			// phi(i-1,j)
	A(u,u) = -1/dh			// phi(i,j)
}

/* Dirichlet boundary on x=0 */
for (j=1; j<=ny; ++j){
	u=(j-1)*nx+1
	A(u,:) = array(0.0,nn)		// clear row
	A(u,u) = 1			// phi(i,j)
}

/* Dirichlet boundary on nodes of plate*/
for (j=box(2,1);j<=box(2,2);++j){
	for (i=box(1,1);i<=box(1,2);++i){
		u=(j-1)*nx+i
		A(u,:) = array(0.0,nn)	// clear row
		A(u,u) = 1		// phi(i,j)
	}
}


/* Initialize */
phi = phi0*array(1.0,nx,ny);	// set initial potential to phi0
np = 0;				// clear number of particles


print,"Solving potential for the first time.";


/* MAIN LOOP */

for (it=1; it<=ts; ++it){		// iterate for ts time steps
	print, "STEP NUMBER =",it,"WITH TOTAL PARTICLES",np

	// Reset field quantities
	den = array(0.0,nx,ny);		// number density
	efx = array(0.0,nx,ny);		// electric field x-component
	efy = array(0.0,nx,ny);		// electric field y-component
	chg = array(0.0,nx,ny);		// charge distribution

	/* 1. CALCULATE CHARGE DENSITY */

	for (p=1; p<=np; ++p){		// loop over particles
		fi = 1+part_x(p,1)/dh;	// real i index of particle's cell
		i = int(floor(fi));	// integral part
		hx = fi-i;		// the remainder

		fj = 1+part_x(p,2)/dh;	// real j index of particle's cell
		j = int(floor(fj));	// integral part
		hy = fj-j;		// the remainder

		// Interpolate charge to nodes
		chg(i,j) = chg(i,j) + (1-hx)*(1-hy);
		chg(i+1,j) = chg(i+1,j) + hx*(1-hy);
		chg(i,j+1) = chg(i,j+1) + (1-hx)*hy;
		chg(i+1,j+1) = chg(i+1,j+1) + hx*hy;
	}

	// Calculate density on the grid
	den = spwt*mp_q*chg/(dh*dh);

	// Apply boundaries, double the density at boundaries since only half volume contributing
	den(1,:) = 2*den(1,:)
	den(nx,:) = 2*den(nx,:)
	den(:,1) = 2*den(:,1)
	den(:,ny) = 2*den(:,ny)

	den = den+1.0e4		// Add density floor for plotting and to help the solver

	/* 2. CALCULATE POTENTIAL */

	phi = eval_2dpot_GS(phi);

	/* 3. CALCULATE ELECTRIC FIELD */

	efx(2:nx-1,:) = phi(1:nx-2,:) - phi(3:nx,:) 	// central difference on internal nodes
	efy(:,2:ny-1) = phi(:,1:ny-2) - phi(:,3:ny)	// central difference on internal nodes
	efx(1,:) = 2*(phi(1,:) - phi(2,:))		// forward difference on x = 0
	efy(:,1) = 2*(phi(:,1) - phi(:,2))		// forward difference on y = 0
	efx(nx,:) = 2*(phi(nx-1,:) - phi(nx,:))		// backward difference on x = Lx
	efy(:,ny) = 2*(phi(:,ny-1) - phi(:,ny))		// backward difference on y = Ly
	efx = efx/(2*dh)
	efy = efy/(2*dh)

	/* 4. GENERATE NEW PARTICLES */

	// Insert particles randomly distributed in y and in the first cell
	part_x(np+1:np+np_insert,1) = random(np_insert,1)*dh	// x position
	part_x(np+1:np+np_insert,2) = random(np_insert,1)*Ly	// y position

	// Sample Maxwellian in x and y, add drift velocity in +x
	part_v(np+1:np+np_insert,1) = v_drift+vth*(random(np_insert,1)+random(np_insert,1)+random(np_insert,1)-1.5)
	part_v(np+1:np+np_insert,2) = 0.5*vth*(random(np_insert,1)+random(np_insert,1)+random(np_insert,1)-1.5)
	np = np+np_insert;		// increments particle counter


	/* 5. MOVE PARTICLES */

	p=1;
	while (p<=np){				// loop over particles
		fi = 1+part_x(p,1)/dh		// i index of particle's cell
		i = int(floor(fi));
		hx = fi-i			// fractional x position in cell

		fj = 1+part_x(p,2)/dh		// j index of particle's cell
		j = int(floor(fj));
		hy = fj-j;			// fractional y position in cell

		// Gather electric field
		E = array(0.0,2)				// makes E=[0,0] (i,j) for this particle
		E += [efx(i,j),efy(i,j)]*(1-hx)*(1-hy)			// (i,j) contribution
		E += [efx(i+1,j),efy(i+1,j)]*(hx)*(1-hy)		// (i+1,j) contribution
		E += [efx(i,j+1),efy(i,j+1)]*(1-hx)*(hy)		// (i,j+1) contribution
		E += [efx(i+1,j+1),efy(i+1,j+1)]*(hx)*(hy)		// (i+1,j+1) contribution

		// Update velocity and position of particle


//if (setB == 0){
//		F = QE*E					// Lorentz force, F=qE
//		a = F/M						// acceleration
//		part_v(p,:) = part_v(p,:)+a*dt			// update velocity
//		part_x(p,:) = part_x(p,:)+(part_v(p,:)*dt);	// update position
//} else if (setB == 1){
		F = QE*E
		ahalf = 0.5*F/M

		part_v(p,:) = part_v(p,:)+ahalf*dt		// this is now v-
		vprime = part_v(p,2)+(part_v(p,2)*(QE/M)*bx*dt)		// this is now v'
		part_v(p,2) = part_v(p,2)-(vprime*(2*(QE/M)*bx/(1+bx^2))*dt)	// this is now v+
		part_v(p,:) = part_v(p,:)+ahalf*dt		// this is updated v
		part_x(p,:) = part_x(p,:)+(part_v(p,:)*dt);	// update position
//}

	// Process boundaries:
		in_box = 0;					// initialize "in box?" flag
		// Reflective boundary on bottom
		if (part_x(p,2)<= 0.0){				// check for y < 0 (2nd index)
			part_x(p,2) = -part_x(p,2);		// move particle back to domain
			part_v(p,2) = -part_v(p,2);		// reverse y velocity
		}
		// Particle inside object?
		if ( (i>=box(1,1) && i<box(1,2)) && (j>=box(2,1) && j<=box(2,2)-1)){
			in_box = 1;
		}
		// Absorbing boundary on left, right, top or if in box (object)
		if ( part_x(p,1)<0 || part_x(p,1) >= Lx || part_x(p,2) >= Ly || in_box == 1){
			part_x(p,:) = part_x(np,:)		// kill particle by replacing it with last particle
			part_v(p,:) = part_v(np,:)	
			np = np-1;
			p = p-1;
		}

		p = p+1;
	}


	if (it%5 == 0){
		makeplots(np)
		pause,300;
	}

	Etot = sum(0.5*M*sqrt(part_v(,1)^2+part_v(,2)^2)^2.0)/(Lx*(2*Ly)*(2*Ly))
	window,3;
if (setB = 0){
	plmk,Etot,it*dt,marker=4,msize=0.2,width=11,color="black"
} else if (setB = 1){
	plmk,Etot,it*dt,marker=4,msize=0.2,width=11,color="green"
}
	xytitles,"Time","Energy"
	pltitle,"Energy diagnostic"


}
