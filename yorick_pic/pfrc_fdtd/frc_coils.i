/* ################################################################
    Yorick implementation of the PFRC-1 antenna
    and explicit Yee method for Maxwell's equations

    A. B. Sefkow, on behalf of the TriForce team
    Version 1.1
    Date: July 9th, 2021
 ################################################################### */

/* ######################################################################## */
// The following just initializes some timing categories to be used
print, "Simulation begun at",timestamp()
elapsed0= elapsed= total= cat1= cat2= cat3= cat4= array(double,3);
cat5= cat6= cat7= cat8= cat9= cat10= cat11= cat12= array(double,3);
timer, elapsed0;
elapsed = elapsed0;
/* ######################################################################## */

  c0 = 299792458.0;     // speed of light, m/s
eps0 = 8.85418782e-12;  // permittivity of vacuum
 mu0 = 1.25663706e-6;   // permeability of vacuum

freq = 14.0e6;      // frequency in Hz of the PFRC-1 antenna
lambda = c0/freq;   // wavelength of the signal
T = freq^-1.0;      // Period

zmin = -0.25;
zmax = +0.25;
nz = 201;
dz = (zmax-zmin)/(nz-1);

xmin = ymin = -0.125;
xmax = ymax = +0.125;
nx = ny = 101;
dx = (xmax-xmin)/(nx-1);
dy = (ymax-ymin)/(ny-1);

SPATIAL_ORDER = 2;
ifdamping = 0;          // flag for damping fields near the boundaries (or not)

// Create 3D arrays for x, y, and z.
x = y = z = array(0.0, nx, ny, nz);
for (xx=1; xx<=nx; ++xx){;
    for(yy=1;yy<=ny;++yy){;
        z(xx,yy,:) = span(zmin, zmax, nz);};};
for (xx=1; xx<=nx; ++xx){;
    for(zz=1; zz<=nz; ++zz){;
        y(xx,:,zz) = span(ymin, ymax, ny);};};
for (yy=1; yy<=ny; ++yy){;
    for(zz=1; zz<=nz; ++zz){;
        x(:,yy,zz) = span(xmin, xmax, nx);};};

include, "frc_plot_funcs.i";
include, "frc_coils_defs.i";

ifcalc = 1;         // set flag to include externally applied B field
if (ifcalc == 1){;
    include, "frc_helmholtz.i";
};

//plotTheCoils;         // this makes the file to plot the antenna in Paraview.
//plotThePhases;        // neither of these are generally needed.

C = 1.0;    // Courant multiplier. <= 6/7 (0.857) for 4th order, 1 for 2nd order

if ( SPATIAL_ORDER == 2 ){;
    dt = C/(c0*sqrt(1.0/dz^2 + 1.0/dy^2 + 1.0/dx^2)); // time step
    dt = 4.0e-12;
} else if ( SPATIAL_ORDER == 4){;
    dt = C*(6./7.)/(c0*sqrt(1.0/dz^2 + 1.0/dy^2 + 1.0/dx^2)); // time step
};
    
// Initialize the three components of both fields
Ex = array(0.0, nx-1, ny, nz);
Ey = array(0.0, nx, ny-1, nz);
Ez = array(0.0, nx, ny, nz-1);
Bx = BxExt = array(0.0, nx, ny-1, nz-1);
By = ByExt = array(0.0, nx-1, ny, nz-1);
Bz = BzExt = array(0.0, nx-1, ny-1, nz);

print, "Setting the external fields...";
for (k=1; k<=nz; ++k){;
    for (j=1; j<=ny; ++j){;
        for (i=1; i<=nx; ++i){;
            xVal = x(i,j,k);
            yVal = y(i,j,k);
            zVal = z(i,j,k);
            rVal = sqrt( xVal^2 + yVal^2 );
            Bri = interp(extBr(:, k), rExt(:, k), rVal);
            rind = int(interp(span(1,nr,nr), rExt(:, k), rVal));
            theTheta = atan(yVal, xVal);
            if ( j<= ny-1 && k<=nz-1 ){;
                BxExt(i,j,k) = Bri*cos(theTheta);};
            if ( i<= nx-1 && k<=nz-1 ){;
                ByExt(i,j,k) = Bri*sin(theTheta);};
            if ( i <= nx-1 && j <=ny-1){;
                BzExt(i,j,k) = extBz(rind, k);};
        };
    };
};

print, "Adding Bz=5 G everywhere (like LSP)";
BzExt += 4.2; // 5.0; LSP does 5 G but something is different (interpolation?)

print, "####################################################";

/* Maxwell's equations:
 dB/dt = - curl(E)
        dBx/dt = - ( dEz/dy - dEy/dz)
        dBy/dt = - ( dEx/dz - dEz/dx)
        dBz/dt = - ( dEy/dx - dEx/dy)
 dE/dt = c^2*curl(B) - c^2*mu0*J
        dEx/dt = c^2 ( dBz/dy - dBy/dz ) - c^2*mu0 * Jx
        dEy/dt = c^2 ( dBx/dz - dBz/dx ) - c^2*mu0 * Jy
        dEz/dt = c^2 ( dBy/dx - dBx/dy ) - c^2*mu0 * Jz       (  J -> E for E/(c^2*mu0) )
 */

// nt is the number of total time steps we will take:
nt = int(25.0e-9/dt);
print, "DT", dt;
/* These are just comments I want to keep, you can ignore them. */
// 5/8 of a period is when the two sets have equal and opposite amplitudes: 5.*T/8 = 4.46429e-8
//nt = int( ((T/4) + (11.0* (T/8.)))/dt);

print, "INITIALIZATION COMPLETE. BEGINNING TIME LOOP.";
interval = 10;              // interval to print time step and make plot
report_interval = 500;      // interval to print timing information to output
timer, elapsed, cat1;       // add elapsed time this timing category

// ********* EXPLICIT YEE ALGORITHM ************
for (t=1;t<=nt;++t){        // Time advance loop
    time = t*dt;            // determine the current time
    if ( t%interval == 0 || t == nt){;
      print, "Explicit time step",t,"and time is",time*1e9,"ns";
    };
    
    /* ######################################################################
                 Advance B field using old B and E fields
      ####################################################################### */
    
    if ( SPATIAL_ORDER == 2 ){;
    Bx(1:nx, 1:ny-1, 1:nz-1) -= dt*( (Ez(1:nx,2:ny,1:nz-1)-Ez(1:nx,1:ny-1,1:nz-1))/dy
                                    -(Ey(1:nx,1:ny-1,2:nz)-Ey(1:nx,1:ny-1,1:nz-1))/dz );
    By(1:nx-1, 1:ny, 1:nz-1) -= dt*( (Ex(1:nx-1,1:ny,2:nz)-Ex(1:nx-1,1:ny,1:nz-1))/dz
                                    -(Ez(2:nx,1:ny,1:nz-1)-Ez(1:nx-1,1:ny,1:nz-1))/dx );
    Bz(1:nx-1, 1:ny-1, 1:nz) -= dt*( (Ey(2:nx,1:ny-1,1:nz)-Ey(1:nx-1,1:ny-1,1:nz))/dx
                                    -(Ex(1:nx-1,2:ny,1:nz)-Ex(1:nx-1,1:ny-1,1:nz))/dy );
    };
    
    /* ######################################################################
       Advance E field (ahead by 1/2 cell) using old E field and new B fields
      ####################################################################### */
    
    if ( SPATIAL_ORDER == 2 ){;
        Ex(1:nx-1, 2:ny-1, 2:nz-1) += c0^2*dt*( (Bz(1:nx-1,2:ny-1,2:nz-1)-Bz(1:nx-1,1:ny-2,2:nz-1))/(dy*epsR(1:nx-1, 2:ny-1, 2:nz-1))
                                    -(By(1:nx-1,2:ny-1,2:nz-1)-By(1:nx-1,2:ny-1,1:nz-2))
                                    /(dz*epsR(1:nx-1, 2:ny-1, 2:nz-1)) );
        Ey(2:nx-1, 1:ny-1, 2:nz-1) += c0^2*dt*( (Bx(2:nx-1,1:ny-1,2:nz-1)-Bx(2:nx-1,1:ny-1,1:nz-2))
                                    /(dz*epsR(2:nx-1, 1:ny-1, 2:nz-1))
                                    -(Bz(2:nx-1,1:ny-1,2:nz-1)-Bz(1:nx-2,1:ny-1,2:nz-1))
                                    /(dx*epsR(2:nx-1, 1:ny-1, 2:nz-1)) );
        Ez(2:nx-1, 2:ny-1, 1:nz-1) += c0^2*dt*( (By(2:nx-1,2:ny-1,1:nz-1)-By(1:nx-2,2:ny-1,1:nz-1))
                                    /(dx*epsR(2:nx-1, 2:ny-1, 1:nz-1))
                                    -(Bx(2:nx-1,2:ny-1,1:nz-1)-Bx(2:nx-1,1:ny-2,1:nz-1))
                                    /(dy*epsR(2:nx-1, 2:ny-1, 1:nz-1)) );
    };
    /* ######################################################################
                    Apply source currents
      ####################################################################### */
    source = -8.100e3 ;

    // Determine the "TB" (Top/Bottom) source value
    if ( time <= T/4){;
        time_dep_TB = 0.0;
    } else {;
        if ( time < (T/4 + 1.0e-9) ){;
            time_dep_TB = ((time-T/4)/1.0e-9) * source * sin(2.0*pi*freq*time );
        } else {;
            time_dep_TB = source * sin(2.0*pi*freq*time );
        };
    };
    // Determine the "LR" (Left/Right) source value
    if ( time < 1.0e-9 ){;
        time_dep_LR = (time/1.0e-9) * source * sin(2.0*pi*freq*time + pi/2);
    } else {;
        time_dep_LR = source * sin(2.0*pi*freq*time + pi/2.0 );
    };
           
    foo = applyBCsTop(Ex, Ey, Ez, Bx, By, Bz);
    foo = applyBCsBottom(Ex, Ey, Ez, Bx, By, Bz);
    foo = applyBCsRight(Ex, Ey, Ez, Bx, By, Bz);
    foo = applyBCsLeft(Ex, Ey, Ez, Bx, By, Bz);
    
    /* ######################################################################
                    Apply internal boundary conditions
      ####################################################################### */

    bar = applyICsTop(Ex, Ey, Ez, Bx, By, Bz);
    bar = applyICsBottom(Ex, Ey, Ez, Bx, By, Bz);
    bar = applyICsRight(Ex, Ey, Ez, Bx, By, Bz);
    bar = applyICsLeft(Ex, Ey, Ez, Bx, By, Bz);

    // Flux conserving rings are grounded:
    Ex *= intBC_Ex; Ey *= intBC_Ey; Ez *= intBC_Ez;
    Ex *= intBC_Bx; By *= intBC_By; Bz *= intBC_Bz;
        
    /* ######################################################################
            Apply damping at boundaries
      ####################################################################### */
   
    if (ifdamping == 1){;
        foo = dampMin(Ex); foo = dampMin(Ey); foo = dampMin(Ez);
        foo = dampMin(Bx); foo = dampMin(By); foo = dampMin(Bz);
        foo = dampMax(Ex); foo = dampMax(Ey); foo = dampMax(Ez);
        foo = dampMax(Bx); foo = dampMax(By); foo = dampMax(Bz);
    };
    timer, elapsed, cat3;     // add elapsed time this timing category
        
    /* ######################################################################
            Make plots and save information to files, if desired
      ####################################################################### */
    
    if ( t%interval == 0 || t == nt){;
        timer, elapsed, cat3;     // add elapsed time this timing category
        dummy= myplot2d(5);
        timer, elapsed, cat2;     // add elapsed time this timing category    };
    };

    if ( time > (T/4.+1.0e-9) && time%(T/32) < dt ){;
  //      dummy = save_fields(t);
    };
    
    /* ######################################################################
            Report timings, if desired
      ####################################################################### */
    
    if ( t%report_interval == 0 || t == nt){;
        timer, elapsed0, total;      // add elapsed time this timing category
        category_other=total-(cat1+cat2+cat3);
        timer_print,                 // Print some timing reports
        "Initialization", cat1,
        "Plotting",cat2,
        "Time loop",cat3,
        "Others...",category_other,"TOTAL", total;
    };
};
