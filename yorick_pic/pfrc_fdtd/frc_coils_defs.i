/* ################################################################
    Implementation of the FRC antenna
    and explicit Yee method for Maxwell's equations

    --> This file contains functions related to defining
        the (1) antenna, (2) faraday cages, and (3) flux-conserving
        rings that make up the PFRC-1 geometry. <--
 
    A. B. Sefkow, on behalf of the TriForce team
    Version 1.1
    Date: July 9th, 2021
 ################################################################### */

print, "LOADED 'frc_coils_defs.i'!"

/* ########  Definitions for FC (flux-conserving) rings ######### */

/*
 FC ring #1
 zloc = -3.425 cm
 dz = 0.5 cm (height)
 outerRad = 5.08 cm
 innerRad = 4.1275 cm

 FC ring #2
 zloc = 2.9025 cm <-- probably wrong! should be 2.925 (to 3.425)
 dz = 0.5 cm (height)
 outerRad = 5.08 cm
 innerRad = 4.1275 cm
 
 FC ring #3
 zloc = 9.275 cm
 dz = 0.5 cm (height)
 outerRad = 5.08 cm
 innerRad = 3.683 cm

 FC ring #4
 zloc = -9.775 cm
 dz = 0.5 cm (height)
 outerRad = 5.08 cm
 innerRad = 3.683 cm

 FC ring #5
 zloc = 16.26 cm
 dz = 0.635 cm (height)
 outerRad = 5.08 cm
 innerRad = 2.73 cm

 FC ring #6
 zloc = -16.76 cm
 dz = 0.635 cm (height)
 outerRad = 5.08 cm
 innerRad = 2.73 cm
*/

// Define the min/max of each coordinate for each of the 28 segments
fc_R = array(0.0, 2, 6);     // (1, :) are the rmin, (2, :) are the rmax
fc_Z = array(0.0, 2, 6);     // (1, :) are the zmin, (2, :) are the zmax

fc_R(1, 1) = 0.043; // 0.041275;
fc_R(2, 1) = 0.050800;
fc_Z(1, 1) = 0.02925;
fc_Z(2, 1) = 0.03425;

fc_R(1, 2) = 0.043; // 0.041275;
fc_R(2, 2) = 0.050800;
fc_Z(1, 2) = -0.03425;
fc_Z(2, 2) = -0.02925;

fc_R(1, 3) = 0.038; // 0.036830;
fc_R(2, 3) = 0.050800;
fc_Z(1, 3) = 0.092750;
fc_Z(2, 3) = 0.097750;

fc_R(1, 4) = 0.038; // 0.036830;
fc_R(2, 4) = 0.050800;
fc_Z(1, 4) = -0.09775;
fc_Z(2, 4) = -0.09275;

fc_R(1, 5) = 0.028; // 0.027300;
fc_R(2, 5) = 0.050800;
fc_Z(1, 5) = 0.162600;
fc_Z(2, 5) = 0.167600;

fc_R(1, 6) = 0.028; // 0.027300;
fc_R(2, 6) = 0.050800;
fc_Z(1, 6) = -0.16760;
fc_Z(2, 6) = -0.16260;

/* Initialize the needed arrays for the internal boundary conditions */
intBC_Ex = array(1.0, nx-1,   ny,   nz);
intBC_Ey = array(1.0,   nx, ny-1,   nz);
intBC_Ez = array(1.0,   nx,   ny, nz-1);
intBC_Bx = array(1.0,   nx, ny-1, nz-1);
intBC_By = array(1.0, nx-1,   ny, nz-1);
intBC_Bz = array(1.0, nx-1, ny-1,   nz);
epsR = array(1.0, nx-1, ny-1, nz);

print, "Setting the flux-conserving rings...";
for (k=1; k<=nz-1; ++k){;
    for (j=1; j<=ny-1; ++j){;
        for (i=1; i<=nx-1; ++i){;
            xVal = x(i,j,k)+dx/2.;
            yVal = y(i,j,k)+dy/2.;
            zVal = z(i,j,k)+dz/2.;
            xVal1 = x(i,j,k);
            yVal1 = y(i,j,k);
            zVal1 = z(i,j,k);
            xVal2 = x(i+1,j,k);
            yVal2 = y(i,j+1,k);
            zVal2 = z(i,j,k+1);

            rVal = sqrt( xVal^2 + yVal^2 );
            rVal1 = sqrt( xVal1^2 + yVal1^2 );
            rVal2 = sqrt( xVal2^2 + yVal2^2 );

            if ( ( rVal >= 0.044 ) && ( rVal <= 0.050 )){;
                epsR(i,j,k) = 2.0;
            };
            
            for (nn =1; nn<=6; ++nn){;
                // There are two different ways (at least) we could determine
                // whether a cell should be flagged to be "within" the geometry:
    /* If center of the cell is within the boundaries of the object:  */
            if ( (rVal >= fc_R(1,nn)) && (rVal <= fc_R(2, nn)) ){;
            if ( (zVal >= fc_Z(1,nn)) && (zVal <= fc_Z(2, nn)) ){;

    /* If the whole cell is fully within the boundaries of the object: */
//           if ( (rVal1 >= fc_R(1,nn)) && (rVal2 <= fc_R(2,nn)) ){;
//           if ( (zVal1 >= fc_Z(1,nn)) && (zVal2 <= fc_Z(2,nn)) ){;
                epsR(i,j,k) = 1.0;

    /* The following are the 12 E field values and 6 B field values that must be set
        for the Yee cell to be an ideal conductor. */
                intBC_Ex(i,j,k) = 0.0;
                intBC_Ex(i,j+1,k) = 0.0;
                intBC_Ex(i,j,k+1) = 0.0;
                intBC_Ex(i,j+1,k+1) = 0.0;

                intBC_Ey(i,j,k) = 0.0;
                intBC_Ey(i+1,j,k) = 0.0;
                intBC_Ey(i,j,k+1) = 0.0;
                intBC_Ey(i+1,j,k+1) = 0.0;

                intBC_Ez(i,j,k) = 0.0;
                intBC_Ez(i+1,j,k) = 0.0;
                intBC_Ez(i,j+1,k) = 0.0;
                intBC_Ez(i+1,j+1,k) = 0.0;

                intBC_Bx(i,j,k) = 0.0
                intBC_Bx(i+1,j,k) = 0.0

                intBC_By(i,j,k) = 0.0
                intBC_By(i,j+1,k) = 0.0

                intBC_Bz(i,j,k) = 0.0
                intBC_Bz(i,j,k+1) = 0.0
            };
            };
            };
        };
    };
};

/* ###################################################################
                    PFRC-1 ANTENNA DEFINITIONS
   #################################################################

   There are 28 total segments (7*4=28) making up the PFRC-1 antenna:

   TOP COIL (at xmax). There are seven segments.
   BOTTOM COIL (at xmin). There are seven segments.
   RIGHT COIL (at ymax). There are seven segments.
   LEFT COIL (at ymin). There are seven segments.
 
    The TOP and BOTTOM have the same time dependence and phase
    The RIGHT and LEFT have the same time dependence and phase
    However, the two sets are 90 degrees (pi/2) out of phase
    relative to each other.
 
####################################################################### */

func getCoilLocation(num, &x_coils, &y_coils, &z_coils){;
    x_coils(1, num) = int(1+round((coil_X(1,num)-xmin)/dx));
    x_coils(2, num) = int(1+round((coil_X(2,num)-xmin)/dx));
    y_coils(1, num) = int(1+round((coil_Y(1,num)-ymin)/dy));
    y_coils(2, num) = int(1+round((coil_Y(2,num)-ymin)/dy));
    z_coils(1, num) = int(1+round((coil_Z(1,num)-zmin)/dz));
    z_coils(2, num) = int(1+round((coil_Z(2,num)-zmin)/dz));
};

func getWireLocation(num, &x_wires, &y_wires, &z_wires){;
    x_wires(1, num) = int(1+round((wire_X(1,num)-xmin)/dx));
    x_wires(2, num) = int(1+round((wire_X(2,num)-xmin)/dx));
    y_wires(1, num) = int(1+round((wire_Y(1,num)-ymin)/dy));
    y_wires(2, num) = int(1+round((wire_Y(2,num)-ymin)/dy));
    z_wires(1, num) = int(1+round((wire_Z(1,num)-zmin)/dz));
    z_wires(2, num) = int(1+round((wire_Z(2,num)-zmin)/dz));
};


// Define the min/max of each coordinate for each of the 28 segments
coil_X = array(0.0, 2, 28);     // (1, :) are the xmin, (2, :) are the xmax
coil_Y = array(0.0, 2, 28);     // (1, :) are the ymin, (2, :) are the ymax
coil_Z = array(0.0, 2, 28);     // (1, :) are the zmin, (2, :) are the zmax

x_coils = array(0, 2, 28);    // (1, num) is min, (2, num) is max
y_coils = array(0, 2, 28);    // (1, num) is min, (2, num) is max
z_coils = array(0, 2, 28);    // (1, num) is min, (2, num) is max

dzShift = 0.0;      // Optional shifting of z origin.

/* ######################################################################
        TOP COIL (at xmax). There are seven segments.
  ####################################################################### */
    
// Segment #1. Time function TB.  Voltage in +Z.
num=1;
coil_X(1, num) = coil_X(2, num) = 0.075;
coil_Y(1, num) = coil_Y(2, num) = 0.065;
coil_Z(1, num) = -0.14+dzShift; coil_Z(2, num) = 0.0+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);
     
// Segment #2. Time function TB.  Voltage in -Z.
num=2;
coil_X(1, num) = coil_X(2, num) = 0.075;
coil_Y(1, num) = coil_Y(2, num) = 0.065;
coil_Z(1, num) = 0.0+dzShift; coil_Z(2, num) = 0.14+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

// Segment #3. Time function TB.  Voltage in -Z.
num=3;
coil_X(1, num) = coil_X(2, num) = 0.075;
coil_Y(1, num) = coil_Y(2, num) = -0.065;
coil_Z(1, num) = -0.14+dzShift; coil_Z(2, num) = 0.0+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

// Segment #4. Time function TB.  Voltage in +Z.
num=4;
coil_X(1, num) = coil_X(2, num) = 0.075;
coil_Y(1, num) = coil_Y(2, num) = -0.065;
coil_Z(1, num) = 0.0+dzShift; coil_Z(2, num) = 0.14+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

// Segment #5. Time function TB.  Voltage in +Y.
num=5;
coil_X(1, num) = coil_X(2, num) = 0.075;
coil_Y(1, num) = -0.065; coil_Y(2, num) = 0.065;
coil_Z(1, num) = coil_Z(2, num) = -0.14+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

// Segment #6. Time function TB.  Voltage doubled in -Y.
num=6;
coil_X(1, num) = coil_X(2, num) = 0.075;
coil_Y(1, num) = -0.065; coil_Y(2, num) = 0.065;
coil_Z(1, num) = coil_Z(2, num) = 0.0+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

// Segment #7. Time function TB.  Voltage in +Y.
num=7;
coil_X(1, num) = coil_X(2, num) = 0.075;
coil_Y(1, num) = -0.065; coil_Y(2, num) = 0.065;
coil_Z(1, num) = coil_Z(2, num) = 0.14+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

/* ######################################################################
        BOTTOM COIL (at xmin). There are seven segments.
  ####################################################################### */
 
// Segment #8. Time function TB.  Voltage in +Z.
num=8;
coil_X(1, num) = coil_X(2, num) = -0.075;
coil_Y(1, num) = coil_Y(2, num) = 0.065;
coil_Z(1, num) = -0.14+dzShift; coil_Z(2, num) = 0.0+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

// Segment #9. Time function TB.  Voltage in -Z.
num=9;
coil_X(1, num) = coil_X(2, num) = -0.075;
coil_Y(1, num) = coil_Y(2, num) = 0.065;
coil_Z(1, num) = 0.0+dzShift; coil_Z(2, num) = 0.14+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

// Segment #10. Time function TB.  Voltage in -Z.
num=10;
coil_X(1, num) = coil_X(2, num) = -0.075;
coil_Y(1, num) = coil_Y(2, num) = -0.065;
coil_Z(1, num) = -0.14+dzShift; coil_Z(2, num) = 0.0+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

// Segment #11. Time function TB.  Voltage in +Z.
num=11;
coil_X(1, num) = coil_X(2, num) = -0.075;
coil_Y(1, num) = coil_Y(2, num) = -0.065;
coil_Z(1, num) = 0.0+dzShift; coil_Z(2, num) = 0.140+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

// Segment #12. Time function TB.  Voltage in +Y.
num=12;
coil_X(1, num) = coil_X(2, num) = -0.075;
coil_Y(1, num) = -0.065; coil_Y(2, num) = 0.065;
coil_Z(1, num) = coil_Z(2, num) = -0.14+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

// Segment #13. Time function TB.  Voltage doubled in -Y.
num=13;
coil_X(1, num) = coil_X(2, num) = -0.075;
coil_Y(1, num) = -0.065; coil_Y(2, num) = 0.065;
coil_Z(1, num) = coil_Z(2, num) = 0.0+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

// Segment #14. Time function TB.  Voltage in +Y.
num=14;
coil_X(1, num) = coil_X(2, num) = -0.075;
coil_Y(1, num) = -0.065; coil_Y(2, num) = 0.065;
coil_Z(1, num) = coil_Z(2, num) = 0.14+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

/* ######################################################################
        RIGHT COIL (at ymax). There are seven segments.
  ####################################################################### */
 
// Segment #15. Time function LR.  Voltage in +Z.
num=15;
coil_X(1, num) = coil_X(2, num) = 0.065;
coil_Y(1, num) = coil_Y(2, num) = 0.075;
coil_Z(1, num) = -0.14+dzShift; coil_Z(2, num) = 0.0+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

// Segment #16. Time function LR.  Voltage in -Z.
num=16;
coil_X(1, num) = coil_X(2, num) = 0.065;
coil_Y(1, num) = coil_Y(2, num) = 0.075;
coil_Z(1, num) = 0.0+dzShift; coil_Z(2, num) = 0.14+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

// Segment #17. Time function LR.  Voltage in -Z.
num=17;
coil_X(1, num) = coil_X(2, num) = -0.065;
coil_Y(1, num) = coil_Y(2, num) = 0.075;
coil_Z(1, num) = -0.14+dzShift; coil_Z(2, num) = 0.0+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

// Segment #18. Time function LR.  Voltage in +Z.
num=18;
coil_X(1, num) = coil_X(2, num) = -0.065;
coil_Y(1, num) = coil_Y(2, num) = 0.075;
coil_Z(1, num) = 0.0+dzShift; coil_Z(2, num) = 0.14+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

// Segment #19. Time function LR.  Voltage in +X.
num=19;
coil_X(1, num) = -0.065; coil_X(2, num) = 0.065;
coil_Y(1, num) = coil_Y(2, num) = 0.075;
coil_Z(1, num) = coil_Z(2, num) = -0.14+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

// Segment #20. Time function LR.  Voltage doubled in -X.
num=20;
coil_X(1, num) = -0.065; coil_X(2, num) = 0.065;
coil_Y(1, num) = coil_Y(2, num) = 0.075;
coil_Z(1, num) = coil_Z(2, num) = 0.0+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

// Segment #21. Time function LR.  Voltage in +X.
num=21;
coil_X(1, num) = -0.065; coil_X(2, num) = 0.065;
coil_Y(1, num) = coil_Y(2, num) = 0.075;
coil_Z(1, num) = coil_Z(2, num) = 0.14+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

/* ######################################################################
        LEFT COIL (at ymin). There are seven segments.
  ####################################################################### */
 
// Segment #22. Time function LR.  Voltage in +Z.
num=22;
coil_X(1, num) = coil_X(2, num) = 0.065;
coil_Y(1, num) = coil_Y(2, num) = -0.075;
coil_Z(1, num) = -0.14+dzShift; coil_Z(2, num) = 0.0+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

// Segment #23. Time function LR.  Voltage in -Z.
num=23;
coil_X(1, num) = coil_X(2, num) = 0.065;
coil_Y(1, num) = coil_Y(2, num) = -0.075;
coil_Z(1, num) = 0.0+dzShift; coil_Z(2, num) = 0.14+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

// Segment #24. Time function LR.  Voltage in -Z.
num=24;
coil_X(1, num) = coil_X(2, num) = -0.065;
coil_Y(1, num) = coil_Y(2, num) = -0.075;
coil_Z(1, num) = -0.14+dzShift; coil_Z(2, num) = 0.0+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

// Segment #25. Time function LR.  Voltage in +Z.
num=25;
coil_X(1, num) = coil_X(2, num) = -0.065;
coil_Y(1, num) = coil_Y(2, num) = -0.075;
coil_Z(1, num) = 0.0+dzShift; coil_Z(2, num) = 0.14+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

// Segment #26. Time function LR.  Voltage in +X.
num=26;
coil_X(1, num) = -0.065; coil_X(2, num) = 0.065;
coil_Y(1, num) = coil_Y(2, num) = -0.075;
coil_Z(1, num) = coil_Z(2, num) = -0.14+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

// Segment #27. Time function LR.  Voltage doubled in -X.
num=27;
coil_X(1, num) = -0.065; coil_X(2, num) = 0.065;
coil_Y(1, num) = coil_Y(2, num) = -0.075;
coil_Z(1, num) = coil_Z(2, num) = 0.0+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

// Segment #28. Time function LR.  Voltage in +X.
num=28;
coil_X(1, num) = -0.065; coil_X(2, num) = 0.065;
coil_Y(1, num) = coil_Y(2, num) = -0.075;
coil_Z(1, num) = coil_Z(2, num) = 0.14+dzShift;
foo = getCoilLocation(num, x_coils, y_coils, z_coils);

 /* This is currently unused code. It will be used for nonuniform dx/dy/dz.
x1_coil = x2_coil = y1_coil = y2_coil = z1_coil = z2_coil = array(0, 28);
for (num=1; num<=28; ++num){
    for (xx=2;xx<=nx-1;++xx){;
        if ( coil_X(1,num) > 0.0 ){
            if ( (x(xx,1,1) <= coil_X(1,num)+1.0e-15) && (x(xx+1,1,1) > coil_X(1,num)+1.0e-15) ){
                print, "x1_ =",xx,"for num=",num;
                x1_coil(num) = xx;
            };
        } else {
            if ( (x(xx,1,1) <= coil_X(1,num)-1.0e-15) && (x(xx+1,1,1) > coil_X(1,num)-1.0e-15) ){
                print, "x1_ =",xx-1,"for num=",num;
                x1_coil(num) = xx-1;
            };
        };
    };
};
*/
    /* ######################################################################
          There are 4 Faraday cages, composed of 4 wires each
    ####################################################################### */
        
    // Define the min/max of each coordinate for each of the 16 segments
    wire_X = array(0.0, 2, 16);     // (1, :) are the xmin, (2, :) are the xmax
    wire_Y = array(0.0, 2, 16);     // (1, :) are the ymin, (2, :) are the ymax
    wire_Z = array(0.0, 2, 16);     // (1, :) are the zmin, (2, :) are the zmax

    x_wires = array(0, 2, 16);    // (1, num) is min, (2, num) is max
    y_wires = array(0, 2, 16);    // (1, num) is min, (2, num) is max
    z_wires = array(0, 2, 16);    // (1, num) is min, (2, num) is max

     /* ######################################################################
           TOP CAGE (at xmax). There are four wires.
     ####################################################################### */

     // Wire #1
     num=1;
     wire_X(1, num) = wire_X(2, num) = 0.0654421;
     wire_Y(1, num) = wire_Y(2, num) = 0.0567176;
     wire_Z(1, num) = -0.14+dzShift; wire_Z(2, num) = 0.14+dzShift;
     foo = getWireLocation(num, x_wires, y_wires, z_wires);

     // Wire #2
     num=2;
     wire_X(1, num) = wire_X(2, num) = 0.0654421;
     wire_Y(1, num) = wire_Y(2, num) = -0.0567176;
     wire_Z(1, num) = -0.14+dzShift; wire_Z(2, num) = 0.14+dzShift;
     foo = getWireLocation(num, x_wires, y_wires, z_wires);

     // Wire #3
     num=3;
     wire_X(1, num) = wire_X(2, num) = 0.0654421;
     wire_Y(1, num) = -0.0567176;
     wire_Y(2, num) = 0.0567176;
     wire_Z(1, num) = wire_Z(2,num) = -0.14+dzShift;
     foo = getWireLocation(num, x_wires, y_wires, z_wires);

     // Wire #4
     num=4;
     wire_X(1, num) = wire_X(2, num) = 0.0654421;
     wire_Y(1, num) = -0.0567176;
     wire_Y(2, num) = 0.0567176;
     wire_Z(1, num) = wire_Z(2,num) = +0.14+dzShift;
     foo = getWireLocation(num, x_wires, y_wires, z_wires);

     /* ######################################################################
           BOTTOM CAGE (at xmin). There are four wires.
     ####################################################################### */
     
      // Wire #5
      num=5;
      wire_X(1, num) = wire_X(2, num) = -0.0654421;
      wire_Y(1, num) = wire_Y(2, num) = 0.0567176;
      wire_Z(1, num) = -0.14+dzShift; wire_Z(2, num) = 0.14+dzShift;
      foo = getWireLocation(num, x_wires, y_wires, z_wires);

      // Wire #6
      num=6;
      wire_X(1, num) = wire_X(2, num) = -0.0654421;
      wire_Y(1, num) = wire_Y(2, num) = -0.0567176;
      wire_Z(1, num) = -0.14+dzShift; wire_Z(2, num) = 0.14+dzShift;
      foo = getWireLocation(num, x_wires, y_wires, z_wires);

      // Wire #7
      num=7;
      wire_X(1, num) = wire_X(2, num) = -0.0654421;
      wire_Y(1, num) = -0.0567176;
      wire_Y(2, num) = 0.0567176;
      wire_Z(1, num) = wire_Z(2,num) = -0.14+dzShift;
      foo = getWireLocation(num, x_wires, y_wires, z_wires);

      // Wire #8
      num=8;
      wire_X(1, num) = wire_X(2, num) = -0.0654421;
      wire_Y(1, num) = -0.0567176;
      wire_Y(2, num) = 0.0567176;
      wire_Z(1, num) = wire_Z(2,num) = +0.14+dzShift;
      foo = getWireLocation(num, x_wires, y_wires, z_wires);

     /* ######################################################################
           RIGHT CAGE (at ymax). There are four wires.
     ####################################################################### */
     
      // Wire #9
      num=9;
      wire_X(1, num) = wire_X(2, num) = 0.0567176;
      wire_Y(1, num) = wire_Y(2, num) = 0.0654421;
      wire_Z(1, num) = -0.14+dzShift; wire_Z(2, num) = 0.14+dzShift;
      foo = getWireLocation(num, x_wires, y_wires, z_wires);

      // Wire #10
      num=10;
      wire_X(1, num) = wire_X(2, num) = -0.0567176;
      wire_Y(1, num) = wire_Y(2, num) = 0.0654421;
      wire_Z(1, num) = -0.14+dzShift; wire_Z(2, num) = 0.14+dzShift;
      foo = getWireLocation(num, x_wires, y_wires, z_wires);

      // Wire #11
      num=11;
      wire_X(1, num) = -0.0567176;
      wire_X(2, num) = 0.0567176;
      wire_Y(1, num) = wire_Y(2, num) = 0.0654421;
      wire_Z(1, num) = wire_Z(2,num) = -0.14+dzShift;
      foo = getWireLocation(num, x_wires, y_wires, z_wires);

      // Wire #12
      num=12;
      wire_X(1, num) = -0.0567176;
      wire_X(2, num) = 0.0567176;
      wire_Y(1, num) = wire_Y(2, num) = 0.0654421;
      wire_Z(1, num) = wire_Z(2,num) = +0.14+dzShift;
      foo = getWireLocation(num, x_wires, y_wires, z_wires);

     /* ######################################################################
           LEFT CAGE (at ymin). There are four wires.
     ####################################################################### */
     
      // Wire #13
      num=13;
      wire_X(1, num) = wire_X(2, num) = 0.0567176;
      wire_Y(1, num) = wire_Y(2, num) = -0.0654421;
      wire_Z(1, num) = -0.14+dzShift; wire_Z(2, num) = 0.14+dzShift;
      foo = getWireLocation(num, x_wires, y_wires, z_wires);

      // Wire #14
      num=14;
      wire_X(1, num) = wire_X(2, num) = -0.0567176;
      wire_Y(1, num) = wire_Y(2, num) = -0.0654421;
      wire_Z(1, num) = -0.14+dzShift; wire_Z(2, num) = 0.14+dzShift;
      foo = getWireLocation(num, x_wires, y_wires, z_wires);

      // Wire #15
      num=15;
      wire_X(1, num) = -0.0567176;
      wire_X(2, num) = 0.0567176;
      wire_Y(1, num) = wire_Y(2, num) = -0.0654421;
      wire_Z(1, num) = wire_Z(2,num) = -0.14+dzShift;
      foo = getWireLocation(num, x_wires, y_wires, z_wires);

      // Wire #16
      num=16;
      wire_X(1, num) = -0.0567176;
      wire_X(2, num) = 0.0567176;
      wire_Y(1, num) = wire_Y(2, num) = -0.0654421;
      wire_Z(1, num) = wire_Z(2,num) = +0.14+dzShift;
      foo = getWireLocation(num, x_wires, y_wires, z_wires);

/* ################################################################### */

// The following are definitions needed for applying damping at the edges
extraForDamping = 0.025;    // The amount of space at the edge where damping happens

xdampL = int(ceil(extraForDamping/dx));
xdampR = nx - int(ceil(extraForDamping/dx));
ydampL = int(ceil(extraForDamping/dy)) ;
ydampR = ny - int(ceil(extraForDamping/dy));
zdampL = int(ceil(extraForDamping/dz)) ;
zdampR = nz - int(ceil(extraForDamping/dz));

zfac = 0.0125;          // This is the damping factor ~exp(-z/zfac)
preEval1 = exp(-1.0*(x(xdampL,1,1)-x(1:xdampL,:,:))/zfac);
preEval2 = exp(-1.0*(y(1,ydampL,1)-y(:,1:ydampL,:))/zfac);
preEval3 = exp(-1.0*(z(1,1,zdampL)-z(:,:,1:zdampL))/zfac);
preEval4 = exp(-1.0*(x(xdampR:nx,:,:)-x(xdampR,1,1))/zfac);
preEval5 = exp(-1.0*(y(:,ydampR:ny,:)-y(1,ydampR,1))/zfac);
preEval6 = exp(-1.0*(z(:,:,zdampR:nz)-z(1,1,zdampR))/zfac);

func dampMin(&Field){;
    theNx = dimsof(Field)(2);
    theNy = dimsof(Field)(3);
    theNz = dimsof(Field)(4);
    Field(1:xdampL, 1:theNy, 1:theNz) *= preEval1(1:xdampL, 1:theNy, 1:theNz);
    Field(1:theNx, 1:ydampL, 1:theNz) *= preEval2(1:theNx, 1:ydampL, 1:theNz);
    Field(1:theNx, 1:theNy, 1:zdampL) *= preEval3(1:theNx, 1:theNy, 1:zdampL);
};

func dampMax(&Field){;
    theNx = dimsof(Field)(2);
    theNy = dimsof(Field)(3);
    theNz = dimsof(Field)(4);
    numX = theNx - xdampR +1;
    numY = theNy - ydampR +1;
    numZ = theNz - zdampR +1;
    Field(xdampR:theNx, 1:theNy, 1:theNz) *= preEval4(1:numX, 1:theNy, 1:theNz);
    Field(1:theNx, ydampR:theNy, 1:theNz) *= preEval5(1:theNx, 1:numY, 1:theNz);
    Field(1:theNx, 1:theNy, zdampR:theNz) *= preEval6(1:theNx, 1:theNy, 1:numZ);
};

/* ################################################################### */

// THE FOLLOWING FUNCTIONS SET INTERNAL E FIELDS FOR THE COILS

func setAntennaEx(dependence, indices, &ExField, &EyField, &EzField, &BxField, &ByField, &BzField){;
    numx1=indices(1); numx2=indices(2);
    numy1=indices(3); numy2=indices(4);
    numz1=indices(5); numz2=indices(6);
    ExField( numx1:numx2,  numy1:numy2+1, numz1:numz2+1) = dependence;
    EyField(numx1:numx2+1,   numy1:numy2, numz1:numz2+1) = 0.0;
    EzField(numx1:numx2+1, numy1:numy2+1,   numz1:numz2) = 0.0;
    BxField(numx1:numx2+1,   numy1:numy2,   numz1:numz2) = 0.0;
    ByField(  numx1:numx2, numy1:numy2+1,   numz1:numz2) = 0.0;
    BzField(  numx1:numx2,   numy1:numy2, numz1:numz2+1) = 0.0;
};

func setAntennaEy(dependence, indices, &ExField, &EyField, &EzField, &BxField, &ByField, &BzField){;
    numx1=indices(1); numx2=indices(2);
    numy1=indices(3); numy2=indices(4);
    numz1=indices(5); numz2=indices(6);
    ExField( numx1:numx2,  numy1:numy2+1, numz1:numz2+1) = 0.0;
    EyField(numx1:numx2+1,   numy1:numy2, numz1:numz2+1) = dependence;
    EzField(numx1:numx2+1, numy1:numy2+1,   numz1:numz2) = 0.0;
    BxField(numx1:numx2+1,   numy1:numy2,   numz1:numz2) = 0.0;
    ByField(  numx1:numx2, numy1:numy2+1,   numz1:numz2) = 0.0;
    BzField(  numx1:numx2,   numy1:numy2, numz1:numz2+1) = 0.0;
};

func setAntennaEz(dependence, indices, &ExField, &EyField, &EzField, &BxField, &ByField, &BzField){;
    numx1=indices(1); numx2=indices(2);
    numy1=indices(3); numy2=indices(4);
    numz1=indices(5); numz2=indices(6);
    ExField( numx1:numx2,  numy1:numy2+1, numz1:numz2+1) = 0.0;
    EyField(numx1:numx2+1,   numy1:numy2, numz1:numz2+1) = 0.0;
    EzField(numx1:numx2+1, numy1:numy2+1,   numz1:numz2) = dependence;
    BxField(numx1:numx2+1,   numy1:numy2,   numz1:numz2) = 0.0;
    ByField(  numx1:numx2, numy1:numy2+1,   numz1:numz2) = 0.0;
    BzField(  numx1:numx2,   numy1:numy2, numz1:numz2+1) = 0.0;
};

func applyBCsTop(&ExField, &EyField, &EzField, &BxField, &ByField, &BzField){;

    /* In the following setAntenna calls, there are extra -1s on just some *_coils terms.
       This is because of the staggered nature of the Yee cell. Please ask me about it.
       Essentially, positive position values set the fields to slightly larger spatial
       values (half a cell in the direction under consideration).
       In addition, negative position values also set the fields to slightly larger
       spatial values, and hence we actually need to set the cell located one step back
       in the direction under consideration. This is how we will end up having symmetry. */

    /* ######################################################################
            TOP COIL (at xmax). There are seven segments.
      ####################################################################### */
    // Segment #1. Time function TB.  Voltage in +Z.
    num = 1;
    foo = setAntennaEz(time_dep_TB*1.0, [x_coils(1,num)  , x_coils(2,num)  ,
                                         y_coils(1,num)  , y_coils(2,num)  ,
                                         z_coils(1,num)  , z_coils(2,num)-1],
                       ExField, EyField, EzField, BxField, ByField, BzField);
    // Segment #2. Time function TB.  Voltage in -Z.
    num = 2;
    foo = setAntennaEz(time_dep_TB*-1.0,[x_coils(1,num)  , x_coils(2,num)  ,
                                         y_coils(1,num)  , y_coils(2,num)  ,
                                         z_coils(1,num)  , z_coils(2,num)-1],
                       ExField, EyField, EzField, BxField, ByField, BzField);
    // Segment #3. Time function TB.  Voltage in -Z.
    num = 3;
    foo = setAntennaEz(time_dep_TB*-1.0,[x_coils(1,num)  , x_coils(2,num)  ,
                                         y_coils(1,num)-1, y_coils(2,num)-1,
                                         z_coils(1,num)  , z_coils(2,num)-1],
                       ExField, EyField, EzField, BxField, ByField, BzField);
    // Segment #4. Time function TB.  Voltage in +Z.
    num = 4;
    foo = setAntennaEz(time_dep_TB*1.0, [x_coils(1,num)  , x_coils(2,num)  ,
                                         y_coils(1,num)-1, y_coils(2,num)-1,
                                         z_coils(1,num)  , z_coils(2,num)-1],
                       ExField, EyField, EzField, BxField, ByField, BzField);
    // Segment #5. Time function TB.  Voltage in +Y.
    num = 5;
    foo = setAntennaEy(time_dep_TB*1.0, [x_coils(1,num)  , x_coils(2,num)  ,
                                         y_coils(1,num)-1, y_coils(2,num)  ,
                                         z_coils(1,num)  , z_coils(2,num)  ],
                       ExField, EyField, EzField, BxField, ByField, BzField);
    // Segment #6. Time function TB.  Voltage doubled in -Y.  Uses two cells!
    num = 6;
    foo = setAntennaEy(time_dep_TB*-2.0,[x_coils(1,num)  , x_coils(2,num)  ,
                                         y_coils(1,num)-1, y_coils(2,num)  ,
                                         z_coils(1,num)-1, z_coils(2,num)  ],
                       ExField, EyField, EzField, BxField, ByField, BzField);
    // Segment #7. Time function TB.  Voltage in +Y.
    num = 7;
    foo = setAntennaEy(time_dep_TB*1.0, [x_coils(1,num)  , x_coils(2,num)  ,
                                         y_coils(1,num)-1, y_coils(2,num)  ,
                                         z_coils(1,num)-1, z_coils(2,num)-1],
                       ExField, EyField, EzField, BxField, ByField, BzField);
};

func applyBCsBottom(&ExField, &EyField, &EzField, &BxField, &ByField, &BzField){;
    
    /* In the following setAntenna calls, there are extra -1s on just some *_coils terms.
       This is because of the staggered nature of the Yee cell. Please ask me about it.
       Essentially, positive position values set the fields to slightly larger spatial
       values (half a cell in the direction under consideration).
       In addition, negative position values also set the fields to slightly larger
       spatial values, and hence we actually need to set the cell located one step back
       in the direction under consideration. This is how we will end up having symmetry. */

    /* ######################################################################
            BOTTOM COIL (at xmin). There are seven segments.
      ####################################################################### */
    // Segment #8. Time function TB.  Voltage in +Z.
    num = 8;
    foo = setAntennaEz(time_dep_TB*1.0 ,[x_coils(1,num)-1, x_coils(2,num)-1,
                                         y_coils(1,num)  , y_coils(2,num)  ,
                                         z_coils(1,num)  , z_coils(2,num)-1],
                       ExField, EyField, EzField, BxField, ByField, BzField);
    // Segment #9. Time function TB.  Voltage in -Z.
    num = 9;
    foo = setAntennaEz(time_dep_TB*-1.0,[x_coils(1,num)-1, x_coils(2,num)-1,
                                         y_coils(1,num)  , y_coils(2,num)  ,
                                         z_coils(1,num)  , z_coils(2,num)-1],
                       ExField, EyField, EzField, BxField, ByField, BzField);
    // Segment #10. Time function TB.  Voltage in -Z.
    num = 10;
    foo = setAntennaEz(time_dep_TB*-1.0,[x_coils(1,num)-1, x_coils(2,num)-1,
                                         y_coils(1,num)-1, y_coils(2,num)-1,
                                         z_coils(1,num)  , z_coils(2,num)-1],
                       ExField, EyField, EzField, BxField, ByField, BzField);
    // Segment #11. Time function TB.  Voltage in +Z.
    num = 11;
    foo = setAntennaEz(time_dep_TB*1.0 ,[x_coils(1,num)-1, x_coils(2,num)-1,
                                         y_coils(1,num)-1, y_coils(2,num)-1,
                                         z_coils(1,num)  , z_coils(2,num)-1],
                       ExField, EyField, EzField, BxField, ByField, BzField);
    // Segment #12. Time function TB.  Voltage in +Y.
    num = 12;
    foo = setAntennaEy(time_dep_TB*1.0 ,[x_coils(1,num)-1, x_coils(2,num)-1,
                                         y_coils(1,num)-1, y_coils(2,num)  ,
                                         z_coils(1,num)  , z_coils(2,num)  ],
                       ExField, EyField, EzField, BxField, ByField, BzField);
    // Segment #13. Time function TB.  Voltage doubled in -Y.  Uses two cells!
    num = 13;
    foo = setAntennaEy(time_dep_TB*-2.0,[x_coils(1,num)-1, x_coils(2,num)-1,
                                         y_coils(1,num)-1, y_coils(2,num)  ,
                                         z_coils(1,num)-1, z_coils(2,num)  ],
                       ExField, EyField, EzField, BxField, ByField, BzField);
    // Segment #14. Time function TB.  Voltage in +Y.
    num = 14;
    foo = setAntennaEy(time_dep_TB*1.0, [x_coils(1,num)-1, x_coils(2,num)-1,
                                         y_coils(1,num)-1, y_coils(2,num)  ,
                                         z_coils(1,num)-1, z_coils(2,num)-1],
                       ExField, EyField, EzField, BxField, ByField, BzField);
};

func applyBCsRight(&ExField, &EyField, &EzField, &BxField, &ByField, &BzField){;
    
    /* In the following setAntenna calls, there are extra -1s on just some *_coils terms.
       This is because of the staggered nature of the Yee cell. Please ask me about it.
       Essentially, positive position values set the fields to slightly larger spatial
       values (half a cell in the direction under consideration).
       In addition, negative position values also set the fields to slightly larger
       spatial values, and hence we actually need to set the cell located one step back
       in the direction under consideration. This is how we will end up having symmetry. */

    /* ######################################################################
            RIGHT COIL (at ymax). There are seven segments.
      ####################################################################### */
    // Segment #15. Time function LR.  Voltage in +Z.
    num = 15;
    foo = setAntennaEz(time_dep_LR*1.0 ,[x_coils(1,num)  , x_coils(2,num)  ,
                                         y_coils(1,num)  , y_coils(2,num)  ,
                                         z_coils(1,num)  , z_coils(2,num)-1],
                       ExField, EyField, EzField, BxField, ByField, BzField);
    // Segment #16. Time function LR.  Voltage in -Z.
    num = 16;
    foo = setAntennaEz(time_dep_LR*-1.0,[x_coils(1,num)  , x_coils(2,num)  ,
                                         y_coils(1,num)  , y_coils(2,num)  ,
                                         z_coils(1,num)  , z_coils(2,num)-1],
                       ExField, EyField, EzField, BxField, ByField, BzField);
    // Segment #17. Time function LR.  Voltage in -Z.
    num = 17;
    foo = setAntennaEz(time_dep_LR*-1.0,[x_coils(1,num)-1, x_coils(2,num)-1,
                                         y_coils(1,num)  , y_coils(2,num)  ,
                                         z_coils(1,num)  , z_coils(2,num)-1],
                       ExField, EyField, EzField, BxField, ByField, BzField);
    // Segment #18. Time function LR.  Voltage in +Z.
    num = 18;
    foo = setAntennaEz(time_dep_LR*1.0 ,[x_coils(1,num)-1, x_coils(2,num)-1,
                                         y_coils(1,num)  , y_coils(2,num)  ,
                                         z_coils(1,num)  , z_coils(2,num)-1],
                       ExField, EyField, EzField, BxField, ByField, BzField);
    // Segment #19. Time function LR.  Voltage in +X.
    num = 19;
    foo = setAntennaEx(time_dep_LR*1.0 ,[x_coils(1,num)-1, x_coils(2,num)  ,
                                         y_coils(1,num)  , y_coils(2,num)  ,
                                         z_coils(1,num)  , z_coils(2,num)  ],
                       ExField, EyField, EzField, BxField, ByField, BzField);
    // Segment #20. Time function LR.  Voltage doubled in -X.  Uses two cells!
    num = 20;
    foo = setAntennaEx(time_dep_LR*-2.0,[x_coils(1,num)-1, x_coils(2,num)  ,
                                         y_coils(1,num)  , y_coils(2,num)  ,
                                         z_coils(1,num)-1, z_coils(2,num)  ],
                       ExField, EyField, EzField, BxField, ByField, BzField);
    // Segment #21. Time function LR.  Voltage in +X.
    num = 21;
    foo = setAntennaEx(time_dep_LR*1.0 ,[x_coils(1,num)-1, x_coils(2,num)  ,
                                         y_coils(1,num)  , y_coils(2,num)  ,
                                         z_coils(1,num)-1, z_coils(2,num)-1],
                       ExField, EyField, EzField, BxField, ByField, BzField);
};

func applyBCsLeft(&ExField, &EyField, &EzField, &BxField, &ByField, &BzField){;
    
    /* In the following setAntenna calls, there are extra -1s on just some *_coils terms.
       This is because of the staggered nature of the Yee cell. Please ask me about it.
       Essentially, positive position values set the fields to slightly larger spatial
       values (half a cell in the direction under consideration).
       In addition, negative position values also set the fields to slightly larger
       spatial values, and hence we actually need to set the cell located one step back
       in the direction under consideration. This is how we will end up having symmetry. */

    /* ######################################################################
            LEFT COIL (at ymin). There are seven segments.
      ####################################################################### */
    // Segment #22. Time function LR.  Voltage in +Z.
    num = 22;
    foo = setAntennaEz(time_dep_LR*1.0 ,[x_coils(1,num)  , x_coils(2,num)  ,
                                         y_coils(1,num)-1, y_coils(2,num)-1,
                                         z_coils(1,num)  , z_coils(2,num)-1],
                       ExField, EyField, EzField, BxField, ByField, BzField);
    // Segment #23. Time function LR.  Voltage in -Z.
    num = 23;
    foo = setAntennaEz(time_dep_LR*-1.0,[x_coils(1,num)  , x_coils(2,num)  ,
                                         y_coils(1,num)-1, y_coils(2,num)-1,
                                         z_coils(1,num)  , z_coils(2,num)-1],
                       ExField, EyField, EzField, BxField, ByField, BzField);
    // Segment #24. Time function LR.  Voltage in -Z.
    num = 24;
    foo = setAntennaEz(time_dep_LR*-1.0,[x_coils(1,num)-1, x_coils(2,num)-1,
                                         y_coils(1,num)-1, y_coils(2,num)-1,
                                         z_coils(1,num)  , z_coils(2,num)-1],
                       ExField, EyField, EzField, BxField, ByField, BzField);
    // Segment #25. Time function LR.  Voltage in +Z.
    num = 25;
    foo = setAntennaEz(time_dep_LR*1.0 ,[x_coils(1,num)-1, x_coils(2,num)-1,
                                         y_coils(1,num)-1, y_coils(2,num)-1,
                                         z_coils(1,num)  , z_coils(2,num)-1],
                       ExField, EyField, EzField, BxField, ByField, BzField);
    // Segment #26. Time function LR.  Voltage in +X.
    num = 26;
    foo = setAntennaEx(time_dep_LR*1.0 ,[x_coils(1,num)-1, x_coils(2,num)  ,
                                         y_coils(1,num)-1, y_coils(2,num)-1,
                                         z_coils(1,num)  , z_coils(2,num)  ],
                       ExField, EyField, EzField, BxField, ByField, BzField);
    // Segment #27. Time function LR.  Voltage doubled in -X.  Uses two cells!
    num = 27;
    foo = setAntennaEx(time_dep_LR*-2.0,[x_coils(1,num)-1, x_coils(2,num)  ,
                                         y_coils(1,num)-1, y_coils(2,num)-1,
                                         z_coils(1,num)-1, z_coils(2,num)  ],
                       ExField, EyField, EzField, BxField, ByField, BzField);
    // Segment #28. Time function LR.  Voltage in +X.
    num = 28;
    foo = setAntennaEx(time_dep_LR*1.0 ,[x_coils(1,num)-1, x_coils(2,num)  ,
                                         y_coils(1,num)-1, y_coils(2,num)-1,
                                         z_coils(1,num)-1, z_coils(2,num)-1],
                       ExField, EyField, EzField, BxField, ByField, BzField);
};

/* ################################################################### */

// THE FOLLOWING FUNCTIONS SET INTERNAL IDEAL CONDUCTING BC FIELDS FOR THE FARADAY CAGE WIRES

func applyICsTop(&ExField, &EyField, &EzField, &BxField, &ByField, &BzField){;
    /* ######################################################################
            TOP CAGE (at xmax). There are four segments.
      ####################################################################### */
    num = 1;
    EzField(x_wires(1,num), y_wires(1,num), z_wires(1,num):z_wires(2,num)) = 0.0;
    num = 2;
    EzField(x_wires(1,num), y_wires(1,num), z_wires(1,num):z_wires(2,num)) = 0.0;
    num = 3;
    EyField(x_wires(1,num), y_wires(1,num):y_wires(2,num), z_wires(1,num)) = 0.0;
    num = 4;
    EyField(x_wires(1,num), y_wires(1,num):y_wires(2,num), z_wires(1,num)) = 0.0;
};

func applyICsBottom(&ExField, &EyField, &EzField, &BxField, &ByField, &BzField){;
    /* ######################################################################
            BOTTOM CAGE (at xmin). There are four segments.
      ####################################################################### */
    num = 5;
    EzField(x_wires(1,num), y_wires(1,num), z_wires(1,num):z_wires(2,num)) = 0.0;
    num = 6;
    EzField(x_wires(1,num), y_wires(1,num), z_wires(1,num):z_wires(2,num)) = 0.0;
    num = 7;
    EyField(x_wires(1,num), y_wires(1,num):y_wires(2,num), z_wires(1,num)) = 0.0;
    num = 8;
    EyField(x_wires(1,num), y_wires(1,num):y_wires(2,num), z_wires(1,num)) = 0.0;
};

func applyICsRight(&ExField, &EyField, &EzField, &BxField, &ByField, &BzField){;
    /* ######################################################################
            RIGHT CAGE (at ymax). There are seven segments.
      ####################################################################### */
   num = 9;
   EzField(x_wires(1,num), y_wires(1,num), z_wires(1,num):z_wires(2,num)) = 0.0;
   num = 10;
   EzField(x_wires(1,num), y_wires(1,num), z_wires(1,num):z_wires(2,num)) = 0.0;
   num = 11;
   ExField(x_wires(1,num):x_wires(2,num), y_wires(1,num), z_wires(1,num)) = 0.0;
   num = 12;
   ExField(x_wires(1,num):x_wires(2,num), y_wires(1,num), z_wires(1,num)) = 0.0;
};

func applyICsLeft(&ExField, &EyField, &EzField, &BxField, &ByField, &BzField){;
    /* ######################################################################
            LEFT CAGE (at ymin). There are seven segments.
      ####################################################################### */
   num = 13;
   EzField(x_wires(1,num), y_wires(1,num), z_wires(1,num):z_wires(2,num)) = 0.0;
   num = 14;
   EzField(x_wires(1,num), y_wires(1,num), z_wires(1,num):z_wires(2,num)) = 0.0;
   num = 15;
   ExField(x_wires(1,num):x_wires(2,num), y_wires(1,num), z_wires(1,num)) = 0.0;
   num = 16;
   ExField(x_wires(1,num):x_wires(2,num), y_wires(1,num), z_wires(1,num)) = 0.0;
};

/* ################################################################### */

func save_fields(timestep){
    s = swrite(format="%05.5hd", timestep);
    filename = sum(["output_data_wExt/pfrc_snapshot_",s,".hdf5"]);
    print, " .... SAVING DATA AS", filename;
    f = h5open(filename, "w");
    h5write, f, "/Coordinates/x/", x;
    h5write, f, "/Coordinates/y/", y;
    h5write, f, "/Coordinates/z/", z;

    Emag = sqrt(Ex(1:0,2:0,2:0)^2.0 + Ey(2:0,1:0,2:0)^2.0 + Ez(2:0,2:0,1:0)^2.0);
    Bmag = sqrt((1.0e4*Bx(2:0,1:0,1:0) + BxExt(2:0,1:0,1:0))^2.0 + (1.0e4*By(1:0,2:0,1:0) + ByExt(1:0,2:0,1:0))^2.0 + (1.0e4*Bz(1:0,1:0,2:0) + BzExt(1:0,1:0,2:0))^2.0);
    
    h5write, f, "/Field_quantities/Ex/", Ex(1:0,2:0,2:0);
    h5write, f, "/Field_quantities/Ey/", Ey(2:0,1:0,2:0);
    h5write, f, "/Field_quantities/Ez/", Ez(2:0,2:0,1:0);
    h5write, f, "/Field_quantities/Emag/", Emag;
    h5write, f, "/Field_quantities/Bx/", 1.0e4*Bx(2:0,1:0,1:0) + BxExt(2:0,1:0,1:0);
    h5write, f, "/Field_quantities/By/", 1.0e4*By(1:0,2:0,1:0) + ByExt(1:0,2:0,1:0);
    h5write, f, "/Field_quantities/Bz/", 1.0e4*Bz(1:0,1:0,2:0) + BzExt(1:0,1:0,2:0);
    h5write, f, "/Field_quantities/Bmag/", Bmag;

    h5close,f;
    h5off;          // not sure why I need this. Without it, hdf5 plugin seems to die
                    // after making 250 files. Something to do with the library.
        
    filename = sum(["/Users/asef/Desktop/pfrc_E_snapshot_",s,".vtk"]);
    f = open(filename, "w");
    write(f,format="%- .26s","# vtk DataFile Version 2.0"); // NO SPACE in front of #, required by Paraview. The format - sign does left justify. Otherwise adds a space and Paraview doesn't work.
    write,f,""; // We need a newline but yorick is quirky
    write,f,"Sample rectilinear grid";
    write,f,"ASCII";
    write,f,"DATASET RECTILINEAR_GRID";
    write,f,"DIMENSIONS",nx-1,ny-1,nz-1;
    write,f,"X_COORDINATES",nx-1,"float";
    write,f,x(2:0,1,1);
    write,f,"Y_COORDINATES",ny-1,"float";
    write,f,y(1,2:0,1);
    write,f,"Z_COORDINATES",nz-1,"float";
    write,f,z(1,1,2:0);
    write,f,"POINT_DATA",(nx-1)*(ny-1)*(nz-1);
    write,f,"SCALARS scalars float";
    write,f,"LOOKUP_TABLE default";
    write,f,Emag;
    write,f,"";
    write,f,"";
    write,f,"VECTORS vectors float";
    write,f,Ex(1:0,2:0,2:0), Ey(2:0,1:0,2:0), Ez(2:0,2:0,1:0);
    close,f;
    
    filename = sum(["/Users/asef/Desktop/pfrc_B_snapshot_",s,".vtk"]);
    f = open(filename, "w");
    write(f,format="%- .26s","# vtk DataFile Version 2.0"); // NO SPACE in front of #, required by Paraview. The format - sign does left justify. Otherwise adds a space and Paraview doesn't work.
    write,f,""; // We need a newline but yorick is quirky
    write,f,"Sample rectilinear grid";
    write,f,"ASCII";
    write,f,"DATASET RECTILINEAR_GRID";
    write,f,"DIMENSIONS",nx-1,ny-1,nz-1;
    write,f,"X_COORDINATES",nx-1,"float";
    write,f,x(2:0,1,1);
    write,f,"Y_COORDINATES",ny-1,"float";
    write,f,y(1,2:0,1);
    write,f,"Z_COORDINATES",nz-1,"float";
    write,f,z(1,1,2:0);
    write,f,"POINT_DATA",(nx-1)*(ny-1)*(nz-1);
    write,f,"SCALARS scalars float";
    write,f,"LOOKUP_TABLE default";
    write,f,Bmag;
    write,f,"";
    write,f,"";
    write,f,"VECTORS vectors float";
    write,f,1.0e4*Bx(2:0,1:0,1:0) + BxExt(2:0,1:0,1:0), 1.0e4*By(1:0,2:0,1:0) + ByExt(1:0,2:0,1:0), 1.0e4*Bz(1:0,1:0,2:0) + BzExt(1:0,1:0,2:0);
    close,f;
};

func save_objects{;
    object1 = array(0, nx, ny, nz);
    object2 = array(0, nx, ny, nz);
    object3 = array(0, nx, ny, nz);
    for (i=1;i<=nx;++i){;
        for (j=1;j<=ny;++j){;
            for (k=1;k<=nz;++k){;
                /* ################################################
                /* The following are for the RF antenna
                 ################################################# */
                if ( i >= x1_01 && i <= x2_01 && j >= y1_01 && j <= y2_01 && k >= z1_01 && k <= z2_01 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_02 && i <= x2_02 && j >= y1_02 && j <= y2_02 && k >= z1_02 && k <= z2_02 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_03 && i <= x2_03 && j >= y1_03 && j <= y2_03 && k >= z1_03 && k <= z2_03 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_04 && i <= x2_04 && j >= y1_04 && j <= y2_04 && k >= z1_04 && k <= z2_04 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_05 && i <= x2_05 && j >= y1_05 && j <= y2_05 && k >= z1_05 && k <= z2_05 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_06 && i <= x2_06 && j >= y1_06 && j <= y2_06 && k >= z1_06 && k <= z2_06 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_07 && i <= x2_07 && j >= y1_07 && j <= y2_07 && k >= z1_07 && k <= z2_07 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_08 && i <= x2_08 && j >= y1_08 && j <= y2_08 && k >= z1_08 && k <= z2_08 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_09 && i <= x2_09 && j >= y1_09 && j <= y2_09 && k >= z1_09 && k <= z2_09 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_10 && i <= x2_10 && j >= y1_10 && j <= y2_10 && k >= z1_10 && k <= z2_10 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_11 && i <= x2_11 && j >= y1_11 && j <= y2_11 && k >= z1_11 && k <= z2_11 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_12 && i <= x2_12 && j >= y1_12 && j <= y2_12 && k >= z1_12 && k <= z2_12 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_13 && i <= x2_13 && j >= y1_13 && j <= y2_13 && k >= z1_13 && k <= z2_13 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_14 && i <= x2_14 && j >= y1_14 && j <= y2_14 && k >= z1_14 && k <= z2_14 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_15 && i <= x2_15 && j >= y1_15 && j <= y2_15 && k >= z1_15 && k <= z2_15 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_16 && i <= x2_16 && j >= y1_16 && j <= y2_16 && k >= z1_16 && k <= z2_16 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_17 && i <= x2_17 && j >= y1_17 && j <= y2_17 && k >= z1_17 && k <= z2_17 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_18 && i <= x2_18 && j >= y1_18 && j <= y2_18 && k >= z1_18 && k <= z2_18 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_19 && i <= x2_19 && j >= y1_19 && j <= y2_19 && k >= z1_19 && k <= z2_19 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_20 && i <= x2_20 && j >= y1_20 && j <= y2_20 && k >= z1_20 && k <= z2_20 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_21 && i <= x2_21 && j >= y1_21 && j <= y2_21 && k >= z1_21 && k <= z2_21 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_22 && i <= x2_22 && j >= y1_22 && j <= y2_22 && k >= z1_22 && k <= z2_22 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_23 && i <= x2_23 && j >= y1_23 && j <= y2_23 && k >= z1_23 && k <= z2_23 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_24 && i <= x2_24 && j >= y1_24 && j <= y2_24 && k >= z1_24 && k <= z2_24 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_25 && i <= x2_25 && j >= y1_25 && j <= y2_25 && k >= z1_25 && k <= z2_25 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_26 && i <= x2_26 && j >= y1_26 && j <= y2_26 && k >= z1_26 && k <= z2_26 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_27 && i <= x2_27 && j >= y1_27 && j <= y2_27 && k >= z1_27 && k <= z2_27 ){;
                    object1(i,j,k) = 1;
                };
                if ( i >= x1_28 && i <= x2_28 && j >= y1_28 && j <= y2_28 && k >= z1_28 && k <= z2_28 ){;
                    object1(i,j,k) = 1;
                };
                /* ################################################
                /* The following are for the faraday cages
                 ################################################# */
                if ( i >= x1_w01 && i <= x2_w01 && j >= y1_w01 && j <= y2_w01 && k >= z1_w01 && k <= z2_w01 ){;
                    object2(i,j,k) = 1;
                };
                if ( i >= x1_w02 && i <= x2_w02 && j >= y1_w02 && j <= y2_w02 && k >= z1_w02 && k <= z2_w02 ){;
                    object2(i,j,k) = 1;
                };
                if ( i >= x1_w03 && i <= x2_w03 && j >= y1_w03 && j <= y2_w03 && k >= z1_w03 && k <= z2_w03 ){;
                    object2(i,j,k) = 1;
                };
                if ( i >= x1_w04 && i <= x2_w04 && j >= y1_w04 && j <= y2_w04 && k >= z1_w04 && k <= z2_w04 ){;
                    object2(i,j,k) = 1;
                };
                if ( i >= x1_w05 && i <= x2_w05 && j >= y1_w05 && j <= y2_w05 && k >= z1_w05 && k <= z2_w05 ){;
                    object2(i,j,k) = 1;
                };
                if ( i >= x1_w06 && i <= x2_w06 && j >= y1_w06 && j <= y2_w06 && k >= z1_w06 && k <= z2_w06 ){;
                    object2(i,j,k) = 1;
                };
                if ( i >= x1_w07 && i <= x2_w07 && j >= y1_w07 && j <= y2_w07 && k >= z1_w07 && k <= z2_w07 ){;
                    object2(i,j,k) = 1;
                };
                if ( i >= x1_w08 && i <= x2_w08 && j >= y1_w08 && j <= y2_w08 && k >= z1_w08 && k <= z2_w08 ){;
                    object2(i,j,k) = 1;
                };
                if ( i >= x1_w09 && i <= x2_w09 && j >= y1_w09 && j <= y2_w09 && k >= z1_w09 && k <= z2_w09 ){;
                    object2(i,j,k) = 1;
                };
                if ( i >= x1_w10 && i <= x2_w10 && j >= y1_w10 && j <= y2_w10 && k >= z1_w10 && k <= z2_w10 ){;
                    object2(i,j,k) = 1;
                };
                if ( i >= x1_w11 && i <= x2_w11 && j >= y1_w11 && j <= y2_w11 && k >= z1_w11 && k <= z2_w11 ){;
                    object2(i,j,k) = 1;
                };
                if ( i >= x1_w12 && i <= x2_w12 && j >= y1_w12 && j <= y2_w12 && k >= z1_w12 && k <= z2_w12 ){;
                    object2(i,j,k) = 1;
                };
                if ( i >= x1_w13 && i <= x2_w13 && j >= y1_w13 && j <= y2_w13 && k >= z1_w13 && k <= z2_w13 ){;
                    object2(i,j,k) = 1;
                };
                if ( i >= x1_w14 && i <= x2_w14 && j >= y1_w14 && j <= y2_w14 && k >= z1_w14 && k <= z2_w14 ){;
                    object2(i,j,k) = 1;
                };
                if ( i >= x1_w15 && i <= x2_w15 && j >= y1_w15 && j <= y2_w15 && k >= z1_w15 && k <= z2_w15 ){;
                    object2(i,j,k) = 1;
                };
                if ( i >= x1_w16 && i <= x2_w16 && j >= y1_w16 && j <= y2_w16 && k >= z1_w16 && k <= z2_w16 ){;
                    object2(i,j,k) = 1;
                };
                /* ################################################
                /* The following are for the flux conserving rings
                 ################################################# */
                xVal = x(i,j,k);
                yVal = y(i,j,k);
                zVal = z(i,j,k);

                // One vs the other below must have to do with how the ring conforms to cartesian
              // for nx = ny = 101
//                rVal = sqrt( (xVal+2.*dx/2.)^2 + (yVal+2.*dy/2.)^2 );
              // for nx = ny = 201
//                rVal = sqrt( (xVal+3.*dx/2.)^2 + (yVal+3.*dy/2.)^2 );
                // for nx = ny = 401
                  rVal = sqrt( (xVal+5.*dx/2.)^2 + (yVal+5.*dy/2.)^2 );

                for (nn =1; nn<=6; ++nn){;
                    if ( (rVal >= fc_R(1,nn)) && (rVal <= fc_R(2, nn)) ){;
                        if ( (zVal >= fc_Z(1,nn)-dz/2.) && (zVal <= fc_Z(2, nn)-dz/2.) ){;
                            object3(i,j,k) = 1;
                        };
                    };
                };
            }; // end of z loop
        }; // end of y loop
    }; // end of x loop

    filename = "pfrc_obj_antenna.hdf5";
    print, " .... SAVING DATA AS", filename;
    f = h5open(filename, "w");
    h5write, f, "/Coordinates/x/", x;
    h5write, f, "/Coordinates/y/", y;
    h5write, f, "/Coordinates/z/", z;
    h5write, f, "/Objects/", object1;
    h5close,f;

    filename = "pfrc_obj_faradaycage.hdf5";
    print, " .... SAVING DATA AS", filename;
    f = h5open(filename, "w");
    h5write, f, "/Coordinates/x/", x;
    h5write, f, "/Coordinates/y/", y;
    h5write, f, "/Coordinates/z/", z;
    h5write, f, "/Objects/", object2;
    h5close,f;

    /*
    filename = "pfrc_obj_fluxrings.hdf5";
    print, " .... SAVING DATA AS", filename;
    f = h5open(filename, "w");
    h5write, f, "/Coordinates/x/", x;
    h5write, f, "/Coordinates/y/", y;
    h5write, f, "/Coordinates/z/", z;
    h5write, f, "/Objects/", object3;
    h5close,f;
    */
    h5off;          // not sure why I need this. Without it, hdf5 plugin seems to die
                    // after making 250 files. Something to do with the library.
    print, "sum(object)=",sum(object);
};
