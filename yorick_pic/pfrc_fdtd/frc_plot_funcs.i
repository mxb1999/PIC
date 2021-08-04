/* ################################################################
    Implementation of the FRC antenna
    and explicit Yee method for Maxwell's equations

    This file contains functions related to plotting.
    Version 1.0
    Date: July 3rd, 2021
 ################################################################### */

print, "LOADED 'frc_plot_funcs.i'!"

/* ################################################################### */

func colorbar(cmin, cmax){;
    plsys, 0;
    pli, span(0,1,200)(-,), .625, .46, .67, .84, legend="";
    plg, [.46,.84,.84,.46],[.67,.67,.625,.625], closed=1, marks=0, color="fg", width=1, type=1, legend="";
    plsys, 1;
    if (!is_void(cmin)){;
        c1=(cmax+cmin)/2.;      c2=(cmax+c1)/2.;
        c3=(cmax+c2)/2.;        c4=(c1+c2)/2.;        c5=(cmin+c1)/2.;
        c6=(cmin+c5)/2.;        c7=(c1+c5)/2.;        loc1=(0.84+0.46)/2.0;
        loc2=(0.84+loc1)/2.0;        loc3=(0.84+loc2)/2.0;        loc4=(loc1+loc2)/2.0;
        loc5=(0.46+loc1)/2.0;        loc6=(0.46+loc5)/2.0;        loc7=(loc1+loc5)/2.0;
        plt, pr1(cmin), .6475, .46, justify="LH";
        plt, pr1(cmax), .6475, .84, justify="LH";
        plt, pr1(c1), 0.6475, loc1, justify="LH";
        plt, pr1(c2), 0.6475, loc2, justify="LH";
        plt, pr1(c3), 0.6475, loc3, justify="LH";
        plt, pr1(c4), 0.6475, loc4, justify="LH";
        plt, pr1(c5), 0.6475, loc5, justify="LH";
        plt, pr1(c6), 0.6475, loc6, justify="LH";
        plt, pr1(c7), 0.6475, loc7, justify="LH";
    };
};

/* ################################################################### */

func plotThePhases{;
    // This simply plots the phases of the two sets of coils.
    nt=1001;
    time = span(0.0, 3*T, nt);
    amplitude1 = array(0.0, nt);    // one for TOP & BOTTOM
    amplitude2 = array(0.0, nt);    // one for RIGHT & LEFT
    for (tt=1;tt<=nt;++tt){;
        amplitude1(tt) = 140.*riseTime1(time(tt))*function3(time(tt));
        amplitude2(tt) = 140.*riseTime7(time(tt))*function2(time(tt));
    };
    window,4;fma;
    plg, amplitude1, 1e9*time, marks=0, width=5, color="black";
    plg, amplitude2, 1e9*time, marks=0, width=5, color="red";
    xytitles,"Time (ns)", "Amplitude (V)";
    pltitle,"Phases of the antennae";
    gridxy,1,1;
};

/* ################################################################### */

func plotTheCoils{;
    /* ########  Plot TOP COIL ########### */
    window,0;fma;
    plg, [coil_Y(2,1), coil_Y(1,1)], [coil_Z(2,1), coil_Z(1,1)], marks=0, width=7, color="black", type="solid";
    plg, [coil_Y(2,3), coil_Y(1,3)], [coil_Z(2,3), coil_Z(1,3)], marks=0, width=7, color="black", type="dash";
    plg, [coil_Y(2,5), coil_Y(1,5)], [coil_Z(2,5), coil_Z(1,5)], marks=0, width=7, color="black", type="solid";

    plg, [coil_Y(2,6), coil_Y(1,6)], [coil_Z(2,6), coil_Z(1,6)], marks=0, width=7, color="blue", type="dash";

    plg, [coil_Y(2,2), coil_Y(1,2)], [coil_Z(2,2), coil_Z(1,2)], marks=0, width=7, color="red", type="dash";
    plg, [coil_Y(2,4), coil_Y(1,4)], [coil_Z(2,4), coil_Z(1,4)], marks=0, width=7, color="red", type="solid";
    plg, [coil_Y(2,7), coil_Y(1,7)], [coil_Z(2,7), coil_Z(1,7)], marks=0, width=7, color="red", type="solid";
    gridxy,1,1;
    xytitles,"Z (m)", "Y (m)";
    pltitle, "COIL #1   (X location = +7.5 cm)";
    limits,-0.15,0.15,-0.15,0.15;

    /* ########  Plot BOTTOM COIL ########### */
    window,1;fma;
    plg, [coil_Y(2,8), coil_Y(1,8)], [coil_Z(2,8), coil_Z(1,8)], marks=0, width=7, color="black", type="solid";
    plg, [coil_Y(2,10), coil_Y(1,10)], [coil_Z(2,10), coil_Z(1,10)], marks=0, width=7, color="black", type="dash";
    plg, [coil_Y(2,12), coil_Y(1,12)], [coil_Z(2,12), coil_Z(1,12)], marks=0, width=7, color="black", type="solid";

    plg, [coil_Y(2,13), coil_Y(1,13)], [coil_Z(2,13), coil_Z(1,13)], marks=0, width=7, color="blue", type="dash";

    plg, [coil_Y(2,9), coil_Y(1,9)], [coil_Z(2,9), coil_Z(1,9)], marks=0, width=7, color="red", type="dash";
    plg, [coil_Y(2,11), coil_Y(1,11)], [coil_Z(2,11), coil_Z(1,11)], marks=0, width=7, color="red", type="solid";
    plg, [coil_Y(2,14), coil_Y(1,14)], [coil_Z(2,14), coil_Z(1,14)], marks=0, width=7, color="red", type="solid";
    gridxy,1,1;
    xytitles,"Z (m)", "Y (m)";
    pltitle, "COIL #2   (X location = -7.5 cm)";
    limits,-0.15,0.15,-0.15,0.15;

    /* ########  Plot RIGHT COIL ########### */
    window,2;fma;
    plg, [coil_X(2,15), coil_X(1,15)], [coil_Z(2,15), coil_Z(1,15)], marks=0, width=7, color="black", type="solid";
    plg, [coil_X(2,17), coil_X(1,17)], [coil_Z(2,17), coil_Z(1,17)], marks=0, width=7, color="black", type="dash";
    plg, [coil_X(2,19), coil_X(1,19)], [coil_Z(2,19), coil_Z(1,19)], marks=0, width=7, color="black", type="solid";

    plg, [coil_X(2,20), coil_X(1,20)], [coil_Z(2,20), coil_Z(1,20)], marks=0, width=7, color="blue", type="dash";

    plg, [coil_X(2,16), coil_X(1,16)], [coil_Z(2,16), coil_Z(1,16)], marks=0, width=7, color="red", type="dash";
    plg, [coil_X(2,18), coil_X(1,18)], [coil_Z(2,18), coil_Z(1,18)], marks=0, width=7, color="red", type="solid";
    plg, [coil_X(2,21), coil_X(1,21)], [coil_Z(2,21), coil_Z(1,21)], marks=0, width=7, color="red", type="solid";
    gridxy,1,1;
    xytitles,"Z (m)", "Y (m)";
    pltitle, "COIL #3   (Y location = +7.5 cm)";
    limits,-0.15,0.15,-0.15,0.15;

    /* ########  Plot LEFT COIL ########### */
    window,3;fma;
    plg, [coil_X(2,22), coil_X(1,22)], [coil_Z(2,22), coil_Z(1,22)], marks=0, width=7, color="black", type="solid";
    plg, [coil_X(2,24), coil_X(1,24)], [coil_Z(2,24), coil_Z(1,24)], marks=0, width=7, color="black", type="dash";
    plg, [coil_X(2,26), coil_X(1,26)], [coil_Z(2,26), coil_Z(1,26)], marks=0, width=7, color="black", type="solid";

    plg, [coil_X(2,27), coil_X(1,27)], [coil_Z(2,27), coil_Z(1,27)], marks=0, width=7, color="blue", type="dash";

    plg, [coil_X(2,23), coil_X(1,23)], [coil_Z(2,23), coil_Z(1,23)], marks=0, width=7, color="red", type="dash";
    plg, [coil_X(2,25), coil_X(1,25)], [coil_Z(2,25), coil_Z(1,25)], marks=0, width=7, color="red", type="solid";
    plg, [coil_X(2,28), coil_X(1,28)], [coil_Z(2,28), coil_Z(1,28)], marks=0, width=7, color="red", type="solid";
    gridxy,1,1;
    xytitles,"Z (m)", "Y (m)";
    pltitle, "COIL #4   (Y location = -7.5 cm)";
    limits,-0.15,0.15,-0.15,0.15;
};

/* ################################################################### */

func myplot2d(pausetime){;
// The flags below set whether to plot the component:
    if_win_06 = 1;  // Ex(y,z)
    if_win_07 = 1;  // Ex(x,z)
    if_win_08 = 1;  // Ex(x,y)
    if_win_09 = 1;  // Ey(y,z)
    if_win_10 = 1;  // Ey(x,z)
    if_win_11 = 1;  // Ey(x,y)
    if_win_12 = 1;  // Ez(y,z)
    if_win_13 = 1;  // Ez(x,z)
    if_win_14 = 1;  // Ez(x,y)
    if_win_15 = 1;  // Bx(y,z)
    if_win_16 = 1;  // Bx(x,z)
    if_win_17 = 1;  // Bx(x,y)
    if_win_18 = 1;  // By(y,z)
    if_win_19 = 1;  // By(x,z)
    if_win_20 = 1;  // By(x,y)
    if_win_21 = 1;  // Bz(y,z)
    if_win_22 = 1;  // Bz(x,z)
    if_win_23 = 1;  // Bz(x,y)
    if_win_24 = 0;  // Ex(z) and Ey(z)
    if_win_25 = 0;  // Etheta(x,y)
    if_win_26 = 0;  // Er(x,y)
    if_win_27 = 0;  // sqrt(Er(x,y)^2 + Etheta(x,y)^2)
    
    
    if ( if_win_06 == 1 ){;
    window,6;fma;
    xx = int((nx+1)/2);
    parameter = Ex(xx, 1:ny-1, 1:nz-1);
    palette,"hot_desaturated.gp";//palette,"idl33.gp"; // palette,"hot_desaturated.gp";
    clo=min(parameter); chi=max(parameter);
    plf, parameter, y(xx,1:ny,1:nz)-dy/2, z(xx,1:ny,1:nz)-dz/2, cmin=clo, cmax=chi;
    plg, [0.15,-0.15],[-0.14,-0.14],marks=0,width=3,type="dash";
    plg, [0.15,-0.15],[ 0.00, 0.00],marks=0,width=3,type="dash";
    plg, [0.15,-0.15],[ 0.14, 0.14],marks=0,width=3,type="dash";
    plg, [0.075,0.075],[-0.20,0.20],marks=0,width=3,type="dash";
    plg, [-0.075,-0.075],[-0.20,0.20],marks=0,width=3,type="dash";
    colorbar, clo, chi;
    xytitles,"Z (m)", "Y (m)";
    title=strtrim(swrite("E_x_ at x=", x(xx,1,1),"at t=",time*1e9,"ns"));
    pltitle, title;
    limits,zmin, zmax, zmin, zmax;
    };
    
    if ( if_win_07 == 1 ){;
    window,7;fma;
    yy = (ny+1)/2;
    parameter = Ex(1:nx-1,yy,1:nz-1);
    palette,"hot_desaturated.gp";//palette,"idl33.gp"; // palette,"hot_desaturated.gp";
    clo=min(parameter); chi=max(parameter);
    plf, parameter, x(1:nx,yy,1:nz), z(1:nx,yy,1:nz)-dz/2, cmin=clo, cmax=chi;
    plg, [0.15,-0.15],[-0.14,-0.14],marks=0,width=3,type="dash";
    plg, [0.15,-0.15],[ 0.00, 0.00],marks=0,width=3,type="dash";
    plg, [0.15,-0.15],[ 0.14, 0.14],marks=0,width=3,type="dash";
    plg, [0.075,0.075],[-0.20,0.20],marks=0,width=3,type="dash";
    plg, [-0.075,-0.075],[-0.20,0.20],marks=0,width=3,type="dash";
    colorbar, clo, chi;
    xytitles,"Z (m)", "X (m)";
    title=strtrim(swrite("E_x_ at y=", y(1,yy,1),"at t=",time*1e9,"ns"));
    pltitle, title;
    limits,zmin, zmax, zmin, zmax;
    };

    if ( if_win_08 == 1 ){;
    window,8;fma;
    zz = int((nz+1)/2);
    parameter = Ex(1:nx-1,1:ny,zz);
    palette,"hot_desaturated.gp";
    clo=min(parameter); chi=max(parameter);
    plf, parameter, y(1:nx-1,1:ny,zz)+dy/2, x(1:nx-1,1:ny,zz)+dx, cmin=clo, cmax=chi;
    theta=span(0.0,2*pi,1000);
    plg, 0.044*sin(theta), 0.044*cos(theta), marks=0, width=3, color="black";
    plg, 0.05*sin(theta), 0.05*cos(theta), marks=0, width=3, color="black";
    plg, [0.075,0.075],[-0.11,0.11],marks=0,width=3,type="dash";
    plg, [-0.075,-0.075],[-0.11,0.11],marks=0,width=3,type="dash";
    plg, [0.11,-0.11],[0.075,0.075],marks=0,width=3,type="dash";
    plg, [0.11,-0.11],[-0.075,-0.075],marks=0,width=3,type="dash";
    theta=span(0.0,2*pi,1000);
    plg, 0.044*sin(theta), 0.044*cos(theta), marks=0, width=3, color="white";
    plg, 0.05*sin(theta), 0.05*cos(theta), marks=0, width=3, color="white";
    plg,[0.1,-0.1],[0.,0.], marks=0,width=1,type="dash", color="white";
    plg,[0.,0.],[-0.1,0.1], marks=0,width=1,type="dash", color="white";

    colorbar, clo, chi;
    xytitles,"X (m)", "Y (m)";
    title=strtrim(swrite("E_x_ at z=", z(1,1,zz),"at t=",time*1e9,"ns"));
    pltitle, title;
    limits//,xmin, xmax, xmin, xmax;
    };

    if ( if_win_09 == 1 ){;
    window,9;fma;
    xx = int((nx+1)/2);
    parameter = Ey(xx, 1:ny-1, 1:nz);
    palette,"hot_desaturated.gp";//palette,"idl33.gp"; // palette,"hot_desaturated.gp";
    clo=min(parameter); chi=max(parameter);
    plf, parameter, y(xx,1:ny-1,1:nz)+dy, z(xx,1:ny-1,1:nz)+dz/2, cmin=clo, cmax=chi;
    plg, [0.15,-0.15],[-0.14,-0.14],marks=0,width=3,type="dash";
    plg, [0.15,-0.15],[ 0.00, 0.00],marks=0,width=3,type="dash";
    plg, [0.15,-0.15],[ 0.14, 0.14],marks=0,width=3,type="dash";
    plg, [0.075,0.075],[-0.20,0.20],marks=0,width=3,type="dash";
    plg, [-0.075,-0.075],[-0.20,0.20],marks=0,width=3,type="dash";
    colorbar, clo, chi;
    xytitles,"Z (m)", "Y (m)";
    title=strtrim(swrite("E_y_ at x=", x(xx,1,1),"at t=",time*1e9,"ns"));
    pltitle, title;
    limits,zmin, zmax, zmin, zmax;
    };

    if ( if_win_10 == 1 ){;
    window,10;fma;
    yy = (ny+1)/2;
    parameter = Ey(1:nx-1,yy,1:nz-1);
    palette,"hot_desaturated.gp";//palette,"idl33.gp"; // palette,"hot_desaturated.gp";
    clo=min(parameter); chi=max(parameter);
    plf, parameter, x(1:nx,yy,1:nz)-dx/2, z(1:nx,yy,1:nz)-dz/2, cmin=clo, cmax=chi;
    plg, [0.15,-0.15],[-0.14,-0.14],marks=0,width=3,type="dash";
    plg, [0.15,-0.15],[ 0.00, 0.00],marks=0,width=3,type="dash";
    plg, [0.15,-0.15],[ 0.14, 0.14],marks=0,width=3,type="dash";
    plg, [0.075,0.075],[-0.20,0.20],marks=0,width=3,type="dash";
    plg, [-0.075,-0.075],[-0.20,0.20],marks=0,width=3,type="dash";
    colorbar, clo, chi;
    xytitles,"Z (m)", "X (m)";
    title=strtrim(swrite("E_y_ at y=", y(1,yy,1),"at t=",time*1e9,"ns"));
    pltitle, title;
    limits,zmin, zmax, zmin, zmax;
    };

    if ( if_win_11 == 1 ){;
    window,11;fma;
    zz = int((nz+1)/2);
    parameter = Ey(1:nx,1:ny-1,zz);
    palette,"hot_desaturated.gp";//palette,"idl33.gp"; // palette,"hot_desaturated.gp";
    clo=min(parameter); chi=max(parameter);
    plf, parameter, y(1:nx,1:ny-1,zz)+dy, x(1:nx,1:ny-1,zz)+dx/2, cmin=clo, cmax=chi;
    plg, [0.075,0.075],[-0.11,0.11],marks=0,width=3,type="dash";
    plg, [-0.075,-0.075],[-0.11,0.11],marks=0,width=3,type="dash";
    plg, [0.11,-0.11],[0.075,0.075],marks=0,width=3,type="dash";
    plg, [0.11,-0.11],[-0.075,-0.075],marks=0,width=3,type="dash";
        theta=span(0.0,2*pi,1000);
        plg, 0.044*sin(theta), 0.044*cos(theta), marks=0, width=3, color="white";
        plg, 0.05*sin(theta), 0.05*cos(theta), marks=0, width=3, color="white";
    plg,[0.1,-0.1],[0.,0.], marks=0,width=1,type="dash", color="white";
    plg,[0.,0.],[-0.1,0.1], marks=0,width=1,type="dash", color="white";
    colorbar, clo, chi;
    xytitles,"X (m)", "Y (m)";
    title=strtrim(swrite("E_y_ at z=", z(1,1,zz),"at t=",time*1e9,"ns"));
    pltitle, title;
    limits//,xmin, xmax, xmin, xmax;
    };

    if ( if_win_12 == 1 ){;
    window,12;fma;
    palette,"hot_desaturated.gp";//palette,"idl33.gp"; // palette,"hot_desaturated.gp";
    xx = int((nx+1)/2);
    parameter = Ez(xx, 1:ny-1, 1:nz-1);
    clo=min(parameter); chi=max(parameter);
    plf, parameter, y(xx,1:ny-1,1:nz-1)+dy/2, z(xx,1:ny-1,1:nz-1)+dz, cmin=clo, cmax=chi;
    plg, [0.15,-0.15],[-0.14,-0.14],marks=0,width=3,type="dash";
    plg, [0.15,-0.15],[ 0.00, 0.00],marks=0,width=3,type="dash";
    plg, [0.15,-0.15],[ 0.14, 0.14],marks=0,width=3,type="dash";
    plg, [0.075,0.075],[-0.20,0.20],marks=0,width=3,type="dash";
    plg, [-0.075,-0.075],[-0.20,0.20],marks=0,width=3,type="dash";
    colorbar, clo, chi;
    xytitles,"Z (m)", "Y (m)";
    title=strtrim(swrite("E_z_ at x=", x(xx,1,1),"at t=",time*1e9,"ns"));
    pltitle, title;
    limits,zmin, zmax, zmin, zmax;
    };

    if ( if_win_13 == 1 ){;
    window,13;fma;
    palette,"hot_desaturated.gp";//palette,"idl33.gp"; // palette,"hot_desaturated.gp";
    yy = (ny+1)/2;
    parameter = Ez(1:nx, yy, 1:nz-1);
    clo=min(parameter); chi=max(parameter);
    plf, parameter, x(1:nx,yy,2:nz)+dx/2, z(1:nx,yy,2:nz), cmin=clo, cmax=chi;
    plg, [0.15,-0.15],[-0.14,-0.14],marks=0,width=3,type="dash";
    plg, [0.15,-0.15],[ 0.00, 0.00],marks=0,width=3,type="dash";
    plg, [0.15,-0.15],[ 0.14, 0.14],marks=0,width=3,type="dash";
    plg, [0.075,0.075],[-0.20,0.20],marks=0,width=3,type="dash";
    plg, [-0.075,-0.075],[-0.20,0.20],marks=0,width=3,type="dash";
    colorbar, clo, chi;
    xytitles,"Z (m)", "X (m)";
    title=strtrim(swrite("E_z_ at y=", y(1,yy,1),"at t=",time*1e9,"ns"));
    pltitle, title;
    limits,zmin, zmax, zmin, zmax;
    };

    if ( if_win_14 == 1 ){;
    window,14;fma;
    palette,"hot_desaturated.gp";//palette,"idl33.gp"; // palette,"hot_desaturated.gp";
    zz = (nz+1)/2 ;//- int((0.14/dz));
    parameter = Ez(1:nx-1, 1:ny-1, zz);
    clo=min(parameter); chi=max(parameter);
    plf, parameter, y(1:nx-1,1:ny-1,zz)+dy, x(1:nx-1,1:ny-1,zz)+dx, cmin=clo, cmax=chi;
    plg, [0.075,0.075],[-0.11,0.11],marks=0,width=3,type="dash";
    plg, [-0.075,-0.075],[-0.11,0.11],marks=0,width=3,type="dash";
    plg, [0.11,-0.11],[0.075,0.075],marks=0,width=3,type="dash";
    plg, [0.11,-0.11],[-0.075,-0.075],marks=0,width=3,type="dash";
        theta=span(0.0,2*pi,1000);
        plg, 0.044*sin(theta), 0.044*cos(theta), marks=0, width=3, color="white";
        plg, 0.05*sin(theta), 0.05*cos(theta), marks=0, width=3, color="white";
    colorbar, clo, chi;
    xytitles,"X (m)", "Y (m)";
    title=strtrim(swrite("E_z_ at z=", z(1,1,zz),"at t=",time*1e9,"ns"));
    pltitle, title;
    limits;
    };

    if ( if_win_15 == 1 ){;
    window,15;fma;
    palette,"hot_desaturated.gp";//palette,"idl33.gp"; // palette,"hot_desaturated.gp";
    xx = (nx+1)/2;
    parameter = 1.0e4*Bx(xx, 1:ny-1, 1:nz-1) + BxExt(xx, 1:ny-1, 1:nz-1);
    clo=min(parameter); chi=max(parameter);
    plf, parameter, y(xx,1:ny,1:nz)+dy/2, z(xx,1:ny,1:nz), cmin=clo, cmax=chi;
    plg, [0.15,-0.15],[-0.14,-0.14],marks=0,width=3,type="dash";
    plg, [0.15,-0.15],[ 0.00, 0.00],marks=0,width=3,type="dash";
    plg, [0.15,-0.15],[ 0.14, 0.14],marks=0,width=3,type="dash";
    plg, [0.075,0.075],[-0.20,0.20],marks=0,width=3,type="dash";
    plg, [-0.075,-0.075],[-0.20,0.20],marks=0,width=3,type="dash";
    colorbar, clo, chi;
    xytitles,"Z (m)", "Y (m)";
    title=strtrim(swrite("B_x_ (G) at x=", x(xx,1,1),"at t=",time*1e9,"ns"));
    pltitle, title;
    limits,zmin, zmax, zmin, zmax;
    };

    if ( if_win_16 == 1 ){;
    window,16;fma;
    palette,"hot_desaturated.gp";//palette,"idl33.gp"; // palette,"hot_desaturated.gp";
    yy = (ny+1)/2;
    parameter = 1.0e4*Bx(1:nx, yy, 1:nz-1) + BxExt(1:nx, yy, 1:nz-1);
    clo=min(parameter); chi=max(parameter);
    plf, parameter, x(1:nx, yy, 1:nz-1)+dx, z(1:nx, yy, 1:nz-1)+dz, cmin=clo, cmax=chi;
    plg, [0.15,-0.15],[-0.14,-0.14],marks=0,width=3,type="dash";
    plg, [0.15,-0.15],[ 0.00, 0.00],marks=0,width=3,type="dash";
    plg, [0.15,-0.15],[ 0.14, 0.14],marks=0,width=3,type="dash";
    plg, [0.075,0.075],[-0.20,0.20],marks=0,width=3,type="dash";
    plg, [-0.075,-0.075],[-0.20,0.20],marks=0,width=3,type="dash";
    colorbar, clo, chi;
    xytitles,"Z (m)", "X (m)";
    title=strtrim(swrite("B_x_ (G) at y=", y(1,yy,1),"at t=",time*1e9,"ns"));
    pltitle, title;
    limits,zmin, zmax, zmin, zmax;
    };

    if ( if_win_17 == 1 ){;
    window,17;fma;
    palette,"hot_desaturated.gp";//palette,"idl33.gp"; // palette,"hot_desaturated.gp";
    zz = (nz+1)/2;
    parameter = 1.0e4*Bx(1:nx, 1:ny-1, zz) + BxExt(1:nx, 1:ny-1, zz);
    clo=min(parameter); chi=max(parameter);
    plf, parameter, y(1:nx, 1:ny-1, zz)+3.0*dy/2, x(1:nx, 1:ny-1, zz)+dx, cmin=clo, cmax=chi;
    plg, [0.075,0.075],[-0.11,0.11],marks=0,width=3,type="dash";
    plg, [-0.075,-0.075],[-0.11,0.11],marks=0,width=3,type="dash";
    plg, [0.11,-0.11],[0.075,0.075],marks=0,width=3,type="dash";
    plg, [0.11,-0.11],[-0.075,-0.075],marks=0,width=3,type="dash";
    theta=span(0.0,2*pi,1000);
    plg, 0.044*sin(theta), 0.044*cos(theta), marks=0, width=3, color="white";
    plg, 0.05*sin(theta), 0.05*cos(theta), marks=0, width=3, color="white";
    colorbar, clo, chi;
    xytitles,"X (m)", "Y (m)";
    title=strtrim(swrite("B_x_ (G) at z=", z(1,1,zz),"at t=",time*1e9,"ns"));
    pltitle, title;
    limits,xmin, xmax, xmin, xmax;
    };
    
    if ( if_win_18 == 1 ){;
    window,18;fma;
    xx = (nx+1)/2;
    parameter = 1.0e4*By(xx, 1:ny, 1:nz-1) + ByExt(xx, 1:ny, 1:nz-1);
    palette,"hot_desaturated.gp";//palette,"idl33.gp"; // palette,"hot_desaturated.gp";
    clo=min(parameter); chi=max(parameter);
    plf, parameter, y(xx,1:ny,1:nz-1)+dy, z(xx,1:ny,1:nz-1)+dz, cmin=clo, cmax=chi;
    plg, [0.15,-0.15],[-0.14,-0.14],marks=0,width=3,type="dash";
    plg, [0.15,-0.15],[ 0.00, 0.00],marks=0,width=3,type="dash";
    plg, [0.15,-0.15],[ 0.14, 0.14],marks=0,width=3,type="dash";
    plg, [0.075,0.075],[-0.20,0.20],marks=0,width=3,type="dash";
    plg, [-0.075,-0.075],[-0.20,0.20],marks=0,width=3,type="dash";
    colorbar, clo, chi;
    xytitles,"Z (m)", "Y (m)";
    title=strtrim(swrite("B_y_ (G) at x=", x(xx,1,1),"at t=",time*1e9,"ns"));
    pltitle, title;
    limits,zmin, zmax, zmin, zmax;
    };

    if ( if_win_19 == 1 ){;
    window,19;fma;
    palette,"hot_desaturated.gp";//palette,"idl33.gp"; // palette,"hot_desaturated.gp";
        yy = (ny+1)/2;
    parameter = 1.0e4*By(1:nx-1, yy, 1:nz-1) + ByExt(1:nx-1, yy, 1:nz-1);
    clo=min(parameter); chi=max(parameter);
    plf, parameter, x(1:nx-1,(ny+1)/2,1:nz-1)+3.0*dx/2, z(1:nx-1,(ny+1)/2,1:nz-1)+dz, cmin=clo, cmax=chi;
    plg, [0.15,-0.15],[-0.14,-0.14],marks=0,width=3,type="dash";
    plg, [0.15,-0.15],[ 0.00, 0.00],marks=0,width=3,type="dash";
    plg, [0.15,-0.15],[ 0.14, 0.14],marks=0,width=3,type="dash";
    plg, [0.075,0.075],[-0.20,0.20],marks=0,width=3,type="dash";
    plg, [-0.075,-0.075],[-0.20,0.20],marks=0,width=3,type="dash";
    colorbar, clo, chi;
    xytitles,"Z (m)", "X (m)";
    title=strtrim(swrite("B_y_ (G) at y=", y(1,yy,1),"at t=",time*1e9,"ns"));
    pltitle, title;
    limits,zmin, zmax, zmin, zmax;
    };

    if ( if_win_20 == 1 ){;
    window,20;fma;
    palette,"hot_desaturated.gp";//palette,"idl33.gp"; // palette,"hot_desaturated.gp";
    zz = (nz+1)/2;
    parameter = 1.0e4*By(1:nx-1, 1:ny, zz) + ByExt(1:nx-1, 1:ny, zz);
    clo=min(parameter); chi=max(parameter);
    plf, parameter, y(1:nx-1,1:ny,zz)+dy, x(1:nx-1,1:ny,zz)+3.0*dx/2, cmin=clo, cmax=chi;
    plg, [0.075,0.075],[-0.11,0.11],marks=0,width=3,type="dash";
    plg, [-0.075,-0.075],[-0.11,0.11],marks=0,width=3,type="dash";
    plg, [0.11,-0.11],[0.075,0.075],marks=0,width=3,type="dash";
    plg, [0.11,-0.11],[-0.075,-0.075],marks=0,width=3,type="dash";
    theta=span(0.0,2*pi,1000);
    plg, 0.044*sin(theta), 0.044*cos(theta), marks=0, width=3, color="white";
    plg, 0.05*sin(theta), 0.05*cos(theta), marks=0, width=3, color="white";
    colorbar, clo, chi;
    xytitles,"X (m)", "Y (m)";
    title=strtrim(swrite("B_y_ (G) at z=", z(1,1,zz),"at t=",time*1e9,"ns"));
    pltitle, title;
    limits,xmin, xmax, xmin, xmax;
    };
    
    if ( if_win_21 == 1 ){;
    window,21;fma;
    xx = (nx+1)/2;
    parameter = 1.0e4*Bz(xx, 1:ny-1, 1:nz) + BzExt(xx, 1:ny-1, 1:nz);
    palette,"hot_desaturated.gp";//palette,"idl33.gp"; // palette,"hot_desaturated.gp";
    clo=min(parameter); chi=max(parameter);
    plf, parameter, y(xx,1:ny-1,1:nz)+3.0*dy/2, z(xx,1:ny-1,1:nz)+dz/2, cmin=clo, cmax=chi;
    plg, [0.15,-0.15],[-0.14,-0.14],marks=0,width=3,type="dash";
    plg, [0.15,-0.15],[ 0.00, 0.00],marks=0,width=3,type="dash";
    plg, [0.15,-0.15],[ 0.14, 0.14],marks=0,width=3,type="dash";
    plg, [0.075,0.075],[-0.20,0.20],marks=0,width=3,type="dash";
    plg, [-0.075,-0.075],[-0.20,0.20],marks=0,width=3,type="dash";
    colorbar, clo, chi;
    xytitles,"Z (m)", "Y (m)";
    title=strtrim(swrite("B_z_ (G) at x=", x(xx,1,1),"at t=",time*1e9,"ns"));
    pltitle, title;
    limits,zmin, zmax, zmin, zmax;
    };

    if ( if_win_22 == 1 ){;
    window,22;fma;
    yy = (ny+1)/2;
    parameter = 1.0e4*Bz(1:nx-1, yy, 1:nz) + BzExt(1:nx-1, yy, 1:nz);
    palette,"hot_desaturated.gp";//palette,"idl33.gp"; // palette,"hot_desaturated.gp";
    clo=min(parameter); chi=max(parameter);
    plf, parameter, x(1:nx-1,yy,1:nz)+3.0*dx/2, z(1:nx-1,yy,1:nz)+dz/2, cmin=clo, cmax=chi;
    plg, [0.15,-0.15],[-0.14,-0.14],marks=0,width=3,type="dash";
    plg, [0.15,-0.15],[ 0.00, 0.00],marks=0,width=3,type="dash";
    plg, [0.15,-0.15],[ 0.14, 0.14],marks=0,width=3,type="dash";
    plg, [0.075,0.075],[-0.20,0.20],marks=0,width=3,type="dash";
    plg, [-0.075,-0.075],[-0.20,0.20],marks=0,width=3,type="dash";
    colorbar, clo, chi;
    xytitles,"Z (m)", "X (m)";
    title=strtrim(swrite("B_z_ (G) at y=", y(1,yy,1),"at t=",time*1e9,"ns"));
    pltitle, title;
    limits,zmin, zmax, zmin, zmax;
    };

    if ( if_win_23 == 1 ){;
    window,23;fma;
    palette,"hot_desaturated.gp";
    zz = (nz+1)/2;
    parameter = 1.0e4*Bz(1:nx-1,1:ny-1,zz) + BzExt(1:nx-1,1:ny-1,zz);
    clo=min(parameter); chi=max(parameter);
    clo=25;chi=55;
    plf, parameter, y(1:nx-1,1:ny-1,zz)+dy/2, x(1:nx-1,1:ny-1,zz)+dx/2, cmin=clo, cmax=chi;
    plfc, parameter, y(1:nx-1,1:ny-1,zz)+dy/2, x(1:nx-1,1:ny-1,zz)+dx/2, levs=span(clo,chi,200);
    plg, [0.075,0.075],[-0.11,0.11],marks=0,width=3,type="dash";
    plg, [-0.075,-0.075],[-0.11,0.11],marks=0,width=3,type="dash";
    plg, [0.11,-0.11],[0.075,0.075],marks=0,width=3,type="dash";
    plg, [0.11,-0.11],[-0.075,-0.075],marks=0,width=3,type="dash";
    theta=span(0.0,2*pi,1000);
    plg, 0.044*sin(theta), 0.044*cos(theta), marks=0, width=3, color="white";
    plg, 0.05*sin(theta), 0.05*cos(theta), marks=0, width=3, color="white";
    colorbar, clo, chi;
    xytitles,"X (m)", "Y (m)";
    title=strtrim(swrite("B_z_ (G) at z=", z(1,1,zz),"at t=",time*1e9,"ns"));
    pltitle, title;
    limits,xmin, xmax, xmin, xmax;
    };
    
    if ( if_win_24 == 1 ){;
    window,24;fma;
    plg, Ex((nx+1)/2,(ny+1)/2,1:nz-1), z((nx+1)/2,(ny+1)/2,1:nz-1), marks=0, width=5, color="green";
    plg, Ex((nx+1)/4-1,(ny+1)/2-1,1:nz-1), z((nx+1)/4-1,(ny+1)/2-1,1:nz-1), marks=0, width=5, color=[0,180,0];
    plg, Ex(3*(nx+1)/4-1,(ny+1)/2-1,1:nz-1), z(3*(nx+1)/4-1,(ny+1)/2-1,1:nz-1), marks=0, width=5, color=[0,110,0];
    xytitles, "Z (m)", "V/m";
    pltitle,"E_x_ (green), E_y_ (blue) vs z at x=y=0";
    gridxy,1,1;
    range,e,e;
    limits,e,e;
    };

    myfac=0.5;
    if ( if_win_25 == 1 ){;
    window,25;fma;
    zz = int((nz+1)/2);
    theta = atan(y(1:nx,1:ny,zz),x(1:nx,1:ny,zz)+1.0e-20);
    parameter = -Ex(5:96,5:96,zz)*sin(theta(5:96,5:96)) + Ey(5:96,5:96,zz)*cos(theta(5:96,5:96));
    palette,"hot_desaturated.gp";
    clo=myfac*min(parameter); chi=myfac*max(parameter);
    chi = myfac*max(parameter); clo = -1*chi;
    plf, parameter, y(5:96,5:96,zz)+dy, x(5:96,5:96,zz)+dx, cmin=clo, cmax=chi;
    plfc, parameter, y(5:96,5:96,zz)+dy, x(5:96,5:96,zz)+dx, levs=span(clo,chi,800);
    plg, [0.075,0.075],[-0.11,0.11],marks=0,width=3,type="dash";
    plg, [-0.075,-0.075],[-0.11,0.11],marks=0,width=3,type="dash";
    plg, [0.11,-0.11],[0.075,0.075],marks=0,width=3,type="dash";
    plg, [0.11,-0.11],[-0.075,-0.075],marks=0,width=3,type="dash";
        theta=span(0.0,2*pi,1000);
        plg, 0.044*sin(theta), 0.044*cos(theta), marks=0, width=3, color="white";
        plg, 0.05*sin(theta), 0.05*cos(theta), marks=0, width=3, color="white";
    colorbar, clo, chi;
    xytitles,"X (m)", "Y (m)";
    title=strtrim(swrite("E_!q_ at z=", z(1,1,zz),"at t=",time*1e9,"ns"));
    pltitle, title;
    limits//,xmin, xmax, xmin, xmax;
    };
    
    if ( if_win_26 == 1 ){;
    window,26;fma;
    zz = int((nz+1)/2);
    theta = atan(y(1:nx,1:ny,zz),x(1:nx,1:ny,zz)+1.0e-20);
    parameter = Ex(5:96,5:96,zz)*cos(theta(5:96,5:96)) + Ey(5:96,5:96,zz)*sin(theta(5:96,5:96));
    palette,"hot_desaturated.gp";
    clo=myfac*min(parameter); chi=myfac*max(parameter);
    chi = myfac*max(parameter); clo = -1*chi;
    plf, parameter, y(5:96,5:96,zz)+dy, x(5:96,5:96,zz)+dx, cmin=clo, cmax=chi;
    plfc, parameter, y(5:96,5:96,zz)+dy, x(5:96,5:96,zz)+dx, levs=span(clo,chi,800);
    plg, [0.075,0.075],[-0.11,0.11],marks=0,width=3,type="dash";
    plg, [-0.075,-0.075],[-0.11,0.11],marks=0,width=3,type="dash";
    plg, [0.11,-0.11],[0.075,0.075],marks=0,width=3,type="dash";
    plg, [0.11,-0.11],[-0.075,-0.075],marks=0,width=3,type="dash";
        theta=span(0.0,2*pi,1000);
        plg, 0.044*sin(theta), 0.044*cos(theta), marks=0, width=3, color="white";
        plg, 0.05*sin(theta), 0.05*cos(theta), marks=0, width=3, color="white";

    colorbar, clo, chi;
    xytitles,"X (m)", "Y (m)";
    title=strtrim(swrite("E_r_ at z=", z(1,1,zz),"at t=",time*1e9,"ns"));
    pltitle, title;
    limits//,xmin, xmax, xmin, xmax;
    };
    
    if ( if_win_27 == 1 ){;
    window,27;fma;
    zz = int((nz+1)/2);
    Er   = Ex(1:nx-1,1:ny,zz)*cos(theta(1:nx-1,1:ny)) + Ey(1:nx-1,1:ny,zz)*sin(theta(1:nx-1,1:ny));
    Eth = -Ex(1:nx-1,1:ny,zz)*sin(theta(1:nx-1,1:ny)) + Ey(1:nx-1,1:ny,zz)*cos(theta(1:nx-1,1:ny));
    parameter = sqrt(Er(5:96,5:96)^2.0 + Eth(5:96,5:96)^2.0);
    palette,"hot_desaturated.gp";
    clo=myfac*min(parameter); chi=myfac*max(parameter);
//    plf, parameter, y(1:nx-1,1:ny,zz)+dx, x(1:nx-1,1:ny,zz)+dy, cmin=clo, cmax=chi;
    plf, parameter, y(5:96,5:96,zz)+dx, x(5:96,5:96,zz)+dy, cmin=clo, cmax=chi;
    plg, [0.075,0.075],[-0.11,0.11],marks=0,width=3,type="dash";
    plg, [-0.075,-0.075],[-0.11,0.11],marks=0,width=3,type="dash";
    plg, [0.11,-0.11],[0.075,0.075],marks=0,width=3,type="dash";
    plg, [0.11,-0.11],[-0.075,-0.075],marks=0,width=3,type="dash";
    colorbar, clo, chi;
    xytitles,"X (m)", "Y (m)";
    title=strtrim(swrite("sqrt(E_r_^2^ + E_!q_^2^) at z=", z(1,1,zz),"at t=",time*1e9,"ns"));
    pltitle, title;
    limits//,xmin, xmax, xmin, xmax;
};
    
    window,29;fma;
    plg, (1.0e4*Bz(51,:,101)+BzExt(51,:,101)),y(51,1:ny-1,101),marks=0,width=5,color="blue";
    limits,-0.06,0.06,35,48;
    gridxy,1,1;
    xytitles,"Y (m)","B_z_ (G)";
    pltitle,"Bz lineout at z=0, x=0"

    window,31;fma;
//        xx = 25-1 // 65;  // x = 5
    xx = 59;  // x = 2
//        yy = 65; // y = 5
    yy = 59; // y = 2
    plg, 1e-5*Ex(xx,yy,:), z(xx,yy,:), marks=0, width=5;
    plg, 1e-5*Ex(xx-1,yy,:), z(xx-1,yy,:), marks=0, width=5, type="dash";
    plg, 1e-5*Ex(xx+1,yy,:), z(xx+1,yy,:), marks=0, width=5, type="dot";
    plg, 1e-5*Ey(xx,yy,:), z(xx,yy,:), marks=0, width=5, color="red";
    plg, 1e-5*Ey(xx-1,yy,:), z(xx-1,yy,:), marks=0, width=5, color="red",type="dash";
    plg, 1e-5*Ey(xx+1,yy,:), z(xx+1,yy,:), marks=0, width=5, color="red",type="dot";
    plg, 1e-5*Ez(xx,yy,1:nz-1), z(xx,yy,1:nz-1), marks=0, width=5, color="blue";
    plg, 1e-5*Ez(xx-1,yy,1:nz-1), z(xx-1,yy,1:nz-1), marks=0, width=5, color="blue",type="dash";
    plg, 1e-5*Ez(xx+1,yy,1:nz-1), z(xx+1,yy,1:nz-1), marks=0, width=5, color="blue",type="dot";
    xytitles,"Z (m)", "kV/cm";
    //       pltitle,"x=5,y=5, E_x_ (bk), E_y_ (rd), E_z_ (bl)";
    pltitle,"x=2,y=2, E_x_ (bk), E_y_ (rd), E_z_ (bl)";
    gridxy,1,1;

    pause, pausetime;
};

/* ################################################################### */
