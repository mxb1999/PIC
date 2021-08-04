/* ################################################################
    Implementation of the FRC antenna
    and explicit Yee method for Maxwell's equations

    --> This file contains functions related to defining
        the magnetic fields for the six PFRC-1
        Helmoltz coils (externally applied fields). <--

    A. B. Sefkow, on behalf of the TriForce team
    Version 1.1
    Date: July 9th, 2021
 ################################################################### */

print, "LOADED 'frc_helmholtz.i'!"

/* ###############  Externally applied solenoidal coils ######################## */

func Pi(n,m){;
    x=span(0.0,pi/2.0,201);
    r=integ(((1-n*sin(x)^2.0)*sqrt(1 - m*sin(x)^2))^-1.0,x,x(0));
    return r;
};

mult1 = mult2 = 100*0.102041;
mult3 = mult4 = 100*0.102041*0.9*(5./2.);
mult5 = mult6 = 100*0.102041*0.57;
curr1 = mult1*1e7; //current
curr2 = mult2*1e7; //current
curr3 = mult3*1e7; //current
curr4 = mult4*1e7; //current
curr5 = mult5*1e7; //current
curr6 = mult6*1e7; //current

epsilon = 1.0e-8; //1.0e-10; Material coefficient in cgs
// a = inner radius
a1 = 0.135 + epsilon;
a2 = 0.135 + epsilon;
a3 = 0.135 + epsilon;
a4 = 0.135 + epsilon;
a5 = 0.020 + epsilon;
a6 = 0.020 + epsilon;

// L = lengths
L1 = 0.80;
L2 = 0.80;
L3 = 0.160;
L4 = 0.160;
L5 = 0.020;
L6 = 0.020;

// zc = zcenters of coils
zc1 = -1.045;
zc2 = 1.045;
zc3 = -0.520;
zc4 = 0.520;
zc5 = -0.430;
zc6 = 0.430;

rmin=0.0; rmax=sqrt(0.125^2 + 0.125^2);
//zmin=-0.225; zmax=0.225;
nr=(nx-1)/2 + 1;
rExt = array(0.0,nr,nz); zExt = array(0.0,nr,nz);
Br = array(0.0,nr,nz); Bz = array(0.0,nr,nz);
Br1 = array(0.0,nr,nz); Bz1 = array(0.0,nr,nz);
Br2 = array(0.0,nr,nz); Bz2 = array(0.0,nr,nz);
Br3 = array(0.0,nr,nz); Bz3 = array(0.0,nr,nz);
Br4 = array(0.0,nr,nz); Bz4 = array(0.0,nr,nz);
Br5 = array(0.0,nr,nz); Bz5 = array(0.0,nr,nz);
Br6 = array(0.0,nr,nz); Bz6 = array(0.0,nr,nz);

for(i=1;i<=nr;++i){;zExt(i,) = span(zmin,zmax,nz);};
for(j=1;j<=nz;++j){;rExt(,j) = span(rmin,rmax,nr);};

// CALCULATE SOLENOID #1
print,"Calculating solenoid #1..."
a=a1; L=L1; zc=zc1;

xip = zExt(,)+(L/2.0)-zc;xim = zExt(,)-(L/2.0)-zc;
kp = sqrt((4.*a*rExt(,))/((a+rExt(,))^2.0+xip^2.0));
km = sqrt((4.*a*rExt(,))/((a+rExt(,))^2.0+xim^2.0));
h = sqrt((4.*a*rExt(,))/((a+rExt(,))^2.0));

for(i=2;i<=nr;++i){;for(j=1;j<=nz;++j){;
    if (rExt(i,j)!=a){;
        Br1(i,j) = (mu0*curr1/(2.0*pi*L))*sqrt(a/(rExt(i,j)+1.0e-10))*(((((kp(i,j)^2.0-2.0)/kp(i,j))*ellip_k(kp(i,j)^2.0))+((2.0/kp(i,j))*ellip_e(kp(i,j)^2.0)))-((((km(i,j)^2.0-2.0)/km(i,j))*ellip_k(km(i,j)^2.0))+((2.0/km(i,j))*ellip_e(km(i,j)^2.0))))
        Bz1(i,j) = (mu0*curr1/(4.0*pi*L*sqrt(a*rExt(i,j))))*((xip(i,j)*kp(i,j)*(ellip_k(kp(i,j)^2.0)+(((a-rExt(i,j))/(a+rExt(i,j)))*Pi(h(i,j)^2.0,kp(i,j)^2.0))))-(xim(i,j)*km(i,j)*(ellip_k(km(i,j)^2.0)+(((a-rExt(i,j))/(a+rExt(i,j)))*Pi(h(i,j)^2.0,km(i,j)^2.0)))))
    };
};};

// CALCULATE SOLENOID #2
print,"Calculating solenoid #2..."
a=a2; L=L2; zc=zc2;

xip = zExt(,)+(L/2.0)-zc;xim = zExt(,)-(L/2.0)-zc;
kp = sqrt((4.*a*rExt(,))/((a+rExt(,))^2.0+xip^2.0));
km = sqrt((4.*a*rExt(,))/((a+rExt(,))^2.0+xim^2.0));
h = sqrt((4.*a*rExt(,))/((a+rExt(,))^2.0));

for(i=2;i<=nr;++i){;for(j=1;j<=nz;++j){;
    if (rExt(i,j)!=a){;
        Br2(i,j) = (mu0*curr1/(2.0*pi*L))*sqrt(a/(rExt(i,j)+1.0e-10))*(((((kp(i,j)^2.0-2.0)/kp(i,j))*ellip_k(kp(i,j)^2.0))+((2.0/kp(i,j))*ellip_e(kp(i,j)^2.0)))-((((km(i,j)^2.0-2.0)/km(i,j))*ellip_k(km(i,j)^2.0))+((2.0/km(i,j))*ellip_e(km(i,j)^2.0))))
        Bz2(i,j) = (mu0*curr1/(4.0*pi*L*sqrt(a*rExt(i,j))))*((xip(i,j)*kp(i,j)*(ellip_k(kp(i,j)^2.0)+(((a-rExt(i,j))/(a+rExt(i,j)))*Pi(h(i,j)^2.0,kp(i,j)^2.0))))-(xim(i,j)*km(i,j)*(ellip_k(km(i,j)^2.0)+(((a-rExt(i,j))/(a+rExt(i,j)))*Pi(h(i,j)^2.0,km(i,j)^2.0)))))
    };
};};

// CALCULATE SOLENOID #3
print,"Calculating solenoid #3..."
a=a3; L=L3; zc=zc3;

xip = zExt(,)+(L/2.0)-zc;xim = zExt(,)-(L/2.0)-zc;
kp = sqrt((4.*a*rExt(,))/((a+rExt(,))^2.0+xip^2.0));
km = sqrt((4.*a*rExt(,))/((a+rExt(,))^2.0+xim^2.0));
h = sqrt((4.*a*rExt(,))/((a+rExt(,))^2.0));

for(i=2;i<=nr;++i){;for(j=1;j<=nz;++j){;
    if (rExt(i,j)!=a){;
        Br3(i,j) = (mu0*curr3/(2.0*pi*L))*sqrt(a/(rExt(i,j)+1.0e-10))*(((((kp(i,j)^2.0-2.0)/kp(i,j))*ellip_k(kp(i,j)^2.0))+((2.0/kp(i,j))*ellip_e(kp(i,j)^2.0)))-((((km(i,j)^2.0-2.0)/km(i,j))*ellip_k(km(i,j)^2.0))+((2.0/km(i,j))*ellip_e(km(i,j)^2.0))))
        Bz3(i,j) = (mu0*curr3/(4.0*pi*L*sqrt(a*rExt(i,j))))*((xip(i,j)*kp(i,j)*(ellip_k(kp(i,j)^2.0)+(((a-rExt(i,j))/(a+rExt(i,j)))*Pi(h(i,j)^2.0,kp(i,j)^2.0))))-(xim(i,j)*km(i,j)*(ellip_k(km(i,j)^2.0)+(((a-rExt(i,j))/(a+rExt(i,j)))*Pi(h(i,j)^2.0,km(i,j)^2.0)))))
    };
};};

// CALCULATE SOLENOID #4
print,"Calculating solenoid #4..."
a=a4; L=L4; zc=zc4;

xip = zExt(,)+(L/2.0)-zc;xim = zExt(,)-(L/2.0)-zc;
kp = sqrt((4.*a*rExt(,))/((a+rExt(,))^2.0+xip^2.0));
km = sqrt((4.*a*rExt(,))/((a+rExt(,))^2.0+xim^2.0));
h = sqrt((4.*a*rExt(,))/((a+rExt(,))^2.0));

for(i=2;i<=nr;++i){;for(j=1;j<=nz;++j){;
    if (rExt(i,j)!=a){;
        Br4(i,j) = (mu0*curr4/(2.0*pi*L))*sqrt(a/(rExt(i,j)+1.0e-10))*(((((kp(i,j)^2.0-2.0)/kp(i,j))*ellip_k(kp(i,j)^2.0))+((2.0/kp(i,j))*ellip_e(kp(i,j)^2.0)))-((((km(i,j)^2.0-2.0)/km(i,j))*ellip_k(km(i,j)^2.0))+((2.0/km(i,j))*ellip_e(km(i,j)^2.0))))
        Bz4(i,j) = (mu0*curr4/(4.0*pi*L*sqrt(a*rExt(i,j))))*((xip(i,j)*kp(i,j)*(ellip_k(kp(i,j)^2.0)+(((a-rExt(i,j))/(a+rExt(i,j)))*Pi(h(i,j)^2.0,kp(i,j)^2.0))))-(xim(i,j)*km(i,j)*(ellip_k(km(i,j)^2.0)+(((a-rExt(i,j))/(a+rExt(i,j)))*Pi(h(i,j)^2.0,km(i,j)^2.0)))))
    };
};};

// CALCULATE SOLENOID #5
print,"Calculating solenoid #5..."
a=a5; L=L5; zc=zc5;

xip = zExt(,)+(L/2.0)-zc;xim = zExt(,)-(L/2.0)-zc;
kp = sqrt((4.*a*rExt(,))/((a+rExt(,))^2.0+xip^2.0));
km = sqrt((4.*a*rExt(,))/((a+rExt(,))^2.0+xim^2.0));
h = sqrt((4.*a*rExt(,))/((a+rExt(,))^2.0));

for(i=2;i<=nr;++i){;for(j=1;j<=nz;++j){;
    if (rExt(i,j)!=a){;
        Br5(i,j) = (mu0*curr5/(2.0*pi*L))*sqrt(a/(rExt(i,j)+1.0e-10))*(((((kp(i,j)^2.0-2.0)/kp(i,j))*ellip_k(kp(i,j)^2.0))+((2.0/kp(i,j))*ellip_e(kp(i,j)^2.0)))-((((km(i,j)^2.0-2.0)/km(i,j))*ellip_k(km(i,j)^2.0))+((2.0/km(i,j))*ellip_e(km(i,j)^2.0))))
        Bz5(i,j) = (mu0*curr5/(4.0*pi*L*sqrt(a*rExt(i,j))))*((xip(i,j)*kp(i,j)*(ellip_k(kp(i,j)^2.0)+(((a-rExt(i,j))/(a+rExt(i,j)))*Pi(h(i,j)^2.0,kp(i,j)^2.0))))-(xim(i,j)*km(i,j)*(ellip_k(km(i,j)^2.0)+(((a-rExt(i,j))/(a+rExt(i,j)))*Pi(h(i,j)^2.0,km(i,j)^2.0)))))
    };
};};

// CALCULATE SOLENOID #6
print,"Calculating solenoid #6..."
a=a6; L=L6; zc=zc6;

xip = zExt(,)+(L/2.0)-zc;xim = zExt(,)-(L/2.0)-zc;
kp = sqrt((4.*a*rExt(,))/((a+rExt(,))^2.0+xip^2.0));
km = sqrt((4.*a*rExt(,))/((a+rExt(,))^2.0+xim^2.0));
h = sqrt((4.*a*rExt(,))/((a+rExt(,))^2.0));

for(i=2;i<=nr;++i){;for(j=1;j<=nz;++j){;
    if (rExt(i,j)!=a){;
        Br6(i,j) = (mu0*curr6/(2.0*pi*L))*sqrt(a/(rExt(i,j)+1.0e-10))*(((((kp(i,j)^2.0-2.0)/kp(i,j))*ellip_k(kp(i,j)^2.0))+((2.0/kp(i,j))*ellip_e(kp(i,j)^2.0)))-((((km(i,j)^2.0-2.0)/km(i,j))*ellip_k(km(i,j)^2.0))+((2.0/km(i,j))*ellip_e(km(i,j)^2.0))))
        Bz6(i,j) = (mu0*curr6/(4.0*pi*L*sqrt(a*rExt(i,j))))*((xip(i,j)*kp(i,j)*(ellip_k(kp(i,j)^2.0)+(((a-rExt(i,j))/(a+rExt(i,j)))*Pi(h(i,j)^2.0,kp(i,j)^2.0))))-(xim(i,j)*km(i,j)*(ellip_k(km(i,j)^2.0)+(((a-rExt(i,j))/(a+rExt(i,j)))*Pi(h(i,j)^2.0,km(i,j)^2.0)))))
    };
};};
 
// Set the on-axis r=0 value to the first r zone's value:
Bz1(1,:)=Bz1(2,:);
Bz2(1,:)=Bz2(2,:);
Bz3(1,:)=Bz3(2,:);
Bz4(1,:)=Bz4(2,:);
Bz5(1,:)=Bz5(2,:);
Bz6(1,:)=Bz6(2,:);

// Add together all
extBr = array(1.0,nr,nz);//Br1 + Br2 + Br3 + Br4 + Br5 + Br6;
extBz = array(1.0,nr,nz);//Bz1 + Bz2 + Bz3 + Bz4 + Bz5 + Bz6;
modB = sqrt(extBr^2.0+extBz^2.0);

// Clear out the memory
Br1 = Br2 = Br3 = Br4 = Br5 = Br6 = [];
Bz1 = Bz2 = Bz3 = Bz4 = Bz5 = Bz6 = [];
