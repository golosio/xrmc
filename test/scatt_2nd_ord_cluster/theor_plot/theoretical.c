#include "xraylib.h"
#include <stdio.h>
#include <math.h>

int main()
{
  double th1_deg = 45;
  int th2_deg;
  double E=50, E1, Erc, Ecc;
  double RhoAl=2.7;
  double SAl = 0.01;
  double PAbsR, PAbsC, PAbsRC, PAbsCC;
  double d1=20, r1=0.2, d2=10, r2=0.1, Omega1, A1, Omega2, A2;

  double th1, th2;
  double PR, PC, P1;
  double PRRa, PRCa, PCRa, PCCa;

  XRayInit();

  th1 = th1_deg*PI/180;


  A1 = PI*r1*r1;
  Omega1 = A1/d1/d1;

  PAbsR = exp(-1.19*2.*(0.0805*CS_Total(1,E) + 0.5999*CS_Total(6,E) +
			  0.3196*CS_Total(8,E)));
  E1 = ComptonEnergy(E, th1);
  PAbsC = exp(-1.19*2.*(0.0805*CS_Total(1,E1) + 0.5999*CS_Total(6,E1) +
			  0.3196*CS_Total(8,E1)));

  PR = RhoAl*SAl*DCSP_Rayl(13, E, th1, PI/2)*Omega1*PAbsR;
  PC = RhoAl*SAl*DCSP_Compt(13, E, th1, PI/2)*Omega1*PAbsC;

  A2 = PI*r2*r2;
  Omega2 = A2/d2/d2;

  for (th2_deg=-90; th2_deg<=90; th2_deg+=1) { 
    if(th2_deg==0) th2=1e-5;
    else th2 = fabs(PI/180*th2_deg);
    Erc = ComptonEnergy(E, th2);
    PAbsRC = exp(-1.19*2.*(0.0805*CS_Total(1,Erc) + 0.5999*CS_Total(6,Erc) +
			   0.3196*CS_Total(8,Erc)));
    Ecc = ComptonEnergy(E1, th2);
    PAbsCC = exp(-1.19*2.*(0.0805*CS_Total(1,Ecc) + 0.5999*CS_Total(6,Ecc) +
			   0.3196*CS_Total(8,Ecc)));
    PRRa = PR*RhoAl*SAl*DCSP_Rayl(13, E, th2, PI/2)*Omega2*PAbsR;
    PRCa = PR*RhoAl*SAl*DCSP_Compt(13, E, th2, PI/2)*Omega2*PAbsRC;
    PCRa = PC*RhoAl*SAl*DCSP_Rayl(13, E1, th2, PI/2)*Omega2*PAbsC;
    PCCa = PC*RhoAl*SAl*DCSP_Compt(13, E1, th2, PI/2)*Omega2*PAbsCC;
    
    P1 = 1e10*(PRRa+PRCa+PCRa+PCCa);
    printf("%d\t%g\n", th2_deg, P1);
  }

  return 0;
}
