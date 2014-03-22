#include <stdio.h>
#include <math.h>
#include "xraylib.h"

int Signal(float s, float y, float z, float *dC, float *dR);

int main(int argc, char *argv[])
{
  XRayInit();
  const int N=100;  
  float R=0.618039/2;
  float s=R/N;
  float y, z, r;
  int iy, iz;
  float dC, dR; 
  float SCompt=0, SRayl=0;

  for (iz=-N; iz<=N; iz++) {
    z=s*iz;
    for (iy=-N; iy<=N; iy++) {
      y=s*iy;
      r = sqrt(y*y+z*z);
      if (r>=R) continue;
      Signal(s, y, z, &dC, &dR);
      SCompt += dC;
      SRayl += dR;
    }
  }
  printf("SCompt: %g\n", SCompt);
  printf("SRayl: %g\n", SRayl);
  printf("STot: %g\n", SCompt+SRayl);

  return 0;
}

int Signal(float s, float y, float z, float *dC, float *dR)
{
  float I0=1e12;
  int NElem=3;
  int Z[]={8, 48, 74};
  float w[]={0.177644, 0.312027, 0.510329};
  float rho=7.9;

  float d=1;
  //float y=0.15;
  //float z=0;
  //float s=0.618039/200;

  float E0 = 28;

  float th, phi;
  float d1, th1, E1, mu0, mu1, l, l1, dOmega;
  float dSCompt, dSRayl;
  float alpha, d2;
  int i;

  d1=sqrt(d*d+y*y+z*z);
  //th=PI/2+atan(y/d);
  th=PI/2-asin(y/d1);
  phi=atan(z/d);
  E1=ComptonEnergy(E0, th);
  //printf("E1: %g\n", E1);

  th1=acos(d/d1);
  //printf("d1: %g\n", d1);
  //printf("th1: %g\n", th1);
  dOmega=s*s*cos(th1)/(d1*d1);
  //printf("dOmega: %g\n", dOmega);

  mu0 = mu1 = dSCompt = dSRayl = 0;
  for(i=0; i<NElem; i++) {
    mu0 += w[i]*CS_Total(Z[i], E0);
    mu1 += w[i]*CS_Total(Z[i], E1);
    dSCompt += w[i]*DCSP_Compt(Z[i], E0, th, phi);
    dSRayl += w[i]*DCSP_Rayl(Z[i], E0, th, phi);
  }
  mu0 *= rho;
  mu1 *= rho;
  dSCompt *= rho*dOmega; 
  dSRayl *= rho*dOmega; 

  //printf("mu0: %g\n", mu0);
  //printf("mu1: %g\n", mu1);
  ///////////////////////////////////////////
  // y = r cos th
  // x = r sin th cos phi
  // z = r sin th sin phi
  // x = 1 + y
  // a = PI/2 - th
  // y = r sin a
  // x = r cos a cos phi
  // z = r cos a sin phi
  // 1 + r sin a = r cos a cos phi
  // r = 1 / (cos a cos phi - sin a)
  alpha = asin(y/d1);
  d2 = 1./(cos(alpha)*cos(phi) - sin(alpha));
  l = 1./(mu0*(1.+d2));
  l1 = 1./(mu0+mu1*d2);

  //l = 1./(2.*mu0);
  //l1 = 1./(mu0+mu1);
  //printf("l: %g\n", l);

  //printf("dSCompt: %g\n", dSCompt);
  //printf("dSRayl: %g\n", dSRayl);

  *dC = dSCompt*l1*I0;
  *dR = dSRayl*l*I0;

  return 0;
}
