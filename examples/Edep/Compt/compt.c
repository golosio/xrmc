#include <stdio.h>
#include <xraylib.h>
#include <math.h>

int main()
{
  int ith, Nth=1000, iZ;
  double th, rho, num, denom, Edep, MeanEdep;
  int Z[] = {1,6,8};
  double w[] = {0.0805, 0.5999, 0.3196}; 
  double E = 50;

  XRayInit();

  num = denom = 0;
  
  for (ith=1; ith<Nth; ith++) {
    th = PI*ith/Nth;
    rho = 0;
    for (iZ=0; iZ<3; iZ++) {
      rho += w[iZ]*sin(th)*DCS_Compt(Z[iZ], E, th);
    }
    Edep = E - ComptonEnergy(E, th);
    num += rho*Edep;
    denom += rho;
  }
  MeanEdep = num/denom;

  printf("%f\n", MeanEdep);

  return 0;
}
