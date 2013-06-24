#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double PI=3.1415927;
double Gauss(double x, double x0, double s)
{
  double A, x1;
  A = 1./sqrt(2.*PI*s*s);
  x1 = x - x0;
  return A*exp(-x1*x1/(2.*s*s));
}

int main()
{
  FILE *fp;
  const int Nx=100, Ny=100, NE=100;
  double dx=0.12, dy = 0.12, Emin=50, Emax=100;
  double x0, y0, x, y, r, E, dE, Ec, Ec0=50, Ec1=100;
  int ix, iy, iE, i;
  double *Image;
  double sx, sy, sE;

  Image = (double*)malloc(Nx*Ny*NE*sizeof(double));
  sx=4.0;
  sy=2.0;
  sE=4;

  y0=(-0.5*Ny + 0.5)*dy;
  x0=(-0.5*Nx + 0.5)*dx;
  dE = (Emax - Emin)/(NE-1);
  i = 0;
  for (iE=0; iE<NE; iE++) {
    E=Emin+dE*iE;
    for (iy=0; iy<Ny; iy++) {
      y = y0 + dy*iy;
      for (ix=0; ix<Nx; ix++) {
	x = x0 + dx*ix;
	r=sqrt(x*x + y*y)/(-x0);
	Ec = Ec0+(Ec1-Ec0)*r;
	Image[i] = Gauss(x, 0, sx)*Gauss(y, 0, sy)*Gauss(E, Ec, sE);
	i++;
	//printf("%d\n", i);
      }
    }
    //printf("%d\n", i);
  }
  fp = fopen("beamscreen_image.dat", "wb");
  fwrite(Image, sizeof(double), Nx*Ny*NE, fp);
  fclose(fp);
  free(Image);

  return 0;
}
