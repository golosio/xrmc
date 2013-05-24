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
  const int Nx=200, Ny=200;
  double dx=0.01, dy = 0.01;
  double x0, y0, x, y;
  int ix, iy, i;
  double *Image;
  double sx, sy;

  Image = (double*)malloc(Nx*Ny*sizeof(double));
  sx=0.4;
  sy=0.2;

  y0=(-0.5*Ny + 0.5)*dy;
  x0=(-0.5*Nx + 0.5)*dx;
  i = 0;
  for (iy=0; iy<Ny; iy++) {
    y = y0 + dy*iy;
    for (ix=0; ix<Nx; ix++) {
      x = x0 + dx*ix;
      Image[i] = Gauss(x, 0, sx)*Gauss(y, 0, sy);
      i++;
      //printf("%d\n", i);
    }
    //printf("%d\n", i);
  }
  fp = fopen("intensityscreen_image.dat", "wb");
  fwrite(Image, sizeof(double), Nx*Ny, fp);
  fclose(fp);
  free(Image);

  return 0;
}
