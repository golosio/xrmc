#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double PI=3.1415927;

int main()
{
  FILE *fp;
  const int Nx=200, Ny=200;
  double dx=0.06, dy = 0.06;
  double x0, y0, x, y, r;
  int ix, iy, i;
  double *Image;
  double R = 0.05*115;

  Image = (double*)malloc(Nx*Ny*sizeof(double));

  y0=(-0.5*Ny + 0.5)*dy;
  x0=(-0.5*Nx + 0.5)*dx;
  i = 0;
  for (iy=0; iy<Ny; iy++) {
    y = y0 + dy*iy;
    for (ix=0; ix<Nx; ix++) {
      x = x0 + dx*ix;
      r=sqrt(x*x + y*y);
      if (r>R) Image[i] = 0;
      else Image[i] = 1;
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
