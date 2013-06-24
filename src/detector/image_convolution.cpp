#include <cmath>
#include "fft.h"
#include "image_convolution.h"

int convolution::BuildFilt(double *filt, int nn, gauss_vect gv,
			   int *ip, double *w, double scale)
{
  if (nn==1) {
    filt[0]=1;
    return 0;
  }

  double sum=0;
  for (int ix=0; ix<nn; ix++) {
    filt[ix] = 0;
    double x = (double)ix;
    if (ix > nn/2-1) x -= (double)nn;
    x *= scale;
    for (unsigned int ig=0; ig<gv.size(); ig++) {
      double x1 = fabs(x - gv[ig].x0);
      if (x1<gv[ig].sigma*4) {
	double k = 0.5/gv[ig].sigma/gv[ig].sigma;
	double val = gv[ig].height*exp(-(k*x1*x1));
	filt[ix] += val;
	sum += val;
      }
    }
  }
  if (sum==0) filt[0] = 1;
  else for (int ix=0; ix<nn; ix++) filt[ix] /= sum;

  rdft(nn, 1, filt, ip, w);

  return 0;
}

int convolution::FiltData(double *datas, double *filt, int nn,
			  int *ip, double *w) {
  rdft(nn, 1, datas, ip, w);
  double fac = 2.0 / nn;
  datas[0] *= fac*filt[0];
  datas[1] *= fac*filt[1];
  double *sp1 = datas+2;
  double *sp2 = filt+2;
  for(int ix=1; ix<nn/2; ix++) {
    double re = sp1[0]*sp2[0] - sp1[1]*sp2[1];
    double im = sp1[0]*sp2[1] + sp1[1]*sp2[0];
    sp1[0] = fac*re;
    sp1[1] = fac*im;
    sp1 += 2;
    sp2 += 2;
  }
  rdft(nn, -1, datas, ip, w);

  return 0;
}

int convolution::PadData(double *datas, int n, int nn)
{
  if (n < nn) {
    double fac = (datas[0] - datas[n-1])/(nn - n + 1);
    for (int i=n; i<nn; i++) {
      datas[i] = datas[n-1] + fac*(i - n + 1);
    }
  }

  return 0;
}


int convolution::XFilt(double *xfilt, double **image, int nx, int ny, int nnx,
		       double *xdata, int *ipx, double *wx)
{
  if (nx>1) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
	xdata[ix] = image[iy][ix];
      }
      PadData(xdata, nx, nnx);
      FiltData(xdata, xfilt, nnx, ipx, wx);
      for (int ix=0; ix<nx; ix++) {
	image[iy][ix] = xdata[ix];
      }
    }
  }

  return 0;
}

int convolution::YFilt(double *yfilt, double **image, int nx, int ny, int nny,
		       double *ydata, int *ipy, double *wy)
{
  if (ny>1) {
    for (int ix=0; ix<nx; ix++) {
      for (int iy=0; iy<ny; iy++) {
	ydata[iy] = image[iy][ix];
      }
      PadData(ydata, ny, nny);
      FiltData(ydata, yfilt, nny, ipy, wy);
      for (int iy=0; iy<ny; iy++) {
	image[iy][ix] = ydata[iy];
      }
    }
  }

  return 0;
}

int convolution::ConvolveFilt(double *filt1, double *filt2, double *out_filt,
			      int nn)
{
  out_filt[0] = filt1[0]*filt2[0];
  out_filt[1] = filt1[1]*filt2[1];
  double *sp1 = filt1+2;
  double *sp2 = filt2+2;
  double *out_sp = out_filt+2;
  for(int ix=1; ix<nn/2; ix++) {
    out_sp[0] = sp1[0]*sp2[0] - sp1[1]*sp2[1];
    out_sp[1] = sp1[0]*sp2[1] + sp1[1]*sp2[0];
    sp1 += 2;
    sp2 += 2;
    out_sp += 2;
  }

  return 0;
}
