#ifndef IMAGE_CONVOLUTION_H
#define IMAGE_CONVOLUTION_H
#include <vector>

struct gauss_par
{
  double height;
  double x0;
  double sigma;
};

typedef std::vector<gauss_par> gauss_vect;

namespace convolution
{
  int XFilt(double *xfilt, double **image, int nx, int ny, int nnx,
	    double *xdata, int *ip, double *w);
  int YFilt(double *yfilt, double **image, int nx, int ny, int nny,
	    double *ydata, int *ip, double *w);
  int BuildFilt(double *filt, int nn, gauss_vect gv, int *ip, double *w,
		double scale);
  
  int FiltData(double *datas, double *filt, int nn, int *ip, double *w);
  
  int PadData(double *datas, int n, int nn);

  int XFilt(double *xfilt, double **image, int nx, int ny, int nnx,
	    double *xdata, int *ipx, double *wx);

  int YFilt(double *yfilt, double **image, int nx, int ny, int nny,
	    double *ydata, int *ipy, double *wy);

  int ConvolveFilt(double *filt1, double *filt2, double *out_filt, int nn);
}

#endif
