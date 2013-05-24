/*
Copyright (C) 2013 Bruno Golosio

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
///////////////////////////////////
//     intensityscreen.cpp       //
//        16/05/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Methods of the class intensityscreen
//
#include <cmath>
#include <iostream>
#include "xrmc_math.h"
#include "intensityscreen.h"
#include "xrmc_exception.h"
#include "xrmc_algo.h"

using namespace std;
using namespace xrmc_algo;
//using namespace arrayNd;

// destructor
intensityscreen::~intensityscreen() {
  if (Image!=NULL) delete[] Image;
  if (CumulImage!=NULL) delete[] CumulImage;
}

// constructor
intensityscreen::intensityscreen(string dev_name) {
  Runnable = false;
  NInputDevices = 0;

  PixelX=NULL;
  Image=NULL;
  CumulImage=NULL;

  NX = NY = N = 0;
  Shape = 0;
  RandomPixelFlag = 1;
  RunningFasterFlag = 0;
  SetDevice(dev_name, "intensityscreen");
}

//////////////////////////////////////////////////////////////////////
// intensityscreen initialization before run method
//////////////////////////////////////////////////////////////////////
int intensityscreen::RunInit()
{
  Init();

  if (Image==NULL || N<=0)
    throw xrmc_exception("intensityscreen image must be loaded before run.\n");

  if (CumulImage!=NULL) delete[] CumulImage;

  // allocate the cumulative image array
  CumulImage = new double [N];

  double sum=0;
  for(int i=0; i<N; i++) {
    CumulImage[i] = sum;
    double val = Image[i];
    if (val<0)
      throw xrmc_exception("intensityscreen pixel contents "
			   "must be nonnegative.\n");
    sum += val;
  }
  if (sum<=0)
      throw xrmc_exception("All intensityscreen pixel contents are null.\n");
    
  for(int i=0; i<N; i++) {
    Image[i] /= sum;
    CumulImage[i] /= sum;
  }

  return 0;
}

// Generates a random pixel on the screen
int intensityscreen::RandomPixel(int &ix, int &iy, randmt_t *rng)
{
  int i;
  double R = Rnd_r(rng);
  
  Locate(R, CumulImage, N, &i);
  ix = i % NX;
  iy = i / NX;
  if (iy>NY)
    throw xrmc_exception("Error in intensityscreen array dimensions.\n");

  return 0;
}

// Generates a random point on the screen
vect3 intensityscreen::RandomPoint(double &w, randmt_t *rng)
{
  int ix, iy;
  double rx, ry;

  RandomPixel(ix, iy, rng);
  rx = Rnd_r(rng) - 0.5;
  ry = Rnd_r(rng) - 0.5;

  // evaluates interpolation weight factor
  w = InterpolWeight(ix, iy, rx, ry);

 // local x, y coordinates of the pixel
  double x = PixelSizeX*(-0.5*NX + 0.5 + ix + rx);
  double y = PixelSizeY*(-0.5*NY + 0.5 + iy + ry);
  vect3 r = X + ui*x + uj*y; // 3d absolute coordinates of the pixel

  return r;
}

// Evaluates the weigth for a given trajectory
int intensityscreen::DirectionWeight(vect3 x0, vect3 u, double &w,
				     randmt_t *rng)
{
  int ix, iy;
  double tx, ty, t;
  if (!Intersect(x0, u, ix, iy, tx, ty, t)) {
    w = 0;
    return 0;
  }
  double dO = dOmega(u*t);
  double num = Image[NX*iy+ix];
  if (num+dO<=num) { // check for division overflow or wrong side of screen
    w=0;
    return 1;
  }

  double w0 = num/dO;
  double rx = tx - 0.5;
  double ry = ty - 0.5;
  w = w0*InterpolWeight(ix, iy, rx, ry);

  return 0;
}


double intensityscreen::InterpolWeight(int ix, int iy, double rx, double ry)
{
  const double rmax = 100;
  int ix1, iy1;
  int i00, i01, i10, i11;
  double c00, c01, c10, c11;
  double w, tx, ty, ux, uy;

  if (InterpolFlag==0) return 1;

  if (rx>=0) {
    tx = rx;
    ix1 = ix+1;
    if (ix1>=NX) ix1=NX-1;
  }
  else {
    tx = -rx;
    ix1 = ix-1;
    if (ix1<0) ix1=0;
  }
  ux = 1. - tx;
    
  if (ry>=0) {
    ty = ry;
    iy1 = iy+1;
    if (iy1>=NY) iy1=NY-1;
  }
  else {
    ty = -ry;
    iy1 = iy-1;
    if (iy1<0) iy1=0;
  }
  uy = 1. - ty;

  i00 = NX*iy + ix;
  i01 = NX*iy + ix1;
  i10 = NX*iy1 + ix;
  i11 = NX*iy1 + ix1;
    
  c00 = Image[i00];
  c01 = Image[i01];
  c10 = Image[i10];
  c11 = Image[i11];

  double norm = c00; 
  if (c00>=rmax*norm) c00 = rmax;
  else c00 /= norm;
  if (c01>=rmax*norm) c01 = rmax;
  else c01 /= norm;
  if (c10>=rmax*norm) c10 = rmax;
  else c10 /= norm;
  if (c11>=rmax*norm) c11 = rmax;
  else c11 /= norm;
    
  w = c00*uy*ux + c01*uy*tx + c10*ty*ux + c11*ty*tx;

  return w;
}


intensityscreen *intensityscreen::Clone(string dev_name) {
  //cout << "Entering intensityscreen::Clone\n";
  intensityscreen *clone = new intensityscreen(dev_name);

  *clone = *this;

  return clone;
}

