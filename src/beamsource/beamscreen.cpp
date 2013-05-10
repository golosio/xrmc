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
//     beamscreen.cpp            //
//        02/05/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Methods of the class detectorarray
//
#include <cmath>
#include <iostream>
#include "xrmc_math.h"
#include "beamscreen.h"
#include "xrmc_exception.h"
#include "xrmc_algo.h"
				    //#include "xrmc_arrayNd.h"
				    //#ifdef _OPENMP
				    //#include <omp.h>
				    //#endif

using namespace std;
using namespace xrmc_algo;
//using namespace arrayNd;

// destructor
beamscreen::~beamscreen() {
  if (Image!=NULL) delete[] Image;
  if (CumulImage!=NULL) delete[] CumulImage;
}

// constructor
beamscreen::beamscreen(string dev_name) {
  Runnable = false;
  NInputDevices = 0;

  PixelX=NULL;
  Image=NULL;
  CumulImage=NULL;
  NX = NY = N = NBins = 0;
  Shape = 0;
  RandomPixelFlag = 1;
  RunningFasterFlag = 0;
  SetDevice(dev_name, "beamscreen");
}

//////////////////////////////////////////////////////////////////////
// beamscreen initialization before run method
//////////////////////////////////////////////////////////////////////
int beamscreen::RunInit()
{
  Init();

  return 0;
}

// Generates a random pixel on the screen
int beamscreen::RandomPixel(int &iE, int &ix, int &iy, int &pol, randmt_t *rng)
{
  int i;
  double R = Rnd_r(rng);
  
  Locate(R, CumulImage, 2*N*NBins, &i);
  iE = i % NBins;
  i /= NBins;
  ix = i % NX;
  i /= NX;
  iy = i % NY;
  pol = i / NY;
  if (pol>1)
    throw xrmc_exception("Error in beam screen array dimensions.\n");

  return 0;
}

// Generates a random point and energy on the screen
vect3 beamscreen::RandomPoint(double &E, int &pol, double &w, randmt_t *rng)
{
  const double rmax = 100;
  int iE, ix, iy, iE1, ix1, iy1, ipol;
  int i000, i001, i010, i011, i100, i101, i110, i111;
  double c000, c001, c010, c011, c100, c101, c110, c111;
  double rE, rx, ry, tE, tx, ty, uE, ux, uy;

  RandomPixel(iE, ix, iy, pol, rng);
  rE = Rnd_r(rng) - 0.5;
  rx = Rnd_r(rng) - 0.5;
  ry = Rnd_r(rng) - 0.5;
  if (InterpolFlag==1) {
    if (rE>=0) {
      tE = rE;
      iE1 = iE+1;
      if (iE1>=NBins) iE1=NBins-1;
    }
    else {
      tE = -rE;
      iE1 = iE-1;
      if (iE1<0) iE1=0;
    }
    uE = 1. - tE;
    
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
    ipol = pol*NY*NX*NBins;
    i000 = ipol + NX*NBins*iy  + NBins*ix   + iE;
    i001 = ipol + NX*NBins*iy  + NBins*ix   + iE1;
    i010 = ipol + NX*NBins*iy  + NBins*ix1  + iE;
    i011 = ipol + NX*NBins*iy  + NBins*ix1  + iE1;
    i100 = ipol + NX*NBins*iy1 + NBins*ix   + iE;
    i101 = ipol + NX*NBins*iy1 + NBins*ix   + iE1;
    i110 = ipol + NX*NBins*iy1 + NBins*ix1  + iE;
    i111 = ipol + NX*NBins*iy1 + NBins*ix1  + iE1;
    
    c000 = Image[i000];
    c001 = Image[i001];
    if (c001>rmax*c000) c001 = rmax;
    else c001 /= c000;
    c010 = Image[i010];
    if (c010>rmax*c000) c010 = rmax;
    else c010 /= c000;
    c011 = Image[i011];
    if (c011>rmax*c000) c011 = rmax;
    else c011 /= c000;
    c100 = Image[i100];
    if (c100>rmax*c000) c100 = rmax;
    else c100 /= c000;
    c101 = Image[i101];
    if (c101>rmax*c000) c101 = rmax;
    else c101 /= c000;
    c110 = Image[i110];
    if (c110>rmax*c000) c110 = rmax;
    else c110 /= c000;
    c111 = Image[i111];
    if (c111>rmax*c000) c111 = rmax;
    else c111 /= c000;
    c000 = 1;
    
    w = c000*uy*ux*uE + c001*uy*ux*tE + c010*uy*tx*uE + c011*uy*tx*tE
      + c100*ty*ux*uE + c101*ty*ux*tE + c110*ty*tx*uE + c111*ty*tx*tE;
  }
  else w = 1;

  // central bin energy
  E = Emin + (Emax - Emin)*(0.5 + iE + rE)/NBins;
 // local x, y coordinates of the pixel
  double x = PixelSizeX*(-0.5*NX + 0.5 + ix + rx);
  double y = PixelSizeY*(-0.5*NY + 0.5 + iy + ry);
  vect3 r = X + ui*x + uj*y; // 3d absolute coordinates of the pixel

  return r;
}
