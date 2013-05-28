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
//     xrmc_screen.cpp           //
//        05/05/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Methods of the class xrmc_screen
//
#include <cmath>
#include <iostream>
#include "xrmc_math.h"
#include "xrmc_screen.h"
#include "xrmc_exception.h"
#include "xrmc_algo.h"
#include "xrmc_arrayNd.h"

				    //#ifdef _OPENMP
				    //#include <omp.h>
				    //#endif

using namespace std;
using namespace xrmc_algo;
using namespace arrayNd;

// destructor
xrmc_screen::~xrmc_screen() {
  if (PixelX!=NULL) delete[] PixelX;
}

//////////////////////////////////////////////////////////////////////
// screen initialization method
//////////////////////////////////////////////////////////////////////
int xrmc_screen::Init()
{
  double x, y;

  OrthoNormal(ui, uj, uk);  // evaluates uj to form a orthonormal basis

  PixelSurf = PixelSizeX*PixelSizeY; // pixel surface
  N = NX*NY; // number of pixels
  if (PixelX!=NULL) delete[] PixelX;
  PixelX = new vect3[N]; // array of pixel coordinates 

  double x0 = -(PixelSizeX*NX)/2 + PixelSizeX/2; // starting (local) coordinates
  double y0 = -(PixelSizeY*NY)/2 + PixelSizeY/2; // on the detector plane
  int i = 0;

  //if (RunningFasterFlag==0) { // columns running faster than rows
    for (int iy=0; iy<NY; iy++) { // loop on detector pixels
      for (int ix=0; ix<NX; ix++) {
	x = x0 + PixelSizeX*ix; // local x coordinate of the pixel
	y = y0 + PixelSizeY*iy; // local y coordinate of the pixel
	PixelX[i] = X + ui*x + uj*y; // 3d absolute coordinates of the pixel
	i++; // increase pixel index
      }
    }
    //}
    /*
  else { // rows running faster than columns
    for (int ix=0; ix<NX; ix++) { // loop on detector pixels
      for (int iy=0; iy<NY; iy++) {
	x = x0 + PixelSizeX*ix; // local x coordinate of the pixel
	y = y0 + PixelSizeY*iy; // local y coordinate of the pixel
	PixelX[i] = X + ui*x + uj*y; // 3d absolute coordinates of the pixel
	i++; // increase pixel index
      }
    }
  }
    */
  
  return 0;
}

//////////////////////////////////////////////////////////////////////
// Generates a random point on the pixel surface
//////////////////////////////////////////////////////////////////////
/* DELETE vect3 xrmc_screen::RandomPointOnPixel(int i)
{
  double x, y, rx, ry;
  vect3 p;

  if (RandomPixelFlag == 1) { // check if random position on pixel is enabled
    do {
      rx=2.*Rnd()-1; // random real numbers between -1 and 1
      ry=2.*Rnd()-1;
    } while (Shape==1 && (rx*rx+ry*ry>1)); // if shape is elliptical check that
                                           // the point is inside the ellipse
    x = PixelSizeX*rx/2; // x,y coordinates of the point on pixel surface
    y = PixelSizeY*ry/2;
  }
  else {
    x = 0;
    y = 0;
  }
  p = PixelX[i] + ui*x + uj*y; // 3d absolute coordinates of the point

  return p;
}
*/

vect3 xrmc_screen::RandomPointOnPixel(int i, randmt_t *rng)
{
  double x, y, rx, ry;
  vect3 p;

  if (RandomPixelFlag == 1) { // check if random position on pixel is enabled
    do {
      rx=2.*Rnd_r(rng)-1; // random real numbers between -1 and 1
      ry=2.*Rnd_r(rng)-1;
    } while (Shape==1 && (rx*rx+ry*ry>1)); // if shape is elliptical check that
                                           // the point is inside the ellipse
    x = PixelSizeX*rx/2; // x,y coordinates of the point on pixel surface
    y = PixelSizeY*ry/2;
  }
  else {
    x = 0;
    y = 0;
  }
  p = PixelX[i] + ui*x + uj*y; // 3d absolute coordinates of the point

  return p;
}

//////////////////////////////////////////////////////////////////////
// evaluates the factor dOmega, related to
// the probability that the last photon trajectory crosses the pixel
//////////////////////////////////////////////////////////////////////
double xrmc_screen::dOmega(vect3 DRp)
{
  double r = DRp.Mod();        // distance from interaction point to the
                               // intersection of the trajectory with the pixel
  double num = -(DRp*uk)*PixelSurf; // numerator
  if (num<=0) return 0;
  double denom = r*r*r; // denominator
  if (num+denom == num) return dOmegaLim; // check for overflow in the ratio
  double dO = num/denom;
  if (dO>dOmegaLim) return dOmegaLim; // threshold to dOmegaLim

  return dO;
}

//////////////////////////////////////////////////////////////////////
// Evaluates the pixel and position of the intersection
// of the straight line x = x0 + u*t with the screen
//////////////////////////////////////////////////////////////////////
bool xrmc_screen::Intersect(vect3 x0, vect3 u, int &ix, int &iy,
			    double &tx, double &ty, double &t)
{
  double num=uk*(X-x0);
  double denom=(uk*u);
  if (num+denom==num) return false;
  t=num/denom;
  if (t<=0) return false;
  double x=ui*(x0-X+u*t)/PixelSizeX;
  double y=uj*(x0-X+u*t)/PixelSizeY;
  if (2.*fabs(x)>NX || 2.*fabs(y)>NY) return false;
  x += (double)NX/2;
  y += (double)NY/2;
  ix=(int)floor(x);
  iy=(int)floor(y);
  tx = x - ix;
  ty = y - iy;

  return true;
}
