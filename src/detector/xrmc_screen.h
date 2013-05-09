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
//     xrmc_screen.h           //
//        05/05/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// detectorarray class definition
//

#ifndef XRMCSCREENH
#define XRMCSCREENH
#include "xrmc_device.h"
				    //#include "xrmc_source.h"
				    //#include "xrmc_arrayNd.h"
#include "randmt.h"

//  xrmc_screen class definition, member variables and functions
class xrmc_screen : public bodydevice
{
 public:
  int NX, NY, N; // number of rows (NY), column (NX) and pixels (NX x NY)
  
  virtual ~xrmc_screen(); // destructor

  // Evaluates the pixel and position of the intersection
  // of the straight line x = x0 + u*t with the screen
  bool Intersect(vect3 x0, vect3 u, int &ix, int &iy, double &tx, double &ty,
		 double &t);
  
 protected:
  double dOmegaLim; // cut on solid angle from interaction point to pixel 
  int Shape; // pixel shape (0: rectangular, 1: elliptical)
  double PixelSizeX, PixelSizeY, PixelSurf; // pixel size and surface (cm2) 
  vect3 *PixelX; // pixel coordinates array
  int RandomPixelFlag; // flag to enable/disable random point on pixel
  int RunningFasterFlag; //columns(0) or rows(1) running faster
  
  int Init(); // screen initialization
  //vect3 RandomPointOnPixel(int i); // Generates a random point
  vect3 RandomPointOnPixel(int i, randmt_t *rng); // on the pixel surface
  double dOmega(vect3 DRp); // evaluates the factor dOmega, related to
  // the probability that the last photon trajectory crosses the pixel  

};

#endif
