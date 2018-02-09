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
//     intensityscreen.h         //
//        16/05/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// intensityscreen class definition
//

#ifndef INTENSITYSCREENH
#define INTENSITYSCREENH
#include "xrmc_device.h"
//#include "xrmc_source.h"
//#include "xrmc_arrayNd.h"
#include "xrmc_screen.h"
#include "randmt.h"

//  intensityscreen class definition, member variables and functions
class intensityscreen : public xrmc_screen
{
 public:
  //int NX, NY, N; // number of rows (NY), column (NX) and pixels (NX x NY)
  double *Image; // image array
  double *CumulImage; // cumulative probability distribution for x, y
  
  ~intensityscreen(); // destructor
  intensityscreen(std::string dev_name); // constructor

  virtual int Load(std::istream &fs); // load intensityscreen from file
  int LoadData(string file_name); // Load screen array contents
  virtual int SetDefault();

  // Generates a random point on the screen
  vect3 RandomPoint(double &w, randmt_t *rng);
  // Evaluates the weigth for a given trajectory
  int DirectionWeight(vect3 x0, vect3 u, double &w, randmt_t *rng);

  // int RunInit(); // initialize before run
  virtual intensityscreen *Clone(std::string dev_name);

 protected:
  //double PixelSizeX, PixelSizeY, PixelSurf; // pixel size and surface (cm2) 
  //vect3 *PixelX; // pixel coordinates array
  int InterpolFlag; // do not use / use interpolation inside pixels/bins (0/1)

  virtual int RunInit(); // beamscreen initialization before run

  // Generates a random point on the pixel surface
  //vect3 RandomPointOnPixel(int i);
  //vect3 RandomPointOnPixel(int i, randmt_t *rng);

  // Generates a random pixel on the screen
  int RandomPixel(int &ix, int &iy, randmt_t *rng);
  /*
  // evaluates the probability weigth for the trajectory x = x0 + u*t
  double WeightTrajectory(vect3 *x0, vect3 *u);
  */
  // evaluates interpolation weight factor
  double InterpolWeight(int ix, int iy, double rx, double ry);

};

#endif
