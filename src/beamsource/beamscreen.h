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
//        beamscreen.h           //
//        02/05/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// nonuniformbeam class definition
//

#ifndef BEAMSCREENH
#define BEAMSCREENH
#include "xrmc_device.h"
//#include "xrmc_source.h"
//#include "xrmc_arrayNd.h"
#include "xrmc_screen.h"
#include "randmt.h"

//  beamscreen class definition, member variables and functions
class beamscreen : public xrmc_screen
{
 public:
  //int NX, NY, N; // number of rows (NY), column (NX) and pixels (NX x NY)
  int NBins; // Num. of energy bins
  double *Image; // image array
  double *CumulImage; // cumulative probability distribution
  double TotalIntensity; // total beam intensity
  
  ~beamscreen(); // destructor
  beamscreen(std::string dev_name); // constructor

  int Load(std::istream &fs); // load beamscreen from file
  int LoadData(string file_name); // Load screen array contents
  int SetDefault();

  // Generates a random point and energy on the screen
  vect3 RandomPoint(double &E, int &pol, double &w, randmt_t *rng);

 protected:
  //double PixelSizeX, PixelSizeY, PixelSurf; // pixel size and surface (cm2) 
  //vect3 *PixelX; // pixel coordinates array
  int PolarizedFlag; // unpolarized/polarized beam flag (0/1)
  int EnergyBinFlag; // energy binning flag (0/1)
  int InterpolFlag; // do not use / use interpolation inside pixels/bins (0/1)
  double Emin, Emax; // minimum and maximum bin energy

  int RunInit(); // beamscreen initialization before run

  // Generates a random point on the pixel surface
  //vect3 RandomPointOnPixel(int i);
  //vect3 RandomPointOnPixel(int i, randmt_t *rng);

  // Generates a random pixel on the screen
  int RandomPixel(int &iE, int &ix, int &iy, int &pol, randmt_t *rng);
  // evaluates the probability weigth for the trajectory x = x0 + u*t
  double WeightTrajectory(vect3 *x0, vect3 *u);
};

#endif
