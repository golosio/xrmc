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
#include "phcdevice.h"

//  beamscreen class definition, member variables and functions
class beamscreen : public xrmc_screen, public phcdevice
{
 private:
  int PhotonNum; // num. of events to be extracted for each interval
  int PhotonIdx; // index of the event extracted for the interval
  int PolIdx; // index of polarization type (0: x, 1: y)
  int EnergyIdx; // index of the energy bin in loop mode

 public:
  //int NX, NY, N; // number of rows (NY), column (NX) and pixels (NX x NY)
  int NBins; // Num. of energy bins
  double *Image; // image array
  double *CumulImage; // cumulative probability distribution for x, y, E
  double *CumulEnergy; // cumulative probability distribution for energy bins
  double *SumEnergyImage; // probability distribution integrated over the energy
  double *BinWeight; // energy bin weight
  double *CumulXY; // cumulative function for position x, y on the beamscreen

  double TotalIntensity; // total beam intensity
  double PhC_rE0; // random number used to set energy in PhC mode

  ~beamscreen(); // destructor
  beamscreen(std::string dev_name); // constructor

  virtual int Load(std::istream &fs); // load beamscreen from file
  int LoadData(string file_name); // Load screen array contents
  virtual int SetDefault();

  // Generates a random point and energy on the screen
  vect3 RandomPoint(double &E, int &pol, double &w, randmt_t *rng);
  // Generates a random energy and polarization for a given trajectory
  bool RandomEnergy(vect3 x0, vect3 u, double &E, int &pol, double &w,
		    randmt_t *rng);
  vect3 RandomXY(double &E, double &w, randmt_t *rng);
  virtual double GetPhC_E0();

  virtual int Begin(); // begin the event loop
  virtual int Next();  // next step of the event loop  
  virtual bool End(); // check if the end of the event loop is reached
  virtual long long EventMulti(); // event multiplicity
  // int RunInit(); // initialize before run
  virtual beamscreen *Clone(std::string dev_name);

 protected:
  //double PixelSizeX, PixelSizeY, PixelSurf; // pixel size and surface (cm2) 
  //vect3 *PixelX; // pixel coordinates array
  int PolarizedFlag; // unpolarized/polarized beam flag (0/1)
  int LoopFlag; // flag for loop on all intervals of the spectrum
  int EnergyBinFlag; // energy binning flag (0/1)
  int InterpolFlag; // do not use / use interpolation inside pixels/bins (0/1)
  double Emin, Emax; // minimum and maximum bin energy

  virtual int RunInit(); // beamscreen initialization before run

  // Generates a random point on the pixel surface
  //vect3 RandomPointOnPixel(int i);
  //vect3 RandomPointOnPixel(int i, randmt_t *rng);

  // Generates a random pixel on the screen
  int RandomPixel(int &iE, int &ix, int &iy, int &pol, randmt_t *rng);
  // evaluates the probability weigth for the trajectory x = x0 + u*t
  double WeightTrajectory(vect3 *x0, vect3 *u);
  // evaluates interpolation weight factor
  double InterpolWeight(int iE, int ix, int iy, int pol,
			double rE, double rx, double ry, int ene_flag);

};

#endif
