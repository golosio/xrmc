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
//     xrmc_detector.h           //
//        31/01/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// detectorarray class definition
//

#ifndef DETECTORH
#define DETECTORH
#include "xrmc_device.h"
#include "xrmc_source.h"
#include "xrmc_arrayNd.h"

//  detectorarray class definition, member variables and functions
class detectorarray : public bodydevice
{
 public:
  int NX, NY, N; // number of rows (NY), column (NX) and pixels (NX x NY)
  double ***Image; // Acquired image array
  double ExpTime; // exposure time
  int PhotonNum; // Multiplicity of simulated events per detector pixel
  int NBins; // Num. of energy bins
  int ModeNum; // Num. of modes (scattering orders)
  
  ~detectorarray() {  // destructor
    if (PixelX!=NULL) delete[] PixelX;
    if (Image!=NULL) arrayNd::free_double_array3d(Image);
  }
  // constructor
  detectorarray(string dev_name) {
    PixelX=NULL;
    Image=NULL;
    NX = NY = N = PhotonNum = NBins = ModeNum = 0;
    xrmc_device(dev_name, "detectorarray");
  }
  int ImportDevice(xrmc_device_map *dev_map); // import device method
  int Load(FILE *fp); // load detector parameters, position, orientation
  int Save(string file_name); // save the acquired image in a file
  int SetDefault(); // set the default values for detector parameters
  int EventMulti(); // event multiplicity
  int Run() {return Acquisition();} // run the acquisition
  int Acquisition(); // run the acquisition
  int Clear(); // clear the detector pixel bin contents

 private:
  basesource *Source; // input device (typically the sample device)
  string SourceName; // name of the input device
  double PixelSizeX, PixelSizeY, PixelSurf; // pixel size and surface (cm2) 
  double dOmegaLim; // cut on solid angle from interaction point to pixel 
  int Shape; // pixel shape (0: rectangular, 1: elliptical)
  vect3 *PixelX; // pixel coordinates array
  int RandomPixelFlag; // flag to enable/disable random point on pixel
  int PoissonFlag; // flag to enable/disable Poisson statistic on counts
  int RoundFlag; // flag to enable/disable round pixel count to integer
  int HeaderFlag; // enable/disable writing header in output file
  int RunningFasterFlag; //columns(0) or rows(1) running faster
  int PixelType; // pixel content type:
                 // 0: fluence,      1: energy fluence,
                 // 2: fluence(E),   3: energy fluence(E)
  int SaturateEmin; // flag to saturate energies lower than Emin
  int  SaturateEmax; // flag to saturate energies greater than Emin
  double Emin, Emax; // minimum and maximum bin energy
  
  int Init(); // detectorarray initialization
  vect3 RandomPointOnPixel(int i); // Generates a random point
                                   // on the pixel surface
  double dOmega(vect3 DRp); // evaluates the factor dOmega, related to
  // the probability that the last photon trajectory crosses the pixel  
  int Poisson(); // generate uncertainty on pixel count using Poisson statistic

};

#endif
