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
#include "xrmc_screen.h"
#include "randmt.h"
#include "xrmc_photon.h"
#include <vector>
#include "image_convolution.h"

//  detectorarray class definition, member variables and functions
class detectorarray : public xrmc_screen
{
 public:
  //int NX, NY, N; // number of rows (NY), column (NX) and pixels (NX x NY)
  double ***Image; // Acquired image array
  double ***ConvolutedImage; // convoluted image array
  double ExpTime; // exposure time
  int PhotonNum; // Multiplicity of simulated events per detector pixel
  int NBins; // Num. of energy bins
  int ModeNum; // Num. of modes (scattering orders)
  double Z12; // object-detector distance used for source size convolution
  std::vector<unsigned long> Seeds;

  virtual ~detectorarray(); // destructor
  detectorarray(std::string dev_name); // constructor

  virtual int CastInputDevices(); // cast input device method
  // method for linking input device
  //int LinkInputDevice(string command, xrmc_device *dev_pt);
  virtual int Load(std::istream &fs); // load detector parameters, position, orientation
 // save the acquired image in a file
  virtual int SaveData(std::string data_name, std::string file_name);
  virtual int SaveData(double ***image, string file_name);

 // save the acquired image in a ascii file
  virtual int SaveAsciiData(double ***image, std::string file_name);
  virtual int SetDefault(); // set the default values for detector parameters
  virtual long long EventMulti(); // event multiplicity
  virtual int Run() {return Acquisition();} // run the acquisition
  virtual int Acquisition(); // run the acquisition
  virtual int Clear(); // clear the detector pixel bin contents

 protected:
  basesource *Source; // input device (typically the sample device)
  //string SourceName; // name of the input device
  ///double PixelSizeX, PixelSizeY, PixelSurf; // pixel size and surface (cm2) 
  //double dOmegaLim; // cut on solid angle from interaction point to pixel 
  //int Shape; // pixel shape (0: rectangular, 1: elliptical)
  //vect3 *PixelX; // pixel coordinates array
  //int RandomPixelFlag; // flag to enable/disable random point on pixel
  int PoissonFlag; // flag to enable/disable Poisson statistic on counts
  int RoundFlag; // flag to enable/disable round pixel count to integer
  int HeaderFlag; // enable/disable writing header in output file
  int AsciiFlag; //  // binary(0) or ascii(1) output file format
  int ConvolveFlag; // generates convoluted image (0/1)
  //int RunningFasterFlag; //columns(0) or rows(1) running faster
  int PixelType; // pixel content type:
                 // 0: fluence,      1: energy fluence,
                 // 2: fluence(E),   3: energy fluence(E)
  int ForceDetectFlag; // Last photon not forced / forced to be detected
  int SaturateEmin; // flag to saturate energies lower than Emin
  int  SaturateEmax; // flag to saturate energies greater than Emin
  double Emin, Emax; // minimum and maximum bin energy
  int NBx, NBy; // size of border (in pixel) used for FFT in convolution
  gauss_vect GaussPSFx, GaussPSFy; // sum-of-gaussians model of detector PSF 
  std::vector<gauss_vect> GaussPSFxBin; // energy-dependent PSF (x component)
  std::vector<gauss_vect> GaussPSFyBin; // energy-dependent PSF (y component)
  gauss_vect GaussSourceX, GaussSourceY; // source size sum-of-gaussians model 
  std::vector<gauss_vect> GaussSourceXBin; // energy-dependent source x size
  std::vector<gauss_vect> GaussSourceYBin; // energy-dependent source y size
  int EfficiencyFlag; // flag for using efficiency before image convolution
  std::vector<double> Efficiency; //  energy-dependent efficiency
  std::vector<int> PrevCompZ;

  virtual int RunInit(); // detectorarray initialization before run
  //vect3 RandomPointOnPixel(int i); // Generates a random point
  //vect3 RandomPointOnPixel(int i, randmt_t *rng);
                                   // on the pixel surface
  //double dOmega(vect3 DRp); // evaluates the factor dOmega, related to
  // the probability that the last photon trajectory crosses the pixel  
  //int Poisson(); // generate uncertainty on pixel count using Poisson statistic
  int Poisson(randmt_t *rng); // generate uncertainty on pixel count using Poisson statistic

// run the forced/unforced acquisition
  int ForcedAcquisition(basesource **SourceClones, photon *PhotonArray,
			randmt_t **rngs);
  int UnforcedAcquisition(basesource **SourceClones, photon *PhotonArray);
  int Convolve();
  bool ParsePSFCommand(istream &fs, string comm);
  int LoadPSFBin(int n, vector<gauss_vect> &PSF_bin, istream &psf_fs);
  int LoadEfficiency(istream &eff_fs);

};

#endif
