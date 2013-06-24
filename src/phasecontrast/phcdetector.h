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

#ifndef PHCDETECTORH
#define PHCDETECTORH
#include "xrmc_device.h"
#include "xrmc_source.h"
#include "xrmc_arrayNd.h"
#include "xrmc_screen.h"
#include "randmt.h"
#include "xrmc_photon.h"
#include "xrmc_detector.h"
#include <vector>

//  phcdetector class definition, member variables and functions
class phcdetector : public detectorarray
{
 public:
  int PhasePhotonNum;
  int ImageNx, ImageNy;
  int NScreenBorderX, NScreenBorderY, NInterpBorderX, NInterpBorderY, NBx, NBy;
  long unsigned int NN1[2];
  double ***Transm, ***Propagator;
  //double **PhaseImage;
  double ***ImageWeight;
  double Csix, Csiy, Csiz;
  double Lx1, Ly1, Sigmax1, Sigmay1;
  vect3 Rv2;
  double R1, R2, R12, Reff, X2, Y2, Z2, Z1;
  double dx1, dy1;
  double L1Coeff, Sigma1Coeff;

  virtual ~phcdetector(); // destructor
  phcdetector(std::string dev_name); // constructor

  virtual int CastInputDevices(); // cast input device method
  // method for linking input device
  //int LinkInputDevice(string command, xrmc_device *dev_pt);
  virtual int Load(std::istream &fs); // load detector parameters, position, orientation
 // save the acquired image in a file
  //virtual int SaveData(std::string data_name, std::string file_name);
 // save the acquired image in a ascii file
  //virtual int SaveAsciiData(std::string data_name, std::string file_name);
  virtual int SetDefault(); // set the default values for detector parameters
  virtual long long EventMulti(); // event multiplicity
  virtual int Clear(); // clear the detector pixel bin contents
  virtual int Run() {return Acquisition();} // run the acquisition
  virtual int Acquisition(); // run the acquisition
  virtual int SaveData(std::string data_name, std::string file_name);
  virtual int CloneBegin(); // begin event loop on clones
  virtual int CloneNext(); // next step on clone event loop
  virtual bool CloneEnd(); // end clone event loop

 protected:
  sample *Sample; // input device (sample device)
  
  virtual int RunInit(); // detectorarray initialization before run

  // run the acquisition
  int PhCAcquisition(randmt_t **rngs);
  int FillBorders(double **transm);
  int CloneAcquisition(int ie, int thread_idx, double **transm, double **propag,
		       randmt_t **rngs);

 private:
  sample **SampleClones;
  photon *PhotonArray;

};

#endif
