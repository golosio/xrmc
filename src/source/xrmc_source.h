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
//     xrmc_source.h             //
//        14/02/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// basesource and source classes definition
//
#ifndef SOURCEH
#define SOURCEH

#include <string>
#include "xrmc_spectrum.h"
#include "xrmc_photon.h"
#include "xrmc_device.h"
#include "xrmc_math.h"
#include "randmt.h"

// generic (virtual) class for sources
class basesource : public bodydevice
{
 public:
  virtual ~basesource() {};
  virtual int ModeNum() {return 1;} // number of modes
  virtual int Out_Photon(photon *Photon) {return 0;} // generate an event
  // generate an event with a photon directed toward the position x1
  virtual int Out_Photon_x1(photon *Photon, vect3 x1) {return 0;}
 // generate an event for a specified mode
  virtual int Out_Photon(photon *Photon, int *ModeIdx)
    {*ModeIdx=0; return Out_Photon(Photon);}
  // generate an event for a specified mode with a photon directed toward
  // the position x1
  virtual int Out_Photon_x1(photon *Photon, vect3 x1, int *ModeIdx)
    {*ModeIdx=0; return Out_Photon_x1(Photon, x1);}  
  virtual basesource *Clone(string dev_name) {return NULL;};
};

//////////////////////////////////////////////////////////////////////
// source class definition
//////////////////////////////////////////////////////////////////////
class source : public basesource
{
 public:
  double Thx, Thy; // Beam divergence (thetax, thetay)
  double Cos2Thx, Sin2Thx, Cos2Thy, Sin2Thy; //related goniometric functions
  double Omega; // source aperture solid angle  
  double Sigmax, Sigmay, Sigmaz; // source size in local coordinate system
  int SizeFlag;
  randmt_t *rng;
  spectrum *Spectrum; // input spectrum device
  string SpectrumName; // name of input spectrum device

  // Constructor
  source(string dev_name) {
	rng = NULL;
  	xrmc_device(dev_name, "source");
  }
  int Load(FILE *fp); // method for loading source parameters from file
 // method for importing input spectrum device
  int ImportDevice(xrmc_device_map *dev_map);
  int Begin(); // begin event loop method
  int Next(); // next event method
  bool End(); // check for end event loop method
  int EventMulti(); // event multiplicity
  int SetDefault(); // set default values for source parameters
  int Out_Photon(photon *Photon); // generate an event
  // generate an event with a photon directed toward the position x1
  int Out_Photon_x1(photon *Photon, vect3 x1);
  // maximum value of polar angle theta for a specified value of phi
  double CosThL(double phi);
  basesource *Clone(string dev_name);

 private:

  // extract the initial direction of a photon produced by the source
  int PhotonDirection(photon *Photon, int pol);
  // build the photon local axis based on its direction and polarization
  int SetPhotonAxes(photon *Photon, int pol);
  // probability per unit solid angle that a photon has direction vr
  double POmega(vect3 vr);
  // maximum value of polar angle theta for specified x, y direction components
  double CosThLxy(double x, double y);

};

#endif
