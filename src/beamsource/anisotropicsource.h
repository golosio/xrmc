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
//     anisotropicsource.h       //
//        16/05/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// anisotropicsource class definition
//
#ifndef ANISOTROPICSOURCEH
#define ANISOTROPICSOURCEH

#include <string>
#include "xrmc_spectrum.h"
#include "xrmc_photon.h"
#include "xrmc_device.h"
#include "xrmc_math.h"
#include "randmt.h"
#include "xrmc_source.h"
#include "intensityscreen.h"

//////////////////////////////////////////////////////////////////////
// anisotropicsource class definition
//////////////////////////////////////////////////////////////////////
class anisotropicsource : public basesource
{
 public:
  double Sigmax, Sigmay, Sigmaz; // source size in local coordinate system
  int SizeFlag;
  spectrum *Spectrum; // input spectrum device
  intensityscreen *IntensityScreen; // input intensityscreen device

  // Constructor
  anisotropicsource(std::string dev_name);

  virtual int Load(istream &fs); // loads source parameters from file
  // method for casting input devices
  virtual int CastInputDevices();
  // method for linking input devices
  //int LinkInputDevice(string command, xrmc_device *dev_pt);
  virtual int Begin(); // begin event loop method
  virtual int Next(); // next event method
  virtual bool End(); // check for end event loop method
  virtual long long EventMulti(); // event multiplicity
  virtual int SetDefault(); // set default values for source parameters
  virtual int Out_Photon(photon *Photon); // generate an event
  // generate an event with a photon directed toward the position x1
  virtual int Out_Photon_x1(photon *Photon, vect3 x1);
  virtual double GetPhC_E0();
  virtual int PhCOn();
  virtual int PhCOff();

  virtual basesource *Clone(string dev_name);
  virtual int RunInit(); // source run initialization method
 // set the random number generator structure
  virtual int SetRng(randmt_t *rng);
  virtual vect3 SourceX() {return X;}
 private:
  // extract the initial direction of a photon produced by the source
  //int PhotonDirection(photon *Photon, int pol);
  // build the photon local axis based on its direction and polarization
  int SetPhotonAxes(photon *Photon, int pol);
  // probability per unit solid angle that a photon has direction vr
  //double POmega(vect3 vr);

};

#endif
