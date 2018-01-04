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
//     anisotropicsource.cpp     //
//        16/05/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Methods of the class anisotropicsource
//

#include <cmath>
#include <string>
#include <iostream>
#include "xrmc_algo.h"
#include "xrmc_math.h"
#include "anisotropicsource.h"
#include "xrmc_exception.h"

using namespace std;
using namespace xrmc_algo;

// Constructor
anisotropicsource::anisotropicsource(string dev_name) {
  Runnable = false;
  NInputDevices=2;
  InputDeviceCommand.push_back("SpectrumName");
  InputDeviceDescription.push_back("Spectrum input device name");
  InputDeviceCommand.push_back("IntensityScreenName");
  InputDeviceDescription.push_back("intensityscreen input device name");
  Rng = NULL;
  InvertScreenAndSourceDirection = 0;

  SetDevice(dev_name, "anisotropicsource");
}

int anisotropicsource::SetRng(randmt_t *rng)
{
  Rng = rng;
  Spectrum->SetRng(rng);
  IntensityScreen->SetRng(rng);

  return 0;
}

int anisotropicsource::RunInit()
{
  PhCFlag = false;
  OrthoNormal(ui, uj, uk);  // evaluates uj to form a orthonormal basis
  
  return 0;
}

//////////////////////////////////////////////////////////////////////
// method for casting input devices
/////////////////////////////////////////////////////////////////////
int anisotropicsource::CastInputDevices()
{
  // cast InputDevice[0] to type spectrum*
  Spectrum = dynamic_cast<spectrum*>(InputDevice[0]);
  if (Spectrum==0)
    throw xrmc_exception(string("Device ") + InputDeviceName[0]
			 + " cannot be casted to type spectrum\n");

  IntensityScreen = dynamic_cast<intensityscreen*>(InputDevice[1]);
  if (IntensityScreen==0)
    throw xrmc_exception(string("Device ") + InputDeviceName[1]
			 + " cannot be casted to type intensityscreen\n");

  return 0;
}


// initialize loop on events
int anisotropicsource::Begin()
{
  return Spectrum->Begin(); // initialize spectrum loop on events
}


// next step of the event loop
int anisotropicsource::Next()
{
  return Spectrum->Next(); // next step of the spectrum event loop  
}

// check if the end of the loop is reached
bool anisotropicsource::End()
{
  return Spectrum->End(); // check if end of spectrum event loop is reached
}

// event multiplicity
long long anisotropicsource::EventMulti()
{
  return Spectrum->EventMulti();
}

//////////////////////////////////////////////////////////////////////
// generate an event with a photon starting from the source
//////////////////////////////////////////////////////////////////////
int anisotropicsource::Out_Photon(photon *Photon)
{
  double E, w0, w1;
  int pol;

  // ask spectrum device to extrace the photon energy and polarization
  Spectrum->ExtractEnergy(&w0, &E, &pol);
 // multiply the event weight by the total beam intensity

  Photon->E = E;
  Photon->x = X; // starting photon position is the source position
  if (SizeFlag != 0 && !PhCFlag) {  // plus gaussian deviations
    //DELETE if (Rng == NULL)
    //  Photon->x += ui*Sigmax*GaussRnd() + uj*Sigmay*GaussRnd()
    //    + uk*Sigmaz*GaussRnd();
    //else
      Photon->x += ui*Sigmax*GaussRnd_r(Rng) + uj*Sigmay*GaussRnd_r(Rng)
        + uk*Sigmaz*GaussRnd_r(Rng);
  }
  // ask intensityscreen device to extract photon endpoint
  vect3 r = IntensityScreen->RandomPoint(w1, Rng);
  // evaluates photon direction
  if (InvertScreenAndSourceDirection) {
    Photon->uk = Photon->x - r;
    Photon->x = r; // photon position in intensityscreen!
  }
  else
    Photon->uk = r - Photon->x;

  Photon->w = w0*w1*Spectrum->TotalIntensity;
  // define local photon axis directions based on direction and polarization
  SetPhotonAxes(Photon, pol);

  return 0;
}

//////////////////////////////////////////////////////////////////////
// define local photon axis directions based on direction and polarization
//////////////////////////////////////////////////////////////////////
int anisotropicsource::SetPhotonAxes(photon *Photon, int pol)
{
  // polarization vector
  if (pol == 0) { // (local) x polarized
    Photon->ui = ui;
  }
  else { // (local) y polarized
    Photon->ui = uj;
  }
  // evaluates uj to form a orthonormal basis
  OrthoNormal(Photon->ui, Photon->uj, Photon->uk);

  return 0;
}


//////////////////////////////////////////////////////////////////////
// Generate an event with a photon starting from the source
// and forced to be directed toward the position x1
//////////////////////////////////////////////////////////////////////
int anisotropicsource::Out_Photon_x1(photon *Photon, vect3 x1)
{
  double E, w0, w1;
  int pol;
 
  Photon->x = X; // starting photon position is the source position
  if (SizeFlag != 0) {  // plus gaussian deviations
    Photon->x += ui*Sigmax*GaussRnd_r(Rng) + uj*Sigmay*GaussRnd_r(Rng)
      + uk*Sigmaz*GaussRnd_r(Rng);
  }
  // ask spectrum device to extrace the photon energy and polarization
  Spectrum->ExtractEnergy(&w0, &E, &pol);
  Photon->E = E;

  vect3 vr = x1 - Photon->x; //relative position
  vr.Normalize(); // normalized direction
  Photon->uk = vr; // update the photon direction

  IntensityScreen->DirectionWeight(Photon->x, vr, w1, Rng, InvertScreenAndSourceDirection);

  // multiply the event weight by the total beam intensity
  Photon->w = w0*w1*Spectrum->TotalIntensity;
  //if (sqrt(vr.Elem[0]*vr.Elem[0]+vr.Elem[2]*vr.Elem[2])<0.12/115) {
  //cout << "w0 " << w0 << endl;
  //cout << "w1 " << w1 << endl;
  //cout << Spectrum->TotalIntensity << endl;
    //}

  // define local photon axis directions based on direction and polarization
  SetPhotonAxes(Photon, pol);

  return 0;
}

basesource *anisotropicsource::Clone(string dev_name) {
	//cout << "Entering anisotropicsource::Clone\n";
	anisotropicsource *clone = new anisotropicsource(dev_name);
	clone->Sigmax = Sigmax;
	clone->Sigmay = Sigmay;
	clone->Sigmaz = Sigmaz;
	clone->SizeFlag = SizeFlag;
	clone->InputDeviceName[0] = InputDeviceName[0];
	clone->InputDeviceName[1] = InputDeviceName[1];
	clone->X = X;
	clone->ui = ui;
	clone->uj = uj;
	clone->uk = uk;

	clone->Spectrum = Spectrum->Clone(InputDeviceName[0]);
	clone->IntensityScreen = IntensityScreen->Clone(InputDeviceName[1]);;
	clone->InvertScreenAndSourceDirection = InvertScreenAndSourceDirection;

	return dynamic_cast<basesource*>(clone);
}

int anisotropicsource::PhCOn()
{
  PhCFlag=true;
  Spectrum->PhCOn();

  return 0;
}

int anisotropicsource::PhCOff()
{
  PhCFlag=false;
  Spectrum->PhCOff();
  
  return 0;
}

double anisotropicsource::GetPhC_E0()
{
  return Spectrum->GetPhC_E0();
}
