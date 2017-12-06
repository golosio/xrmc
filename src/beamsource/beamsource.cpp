/*
Copyright (C) 2017 Bruno Golosio

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
//     beamsource.cpp            //
//        06/12/2017             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Methods of the class beamsource
//

#include <cmath>
#include <string>
#include <iostream>
#include "xrmc_algo.h"
#include "xrmc_math.h"
#include "beamsource.h"
#include "xrmc_exception.h"

using namespace std;
using namespace xrmc_algo;

// Constructor
beamsource::beamsource(string dev_name) {
  Runnable = false;
  NInputDevices=1;
  InputDeviceCommand.push_back("BeamScreenName");
  InputDeviceDescription.push_back("beamscreen input device name");
  Rng = NULL;

  SetDevice(dev_name, "beamsource");
}

//////////////////////////////////////////////////////////////////////
// beamsource run initialization method
//////////////////////////////////////////////////////////////////////
int beamsource::RunInit()
{
  PhCFlag = false;  
  OrthoNormal(ui, uj, uk);  // evaluates uj to form a orthonormal basis
  
  return 0;
}

//////////////////////////////////////////////////////////////////////
// method for casting input device to type beamscreen*
/////////////////////////////////////////////////////////////////////
int beamsource::CastInputDevices()
{
  // cast InputDevice[0] to type beamscreen*
  BeamScreen = dynamic_cast<beamscreen*>(InputDevice[0]);
  if (BeamScreen==0)
    throw xrmc_exception(string("Device ") + InputDeviceName[0]
			 + " cannot be casted to type beamscreen\n");

  return 0;
}

//////////////////////////////////////////////////////////////////////
// generate an event with a photon starting from the source
//////////////////////////////////////////////////////////////////////
int beamsource::Out_Photon(photon *Photon)
{
  double E, w;
  int pol;
  vect3 r;

  // ask beamscreen device to extract photon endpoint, energy and polarization
  r = BeamScreen->RandomPoint(E, pol, w, Rng);

  // multiply the event weight by the total beam intensity
  Photon->w = w*BeamScreen->TotalIntensity;
  Photon->E = E;
  Photon->x = X; // starting photon position is the source position
  if (SizeFlag != 0 && !PhCFlag) {  // plus gaussian deviations
    //if (Rng == NULL)
    //  Photon->x += ui*Sigmax*GaussRnd() + uj*Sigmay*GaussRnd()
    //    + uk*Sigmaz*GaussRnd();
    //else
      Photon->x += ui*Sigmax*GaussRnd_r(Rng) + uj*Sigmay*GaussRnd_r(Rng)
        + uk*Sigmaz*GaussRnd_r(Rng);
  }

  // evaluates photon direction
  Photon->uk = r - Photon->x; /////////////////////////////////////////////
  ///////////////////////// change to give also the possibility to use r - X

  // define local photon axis directions based on direction and polarization
  SetPhotonAxes(Photon, pol);

  //double x = Photon->uk.Elem[0]*115;
  //double y = Photon->uk.Elem[2]*115;
  //cout << "STORE " << x << " " << y << " " << Photon->E << " " << w << endl;

  return 0;
}

//////////////////////////////////////////////////////////////////////
// define local photon axis directions based on direction and polarization
//////////////////////////////////////////////////////////////////////
int beamsource::SetPhotonAxes(photon *Photon, int pol)
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
int beamsource::Out_Photon_x1(photon *Photon, vect3 x1)
{
  double E, w;
  int pol;

  Photon->x = X; // starting photon position is the source position
  if (SizeFlag != 0) {  // plus gaussian deviations
      Photon->x += ui*Sigmax*GaussRnd_r(Rng) + uj*Sigmay*GaussRnd_r(Rng)
        + uk*Sigmaz*GaussRnd_r(Rng);
  }
  vect3 vr = x1 - Photon->x; //relative position
  double r = vr.Mod();
  if (r<Rlim) r = Rlim;
  vr.Normalize(); // normalized direction
  Photon->uk = vr; // update the photon direction

  if (BeamScreen->RandomEnergy(Photon->x, vr, E, pol, w, Rng)) {
    Photon->E = E;
    // multiply the event weight by the total beam intensity
    Photon->w = w*BeamScreen->TotalIntensity/(r*r);
    // define local photon axis based on direction and polarization
    SetPhotonAxes(Photon, pol);
  }
  else {
    Photon->w = 0;
    pol = 0;
    SetPhotonAxes(Photon, pol);
  }    

  return 0;
}

// initialize loop on events
int beamsource::Begin()
{
  return BeamScreen->Begin(); // initialize spectrum loop on events
}


// next step of the event loop
int beamsource::Next()
{
  return BeamScreen->Next(); // next step of the spectrum event loop  
}

// check if the end of the loop is reached
bool beamsource::End()
{
  return BeamScreen->End(); // check if end of spectrum event loop is reached
}

// event multiplicity
long long beamsource::EventMulti()
{
  return BeamScreen->EventMulti();
}

basesource *beamsource::Clone(string dev_name) {
	//cout << "Entering beamsource::Clone\n";
	beamsource *clone = new beamsource(dev_name);
	clone->Sigmax = Sigmax;
	clone->Sigmay = Sigmay;
	clone->Sigmaz = Sigmaz;
	clone->SizeFlag = SizeFlag;
	clone->InputDeviceName[0] = InputDeviceName[0];
	clone->X = X;
	clone->ui = ui;
	clone->uj = uj;
	clone->uk = uk;
	clone->BeamScreen = BeamScreen->Clone(InputDeviceName[0]);;

	return dynamic_cast<basesource*>(clone);
}

int beamsource::SetRng(randmt_t *rng)
{
  Rng = rng;
  BeamScreen->SetRng(rng);

  return 0;
}

int beamsource::PhCOn()
{
  PhCFlag=true;
  BeamScreen->PhCOn();

  return 0;
}

int beamsource::PhCOff()
{
  PhCFlag=false;
  BeamScreen->PhCOff();
  
  return 0;
}

double beamsource::GetPhC_E0()
{
  return BeamScreen->GetPhC_E0();
}

