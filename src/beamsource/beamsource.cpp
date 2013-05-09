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
//     beamsource.cpp            //
//        02/05/2013             //
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
  
  OrthoNormal(ui, uj, uk);  // evaluates uj to form a orthonormal basis
  
  return 0;
}

//////////////////////////////////////////////////////////////////////
// method for casting input device to type spectrum
/////////////////////////////////////////////////////////////////////
int beamsource::CastInputDevices()
{
  // cast InputDevice[0] to type spectrum*
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

  // ask beamscreen device to extract photon endpoint, energy and polarization
  vect3 r = BeamScreen->RandomPoint(E, pol, w, Rng);
 // multiply the event weight by the total beam intensity
  Photon->w = w*BeamScreen->TotalIntensity;
  Photon->E = E;
  Photon->x = X; // starting photon position is the source position
  if (SizeFlag != 0) {  // plus gaussian deviations
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

  double x = Photon->uk.Elem[0]*115;
  double y = Photon->uk.Elem[2]*115;

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
// probability that a photon produced by the source has the direction
// of point_r per unit solid angle (uniform distribution)
//////////////////////////////////////////////////////////////////////
double beamsource::POmega(vect3 point_r)
{
  /*
  double x, y, z, r;
  double cos_th, cos_th_l;

  if (Omega==0) return 0;
  r = point_r.Mod(); // vector module
  if (r==0) return 0;
  // vector components in the local source coordinate system
  x = point_r*ui;
  y = point_r*uj;
  z = point_r*uk;
  cos_th_l = CosThLxy(x, y); // maximum value of theta
  cos_th = z / r;            // actual value of theta
 // check that theta < maximum 
  if (cos_th < cos_th_l) return 0; // if not, weight is zero
  else return 1./Omega; // otherwise weight is 1 / Omega (uniform distribution)
  */
  return 1;
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
    //if (Rng == NULL)
    //  Photon->x += ui*Sigmax*GaussRnd() + uj*Sigmay*GaussRnd()
    //    + uk*Sigmaz*GaussRnd();
    //else
      Photon->x += ui*Sigmax*GaussRnd_r(Rng) + uj*Sigmay*GaussRnd_r(Rng)
        + uk*Sigmaz*GaussRnd_r(Rng);
  }
  // ask spectrum device to extrace the photon energy and polarization
  //Spectrum->ExtractEnergy(&w, &E, &pol);
  E=10;
  w=1;
  pol=0;
  Photon->E = E;
  Photon->uk = x1 - Photon->x; // photon direction
  Photon->uk.Normalize();

  // multiply the event weight by the total beam intensity
  // and by the probability that it has the direction uk per unit solid angle
  Photon->w = w*BeamScreen->TotalIntensity*POmega(Photon->uk);

  // define local photon axis directions based on direction and polarization
  SetPhotonAxes(Photon, pol);

  return 0;
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
	clone->BeamScreen = BeamScreen;

	return dynamic_cast<basesource*>(clone);
}

int beamsource::SetRng(randmt_t *rng)
{
  Rng = rng;
  BeamScreen->SetRng(rng);

  return 0;
}
