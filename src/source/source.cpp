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
//     source.cpp                //
//        14/02/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Methods of the class source
//

#include <cmath>
#include <string>
#include <iostream>
#include "xrmc_algo.h"
#include "xrmc_math.h"
#include "xrmc_source.h"
#include "xrmc_exception.h"

using namespace std;
using namespace xrmc_algo;

// Constructor
source::source(string dev_name) {
  Runnable = false;
  NInputDevices=1;
  InputDeviceCommand.push_back("SpectrumName");
  InputDeviceDescription.push_back("Spectrum input device name");

  SetDevice(dev_name, "source");
}

source *costhl_Source; // temporary object needed to define the following
                       // function costhl

//////////////////////////////////////////////////////////////////////
// Function passed as input to the integration algorithm
// used to evaluate the solid angle 
//////////////////////////////////////////////////////////////////////
double one_min_costhl(double phi)
{
  return -costhl_Source->CosThL(phi) + 1.;
}


//////////////////////////////////////////////////////////////////////
// source run initialization method
//////////////////////////////////////////////////////////////////////
int source::RunInit()
{
  costhl_Source = this;
  Cos2Thx = cos(Thx)*cos(Thx);
  Sin2Thx = sin(Thx)*sin(Thx);
  Cos2Thy = cos(Thy)*cos(Thy);
  Sin2Thy = sin(Thy)*sin(Thy);
  if (Thx==0 || Thy==0) Omega = 0;
  else {
    Omega = Integrate(one_min_costhl, 0, 2*PI); // evaluate the solid angle
  }
  cout << "Omega: " << Omega << "\n";
  
  OrthoNormal(ui, uj, uk);  // evaluates uj to form a orthonormal basis
  
  return 0;
}

//////////////////////////////////////////////////////////////////////
// method for casting input device to type spectrum
/////////////////////////////////////////////////////////////////////
int source::CastInputDevices()
{
  // cast InputDevice[0] to type spectrum*
  Spectrum = dynamic_cast<spectrum*>(InputDevice[0]);
  if (Spectrum==0)
    throw xrmc_exception(string("Device ") + InputDeviceName[0]
			 + " cannot be casted to type spectrum\n");

  return 0;
}

//////////////////////////////////////////////////////////////////////
// method for linking input device
/////////////////////////////////////////////////////////////////////
/*
int source::LinkInputDevice(string command, xrmc_device *dev_pt)
{
  if (command=="SpectrumName") {
    // cast it to type spectrum*
    Spectrum = dynamic_cast<spectrum*>(dev_pt);
    if (Spectrum==0)
      throw xrmc_exception(string("Device cannot be casted to type"
				  " spectrum\n"));
  }
  else
    throw xrmc_exception(string("Unrecognized command: ") + command + "\n");

  return 0;
}
*/

// initialize loop on events
int source::Begin()
{
  return Spectrum->Begin(); // initialize spectrum loop on events
}


// next step of the event loop
int source::Next()
{
  return Spectrum->Next(); // next step of the spectrum event loop  
}

// check if the end of the loop is reached
bool source::End()
{
  return Spectrum->End(); // check if end of spectrum event loop is reached
}

// event multiplicity
int source::EventMulti()
{
  return Spectrum->EventMulti();
}

//////////////////////////////////////////////////////////////////////
// generate an event with a photon starting from the source
//////////////////////////////////////////////////////////////////////
int source::Out_Photon(photon *Photon)
{
  double E, w;
  int pol;

  // ask spectrum device to extrace the photon energy and polarization
  Spectrum->ExtractEnergy(&w, &E, &pol);
 // multiply the event weight by the total beam intensity
  Photon->w = w*Spectrum->TotalIntensity;
  Photon->E = E;
  Photon->x = X; // starting photon position is the source position
  if (SizeFlag != 0) {  // plus gaussian deviations
    if (rng == NULL)
      Photon->x += ui*Sigmax*GaussRnd() + uj*Sigmay*GaussRnd()
        + uk*Sigmaz*GaussRnd();
    else
      Photon->x += ui*Sigmax*GaussRnd_r(rng) + uj*Sigmay*GaussRnd_r(rng)
        + uk*Sigmaz*GaussRnd_r(rng);
  }
  // call source method for extracting photon initial direction
  PhotonDirection(Photon, pol);

  return 0;
}

//////////////////////////////////////////////////////////////////////
// method for extracting photon initial direction
//////////////////////////////////////////////////////////////////////
int source::PhotonDirection(photon *Photon, int pol)
{
  double phi, x, y, z;
  double cos_phi, sin_phi, cos_theta, sin_theta, cos_theta_lim;
  
  if (Thy==0) {
    sin_phi = 0;
    cos_phi = 1;
    double theta = (2.*(rng == NULL ? Rnd() : Rnd_r(rng))-1)*Thx;
    sin_theta=sin(theta);
    cos_theta=cos(theta);
  }
  else if (Thx==0) {
    sin_phi = 1;
    cos_phi = 0;
    double theta = (2.*(rng == NULL ? Rnd() : Rnd_r(rng))-1)*Thy;
    sin_theta=sin(theta);
    cos_theta=cos(theta);
  }
  else {
    phi = 2*PI*(rng == NULL ? Rnd() : Rnd_r(rng)); // random azimuthal angle phi (0,2*PI)
    cos_phi = cos(phi);
    sin_phi = sin(phi);
    cos_theta_lim = CosThL(phi); //maximum value of theta for this value of phi 
    Photon->w *= 2*PI*(1. - cos_theta_lim)/Omega;
    cos_theta = 1. - (rng == NULL ? Rnd() : Rnd_r(rng))*(1. - cos_theta_lim); // cosine of polar angle theta
    sin_theta = sqrt(1 - cos_theta*cos_theta); // sine of theta
  }
  x = sin_theta*cos_phi; // components of the photon direction
  y = sin_theta*sin_phi; // in the local coordinate system ui, uj, uk
  z = cos_theta;
  
 // photon direction in the absolute coordinate system
  Photon->uk = ui*x + uj*y + uk*z;

  // define local photon axis directions based on direction and polarization
  SetPhotonAxes(Photon, pol);

  return 0;
}

//////////////////////////////////////////////////////////////////////
// define local photon axis directions based on direction and polarization
//////////////////////////////////////////////////////////////////////
int source::SetPhotonAxes(photon *Photon, int pol)
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
double source::POmega(vect3 point_r)
{
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
}

//////////////////////////////////////////////////////////////////////
// Evaluate the minimum value of the cosine of the polar angle theta
// for a given value of the azimuthal angle phi
//////////////////////////////////////////////////////////////////////
double source::CosThL(double phi)
{
  double Cos2Phi, Sin2Phi, t1, t2;

  Cos2Phi = cos(phi);
  Cos2Phi *= Cos2Phi;
  Sin2Phi = sin(phi);
  Sin2Phi *= Sin2Phi;

  t1 = Cos2Thx*Cos2Phi;
  if ((t1 + Sin2Thx) == t1) return 1;
  t1 /= Sin2Thx;

  t2 = Cos2Thy*Sin2Phi;
  if ((t2 + Sin2Thy) == t2) return 1;
  t2 /= Sin2Thy;

  return sqrt((t1 + t2) / (t1 + t2 + 1.));
}

//////////////////////////////////////////////////////////////////////
// Evaluate the minimum value of the cosine of the polar angle theta
// for a given value of the direction x and y components
//////////////////////////////////////////////////////////////////////
double source::CosThLxy(double x, double y)
{
  double x2, y2, r2;
  double Cos2Phi, Sin2Phi, t1, t2;

  x2 = x*x;
  y2 = y*y;
  r2 = x2 + y2;
  if (r2==0) return 0;
  Cos2Phi = x2/r2;
  Sin2Phi = y2/r2;

  t1 = Cos2Thx*Cos2Phi;
  if ((t1 + Sin2Thx) == t1) return 1;
  t1 /= Sin2Thx;

  t2 = Cos2Thy*Sin2Phi;
  if ((t2 + Sin2Thy) == t2) return 1;
  t2 /= Sin2Thy;

  return sqrt((t1 + t2) / (t1 + t2 + 1.));
}

//////////////////////////////////////////////////////////////////////
// Generate an event with a photon starting from the source
// and forced to be directed toward the position x1
//////////////////////////////////////////////////////////////////////
int source::Out_Photon_x1(photon *Photon, vect3 x1)
{
  double E, w;
  int pol;
 
  Photon->x = X; // starting photon position is the source position
  if (SizeFlag != 0) {  // plus gaussian deviations
    if (rng == NULL)
      Photon->x += ui*Sigmax*GaussRnd() + uj*Sigmay*GaussRnd()
        + uk*Sigmaz*GaussRnd();
    else
      Photon->x += ui*Sigmax*GaussRnd_r(rng) + uj*Sigmay*GaussRnd_r(rng)
        + uk*Sigmaz*GaussRnd_r(rng);
  }
  // ask spectrum device to extrace the photon energy and polarization
  Spectrum->ExtractEnergy(&w, &E, &pol);
  Photon->E = E;

  Photon->uk = x1 - Photon->x; // photon direction
  Photon->uk.Normalize();

  // multiply the event weight by the total beam intensity
  // and by the probability that it has the direction uk per unit solid angle
  Photon->w = w*Spectrum->TotalIntensity*POmega(Photon->uk);

  // define local photon axis directions based on direction and polarization
  SetPhotonAxes(Photon, pol);

  return 0;
}

basesource *source::Clone(string dev_name) {
	//cout << "Entering source::Clone\n";
	source *clone = new source(dev_name);
	clone->Thx = Thx;
	clone->Thy = Thy;
	clone->Cos2Thx = Cos2Thx;
	clone->Sin2Thx = Sin2Thx;
	clone->Cos2Thy = Cos2Thy;
	clone->Sin2Thy = Sin2Thy;
	clone->Omega = Omega;
	clone->Sigmax = Sigmax;
	clone->Sigmay = Sigmay;
	clone->Sigmaz = Sigmaz;
	clone->SizeFlag = SizeFlag;
	clone->InputDeviceName[0] = InputDeviceName[0];
	clone->X = X;
	clone->ui = ui;
	clone->uj = uj;
	clone->uk = uk;
	//clone Spectrum
	clone->Spectrum = Spectrum->Clone(InputDeviceName[0]);


	return dynamic_cast<basesource*>(clone);
}

