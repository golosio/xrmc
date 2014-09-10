/*
Copyright (C) 2013 Tom Schoonjans and Bruno Golosio 

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

#include <iostream>
#include <string>
#include "xrmc.h"
#include "xrmc_gettoken.h"
#include "xrmc_exception.h"
#include "xrmc_math.h"
#include "xrmc_detector.h"
#include "xrmc_detectorconvolute.h"


using namespace std;
using namespace gettoken;

int detectorconvolute::Load(istream &fs) {
  int i;
  vect3 x0, u;
  matr3 R;
  double theta;
  string comm="";
  //char s[MAXSTRLEN];
  string s;
  cout << "XMI-MSIM detector convolution parameters\n";

  // get a command/variable name from input file
  while (GetToken(fs, comm)) {
    // parse the command and decide what to do
    //
    // check if it's a command for setting an input device name
    if (ParseInputDeviceCommand(fs, comm)) continue;
    else if (comm == "CrystalPhase") {
      GetToken(fs, s);
      CrystalPhaseName = s;
      cout << "CrystalPhase name: " << CrystalPhaseName << "\n";
    }
    else if (comm == "WindowPhase") {
      GetToken(fs, s);
      WindowPhaseName = s;
      cout << "WindowPhase name: " << WindowPhaseName << "\n";
    }
    else if (comm == "PulseWidth") {
      GetDoubleToken(fs, &(xd->pulse_width));
      cout << "PulseWidth: " << xd->pulse_width << "\n";
    }
    else if (comm == "FanoFactor") {
      GetDoubleToken(fs, &(xd->fano));
      cout << "Fano factor: " << xd->fano << "\n";
    }
    else if (comm == "Noise") {
      GetDoubleToken(fs, &(xd->noise));
      cout << "Detector noise: " << xd->noise << "\n";
    }
    else if (comm == "CrystalThickness") {
      GetDoubleToken(fs, &CrystalThickness);
      cout << "Crystal thickness: " << CrystalThickness << "\n";
    }
    else if (comm == "WindowThickness") {
      GetDoubleToken(fs, &WindowThickness);
      cout << "Window thickness: " << WindowThickness << "\n";
    }
    /*else if(comm=="SourceName") { // set the source input device name
      GetToken(fs, s);
      InputDeviceName[0] = s;
      cout << "Source input device name: " << SourceName << "\n";
    }*/
    else if(comm=="NPixels") { // set the number of pixels (rows and columns)
      GetIntToken(fs, &NX);
      GetIntToken(fs, &NY);
      cout << "Pixel number (NX x NY): " << NX << " x " << NY << "\n";
    }
    else if(comm=="PixelSize") { // set the pixel size in x and y directions
      GetDoubleToken(fs, &PixelSizeX);
      GetDoubleToken(fs, &PixelSizeY);
      cout << "Pixel size (Sx x Sy) (cm): " << PixelSizeX << " x "
	   <<  PixelSizeY << "\n";
    }
    else if(comm=="Shape") { // detector elements shape
      // (rectangular or elliptical)
      GetIntToken(fs, &Shape);
      cout << "Detector elements shape (0 rectangular/1 elliptical): "
	   << Shape << "\n";
    }
    else if(comm=="dOmegaLim") { // set the cut on solid angle
      GetDoubleToken(fs, &dOmegaLim);
      if (dOmegaLim==0) dOmegaLim=2*PI; // if 0 set it to the default value
      cout << "Cut on solid angle from interaction point to a single pixel: "
	   << dOmegaLim << "\n";
    }
    else if(comm=="X") { // set the detector center coordinates
      cout << "Detector array center position :\t"; 
      for (i=0; i<3; i++) {
	GetDoubleToken(fs, &X.Elem[i]);
	//cout << X[i] << "\t";
      }
      cout << X << endl;
      //cout << "\n";
    }
    else if(comm=="uk") { // direction of the normal to the detector surface uk
      cout << "Detector orientation :\n"; 
      cout << "\tuk vector:\t"; 
      for (i=0; i<3; i++) {
	GetDoubleToken(fs, &uk.Elem[i]);
	//cout << uk[i] << "\t";
      }
      cout << uk << endl;
      //cout << "\n";
    }	
    else if(comm=="ui") { // direction of ui, parallel to the detector rows
      cout << "Detector orientation :\n"; 
      cout << "\tui vector:\t"; 
      for (i=0; i<3; i++) {
	GetDoubleToken(fs, &ui.Elem[i]);
	//cout << ui[i] << "\t";
      }
      cout << ui << endl;
      //cout << "\n";
    }	
    else if(comm=="ExpTime") { // set the Exposure Time
      GetDoubleToken(fs, &ExpTime);
      cout << "Exposure Time (s): " << ExpTime << "\n";
    }	
    else if(comm=="PhotonNum") { // Multiplicity of simulated events per pixel
      GetIntToken(fs, &PhotonNum);
      cout << "Multiplicity of simulated events per pixel: "
	   << PhotonNum << "\n";
    }
    else if(comm=="RandomPixelFlag") { // enable/disable random point on pixels
      GetIntToken(fs, &RandomPixelFlag);
      cout << "Enable random point on pixels (0/1): " << RandomPixelFlag
	   << "\n";
    }
    else if(comm=="PoissonFlag") { // flag to enable/disable Poisson statistic
      GetIntToken(fs, &PoissonFlag);
      cout << "Enable Poisson statistic on pixel counts (0/1): "
	   << PoissonFlag << "\n";
    }
    else if(comm=="RoundFlag") { // enable/disable round counts to integer
      GetIntToken(fs, &RoundFlag);
      cout << "Round pixel counts to integer (0/1): " << RoundFlag << "\n";
    }
    else if(comm=="HeaderFlag") { // flag to use or not header in output file
      GetIntToken(fs, &HeaderFlag);
      cout << "Use header in output file (0/1): " << HeaderFlag << "\n"; 
    }
    /*
    else if(comm=="RunningFasterFlag") { //columns(0) or rows(1) running faster
      GetIntToken(fs, &RunningFasterFlag);
      if (RunningFasterFlag==0) 
	cout << "Columns running faster than rows in output file\n"; 
      else
	cout << "Rows running faster than columns in output file\n"; 
    }
    */
    else if(comm=="PixelType") { // set the pixel content type
      GetIntToken(fs, &PixelType);
      cout << "Pixel content type: " << PixelType << "\n"; 
      cout << "0: fluence, 1: energy fluence, 2: fluence(E), " 
	"3: energy fluence(E)\n";
      if (PixelType != 2) 
	throw xrmc_exception("when using a detectorconvolute device, the PixelType must be 2!"); 
    }
    else if(comm=="ForceDetectFlag") {     // photon not forced/forced 
      GetIntToken(fs, &ForceDetectFlag); //   to be detected
      if (ForceDetectFlag==0) 
	cout << "Photon not forced to be detected\n"; 
      else
	cout << "Photon forced to be detected\n"; 
    }
    else if(comm=="Emin") { // set the minimum bin energy
      GetDoubleToken(fs, &Emin);
      cout << "\tEmin: " << Emin << "\n"; 
    }	
    else if(comm=="Emax") { // set the maximum bin energy
      GetDoubleToken(fs, &Emax);
      cout << "\tEmax: " << Emax << "\n"; 
    }	
    else if(comm=="NBins") { // set the number of energy bins
      GetIntToken(fs, &NBins);
      cout << "\tNBins: " << NBins << "\n"; 
    }	
    else if(comm=="SaturateEmin") { // flag to saturate energies < Emin
      GetIntToken(fs, &SaturateEmin);
      cout << "\tSaturate energies lower than Emin (0/1):"
	   << SaturateEmin << "\n"; 
    }
    else if(comm=="SaturateEmax") { // flag to saturate energies > Emax
      GetIntToken(fs, &SaturateEmax);
      cout << "\tSaturate energies greater than Emax (0/1):"
	   << SaturateEmax << "\n";
    }
    else if(comm=="Rotate") { // detector rotation
      cout << "Detector rotation :\n"; 
      cout << "\tPoint on rotation axis x0:\t";
      for (int i=0; i<3; i++) {
	GetDoubleToken(fs, &x0.Elem[i]); // translation components
      }
      cout << x0 << endl;
      cout << "\tRotation axis direction u:\t"; 
      for (int i=0; i<3; i++) {
	GetDoubleToken(fs, &u.Elem[i]); // get rotation axis direction
      }
      u.Normalize(); // normalize axis direction
      cout << u << endl;
      cout << "\tRotation angle theta (degrees): ";
      GetDoubleToken(fs, &theta); // get angle in degrees
      cout << theta << endl;
      R = matr3::RotMatr(u, -theta); // build rotation matrix
      X -= x0; // translate detector position: rotation axis-> origin
      X = R*X; // rotate detector position
      uk = R*uk; // rotate detector uk direction
      ui = R*ui; // rotate detector ui direction
      X += x0; // translate back detector position
    }
    else if(comm=="End") {
      break;
    }
    else if(comm=="") {
      cout << "Empty string\n";
    }
    else {
      throw xrmc_exception("syntax error in detectorarray input file"); 
    }

  }
  OrthoNormal(ui, uj, uk);  // evaluates uj to form a orthonormal basis
  return 0;
}

int detectorconvolute::SetDefault() {
	detectorarray::SetDefault();
	InputDeviceName[1] = "Composition";
	CrystalPhaseName = "Crystal";
	WindowPhaseName = "Window";
	xd->detector_type = XMI_DETECTOR_SI_SDD;
	xd->pulse_width = 0.0;
	xd->fano = 0.12;
	xd->noise = 0.1;
	xd->n_crystal_layers = 1;
	xd->crystal_layers = (struct xmi_layer *) malloc(sizeof(struct xmi_layer)); 
	xd->crystal_layers[0].n_elements = 1;
	xd->crystal_layers[0].Z = (int *) malloc(sizeof(int)); 
	xd->crystal_layers[0].weight = (double *) malloc(sizeof(double)); 
	xd->crystal_layers[0].Z[0] = 14;
	xd->crystal_layers[0].weight[0] = 1.0;
	xd->crystal_layers[0].thickness = 0.5;
	xd->crystal_layers[0].density = 2.33;

	det_absorber = (struct xmi_layer *) malloc(sizeof(struct xmi_layer));
	det_absorber->n_elements = 1;
	det_absorber->Z = (int *) malloc(sizeof(int));
	det_absorber->weight = (double *) malloc(sizeof(double));
	*(det_absorber->Z) = 4;
	*(det_absorber->weight) = 1.0;
	det_absorber->thickness = 2E-3;
	det_absorber->density = 1.85;

	CrystalThickness = 0.5;
	WindowThickness = 2E-3;

	return 0;
}
