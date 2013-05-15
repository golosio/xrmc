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
//     loadbeamscreen.cpp        //
//        02/05/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Load beamscreen parameters, position, orientation and size
//

#include <iostream>
#include <string>
#include "xrmc.h"
#include "xrmc_gettoken.h"
#include "xrmc_exception.h"
#include "xrmc_math.h"
#include "beamscreen.h"

using namespace std;
using namespace gettoken;

//////////////////////////////////////////////////////////////////////
int beamscreen::Load(istream &fs)
  // Loads beamscreen parameters, position, orientation, size
{
  int i;
  vect3 x0, u;
  matr3 R;
  double theta;
  string comm="", image_file;

  cout << "Beam Screen parameters, position, orientation and size file\n";

 // get a command/variable name from input file
  while (GetToken(fs, comm)) {
    // parse the command and decide what to do
    //
    if(comm=="NPixels") { // set the number of pixels (rows and columns)
      GetIntToken(fs, &NX);
      GetIntToken(fs, &NY);
      cout << "Pixel number (NX x NY): " << NX << " x " << NY << "\n";
      N = NX*NY;
      if (Image!=NULL) {
	delete[] Image;
	Image = NULL;
      }
    }
    else if(comm=="PixelSize") { // set the pixel size in x and y directions
      GetDoubleToken(fs, &PixelSizeX);
      GetDoubleToken(fs, &PixelSizeY);
      cout << "Pixel size (Sx x Sy) (cm): " << PixelSizeX << " x "
	   <<  PixelSizeY << "\n";
    }
    else if(comm=="X") { // set the screen center coordinates
      cout << "Screen center position :\t"; 
      for (i=0; i<3; i++) {
	GetDoubleToken(fs, &X.Elem[i]);
      }
      cout << X << endl;
    }
    else if(comm=="uk") { // direction of the normal to the screen surface uk
      cout << "Screen orientation :\n"; 
      cout << "\tuk vector:\t"; 
      for (i=0; i<3; i++) {
	GetDoubleToken(fs, &uk.Elem[i]);
      }
      cout << uk << endl;
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
    else if(comm=="PolarizedFlag") { //  unpolarized/polarized beam flag (0/1)
      GetIntToken(fs, &PolarizedFlag);
      cout << "Unpolarized/polarized beam (0/1): " << PolarizedFlag << "\n"; 
    }
    else if(comm=="EnergyBinFlag") { //  energy binning flag (0/1)
      GetIntToken(fs, &EnergyBinFlag);
      cout << "Energy Binning (0/1): " << EnergyBinFlag << "\n"; 
    }
    else if(comm=="InterpolFlag") { //  do not use / use interpolation (0/1)
      GetIntToken(fs, &InterpolFlag);
      cout << "Do not use / use interpolation inside pixels/bins (0/1): "
	   << InterpolFlag << "\n"; 
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
      if (Image!=NULL) {
	delete[] Image;
	Image = NULL;
      }
    }
    else if(comm=="TotalIntensity") { // set the beam total intensity
      GetDoubleToken(fs, &TotalIntensity);
      cout << "\tTotalIntensity: " << TotalIntensity << "\n"; 
    }	

    else if (comm=="ImageFile") { // read screen image from a file
      GetToken(fs, image_file); // image file name
      LoadData(image_file);
    }
    else if(comm=="LoopFlag") { // flag for loop mode
      GetIntToken(fs, &LoopFlag);
      if (LoopFlag==0)
	cout << "Extract random energies on the whole spectrum\n";
      else if (LoopFlag==1)
	cout << "Loop on all energy bins and polarization\n";
      else throw xrmc_exception("Wrong loop flag.\n");
    }
    else if(comm=="PhotonNum") {
      // Multiplicity of events for each energy bin 
      GetIntToken(fs, &PhotonNum);
      cout <<"Multiplicity of events: " 
	   << PhotonNum << "\n";
    }
	
    else if(comm=="Rotate") { // screen rotation
      cout << "Screen rotation :\n"; 
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
      X -= x0; // translate screen position: rotation axis-> origin
      X = R*X; // rotate screen position
      uk = R*uk; // rotate screen uk direction
      ui = R*ui; // rotate screen ui direction
      X += x0; // translate back screen position
    }
    else if(comm=="End") {
      break;
    }
    else if(comm=="") {
      cout << "Empy string\n";
    }
    else {
      throw xrmc_exception("syntax error in beamscreen input file"); 
    }
  }

  return 0;
}

int beamscreen::LoadData(string file_name)
  // Load screen array contents
{
  FILE *fp;

  cout << "Loading beam screen data\n";

  if (N*NBins<=0)
    throw xrmc_exception("Screen must be dimensioned before loading.\n");
  if (Image!=NULL) delete[] Image;
  Image = new double[2*N*NBins]; // allocate the image array

  cout << "Input File: " << file_name << "\n";
  if ((fp = fopen(file_name.c_str(),"rb")) == NULL)
    throw xrmc_exception("Input file cannot be open.\n");

  // check if beam is polarized
  if (PolarizedFlag==1) {
    if (fread(Image, sizeof(double), 2*N*NBins, fp)
	!= (unsigned int)2*N*NBins)
      throw xrmc_exception("Error reading input file.\n");
  }
  else { // if not build two polarization components with equal intensity
    // x-polarized component
    if (fread(Image, sizeof(double), N*NBins, fp) != (unsigned int)N*NBins)
      throw xrmc_exception("Error reading input file.\n");
    for(int i=0; i<N*NBins; i++) {
      Image[N*NBins+i] = Image[i]; // y-polarized component
    }
  }
  fclose(fp);

  double sum=0;
  for(int i=0; i<2*N*NBins; i++) {
    double val = Image[i];
    if (val<0)
      throw xrmc_exception("Image bin contents must be nonnegative.\n");
    sum += val;
  }
  for(int i=0; i<2*N*NBins; i++) {
    Image[i] /= sum;
  }

  return 0;
}


//////////////////////////////////////////////////////////////////////
int beamscreen::SetDefault()
// set default values for beamscreen parameters
{
  NX = NY = 0; // number of rows and columns
  Image = NULL;
  CumulImage = NULL;
  CumulEnergy = NULL;
  PixelSizeX = PixelSizeY = 1; // Pixel size (Sx x Sy) = 1x1 cm2
  X.Set(0, 50, 0); // Screen center position 
  // Screen orientation :
  uk.Set(0,-1,0); // uk has direction opposite to y axis
  ui.Set(1,0,0);// ui has the same direction as the x axis
  OrthoNormal(ui, uj, uk); // evaluates uj to form a orthonormal basis

  PolarizedFlag = 0;  // unpolarized/polarized beam flag (0/1)
  EnergyBinFlag = 1;  // energy binning flag (0/1)
  InterpolFlag = 1; // use interpolation inside pixels/bins (0/1)
  NBins = 1; // 1 energy bin
  Emin = 0; // minimum bin energy is 0 keV
  Emax = 100; // maximum bin energy is 100 keV
  TotalIntensity = 1; // beam total intensity 
  dOmegaLim=2*PI; //Cut on solid angle from interaction point to a single pixel
  LoopFlag = 0; // flag for loop mode
  PhotonNum = 1; // Multiplicity of events

  return 0;
}
