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
//  loadintensityscreen.cpp      //
//        16/05/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Load intensityscreen parameters, position, orientation and size
//

#include <iostream>
#include <string>
#include "xrmc.h"
#include "xrmc_gettoken.h"
#include "xrmc_exception.h"
#include "xrmc_math.h"
#include "intensityscreen.h"

using namespace std;
using namespace gettoken;

//////////////////////////////////////////////////////////////////////
int intensityscreen::Load(istream &fs)
  // Loads intensityscreen parameters, position, orientation, size
{
  int i;
  vect3 x0, u;
  matr3 R;
  double theta;
  string comm="", image_file;

  cout << "intensityscreen parameters, position, orientation and size file\n";

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
    else if(comm=="ui") { // direction of ui, parallel to the screen rows
      cout << "Screen orientation :\n"; 
      cout << "\tui vector:\t"; 
      for (i=0; i<3; i++) {
	GetDoubleToken(fs, &ui.Elem[i]);
	//cout << ui[i] << "\t";
      }
      cout << ui << endl;
      //cout << "\n";
    }	
    else if(comm=="InterpolFlag") { //  do not use / use interpolation (0/1)
      GetIntToken(fs, &InterpolFlag);
      cout << "Do not use / use interpolation inside pixels/bins (0/1): "
	   << InterpolFlag << "\n"; 
    }
    else if (comm=="ImageFile") { // read screen image from a file
      GetToken(fs, image_file); // image file name
      LoadData(image_file);
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
      throw xrmc_exception("syntax error in intensityscreen input file"); 
    }
  }

  return 0;
}

int intensityscreen::LoadData(string file_name)
  // Load screen array contents
{
  FILE *fp;

  cout << "Loading intensityscreen data\n";

  if (N<=0)
    throw xrmc_exception("Screen must be dimensioned before loading.\n");
  if (Image!=NULL) delete[] Image;
  Image = new double[N]; // allocate the image array

  cout << "Input File: " << file_name << "\n";
  if ((fp = fopen(file_name.c_str(),"rb")) == NULL)
    throw xrmc_exception("Input file cannot be open.\n");

  if (fread(Image, sizeof(double), N, fp) != (unsigned int)N)
    throw xrmc_exception("Error reading input file.\n");

  fclose(fp);

  double sum=0;
  for(int i=0; i<N; i++) {
    double val = Image[i];
    if (val<0)
      throw xrmc_exception("Image bin contents must be nonnegative.\n");
    sum += val;
  }
  for(int i=0; i<N; i++) {
    Image[i] /= sum;
  }
  
  return 0;
}


//////////////////////////////////////////////////////////////////////
int intensityscreen::SetDefault()
// set default values for intensityscreen parameters
{
  NX = NY = 0; // number of rows and columns
  Image = NULL;
  CumulImage = NULL;
  PixelSizeX = PixelSizeY = 1; // Pixel size (Sx x Sy) = 1x1 cm2
  X.Set(0, 50, 0); // Screen center position 
  // Screen orientation :
  uk.Set(0,-1,0); // uk has direction opposite to y axis
  ui.Set(1,0,0);// ui has the same direction as the x axis
  OrthoNormal(ui, uj, uk); // evaluates uj to form a orthonormal basis

  InterpolFlag = 1; // use interpolation inside pixels/bins (0/1)
  dOmegaLim=2*PI; //Cut on solid angle from interaction point to a single pixel

  return 0;
}
