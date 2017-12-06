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
//     loaddetector3d.cpp        //
//        19/10/2017             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Load 3d detector parameters, position, orientation and size
//

#include <iostream>
#include <string>
#include "xrmc.h"
#include "xrmc_gettoken.h"
#include "xrmc_exception.h"
#include "xrmc_math.h"
#include "xrmc_loaddetector3d.h"

using namespace std;
using namespace gettoken;

//////////////////////////////////////////////////////////////////////
int detectorarray3d::Load(istream &fs)
  // Loads 3d detector array position/orientation/size
{
  int i;
  vect3 x0, u;
  matr3 R;
  double theta;
  string comm="";

  cout << "3D Detector Array parameters, position, orientation and size file\n";

 // get a command/variable name from input file
  while (GetToken(fs, comm)) {
    // parse the command and decide what to do
    //
    // check if it's a command for setting an input device name
    if (ParseInputDeviceCommand(fs, comm)) continue;
    else if(comm=="NVoxels") { // set the number of voxels (NX,NY,NZ)
      GetIntToken(fs, &NX);
      GetIntToken(fs, &NY);
      GetIntToken(fs, &NZ);
      cout << "Voxel number (NX x NY X NZ): " << NX << " x " << NY << " x "
	   << NZ << "\n";
    }
    else if(comm=="VoxelSize") { // set the voxel size in x,y,z directions
      GetDoubleToken(fs, &VoxelSizeX);
      GetDoubleToken(fs, &VoxelSizeY);
      GetDoubleToken(fs, &VoxelSizeZ);

      cout << "Voxel size (Sx x Sy x Sz) (cm): " << VoxelSizeX << " x "
	   <<  VoxelSizeY << " x " <<  VoxelSizeZ << "\n";
    }
    else if(comm=="Shape") { // detector elements shape
                             // (parallelepiped, sphere or cylinder(x,y,z))
      GetIntToken(fs, &Shape);
      cout << "Detector elements shape (0 parallelepiped, "
	"1 sphere, 2,3,4 cylinder (x,y,z axis)): "
	   << Shape << "\n";
    }
    else if(comm=="X") { // set the detector center coordinates
      cout << "3D detector array center position :\t"; 
      for (i=0; i<3; i++) {
	GetDoubleToken(fs, &X.Elem[i]);
	//cout << X[i] << "\t";
      }
      cout << X << endl;
      //cout << "\n";
    }
    else if(comm=="uk") { // direction of the normal to the detector surface uk
      cout << "3d detector orientation :\n"; 
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
    else if(comm=="PhotonNum") { // Multiplicity of simulated events per voxel
      GetIntToken(fs, &PhotonNum);
      cout << "Multiplicity of simulated events per voxel: "
	   << PhotonNum << "\n";
    }
    else if(comm=="RandomVoxelFlag") { // enable/disable random point in voxels
      GetIntToken(fs, &RandomVoxelFlag);
      cout << "Enable random point in voxels (0/1): " << RandomVoxelFlag
	   << "\n";
    }
    else if(comm=="PoissonFlag") { // flag to enable/disable Poisson statistic
      GetIntToken(fs, &PoissonFlag);
      cout << "Enable Poisson statistic on voxel counts (0/1): "
	   << PoissonFlag << "\n";
    }
    else if(comm=="RoundFlag") { // enable/disable round counts to integer
      GetIntToken(fs, &RoundFlag);
      cout << "Round voxel counts to integer (0/1): " << RoundFlag << "\n";
    }
    else if(comm=="HeaderFlag") { // flag to use or not header in output file
      GetIntToken(fs, &HeaderFlag);
      cout << "Use header in output file (0/1): " << HeaderFlag << "\n"; 
    }
    else if(comm=="AsciiFlag") { // binary(0) or ascii(1) output file format
      GetIntToken(fs, &AsciiFlag);
      cout << "Binary(0) or ascii(1) output file format: " << AsciiFlag
	   << "\n"; 
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
    else if(comm=="VoxelType") { // set the voxel content type
      GetIntToken(fs, &VoxelType);
      cout << "Voxel content type: " << VoxelType << "\n"; 
      cout << "0: fluence, 1: energy fluence, 2: fluence(E), " 
	"3: energy fluence(E)\n";
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
    else if(comm=="Seeds") { // seeds for random number generation
      cout << "Seeds for random number generation\n";
      int NS;
      long l;
      GetIntToken(fs, &NS);
      cout << "\tNumber of seeds: " << NS << "\n";
      Seeds.clear();
      for (int i=0; i<NS; i++) {
	GetLongToken(fs, &l);
	cout << "\tThread n. " << i << " , Seed: " << l << "\n";
	Seeds.push_back(l);
      }	
    }
    else if(comm=="Rotate") { // detector rotation
      cout << "3D detector rotation :\n"; 
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
      cout << "Empy string\n";
    }
    else {
      throw xrmc_exception("syntax error in detectorarray3d input file"); 
    }
  }

  return 0;
}

int detectorarray3d::SaveData(string data_name, string file_name)
  // Save detector array contents
{

  cout << "Saving data: " << data_name << "\n";
  if (data_name!=SaveDataName[0])
    throw xrmc_exception
      (string("Error: detectorarray3d device can only save data of type ")
       + SaveDataName[0] + "\n");  
  else
    SaveData(Image, file_name) ; 

  return 0;
}

int detectorarray3d::SaveData(double *image, string file_name)
  // Save 3d detector array contents
{
  FILE *fp;

  if (image==NULL)
    throw xrmc_exception("Data array was not created.\n");

  cout << "Output File: " << file_name << "\n";

  if (AsciiFlag!=0) return SaveAsciiData(image, file_name);

  if ((fp = fopen(file_name.c_str(),"wb")) == NULL)
    throw xrmc_exception("Output file cannot be open.\n");

  if (HeaderFlag==1) {
    fwrite(&ModeNum, sizeof(int), 1, fp);
    fwrite(&NX, sizeof(int), 1, fp);
    fwrite(&NY, sizeof(int), 1, fp);
    fwrite(&NZ, sizeof(int), 1, fp);
    fwrite(&VoxelSizeX, sizeof(double), 1, fp);
    fwrite(&VoxelSizeY, sizeof(double), 1, fp);
    fwrite(&VoxelSizeZ, sizeof(double), 1, fp);
    fwrite(&ExpTime, sizeof(double), 1, fp);
    fwrite(&VoxelType, sizeof(int), 1, fp);
    fwrite(&NBins, sizeof(int), 1, fp);
    fwrite(&Emin, sizeof(double), 1, fp);
    fwrite(&Emax, sizeof(double), 1, fp);
  }

  fwrite(image, sizeof(double), ModeNum*NX*NY*NZ*NBins, fp);
    
  fclose(fp);

  return 0;
}

int detectorarray3d::SaveAsciiData(double *image, string file_name)
  // Save 3d detector array contents
{
  ofstream ofs(file_name.c_str());
  if (!ofs)
    throw xrmc_exception("Output ascii file cannot be open.\n");

  if (HeaderFlag==1) {
    ofs << ModeNum << endl;
    ofs << NX << endl;
    ofs << NY << endl;
    ofs << NZ << endl;
    ofs << VoxelSizeX << endl;
    ofs << VoxelSizeY << endl;
    ofs << VoxelSizeZ << endl;
    ofs << ExpTime << endl;
    ofs << VoxelType << endl;
    ofs << NBins << endl;
    ofs << Emin << endl;
    ofs << Emax << endl;
  }
  for (int i=0; i<ModeNum*NX*NY*NZ*NBins; i++) {
    ofs << image[i] << endl;
  }
  ofs.close();

  return 0;
}

//////////////////////////////////////////////////////////////////////
int detectorarray3d::SetDefault()
// set default values for detectorarray3d parameters
{
  InputDeviceName[0] = "Sample"; // input device name
  NX = NY = NZ = 1; // number of voxels

  VoxelSizeX = VoxelSizeY = VoxelSizeZ = 1; // Voxel size (Sx,Sy,Sz in cm)
  Shape = 0; // parallelepiped voxel shape

  X.Set(0, 50, 0); // Detector center position 
  // Detector array orientation :
  uk.Set(0,-1,0); // uk has direction opposite to y axis
  ui.Set(1,0,0);// ui has the same direction as the x axis
  OrthoNormal(ui, uj, uk); // evaluates uj to form a orthonormal basis

  ExpTime = 1; // 1 sec. Exposure Time
  PhotonNum = 1; // Num. of simulated photons per detector voxel
  RandomVoxelFlag = 1; // Enable random point in voxels 
  PoissonFlag = 0; // Poisson statistic on voxel counts disabled
  RoundFlag = 0;   // voxel counts round to integer disabled
  HeaderFlag = 0; // header in output file disabled
  AsciiFlag = 0; // Default output format is binary 
  //RunningFasterFlag = 0; // columns running faster than rows in output file
  ForceDetectFlag = 1; //Photon forced to be detected 
  VoxelType = 0; // Voxel content type = 0: fluence
  NBins = 1; // 1 energy bin
  Emin = 0; // minimum bin energy is 0 keV
  Emax = 100; // maximum bin energy is 100 keV
  SaturateEmin = 0; // do not Saturate energies lower than Emin
  SaturateEmax = 0; // do not Saturate energies greater than Emax

  return 0;
}


