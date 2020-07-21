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
//     loaddetector.cpp          //
//        31/01/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Load detector parameters, position, orientation and size
//

#include <iostream>
#include <string>
#include "xrmc.h"
#include "xrmc_gettoken.h"
#include "xrmc_exception.h"
#include "xrmc_math.h"
#include "xrmc_loaddetector.h"
#include "image_convolution.h"

using namespace std;
using namespace gettoken;

//////////////////////////////////////////////////////////////////////
int detectorarray::Load(istream &fs)
  // Loads detector array position/orientation/size
{
  int i;
  vect3 x0, u;
  matr3 R;
  double theta;
  string comm="";

  cout << "Detector Array parameters, position, orientation and size file\n";

 // get a command/variable name from input file
  while (GetToken(fs, comm)) {
    // parse the command and decide what to do
    //
    // check if it's a command for setting an input device name
    if (ParseInputDeviceCommand(fs, comm)) continue;
    else if (ParsePSFCommand(fs, comm)) continue;
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
    else if(comm=="Z12") { // Distance between object plane and detector
      GetDoubleToken(fs, &Z12);
      cout <<  "Distance between object plane and detector: " <<  Z12 << endl;
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
    else if(comm=="AsciiFlag") { // binary(0) or ascii(1) output file format
      GetIntToken(fs, &AsciiFlag);
      cout << "Binary(0) or ascii(1) output file format: " << AsciiFlag
	   << "\n"; 
    }
    else if(comm=="MhdFlag") { // whether generate(1) or not(0) *.mhd file for raw image
      GetIntToken(fs, &MhdFlag);
      cout << "Whether generate(1) or not(0) *.mhd file for raw image: " << MhdFlag
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
    else if(comm=="PixelType") { // set the pixel content type
      GetIntToken(fs, &PixelType);
      cout << "Pixel content type: " << PixelType << "\n"; 
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
      cout << "Empy string\n";
    }
    else {
      throw xrmc_exception("syntax error in detectorarray input file"); 
    }
  }

  return 0;
}

int detectorarray::SaveData(string data_name, string file_name)
  // Save detector array contents
{

  cout << "Saving data: " << data_name << "\n";
  if (data_name!=SaveDataName[0] && data_name!=SaveDataName[1])
    throw xrmc_exception
      (string("Error: detectorarray device can only save data of type ")
       + SaveDataName[0] + " or " + SaveDataName[1] + "\n");  
  else if (data_name==SaveDataName[0]) {
    SaveData(Image, file_name) ; 
  }
  else {
    SaveData(ConvolutedImage, file_name) ;
  }
  if (MhdFlag)
  {
    // generate meta header file for a single raw image
    SaveImageMetaHeaderData(file_name) ;
  }

  return 0;
}

int detectorarray::SaveData(double ***image, string file_name)
  // Save detector array contents
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
    fwrite(&PixelSizeX, sizeof(double), 1, fp);
    fwrite(&PixelSizeY, sizeof(double), 1, fp);
    fwrite(&ExpTime, sizeof(double), 1, fp);
    fwrite(&PixelType, sizeof(int), 1, fp);
    fwrite(&NBins, sizeof(int), 1, fp);
    fwrite(&Emin, sizeof(double), 1, fp);
    fwrite(&Emax, sizeof(double), 1, fp);
  }

  fwrite(image[0][0], sizeof(double), ModeNum*NX*NY*NBins, fp);
    
  fclose(fp);

  return 0;
}

int detectorarray::SaveImageMetaHeaderData(string file_name)
{
  if (AsciiFlag != 0 || Shape != 0 || NBins != 1 || ModeNum != 1)
  {
    return 0;
  }

  string mhd_name;
  size_t pos = file_name.find_last_of('.');
  if (pos < file_name.size() && file_name.substr(pos).size() == 4)
  {
    mhd_name = file_name.substr( 0, pos) + ".mhd";
  }
  else if (pos < file_name.size())
  {
    mhd_name = file_name + ".mhd";
  }

  if (mhd_name.empty())
  {
    return 0;
  }

  ofstream ofs(mhd_name.c_str());
  if (!ofs)
    throw xrmc_exception("Output image meta header file cannot be open.\n");

  ofs << "NDims = " << 2 << endl;
  ofs << "DimSize = " << NX << ' ' << NY << endl;
  ofs << "ElementSpacing = " << PixelSizeX << ' ' << PixelSizeY << endl;
  ofs << "Position = 0 0" << endl;
  ofs << "BinaryData = True" << endl;
  ofs << "ElementByteOrderMSB = False" << endl;
  ofs << "ElementType = MET_DOUBLE" << endl;
  if (HeaderFlag == 1) {
    ofs << "HeaderSize = 60" << endl;
  }
  ofs << "ElementDataFile = " << file_name << endl;

  ofs.close();

  return 0;
}

int detectorarray::SaveAsciiData(double ***image, string file_name)
  // Save detector array contents
{
  ofstream ofs(file_name.c_str());
  if (!ofs)
    throw xrmc_exception("Output ascii file cannot be open.\n");

  if (HeaderFlag==1) {
    ofs << ModeNum << endl;
    ofs << NX << endl;
    ofs << NY << endl;
    ofs << PixelSizeX << endl;
    ofs << PixelSizeY << endl;
    ofs << ExpTime << endl;
    ofs << PixelType << endl;
    ofs << NBins << endl;
    ofs << Emin << endl;
    ofs << Emax << endl;
  }
  for (int i=0; i<ModeNum*NX*NY*NBins; i++) {
    ofs << *(image[0][0]+i) << endl;
  }
  ofs.close();

  return 0;
}

//////////////////////////////////////////////////////////////////////
int detectorarray::SetDefault()
// set default values for detectorarray parameters
{
  InputDeviceName[0] = "Sample"; // input device name
  NX = NY = 1; // number of rows and columns

  PixelSizeX = PixelSizeY = 1; // Pixel size (Sx x Sy) = 1x1 cm2
  Shape = 0; // rectangular pixel shape
  dOmegaLim=2*PI; //Cut on solid angle from interaction point to a single pixel
  X.Set(0, 50, 0); // Detector center position 
  // Detector array orientation :
  uk.Set(0,-1,0); // uk has direction opposite to y axis
  ui.Set(1,0,0);// ui has the same direction as the x axis
  OrthoNormal(ui, uj, uk); // evaluates uj to form a orthonormal basis

  ExpTime = 1; // 1 sec. Exposure Time
  PhotonNum = 1; // Num. of simulated photons per detector pixel
  Z12 = 0; // object-detector distance used for source convolution
  RandomPixelFlag = 1; // Enable random point on pixels 
  PoissonFlag = 0; // Poisson statistic on pixel counts disabled
  RoundFlag = 0;   // pixel counts round to integer disabled
  HeaderFlag = 0; // header in output file disabled
  AsciiFlag = 0; // Default output format is binary 
  MhdFlag = 0; // do not generate *.mhd file for a raw image
  //RunningFasterFlag = 0; // columns running faster than rows in output file
  ForceDetectFlag = 1; //Photon forced to be detected 
  ConvolveFlag = 0; // do not generate convolved image
  PixelType = 0; // Pixel content type = 0: fluence
  NBins = 1; // 1 energy bin
  Emin = 0; // minimum bin energy is 0 keV
  Emax = 100; // maximum bin energy is 100 keV
  SaturateEmin = 0; // do not Saturate energies lower than Emin
  SaturateEmax = 0; // do not Saturate energies greater than Emax
  NBx = NBy = 20; // size of borders (in pixels) used for FFT convolution
  EfficiencyFlag = 0; //EfficiencyFLag off by default

  return 0;
}


int detectorarray::LoadPSFBin(int n, vector<gauss_vect> &PSF_bin,
			      istream &psf_fs)
{
  cout << "Gaussian function parameters :\n";
  gauss_par gp;
  gauss_vect gv;

  if (NBins<=0)
    throw xrmc_exception("Number of energy bins (NBins) undefined.\n");
    
  PSF_bin.clear();    
  for (int ibin=0; ibin<NBins; ibin++) {
    gv.clear();
    for (int ig=0; ig<n; ig++) {
      GetDoubleToken(psf_fs, &gp.height);
      GetDoubleToken(psf_fs, &gp.x0);
      GetDoubleToken(psf_fs, &gp.sigma);
      cout << "\t" << gp.height <<
	"\t" << gp.x0 << "\t" << gp.sigma << "\t"; 
      gv.push_back(gp);
    }
    cout << endl;
    PSF_bin.push_back(gv);    
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////
// parse file for PSF, source size and efficiency commands
//////////////////////////////////////////////////////////////////////
bool detectorarray::ParsePSFCommand(istream &fs, string comm)
{
  if(comm=="EfficiencyFlag") { // flag for using efficiency before convolution
    GetIntToken(fs, &EfficiencyFlag);
    cout << "Use efficiency before image convolution (0/1): "
	 << ConvolveFlag << "\n"; 
  }
  //  load from file energy-dependent efficiency
  else if(comm=="EfficiencyFile") {
    cout << "Energy dependent detector efficiency\n";
    Efficiency.clear();
    string eff_file;
    GetToken(fs, eff_file); // efficiency file name
    cout << "\tefficiency file name: " << eff_file << "\n";
    ifstream eff_fs(eff_file.c_str());
    if (!eff_fs)
      throw xrmc_exception("Efficiency file can not be opened.");
    LoadEfficiency(eff_fs);
    eff_fs.close();
  }
  else if(comm=="ConvolveFlag") {  // generates convoluted image (0/1)
      GetIntToken(fs, &ConvolveFlag);
      cout << "Generates convoluted image (0/1): " << ConvolveFlag << "\n"; 
    }
    else if(comm=="GaussPSFx") { // sum-of-gaussians model of detector PSF
    cout << "Gaussian model of detector PSF (x component)\n";
    int n;
    GetIntToken(fs, &n);
    GaussPSFx.clear();
    if (n==0) cout << "\t\tDisabled\n";
    else {
      GaussPSFxBin.clear();
      cout << "\tSuperposition of " << n << " Gaussian functions\n";
      gauss_par gp;
      for (int i=0; i< n; i++) {
	GetDoubleToken(fs, &gp.height);
	GetDoubleToken(fs, &gp.x0);
	GetDoubleToken(fs, &gp.sigma);
	cout << "\tGaussian function n. " << i << ": height " << gp.height <<
	  ", x0 " << gp.x0 << ", sigma " << gp.sigma << endl; 
	GaussPSFx.push_back(gp);
      }
    }
  }
  else if(comm=="GaussPSFy") { // sum-of-gaussians model of detector PSF
    cout << "Gaussian model of detector PSF (y component)\n";
    int n;
    GetIntToken(fs, &n);
    GaussPSFy.clear();
    if (n==0) cout << "\t\tDisabled\n";
    else {
      GaussPSFyBin.clear();
      cout << "\tSuperposition of " << n << " Gaussian functions\n";
      gauss_par gp;
      for (int i=0; i< n; i++) {
	GetDoubleToken(fs, &gp.height);
	GetDoubleToken(fs, &gp.x0);
	GetDoubleToken(fs, &gp.sigma);
	cout << "\tGaussian function n. " << i << ": height " << gp.height <<
	  ", y0 " << gp.x0 << ", sigma " << gp.sigma << endl; 
	GaussPSFy.push_back(gp);
      }
    }
  }
  // energy bin dependent sum-of-gaussians model of detector PSF
  // load from file
  else if(comm=="GaussPSFxBinFile") {
    cout << "Energy dependent Gaussian model of detector PSF (x component)\n";
    int n;
    GetIntToken(fs, &n);
    GaussPSFxBin.clear();
    if (n==0) cout << "\tDisabled\n";
    else {
      GaussPSFx.clear();
      cout << "\tSuperposition of " << n << " Gaussian functions\n";
      string psf_file;
      GetToken(fs, psf_file); // psf file name
      cout << "\tPSF file name: " << psf_file << "\n";
      ifstream psf_fs(psf_file.c_str());
      if (!psf_fs)
	throw xrmc_exception("PSF file can not be opened.");
      LoadPSFBin(n, GaussPSFxBin, psf_fs);
      psf_fs.close();
    }
  }
  else if(comm=="GaussPSFyBinFile") {
    cout << "Energy dependent Gaussian model of detector PSF (y component)\n";
    int n;
    GetIntToken(fs, &n);
    GaussPSFyBin.clear();
    if (n==0) cout << "\tDisabled\n";
    else {
      GaussPSFy.clear();
      cout << "\tSuperposition of " << n << " Gaussian functions\n";
      string psf_file;
      GetToken(fs, psf_file); // psf file name
      cout << "\tPSF file name: " << psf_file << "\n";
      ifstream psf_fs(psf_file.c_str());
      if (!psf_fs)
	throw xrmc_exception("PSF file can not be opened.");
      LoadPSFBin(n, GaussPSFyBin, psf_fs);
      psf_fs.close();
    }
  }
  
  else if(comm=="GaussSourceX") { // sum-of-gaussians model of source x size
    cout << "Gaussian model of source size (x component)\n";
    int n;
    GetIntToken(fs, &n);
    GaussSourceX.clear();
    if (n==0) cout << "\t\tDisabled\n";
    else {
      GaussSourceXBin.clear();
      cout << "\tSuperposition of " << n << " Gaussian functions\n";
      gauss_par gp;
      for (int i=0; i< n; i++) {
	GetDoubleToken(fs, &gp.height);
	GetDoubleToken(fs, &gp.x0);
	GetDoubleToken(fs, &gp.sigma);
	cout << "\tGaussian function n. " << i << ": height " << gp.height <<
	  ", x0 " << gp.x0 << ", sigma " << gp.sigma << endl; 
	GaussSourceX.push_back(gp);
      }
    }
  }
  else if(comm=="GaussSourceY") { // sum-of-gaussians model of source y size
    cout << "Gaussian model of source size (y component)\n";
    int n;
    GetIntToken(fs, &n);
    GaussSourceY.clear();
    if (n==0) cout << "\t\tDisabled\n";
    else {
      GaussSourceYBin.clear();
      cout << "\tSuperposition of " << n << " Gaussian functions\n";
      gauss_par gp;
      for (int i=0; i< n; i++) {
	GetDoubleToken(fs, &gp.height);
	GetDoubleToken(fs, &gp.x0);
	GetDoubleToken(fs, &gp.sigma);
	cout << "\tGaussian function n. " << i << ": height " << gp.height <<
	  ", y0 " << gp.x0 << ", sigma " << gp.sigma << endl; 
	GaussSourceY.push_back(gp);
      }
    }
  }
  // energy bin dependent sum-of-gaussians model of source size
  // load from file
  else if(comm=="GaussSourceXBinFile") {
    cout << "Energy dependent Gaussian model of source size (x component)\n";
    int n;
    GetIntToken(fs, &n);
    GaussSourceXBin.clear();
    if (n==0) cout << "\tDisabled\n";
    else {
      GaussSourceX.clear();
      cout << "\tSuperposition of " << n << " Gaussian functions\n";
      string psf_file;
      GetToken(fs, psf_file); // psf file name
      cout << "\tPSF file name: " << psf_file << "\n";
      ifstream psf_fs(psf_file.c_str());
      if (!psf_fs)
	throw xrmc_exception("PSF file can not be opened.");
      LoadPSFBin(n, GaussSourceXBin, psf_fs);
      psf_fs.close();
    }
  }
  else if(comm=="GaussSourceYBinFile") {
    cout << "Energy dependent Gaussian model of source size (y component)\n";
    int n;
    GetIntToken(fs, &n);
    GaussSourceYBin.clear();
    if (n==0) cout << "\tDisabled\n";
    else {
      GaussSourceY.clear();
      cout << "\tSuperposition of " << n << " Gaussian functions\n";
      string psf_file;
      GetToken(fs, psf_file); // psf file name
      cout << "\tPSF file name: " << psf_file << "\n";
      ifstream psf_fs(psf_file.c_str());
      if (!psf_fs)
	throw xrmc_exception("PSF file can not be opened.");
      LoadPSFBin(n, GaussSourceYBin, psf_fs);
      psf_fs.close();
    }
  }
  else return false;

  return true;
}

//////////////////////////////////////////////////////////////////////
// Load energy dependent efficiency
//////////////////////////////////////////////////////////////////////
int detectorarray::LoadEfficiency(istream &eff_fs)
{
  if (NBins<=0)
    throw xrmc_exception("Number of energy bins (NBins) undefined.\n");
    
  Efficiency.clear();    
  for (int ibin=0; ibin<NBins; ibin++) {
    double eff;
    GetDoubleToken(eff_fs, &eff);
    Efficiency.push_back(eff);
  }

  return 0;
}
