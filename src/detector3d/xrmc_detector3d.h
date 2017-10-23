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
//     xrmc_detector3d.h           //
//        19/10/2017             //
//   author : Bruno Golosio      //
///////////////////////////////////
// detectorarray3d class definition
//

#ifndef DETECTOR3DH
#define DETECTOR3DH
#include "xrmc_device.h"
#include "xrmc_source.h"
#include "xrmc_arrayNd.h"
#include "randmt.h"
#include "xrmc_photon.h"
#include <vector>

//  detectorarray class definition, member variables and functions
class detectorarray3d : public bodydevice
{
 public:
  int NX, NY, NZ, N; // 3d array size and number of voxels (NX x NY x NZ)
  double *Image; // Acquired image array
  double ExpTime; // exposure time
  int PhotonNum; // Multiplicity of simulated events per detector voxel
  int NBins; // Num. of energy bins
  int ModeNum; // Num. of modes (scattering orders)
  std::vector<unsigned long> Seeds;

  virtual ~detectorarray3d(); // destructor
  detectorarray3d(std::string dev_name); // constructor

  virtual int CastInputDevices(); // cast input device method
  // method for linking input device
  //int LinkInputDevice(string command, xrmc_device *dev_pt);
  virtual int Load(std::istream &fs); // load detector parameters, position, orientation
 // save the acquired image in a file
  virtual int SaveData(std::string data_name, std::string file_name);
  virtual int SaveData(double *image, string file_name);

 // save the acquired image in a ascii file
  virtual int SaveAsciiData(double *image, std::string file_name);
  virtual int SetDefault(); // set the default values for detector parameters
  virtual long long EventMulti(); // event multiplicity
  virtual int Run() {return Acquisition();} // run the acquisition
  virtual int Acquisition(); // run the acquisition
  virtual int Clear(); // clear the detector voxel bin contents

 protected:
  basesource *Source; // input device (typically the sample device)
  //string SourceName; // name of the input device
  double VoxelSizeX, VoxelSizeY, VoxelSizeZ; // voxel size (cm)
  double VoxelVol; // voxel volume (cm3) 
  double fGLim; // cut on geometric factor
  int Shape; // Detector elements shape (0 parallelepiped, 1 sphere,
             // 2,3,4 cylinder (x,y,z axis)): 
  vect3 *VoxelX; // voxel coordinates array
  int RandomVoxelFlag; // flag to enable/disable random point on voxel
  int PoissonFlag; // flag to enable/disable Poisson statistic on counts
  int RoundFlag; // flag to enable/disable round voxel count to integer
  int HeaderFlag; // enable/disable writing header in output file
  int AsciiFlag; //  // binary(0) or ascii(1) output file format
  //int RunningFasterFlag; //columns(0) or rows(1) running faster
  int VoxelType; // voxel content type:
                 // 0: fluence,      1: energy fluence,
                 // 2: fluence(E),   3: energy fluence(E)
  int ForceDetectFlag; // Last photon not forced / forced to be detected
  int SaturateEmin; // flag to saturate energies lower than Emin
  int  SaturateEmax; // flag to saturate energies greater than Emin
  double Emin, Emax; // minimum and maximum bin energy
  std::vector<int> PrevCompZ;

  virtual int RunInit(); // detectorarray3d initialization before run
  int Init(); // detectorarray3d initialization before run
  //vect3 RandomPointInVoxel(int i); // Generates a random point
  vect3 RandomPointInVoxel(int i, randmt_t *rng);
                                   // on the pixel surface
  double fG(vect3 DRp); // evaluates the geometric factor fG, related to
  // the probability that the last photon trajectory crosses the voxel  
  //int Poisson(); // generate uncertainty on voxel count using Poisson statistic
  int Poisson(randmt_t *rng); // generate uncertainty on voxel count using Poisson statistic

// run the forced/unforced acquisition
  int ForcedAcquisition(basesource **SourceClones, photon *PhotonArray,
			randmt_t **rngs);
  int UnforcedAcquisition(basesource **SourceClones, photon *PhotonArray);

};

#endif
