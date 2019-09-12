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
//     detector3d.cpp            //
//        21/12/2017             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Methods of the class detectorarray3d
//
#include <cmath>
#include <cstring>
#include <iostream>
#include "xrmc_math.h"
#include "xrmc_detector3d.h"
#include "xrmc_photon.h"
#include "xrmc_exception.h"
#include "xrmc_algo.h"
#include "xrmc_arrayNd.h"
#include "xrmc_sample.h"

#ifdef TIME_PERF
#include <ctime>
double outph_time, outph_time0, edep_time, soutph_time, mu_time0, scatter_time,
				       phist_time, mu_time1, mux1_time,
				       psw_time, soutphx1_time, phist_time0;
#endif

//#undef _OPENMP
#ifdef _OPENMP
#include <omp.h>
#define THREAD_MAXNUM omp_get_max_threads()
#define THREAD_IDX omp_get_thread_num()
#else
#define THREAD_MAXNUM 1
#define THREAD_IDX 0
#endif

using namespace std;
using namespace xrmc_algo;
using namespace arrayNd;

// destructor
detectorarray3d::~detectorarray3d() {
  if (VoxelX!=NULL) delete[] VoxelX; // should becalled on base class
  if (Image!=NULL) delete[] Image;
}

// constructor
detectorarray3d::detectorarray3d(string dev_name) {
  Runnable = true;
  SaveDataName.push_back("Image");

  NInputDevices = 1;
  InputDeviceCommand.push_back("SourceName");
  InputDeviceDescription.push_back("Source input device name");

  VoxelX=NULL;
  Image=NULL;
  NX = NY = NZ = N = PhotonNum = NBins = ModeNum = 0;
  SetDevice(dev_name, "detectorarray3d");
}

//////////////////////////////////////////////////////////////////////
// detectorarray3d initialization before run method
//////////////////////////////////////////////////////////////////////
int detectorarray3d::RunInit()
{
  Init();

  ModeNum = Source->ModeNum(); // number of modes (scattering orders)
  if (Image!=NULL) delete[] Image;
  Image = new double[ModeNum*NBins*NX*NY*NZ]; // allocate the image array

  return 0;
}

//////////////////////////////////////////////////////////////////////
// run the acquisition
//////////////////////////////////////////////////////////////////////
int detectorarray3d::Acquisition()
{
  vect3 Rp, DRp;

  Clear(); // clear the detector bin contents

  //init_randmt_auto_r(Rng);
  //Photon.Rng = Rng;

  //initialize rngs
  randmt_t **rngs;
  rngs = new randmt_t*[THREAD_MAXNUM];
  basesource **SourceClones = new basesource*[THREAD_MAXNUM];
  photon *PhotonArray = new photon[THREAD_MAXNUM];

  if (PrevCompZ != dynamic_cast<sample*>(Source)->compZ) {
  	//new elements found -> recalculate
	//this could be done in a better way, for example by recalculating only those elements
	//that were not present before
  	dynamic_cast<sample*>(Source)->DopplerClear();
  	dynamic_cast<sample*>(Source)->DopplerInit();
	PrevCompZ = dynamic_cast<sample*>(Source)->compZ;
  }

  cout << "Maximum number of threads: "  << THREAD_MAXNUM << endl;
  if (THREAD_MAXNUM==1) {
    SourceClones[0] = Source;
    rngs[0] = new_randmt();
    if (Seeds.size()>0 && Seeds[0]!=0) init_randmt_r(rngs[0], Seeds[0]);
    else init_randmt_auto_r(rngs[0]);
    SourceClones[0]->SetRng(rngs[0]);
    PhotonArray[0].Rng = rngs[0];
  }
  else {
    for (int i=0; i<THREAD_MAXNUM; i++) {
      SourceClones[i] = Source->Clone(InputDeviceName[0]);
      rngs[i] = new_randmt();
      if (rngs[i] == NULL)
	throw xrmc_exception(string("Could not allocate new RNG\n"));
      if (i<(int)Seeds.size() && Seeds[i]!=0) {
	init_randmt_r(rngs[i], Seeds[i]);
	cout << "\tThread " << i << " , seed " << Seeds[i] << endl;
      } 
      else init_randmt_auto_r(rngs[i]);
      SourceClones[i]->SetRng(rngs[i]);
      PhotonArray[i].Rng = rngs[i];
    }
  }
  // Last photon not forced / forced to be detected
  if(ForceDetectFlag==0) {
    UnforcedAcquisition(SourceClones, PhotonArray);
    //throw xrmc_exception(string("Unforced acquisition not yet implemented for 3d detector\n"));

  }
  else ForcedAcquisition(SourceClones, PhotonArray, rngs);
    
  // generate uncertainty on voxel count using Poisson statistic
  if (PoissonFlag == 1) {
    Poisson(rngs[0]);
  }

  //free rngs
  for (int i=0; i<THREAD_MAXNUM; i++) {
    free_randmt(rngs[i]);
    if (THREAD_MAXNUM > 1)
      delete SourceClones[i];
  }
  delete [] rngs;

  delete [] SourceClones;
  delete [] PhotonArray;

  return 0;
}

//////////////////////////////////////////////////////////////////////
// run the forced acquisition
//////////////////////////////////////////////////////////////////////
int detectorarray3d::ForcedAcquisition(basesource **SourceClones,
				     photon *PhotonArray, randmt_t **rngs)
{
  int ibin, mode_idx;
  vect3 Rp;
  double signal;
  double mu_x1, Edep;
  const int ProgrUpdate=100000;
  long long event_idx=0, event_tot=EventMulti();
  if (event_tot>ProgrUpdate) {
    printf("\nProgress 0 %%\r");
    fflush(stdout);
  }

#ifdef TIME_PERF
  outph_time=outph_time0=edep_time=soutph_time=mu_time0=scatter_time=phist_time
    =mu_time1=mux1_time=psw_time=soutphx1_time=phist_time0=0;
#endif

  int ivox, ix, iy, iz, iph;
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(Rp, mode_idx, signal, mu_x1, Edep, ibin, ivox, ix, iy, iz, iph) collapse(4)
#endif
  for (iz=0; iz<NZ; iz++) { // loop on detector voxels
    for (iy=0; iy<NY; iy++) { // loop on detector voxels
      for (ix=0; ix<NX; ix++) { // loop on detector voxels
	for (iph=0; iph<PhotonNum; iph++) { // loop on event multiplicity
	  ivox=iz*NY*NX+iy*NX+ix;
	  //cout << ivox << " " << ivox << endl;
	  // loop on events generated by the source device
	  for(SourceClones[THREAD_IDX]->Begin();
	      !SourceClones[THREAD_IDX]->End();
	      SourceClones[THREAD_IDX]->Next()) {
	    // extract random point in voxel
	    Rp = RandomPointInVoxel(ivox, rngs[THREAD_IDX]);
#ifdef TIME_PERF
	    clock_t cpu_time = clock();
#endif
	    vect3 dummy_prev_x;
	    // get a photon from the source forcing it to terminate on the voxel
	    // volume:
	    SourceClones[THREAD_IDX]->Out_Photon_x1(&PhotonArray[THREAD_IDX],
						    Rp, &mode_idx, &mu_x1,
						    &Edep, &dummy_prev_x);
#ifdef TIME_PERF
	    outph_time += (double)(clock() - cpu_time)/CLOCKS_PER_SEC;
#endif
	    if (PhotonArray[THREAD_IDX].w != 0) { // check that weight is not 0
	      // if voxel is not parallelepiped, correct the weight
	      if (Shape==1) PhotonArray[THREAD_IDX].w *= PI/6; // sphere
	      else if (Shape>1) PhotonArray[THREAD_IDX].w *= PI/4; // cylinder
	      // evaluate the signal associated to the single event
	      signal = PhotonArray[THREAD_IDX].w*mu_x1*VoxelVol
		*ExpTime/PhotonNum;
	      // Depending on voxel content type, multiply it by the energy
	      if (VoxelType == 1 || VoxelType == 3)
		signal *= Edep;

	      ibin = 0;
	      // check if energy binning is used
	      if ((VoxelType==2 || VoxelType==3) && NBins > 1) {
		// evaluate the energy bin
		ibin = ((int)trunc((Edep-Emin)
				  /(Emax-Emin)*NBins));
		if (ibin >= NBins) { // check for saturation
		  if (SaturateEmax == 1) ibin = NBins - 1;
		  else continue;
		}
		if (ibin <0) {
		  if (SaturateEmin == 1) ibin = 0;
		  else continue;
		}
	      }
	      // add the signal to the voxel bin
#ifdef _OPENMP
#pragma omp atomic
#endif
	      Image[(mode_idx*NBins+ibin)*NZ*NY*NX + ivox] += signal;
	    }
#ifdef _OPENMP
#pragma omp atomic
#endif
	    event_idx++;
	  
	    //if (THREAD_IDX==0) {
	    if (event_idx%ProgrUpdate==0 && event_tot>ProgrUpdate) {
	      printf("Progress %.3f %%\r", 100.*event_idx/event_tot);
	      fflush(stdout);
	      //cout <<outph_time << endl;
	    }
	    //}
	  }
	}
      }
    }
  }
  //end of big for loop
  if (event_tot>ProgrUpdate) {
    printf("\rProgress 100.000 %%\n");
    //cout <<outph_time << endl;
  }
#ifdef TIME_PERF
  cout << "Sample->Out_Photon_x1\n"
    "from detectorarray3d::ForcedAcquisition(...)\n"
       << outph_time << "\n\n";
  cout << "Out_Photon_x1\n"
    "from sample::Out_Photon_x1(photon *Photon, vect3 x1, "
    "double *mu_x1, double *Edep)\n" << outph_time0 << "\n\n";
  cout << "Photon->EnergyDeposition\n"
    "from sample::Out_Photon_x1(photon *Photon, "
    "vect3 x1, double *mu_x1, double *Edep)\n" << edep_time << "\n\n";
  cout << "Source->Out_Photon\n"
    "from sample::Out_Photon(photon *Photon)\n"
       << soutph_time << "\n\n";
  cout << "Comp->Mu(Photon->E)\nfrom sample::Out_Photon(photon *Photon) "
    "(with ScattorderIdx==0)\n"
    "from sample::Out_Photon_x1(photon *Photon, vect3 x1, int scatt_order) "
    "(with ScattorderIdx==0)\n"
    "from PhotonHistory(...)\n" << mu_time0 << "\n\n";
  cout << "Comp->Mu(Photon->E)\nfrom sample::Out_Photon(photon *Photon) "
    "(with ScattorderIdx>0)\n"
    "from sample::Out_Photon_x1(photon *Photon, vect3 x1, int scatt_order) "
    "(with ScattorderIdx>0)\n"
       << mu_time1 << "\n\n";
  cout << "PhotonHistory(...)\nfrom sample::Out_Photon(photon *Photon)\n"
    "from Out_Photon_x1(photon *Photon, vect3 x1, int scatt_order)\n"
       << phist_time << "\n\n";
  cout << "Photon->Scatter(...)\nfrom sample::Out_Photon(photon *Photon)\n"
    "from sample::Out_Photon_x1(photon *Photon, vect3 x1, int scatt_order)\n"
       << scatter_time << "\n\n";
  cout << "Photon->Mu_x1\nfrom "
    "sample::Out_Photon_x1(photon *Photon, vect3 x1, int scatt_order)\n"
       << mux1_time << "\n\n";
  cout << "PhotonSurvivalWeight\nfrom "
    "sample::Out_Photon_x1(photon *Photon, vect3 x1, int scatt_order)\n"
       << psw_time << "\n\n";
  cout << "Source->Out_Photon_x1\nfrom "
    "sample::Out_Photon_x1(photon *Photon, vect3 x1, int scatt_order)\n"
       << soutphx1_time << "\n\n";
  cout << "PhotonHistory main loop (for (int is=1; is<=scatt_order; is++)\n"
       << phist_time0 << "\n\n";
#endif

  return 0;
}


//////////////////////////////////////////////////////////////////////
// generate uncertainty on voxel count using Poisson statistic
//////////////////////////////////////////////////////////////////////
int detectorarray3d::Poisson(randmt_t *rng)
{
  double signal, sigma;
  // loop on modes (scattering orders)
  for (int mode_idx=0; mode_idx<ModeNum; mode_idx++) {
    for (int iz=0; iz<NZ; iz++) { // loop on detector voxels
      for (int iy=0; iy<NY; iy++) {
	for (int ix=0; ix<NX; ix++) {
	  int ivox=iz*NY*NX+iy*NX+ix;
	  for (int ibin=0; ibin<NBins; ibin++) { // loop on energy bins
	    signal = Image[(mode_idx*NBins+ibin)*NZ*NY*NX + ivox];
	    sigma = sqrt(signal);          // approximate poisson deviation 
	    signal +=sigma*GaussRnd_r(rng);// using a Gaussian distribution
	    signal = signal>0 ? signal : 0; // threshold to 0
	    if (RoundFlag == 1) signal = round(signal); // round to integer
	    Image[(mode_idx*NBins+ibin)*NZ*NY*NX + ivox] = signal;
	  }
	}
      }
    }
  }
  
  return 0;
}
//////////////////////////////////////////////////////////////////////
// clear the detector voxel bin contents
//////////////////////////////////////////////////////////////////////
int detectorarray3d::Clear()
{
  for (int mode_idx=0; mode_idx<ModeNum; mode_idx++) {
    for (int iz=0; iz<NZ; iz++) { // loop on detector voxels
      for (int iy=0; iy<NY; iy++) {
	for (int ix=0; ix<NX; ix++) {
	  int ivox=iz*NY*NX+iy*NX+ix;
	  for (int ibin=0; ibin<NBins; ibin++) {
	    Image[(mode_idx*NBins+ibin)*NZ*NY*NX + ivox] = 0;
	  }
	}
      }
    }
  }
  
  return 0;
}

// event multiplicity
long long detectorarray3d::EventMulti()
{
  if(ForceDetectFlag==1) return N*PhotonNum*Source->EventMulti();
  else return PhotonNum*Source->EventMulti();
}

//////////////////////////////////////////////////////////////////////
// method for casting input device to type basesource
/////////////////////////////////////////////////////////////////////
int detectorarray3d::CastInputDevices()

{
  // cast it to type basesource*
  Source = dynamic_cast<basesource*>(InputDevice[0]);
  if (Source==0)
    throw xrmc_exception(string("Device ") + InputDeviceName[0]
			 + " cannot be casted to type basesource\n");

  return 0;
}


//////////////////////////////////////////////////////////////////////
// detector3d initialization method
//////////////////////////////////////////////////////////////////////
int detectorarray3d::Init()
{
  double x, y, z;

  OrthoNormal(ui, uj, uk);  // evaluates uj to form a orthonormal basis

  VoxelVol = VoxelSizeX*VoxelSizeY*VoxelSizeZ; // voxel volume
  N = NX*NY*NZ; // number of voxels
  if (VoxelX!=NULL) delete[] VoxelX;
  VoxelX = new vect3[N]; // array of voxel coordinates 

  double x0 = -(VoxelSizeX*NX)/2 + VoxelSizeX/2; // starting (local) coordinates
  double y0 = -(VoxelSizeY*NY)/2 + VoxelSizeY/2; // on the detector volume
  double z0 = -(VoxelSizeZ*NZ)/2 + VoxelSizeZ/2; //
  int i = 0;

  //if (RunningFasterFlag==0) { // columns running faster than rows
  for (int iz=0; iz<NZ; iz++) { // loop on detector voxels
    for (int iy=0; iy<NY; iy++) {
      for (int ix=0; ix<NX; ix++) {
	x = x0 + VoxelSizeX*ix; // local x coordinate of the voxel
	y = y0 + VoxelSizeY*iy; // local y coordinate of the voxel
	z = z0 + VoxelSizeZ*iz; // local z coordinate of the voxel
	VoxelX[i] = X + ui*x + uj*y + uk*z; // voxel 3d absolute coordinates

	i++; // increase voxel index
      }
    }
  }
  //}
  
  return 0;
}

//////////////////////////////////////////////////////////////////////
// Check if point is inside the detector
//////////////////////////////////////////////////////////////////////
bool detectorarray3d::inVoxel(vect3 x, int *ix, int *iy, int *iz,
	    double *tx, double *ty, double *tz)
{
  double x0, y0, z0, x1, y1, z1;

  // relative coordinates in voxel sides units
  x1=(x-X)*ui/VoxelSizeX + 0.5*NX;
  y1=(x-X)*uj/VoxelSizeY + 0.5*NY;
  z1=(x-X)*uk/VoxelSizeZ + 0.5*NZ;
  
  x0 = floor(x1);
  y0 = floor(y1);
  z0 = floor(z1);

  *ix = (int)x0;
  *iy = (int)y0;
  *iz = (int)z0;

  *tx = x1 - x0;
  *ty = y1 - y0;
  *tz = z1 - z0;
  
  if (x1<0 || y1<0 || z1<0 || *ix>=NX || *iy>=NY || *iz>=NZ) return false;

  return true;
}
  
//////////////////////////////////////////////////////////////////////
// Generates a random point in the voxel volume
//////////////////////////////////////////////////////////////////////
vect3 detectorarray3d::RandomPointInVoxel(int i, randmt_t *rng)
{
  double x, y, z, rx, ry, rz;
  vect3 p;

  if (RandomVoxelFlag == 1) { // check if random position in voxel is enabled
    do {
      rx=2.*Rnd_r(rng)-1; // random real numbers between -1 and 1
      ry=2.*Rnd_r(rng)-1;
      rz=2.*Rnd_r(rng)-1;
    } while ( // if shape is not parallelepiped check that the point is inside
	     ( Shape==1 && (rx*rx+ry*ry+rz*rz>1) ) || //sphere
	     ( Shape==2 && (ry*ry+rz*rz>1) ) || // cylinder, x axis
	     ( Shape==3 && (rx*rx+rz*rz>1) ) || // cylinder, y axis
	     ( Shape==4 && (rx*rx+ry*ry>1) ));    // cylinder, z axis

    x = VoxelSizeX*rx/2; // x,y,z coordinates of the point in voxel volume
    y = VoxelSizeY*ry/2;
    z = VoxelSizeZ*rz/2;
  }
  else {
    x = 0;
    y = 0;
    z = 0;
  }
  p = VoxelX[i] + ui*x + uj*y + uk*z; // 3d absolute coordinates of the point

  return p;
}
//////////////////////////////////////////////////////////////////////
// run the unforced acquisition
//////////////////////////////////////////////////////////////////////
int detectorarray3d::UnforcedAcquisition(basesource **SourceClones,
				       photon *PhotonArray)
{
  int ibin, mode_idx;
  double signal;
  double mu_x1, Edep;
  const int ProgrUpdate=100000;
  long long event_idx=0, event_tot=EventMulti();
  if (event_tot>ProgrUpdate) {
    printf("\nProgress 0 %%\r");
    fflush(stdout);
  }

#ifdef TIME_PERF
  outph_time=0;
  // invoxel_time=0; da aggiungere!!!
#endif

  int ix, iy, iz, ivox;
  double tx, ty, tz;
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(mode_idx, signal, mu_x1, Edep, ibin, ix, iy, iz, ivox, tx, ty, tz) collapse(1)
#endif
  for (int iph=0; iph<PhotonNum; iph++) { // loop on event multiplicity
    // loop on events generated by the source device
    for (SourceClones[THREAD_IDX]->Begin(); !SourceClones[THREAD_IDX]->End();
	 SourceClones[THREAD_IDX]->Next()) {
#ifdef TIME_PERF
      clock_t cpu_time = clock();
#endif
      SourceClones[THREAD_IDX]->Out_Photon(&PhotonArray[THREAD_IDX],
					   &mode_idx, &Edep);
#ifdef TIME_PERF
      outph_time += (double)(clock() - cpu_time)/CLOCKS_PER_SEC;
#endif
#ifdef TIME_PERF
      cpu_time = clock();
#endif
      
      /*
      // TEMP REMOVE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (THREAD_IDX==0 && PhotonArray[THREAD_IDX].x.Elem[1]>-64.0) {
	cout << PhotonArray[THREAD_IDX].x << endl;
	cout << "w: " << std::scientific << PhotonArray[THREAD_IDX].w << endl;
	cout.unsetf(ios_base::floatfield);
	cout << inVoxel(PhotonArray[THREAD_IDX].x, &ix, &iy, &iz,
			&tx, &ty, &tz) << endl;
	cout << ix << " " << iy << " " << iz << endl;
	
      }

      //
      */
      
      if (PhotonArray[THREAD_IDX].w>0 &&
	  inVoxel(PhotonArray[THREAD_IDX].x, &ix, &iy, &iz, &tx, &ty, &tz)) {
	// aggiungere controlli su shape!!!!
	//#ifdef TIME_PERF da attivare!!!!
	// da attivare!!!!
	//invoxel_time += (double)(clock() - cpu_time)/CLOCKS_PER_SEC;
	//#endif
	
	signal = PhotonArray[THREAD_IDX].w*ExpTime/PhotonNum;
	// Depending on voxel content type, multiply it by the energy deposit
	if (VoxelType == 1 || VoxelType == 3) signal *= Edep;
	ibin = 0;
	// check if energy binning is used
	if ((VoxelType==2 || VoxelType==3) && NBins > 1) {
	  // evaluate the energy bin
	  ibin = ((int)trunc((Edep-Emin)
			     /(Emax-Emin)*NBins));
	  if (ibin >= NBins) { // check for saturation
	    if (SaturateEmax == 1) ibin = NBins - 1;
	    else continue;
	  }
	  if (ibin <0) {
	    if (SaturateEmin == 1) ibin = 0;
	    else continue;
	  }
	}
	// add the signal to the voxel bin
	ivox=iz*NY*NX+iy*NX+ix;	
#ifdef _OPENMP
#pragma omp atomic
#endif
	Image[(mode_idx*NBins+ibin)*NZ*NY*NX + ivox] += signal;
      }
#ifdef _OPENMP
#pragma omp atomic
#endif
      event_idx++;
	  
      //if (THREAD_IDX==0) {
      if (event_idx%ProgrUpdate==0 && event_tot>ProgrUpdate) {
	printf("Progress %.3f %%\r", 100.*event_idx/event_tot);
	fflush(stdout);
	//cout <<outph_time << endl;
      }
      //}
    }
  }
  //end of big for loop

  if (event_tot>ProgrUpdate) {
    printf("\rProgress 100.000 %%\n");
    //cout <<outph_time << endl;
  }

#ifdef TIME_PERF
  cout << "Sample->Out_Photon\n"
    "from detectorarray3d::UnforcedAcquisition(...)\n"
       << outph_time << "\n\n";
#endif

  return 0;
}

