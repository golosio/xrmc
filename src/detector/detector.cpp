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
//     detector.cpp              //
//        05/12/2017             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Methods of the class detectorarray
//
#include <cmath>
#include <cstring>
#include <iostream>
#include "xrmc_math.h"
#include "xrmc_detector.h"
#include "xrmc_photon.h"
#include "xrmc_exception.h"
#include "xrmc_algo.h"
#include "xrmc_arrayNd.h"
#include "xrmc_sample.h"
#include "image_convolution.h"
				    //#include "fft.h"
#ifdef TIME_PERF
#include <ctime>
extern double outph_time, soutph_time, mu_time0, scatter_time,
				       phist_time, mu_time1, mux1_time,
				       psw_time, soutphx1_time, phist_time0;
double inters_time;
#endif

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
detectorarray::~detectorarray() {
  //if (PixelX!=NULL) delete[] PixelX; // should becalled on base class
  if (Image!=NULL) arrayNd::free_double_array3d(Image);
  if (ConvolutedImage!=NULL) arrayNd::free_double_array3d(ConvolutedImage);
}

// constructor
detectorarray::detectorarray(string dev_name) {
  Runnable = true;
  SaveDataName.push_back("Image");
  SaveDataName.push_back("ConvolutedImage");

  NInputDevices = 1;
  InputDeviceCommand.push_back("SourceName");
  InputDeviceDescription.push_back("Source input device name");

  PixelX=NULL;
  Image=NULL;
  ConvolutedImage=NULL;
  NX = NY = N = PhotonNum = NBins = ModeNum = 0;
  Z12 = 0;
  SetDevice(dev_name, "detectorarray");
}

//////////////////////////////////////////////////////////////////////
// detectorarray initialization before run method
//////////////////////////////////////////////////////////////////////
int detectorarray::RunInit()
{
  Init();

  ModeNum = Source->ModeNum(); // number of modes (scattering orders)
  if (Image!=NULL) free_double_array3d(Image);
  Image = double_array3d(ModeNum*NBins, NY, NX); // allocate the image array

  if (ConvolutedImage!=NULL) free_double_array3d(ConvolutedImage);
  if (ConvolveFlag!=0)
    ConvolutedImage = double_array3d(ModeNum*NBins, NY, NX);	

  return 0;
}

//////////////////////////////////////////////////////////////////////
// run the acquisition
//////////////////////////////////////////////////////////////////////
int detectorarray::Acquisition()
{

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
  if(ForceDetectFlag==0) UnforcedAcquisition(SourceClones, PhotonArray);
  else ForcedAcquisition(SourceClones, PhotonArray, rngs);
    
  // generate uncertainty on pixel count using Poisson statistic
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

  if (ConvolveFlag!=0) Convolve();

  return 0;
}

//////////////////////////////////////////////////////////////////////
// run the forced acquisition
//////////////////////////////////////////////////////////////////////
int detectorarray::ForcedAcquisition(basesource **SourceClones,
				     photon *PhotonArray, randmt_t **rngs)
{
  int bin, mode_idx;
  vect3 Rp;
  double Pgeom, signal;
  const int ProgrUpdate=100000;
  long long event_idx=0, event_tot=EventMulti();
  if (event_tot>ProgrUpdate) {
    printf("\nProgress 0 %%\r");
    fflush(stdout);
  }

#ifdef TIME_PERF
  outph_time=soutph_time=mu_time0=scatter_time=phist_time
    =mu_time1=mux1_time=psw_time=soutphx1_time=phist_time0=0;
#endif

  int ipix, ix, iy, iph;
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(Rp, mode_idx, Pgeom, signal, bin, ipix, ix, iy, iph) collapse(3)
#endif
  for (iy=0; iy<NY; iy++) { // loop on detector pixels
    for (ix=0; ix<NX; ix++) { // loop on detector pixels
      for (iph=0; iph<PhotonNum; iph++) { // loop on event multiplicity
	ipix=iy*NX+ix;
	//cout << ipix << " " << iph << endl;
	// loop on events generated by the source device
	for(SourceClones[THREAD_IDX]->Begin(); !SourceClones[THREAD_IDX]->End();
	    SourceClones[THREAD_IDX]->Next()) {
	  // extract random point on pixel
	  Rp = RandomPointOnPixel(ipix, rngs[THREAD_IDX]);
	  // get a photon from the source forcing it to terminate on the pixel
	  // surface: 
#ifdef TIME_PERF
	    clock_t cpu_time = clock();
#endif
	  SourceClones[THREAD_IDX]->Out_Photon_x1(&PhotonArray[THREAD_IDX], Rp,
						  &mode_idx);
#ifdef TIME_PERF
	    outph_time += (double)(clock() - cpu_time)/CLOCKS_PER_SEC;
#endif
	  if (PhotonArray[THREAD_IDX].w != 0) { // check that weight is not 0
	    // if pixel is elliptical, correct the weight
	    if (Shape==1) PhotonArray[THREAD_IDX].w *= PI/4;
	    
            // evaluates the geometric factor related to
	    // the probability that the last photon trajectory crosses the pixel
	    Pgeom = -(PhotonArray[THREAD_IDX].uk*uk)*PixelSurf;
	    signal = PhotonArray[THREAD_IDX].w*Pgeom*ExpTime/PhotonNum;
	    // evaluate the signal associated to the single event

	    // Depending on pixel content type, multiply it by the energy
	    if (PixelType == 1 || PixelType == 3)
	      signal *= PhotonArray[THREAD_IDX].E;
	    bin = 0;
	    // check if energy binning is used
	    if ((PixelType==2 || PixelType==3) && NBins > 1) {
	      // evaluate the energy bin
	      bin = ((int)trunc((PhotonArray[THREAD_IDX].E-Emin)
				/(Emax-Emin)*NBins));
	      if (bin >= NBins) { // check for saturation
		if (SaturateEmax == 1) bin = NBins - 1;
		else continue;
	      }
	      if (bin <0) {
		if (SaturateEmin == 1) bin = 0;
		else continue;
	      }
	    }
	    // add the signal to the pixel bin
#ifdef _OPENMP
#pragma omp atomic
#endif
	    Image[mode_idx*NBins+bin][iy][ix] += signal;
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
  //end of big for loop
  if (event_tot>ProgrUpdate) {
    printf("\rProgress 100.000 %%\n");
    //cout <<outph_time << endl;
  }

#ifdef TIME_PERF
  cout << "Sample->Out_Photon_x1\n"
    "from detectorarray::ForcedAcquisition(...)\n"
       << outph_time << "\n\n";
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
// run the unforced acquisition
//////////////////////////////////////////////////////////////////////
int detectorarray::UnforcedAcquisition(basesource **SourceClones,
				       photon *PhotonArray)
{
  int bin, mode_idx;
  double signal;
  const int ProgrUpdate=100000;
  long long event_idx=0, event_tot=EventMulti();
  if (event_tot>ProgrUpdate) {
    printf("\nProgress 0 %%\r");
    fflush(stdout);
  }

#ifdef TIME_PERF
  outph_time=soutph_time=mu_time0=scatter_time=phist_time
    =mu_time1=psw_time=phist_time0=inters_time=0;
#endif

  int ix, iy;
  double tx, ty, t;
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(mode_idx, signal, bin, ix, iy, tx, ty, t) collapse(1)
#endif
  for (int iph=0; iph<PhotonNum; iph++) { // loop on event multiplicity
    // loop on events generated by the source device
    for (SourceClones[THREAD_IDX]->Begin(); !SourceClones[THREAD_IDX]->End();
	 SourceClones[THREAD_IDX]->Next()) {
#ifdef TIME_PERF
	    clock_t cpu_time = clock();
#endif
      SourceClones[THREAD_IDX]->Out_Photon(&PhotonArray[THREAD_IDX],
					   &mode_idx);
#ifdef TIME_PERF
	    outph_time += (double)(clock() - cpu_time)/CLOCKS_PER_SEC;
#endif
#ifdef TIME_PERF
	    cpu_time = clock();
#endif
      if (PhotonArray[THREAD_IDX].w>0 && Intersect(PhotonArray[THREAD_IDX].x,
	    PhotonArray[THREAD_IDX].uk, ix, iy, tx, ty, t)) {
#ifdef TIME_PERF
	    inters_time += (double)(clock() - cpu_time)/CLOCKS_PER_SEC;
#endif

#ifdef TIME_PERF
	    cpu_time = clock();
#endif
	SourceClones[THREAD_IDX]->PhotonSurvivalWeight
	  (&PhotonArray[THREAD_IDX], t);
#ifdef TIME_PERF
	    psw_time += (double)(clock() - cpu_time)/CLOCKS_PER_SEC;
#endif
	signal = PhotonArray[THREAD_IDX].w*ExpTime/PhotonNum;
	// Depending on pixel content type, multiply it by the energy
	if (PixelType == 1 || PixelType == 3)
	  signal *= PhotonArray[THREAD_IDX].E;
	bin = 0;
	// check if energy binning is used
	if ((PixelType==2 || PixelType==3) && NBins > 1) {
	  // evaluate the energy bin
	  bin=((int)trunc((PhotonArray[THREAD_IDX].E-Emin)/(Emax-Emin)*NBins));
	  if (bin >= NBins) { // check for saturation
	    if (SaturateEmax == 1) bin = NBins - 1;
	    else continue;
	  }
	  if (bin <=0) {
	    if (SaturateEmin == 1) bin = 0;
	    else continue;
	  }
	}
	// add the signal to the pixel bin

#ifdef _OPENMP
#pragma omp atomic
#endif
	Image[mode_idx*NBins+bin][iy][ix] += signal;
      }
#ifdef _OPENMP
#pragma omp atomic
#endif
      event_idx++;
      
      if (THREAD_IDX==0) {
	if (event_idx%ProgrUpdate==0 && event_tot>ProgrUpdate) {
	  printf("Progress %.3f %%\r", 100.*event_idx/event_tot);
	  fflush(stdout);
	  //cout <<outph_time << endl;
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
  cout << "Source->Out_Photon\n"
    "from detectorarray::UnforcedAcquisition(...)\n"
       << outph_time << "\n\n";
  cout << "Intersect\n"
    "from detectorarray::UnforcedAcquisition(...)\n"
       << inters_time << "\n\n";
  cout << "Source->Out_Photon\n"
    "from sample::Out_Photon(photon *Photon)\n"
       << soutph_time << "\n\n";
  cout << "Comp->Mu(Photon->E)\nfrom sample::Out_Photon(photon *Photon) "
    "(with ScattorderIdx==0)\n"
    "from PhotonHistory(...)\n" << mu_time0 << "\n\n";
  cout << "Comp->Mu(Photon->E)\nfrom sample::Out_Photon(photon *Photon) "
    "(with ScattorderIdx>0)\n"
       << mu_time1 << "\n\n";
  cout << "PhotonHistory(...)\nfrom sample::Out_Photon(photon *Photon)\n"
       << phist_time << "\n\n";
  cout << "Photon->Scatter(...)\nfrom sample::Out_Photon(photon *Photon)\n"
       << scatter_time << "\n\n";
  cout << "PhotonSurvivalWeight\nfrom "
    "sample::Out_Photon_x1(photon *Photon, vect3 x1, int scatt_order)\n"
    "from detectorarray::UnforcedAcquisition(...)\n"
       << psw_time << "\n\n";
  cout << "PhotonHistory main loop (for (int is=1; is<=scatt_order; is++)\n"
       << phist_time0 << "\n\n";
#endif

  return 0;
}

//////////////////////////////////////////////////////////////////////
// generate uncertainty on pixel count using Poisson statistic
//////////////////////////////////////////////////////////////////////
int detectorarray::Poisson(randmt_t *rng)
{
  double signal, sigma;
  // loop on modes (scattering orders)
  for (int mode_idx=0; mode_idx<ModeNum; mode_idx++) {
    for (int iy=0; iy<NY; iy++) { // loop on detector pixels
      for (int ix=0; ix<NX; ix++) { // loop on detector pixels
	for (int ibin=0; ibin<NBins; ibin++) { // loop on energy bins
	  signal = Image[mode_idx*NBins+ibin][iy][ix];
	  sigma = sqrt(signal);          // approximate poisson deviation 
	  signal +=sigma*GaussRnd_r(rng);// using a random Gaussian distribution
	  signal = signal>0 ? signal : 0; // threshold to 0
	  if (RoundFlag == 1) signal = round(signal); // round to integer
	  Image[mode_idx*NBins+ibin][iy][ix] = signal;
	}
      }
    }
  }

  return 0;
}
//////////////////////////////////////////////////////////////////////
// clear the detector pixel bin contents
//////////////////////////////////////////////////////////////////////
int detectorarray::Clear()
{
  for (int mode_idx=0; mode_idx<ModeNum; mode_idx++) {
    for (int iy=0; iy<NY; iy++) { // loop on detector pixels
      for (int ix=0; ix<NX; ix++) { // loop on detector pixels
	for (int ibin=0; ibin<NBins; ibin++) {
	  Image[mode_idx*NBins+ibin][iy][ix] = 0;
	}
      }
    }
  }
  
  return 0;
}

// event multiplicity
long long detectorarray::EventMulti()
{
  if(ForceDetectFlag==1) return N*PhotonNum*Source->EventMulti();
  else return PhotonNum*Source->EventMulti();
}

//////////////////////////////////////////////////////////////////////
// method for casting input device to type basesource
/////////////////////////////////////////////////////////////////////
int detectorarray::CastInputDevices()

{
  // cast it to type basesource*
  Source = dynamic_cast<basesource*>(InputDevice[0]);
  if (Source==0)
    throw xrmc_exception(string("Device ") + InputDeviceName[0]
			 + " cannot be casted to type basesource\n");

  return 0;
}


//////////////////////////////////////////////////////////////////////
// convolve the image
//////////////////////////////////////////////////////////////////////
int detectorarray::Convolve()
{
  for (int ibin=0; ibin<NBins; ibin++) {
    double eff=1;
    if (EfficiencyFlag && (int)Efficiency.size()==NBins)
      eff = Efficiency[ibin];
    // can evaluate efficiency vs bin energy here
    for (int mode_idx=0; mode_idx<ModeNum; mode_idx++) {
      for (int iy=0; iy<NY; iy++) { // loop on detector pixels
	for (int ix=0; ix<NX; ix++) { // loop on detector pixels
	  ConvolutedImage[mode_idx*NBins+ibin][iy][ix] =
	    eff*Image[mode_idx*NBins+ibin][iy][ix]; // multiply by efficiency
	}
      }
    }
  }

  double Z1 = 0;
  if (GaussSourceX.size()>0 || GaussSourceXBin.size()>0 ||
      GaussSourceY.size()>0 || GaussSourceYBin.size()>0) {
    if (Z12<=0)
      throw xrmc_exception("Object-detector distance (Z12) must be defined "
			   " for using source size convolution\n");
    vect3 v = X - Source->SourceX();
     double Z2 = -(v*uk);
    Z1 = Z2 - Z12;
    if (Z1<=0 || (Z12+Z1==Z12))
      throw xrmc_exception("Error: negative or too small source-object "
			   "distance\n");
    if (Z1+Z12==Z1)
      throw xrmc_exception("Error: too small object-detector distance\n");
  }

  int nnx, nny;

  int n = NX + NBx;
  for(nnx=1; nnx<n; nnx<<=1);
  double *psf_xfilt = new double[nnx];
  double *src_xfilt = new double[nnx];
  double *xfilt = new double[nnx];
  double *xdata = new double[nnx];
  int *ipx = new int[nnx];
  ipx[0] = 0;
  double *wx = new double[nnx];
  
  bool psf_x_flag = false;
  if (GaussPSFx.size()>0) {
    convolution::BuildFilt(psf_xfilt, nnx, GaussPSFx, ipx, wx, PixelSizeX);
    psf_x_flag = true;
  }
  bool src_x_flag = false;
  if (GaussSourceX.size()>0) {
    convolution::BuildFilt(src_xfilt, nnx, GaussSourceX, ipx, wx,
			   Z1/Z12*PixelSizeX);
    src_x_flag = true;
  }

  n = NY + NBy;
  for(nny=1; nny<n; nny<<=1);
  double *psf_yfilt = new double[nny];
  double *src_yfilt = new double[nny];
  double *yfilt = new double[nny];
  double *ydata = new double[nny];
  int *ipy = new int[nny];
  ipy[0] = 0;
  double *wy = new double[nny];

  bool psf_y_flag = false;
  if (GaussPSFy.size()>0) {
    convolution::BuildFilt(psf_yfilt, nny, GaussPSFy, ipy, wy, PixelSizeY);
    psf_y_flag = true;
  }
  bool src_y_flag = false;
  if (GaussSourceY.size()>0) {
    convolution::BuildFilt(src_yfilt, nny, GaussSourceY, ipy, wy,
			   Z1/Z12*PixelSizeY);
    src_y_flag = true;
  }

  for (int ibin=0; ibin<NBins; ibin++) {
    if ((int)GaussPSFxBin.size()==NBins && GaussPSFxBin[ibin].size()>0) {
      convolution::BuildFilt(psf_xfilt, nnx, GaussPSFxBin[ibin], ipx, wx,
			     PixelSizeX);
      psf_x_flag = true;
    }
    if ((int)GaussSourceXBin.size()==NBins && GaussSourceXBin[ibin].size()>0) {
      convolution::BuildFilt(src_xfilt, nnx, GaussSourceXBin[ibin], ipx, wx,
			     Z1/Z12*PixelSizeX);
      src_x_flag = true;
    }
    if (psf_x_flag && src_x_flag)
      convolution::ConvolveFilt(psf_xfilt, src_xfilt, xfilt, nnx);
    else if (psf_x_flag)
      memcpy(xfilt, psf_xfilt, sizeof(double)*nnx);
    else if (src_x_flag)
      memcpy(xfilt, src_xfilt, sizeof(double)*nnx);
      
    if ((int)GaussPSFyBin.size()==NBins && GaussPSFyBin[ibin].size()>0) {
      convolution::BuildFilt(psf_yfilt, nny, GaussPSFyBin[ibin], ipy, wy,
			     PixelSizeY);
      psf_y_flag = true;
    } 
    if ((int)GaussSourceYBin.size()==NBins && GaussSourceYBin[ibin].size()>0) {
      convolution::BuildFilt(src_yfilt, nny, GaussSourceYBin[ibin], ipy, wy,
			     Z1/Z12*PixelSizeY);
      src_y_flag = true;
    }
    if (psf_y_flag && src_y_flag)
      convolution::ConvolveFilt(psf_yfilt, src_yfilt, yfilt, nny);
    else if (psf_y_flag)
      memcpy(yfilt, psf_yfilt, sizeof(double)*nny);
    else if (src_y_flag)
      memcpy(yfilt, src_yfilt, sizeof(double)*nny);
    
    for (int mode_idx=0; mode_idx<ModeNum; mode_idx++) {
      if (psf_x_flag || src_x_flag)
	convolution::XFilt(xfilt, ConvolutedImage[mode_idx*NBins+ibin], NX, NY,
			   nnx, xdata, ipx, wx);
      if (psf_y_flag || src_y_flag)
	convolution::YFilt(yfilt, ConvolutedImage[mode_idx*NBins+ibin], NX, NY,
			   nny, ydata, ipy, wy);
    }
  }

  delete[] xfilt;
  delete[] xdata;
  delete[] ipx;
  delete[] wx;

  delete[] yfilt;
  delete[] ydata;
  delete[] ipy;
  delete[] wy;

  return 0;
}

