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
//     phcdetector.cpp           //
//        31/01/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Methods of the class phcdetector
//
#include <cmath>
#include <iostream>
#include <complex>
#include "xrmc_math.h"
#include "phcdetector.h"
#include "xrmc_photon.h"
#include "xrmc_exception.h"
#include "xrmc_algo.h"
//#include "xrmc_arrayNd.h"
#include "xrmc_sample.h"
#include "alloc.h"
#include "fft.h"
#ifdef _OPENMP
#include <omp.h>
#define THREAD_MAXNUM omp_get_max_threads()
#define THREAD_IDX omp_get_thread_num()
#else
#define THREAD_MAXNUM 1
#define THREAD_IDX 0
#endif
#define MAX(x,y) ((x) > (y) ? (x) : (y))

using namespace std;
using namespace xrmc_algo;
using namespace arrayNd;

// destructor
phcdetector::~phcdetector() {
  //if (PixelX!=NULL) delete[] PixelX; // should be called on base class
  //if (Image!=NULL) arrayNd::free_double_array3d(Image); // as above
  if (Transm!=NULL) free_3d_double(Transm);
  if (Propagator!=NULL) free_3d_double(Propagator);
  if (ImageWeight!=NULL) free_3d_double(ImageWeight);
}

// constructor
phcdetector::phcdetector(string dev_name) : detectorarray(dev_name) {
  Runnable = true;
  //SaveDataName.push_back("Image"); done by base class detectorarray
  NInputDevices = 1;
  InputDeviceCommand.clear();
  InputDeviceDescription.clear();
  InputDeviceCommand.push_back("SampleName");
  InputDeviceDescription.push_back("Sample input device name");

  PixelX=NULL;
  Image=NULL;
  Transm=NULL;
  Propagator=NULL;
  ImageWeight=NULL;

  NX = NY = N = PhotonNum = NBins = 0;
  ModeNum = 1;
  SetDevice(dev_name, "phcdetector");
}


//////////////////////////////////////////////////////////////////////
// run the acquisition
//////////////////////////////////////////////////////////////////////
int phcdetector::Acquisition()
{
  vect3 Rp, DRp;

  Clear(); // clear the detector bin contents

  //init_randmt_auto_r(Rng);
  //Photon.Rng = Rng;

  //initialize rngs
  randmt_t **rngs;
  rngs = new randmt_t*[THREAD_MAXNUM];
  SampleClones = new sample*[THREAD_MAXNUM];
  PhotonArray = new photon[THREAD_MAXNUM];

  cout << "Maximum number of threads: "  << THREAD_MAXNUM << endl;
  if (THREAD_MAXNUM==1) {
    SampleClones[0] = Sample;
    rngs[0] = new_randmt();
    if (Seeds.size()>0 && Seeds[0]!=0) init_randmt_r(rngs[0], Seeds[0]);
    else init_randmt_auto_r(rngs[0]);
    SampleClones[0]->SetRng(rngs[0]);
    PhotonArray[0].Rng = rngs[0];
  }
  else {
    for (int i=0; i<THREAD_MAXNUM; i++) {
      SampleClones[i] =dynamic_cast<sample*>(Sample->Clone(InputDeviceName[0]));
      if (SampleClones[i]==0)
	throw xrmc_exception(string("Device ") + InputDeviceName[0] +
			     " cannot be casted to type sample\n");
      rngs[i] = new_randmt();
      if (rngs[i] == NULL)
	throw xrmc_exception(string("Could not allocate new RNG\n"));
      if (i<(int)Seeds.size() && Seeds[i]!=0) {
	init_randmt_r(rngs[i], Seeds[i]);
	cout << "\tThread " << i << " , seed " << Seeds[i] << endl;
      } 
      else init_randmt_auto_r(rngs[i]);
      SampleClones[i]->SetRng(rngs[i]);
      PhotonArray[i].Rng = rngs[i];
    }
  }

  PhCAcquisition(rngs);
    
  // generate uncertainty on pixel count using Poisson statistic
  if (PoissonFlag == 1) {
    Poisson(rngs[0]);
  }

  //free rngs
  for (int i=0; i<THREAD_MAXNUM; i++) {
    free_randmt(rngs[i]);
  }
  delete [] rngs;

  delete [] SampleClones;
  delete [] PhotonArray;

  return 0;
}

//////////////////////////////////////////////////////////////////////
// run the forced acquisition
//////////////////////////////////////////////////////////////////////
int phcdetector::PhCAcquisition(randmt_t **rngs)
{
  //const double KEV2ANGST=12.398;
  //complex<double> I = complex<double>(0,1);
  //int bin;
  //vect3 Rp, DRp;
  //double Pgeom, signal;

  Sample->PhCOn();
  for (int i=0; i<THREAD_MAXNUM; i++) {
    SampleClones[i]->PhCOn();
  }

  const int ProgrUpdate=1;
  long long event_idx=0, event_tot=EventMulti();
  if (event_tot>ProgrUpdate) {
    printf("\nProgress 0 %%\r");
    fflush(stdout);
  }
  
  //for (bin=0; bin<NBins; bin++) {
  //  for (int ipix=0; ipix<NX*NY; ipix++) {
  //    ImageWeight[ipix][bin] = 0;
  //  }
  //}
  
  
  //cout << "Sample->EventMulti() " << Sample->EventMulti() << endl;

  #ifdef _OPENMP
  #pragma omp parallel for default(shared) collapse(2)
  #endif
  for (int iph=0; iph<PhotonNum; iph++) { // loop on event multiplicity
    for (int ie=0; ie<Sample->EventMulti(); ie++) {
      //for(CloneBegin(); !CloneEnd(); CloneNext()) {//clone (energy) event loop
      CloneAcquisition(ie, THREAD_IDX, Transm[THREAD_IDX],
		       Propagator[THREAD_IDX], rngs);
      /*
      for (int iy=0; iy<NY; iy++) {
	for (int ix=0; ix<NX; ix++) {
	  int ipix = NX*iy + ix;
	  int i1 = 2*ipix;
	  // transm = complex<double>(Propagator[i1],Propagator[i1+1]);
	  complex<double> transm = complex<double>(Transm[i1],Transm[i1+1])
	    * Csiz / (I*Lambda0*R1*R12);	  
	  double signal = norm(transm)*ExpTime*PixelSizeX
	    *PixelSizeY/PhotonNum;
	  int ix1 = ix - (NX - ImageNx)/2;
	  int iy1 = iy - (NY - ImageNy)/2;
	  
	  //#ifdef _OPENMP
	  //#pragma omp atomic
	  //#endif
	  if (ix1>=0 && ix1<ImageNx && iy1>=0 && iy1<ImageNy)
	    Image[0][ImageNx*iy1+ix1][bin] += signal
	      *ImageWeight[THREAD_IDX][ipix][bin];
	  //Image[0][NX*iy+ix][bin] += signal;
	}
      }
      */
      #ifdef _OPENMP
      #pragma omp atomic
      #endif
      event_idx++;
	
      //if (THREAD_IDX==0) {
      if (event_idx%ProgrUpdate==0 && event_tot>ProgrUpdate) {
	  printf("Progress %.3f %%\r", 100.*event_idx/event_tot);
	  fflush(stdout);
      }
	//}
      
      
    }
  }
  
  //end of big for loop
  if (event_tot>ProgrUpdate) {
    printf("\rProgress 100.000 %%\n");
  }
  
  return 0;
}


//////////////////////////////////////////////////////////////////////
// phcdetector initialization before run method
//////////////////////////////////////////////////////////////////////
int phcdetector::RunInit()
{
  int n = ImageNx + 2*(NScreenBorderX+NInterpBorderX);
  for(NX=1; NX<n; NX<<=1);
  NBx = NInterpBorderX + (NX - n)/2;
  n = ImageNy + 2*(NScreenBorderY+NInterpBorderY);
  for(NY=1; NY<n; NY<<=1);
  NBy = NInterpBorderY + (NY - n)/2;
  NN1[0] = NY;
  NN1[1] = NX;
  cout << "NN1[0] " << NN1[0] << endl;
  cout << "NN1[1] " << NN1[1] << endl;

  Init();
  ModeNum = 1; // number of modes
  if (Image!=NULL) free_3d_double(Image);
  // allocate the image array

  // TMP: RESTORE THIS
  cout << "ImageNx " << ImageNx << endl;
  cout << "ImageNy " << ImageNy << endl;
  Image = alloc_3d_double(ModeNum*NBins, ImageNy, ImageNx);
  //Image = double_array3d(ModeNum, N, NBins);
  //if (PhaseImage!=NULL) free_double_array2d(PhaseImage);
  // allocate the phase image array
  //PhaseImage = double_array2d(N, NBins);
  if (ImageWeight!=NULL) free_3d_double(ImageWeight);
  // allocate the image weights array
  ImageWeight = alloc_3d_double(THREAD_MAXNUM*NBins, NY, NX);
  
  Rv2 = X - Sample->Source->X;

  R2 = Rv2.Mod();
  // TO DO: CHECK ON DISTANCES, DIV BY ZERO,....
  X2 = Rv2*ui;
  Y2 = Rv2*uj;
  Z2 = -(Rv2*uk);
  Csix = X2/R2;
  Csiy = Y2/R2;
  Csiz = Z2/R2;
  Z1 = Z2 - Z12;
  R12 = R2*Z12/Z2;
  R1 = R2 - R12;
  Reff = R12*R1/R2;  
  dx1 = PixelSizeX*R1/R2;
  dy1 = PixelSizeY*R1/R2;

  Lx1 = L1Coeff*Reff/PixelSizeX/2;
  Ly1 = L1Coeff*Reff/PixelSizeY/2;
  Sigmax1 = Sigma1Coeff*L1Coeff*Reff/PixelSizeX/2;
  Sigmay1 = Sigma1Coeff*L1Coeff*Reff/PixelSizeY/2;

  if (Transm!=NULL) free_3d_double(Transm);
  Transm = alloc_3d_double(THREAD_MAXNUM, NY, 2*NX);
  if (Propagator!=NULL) free_3d_double(Propagator);
  Propagator = alloc_3d_double(THREAD_MAXNUM, NY, 2*NX);

  return 0;
}

// event multiplicity
long long phcdetector::EventMulti()
{
  return PhotonNum*Sample->EventMulti();
}

//////////////////////////////////////////////////////////////////////
// method for casting input device to type sample
/////////////////////////////////////////////////////////////////////
int phcdetector::CastInputDevices()
{
  // cast it to type sample*
  Sample = dynamic_cast<sample*>(InputDevice[0]);
  if (Sample==0)
    throw xrmc_exception(string("Device ") + InputDeviceName[0]
			 + " cannot be casted to type sample\n");

  return 0;
}

//////////////////////////////////////////////////////////////////////
// method for filling interpolation borders
/////////////////////////////////////////////////////////////////////
int phcdetector::FillBorders(double **transm)
{
  double t, u, re0, re1, re, im0, im1, im;

  // fill left and right borders; loop on rows
  for (int iy=NBy; iy<NY-NBy; iy++) { // loop on rows
    int ix = NBx;
    //int i1 = 2*(NX*iy+ix); // interpolation index
    re0 = 1;               // interpolation left real value
    re1 = transm[iy][2*ix];      // interpolation right real value
    im0 = 0;               // interpolation left imaginary value
    im1 = transm[iy][2*ix+1];    // interpolation right imaginary value
    for (int ix=0; ix<NBx; ix++) { // left border, loop on columns
      t = (double)ix/NBx;
      u = 1. - t;
      re = re0*u + re1*t;
      im = im0*u + im1*t;
      //int ipix = NX*iy + ix;
      //int i1 = 2*ipix;
      transm[iy][2*ix] = re;
      transm[iy][2*ix+1] = im;
      //Image[0][ipix][bin] = im;
    }
  }
  for (int iy=0; iy<NY; iy++) { // right border, loop on rows
    int ix = NX-NBx-1;
    //int i1 = 2*(NX*iy+ix); // interpolation index
    re0 = transm[iy][2*ix];      // interpolation left real value
    re1 = 1;               // interpolation right real value
    im0 = transm[iy][2*ix+1];    // interpolation left imaginary value
    im1 = 0;               // interpolation right imaginary value
    for (int ix=NX-NBx; ix<NX; ix++) { // top border, loop on rows
      t = (double)(ix-(NX-NBx))/NBx;
      u = 1. - t;
      re = re0*u + re1*t;
      im = im0*u + im1*t;
      //int ipix = NX*iy + ix;
      //int i1 = 2*ipix;
      transm[iy][2*ix] = re;
      transm[iy][2*ix+1] = im;
      //Image[0][ipix][bin] = im;
    }
  }

  // fill top and bottom borders; loop on columns
  for (int ix=0; ix<NX; ix++) { // top border loop on columns
    int iy = NBy;
    //int i1 = 2*(NX*iy+ix); // interpolation index
    re0 = 1;               // interpolation top real value
    re1 = transm[iy][2*ix];      // interpolation bottom real value
    im0 = 0;               // interpolation top imaginary value
    im1 = transm[iy][2*ix+1];    // interpolation bottom imaginary value
    for (int iy=0; iy<NBy; iy++) { // top border, loop on rows
      t = (double)iy/NBy;
      u = 1. - t;
      re = re0*u + re1*t;
      im = im0*u + im1*t;
      //int ipix = NX*iy + ix;
      //int i1 = 2*ipix;
      transm[iy][2*ix] = re;
      transm[iy][2*ix+1] = im;
      //Image[0][ipix][bin] = im;
    }
  }
  for (int ix=0; ix<NX; ix++) { // bottom border, loop on columns
    int iy = NY-NBy-1;
    //int i1 = 2*(NX*iy+ix); // interpolation index
    re0 = transm[iy][2*ix];      // interpolation top real value
    re1 = 1;               // interpolation bottom real value
    im0 = transm[iy][2*ix+1];    // interpolation top imaginary value
    im1 = 0;               // interpolation bottom imaginary value
    for (int iy=NY-NBy; iy<NY; iy++) { // top border, loop on rows
      t = (double)(iy-(NY-NBy))/NBy;
      u = 1. - t;
      re = re0*u + re1*t;
      im = im0*u + im1*t;
      //int ipix = NX*iy + ix;
      //int i1 = 2*ipix;
      transm[iy][2*ix] = re;
      transm[iy][2*ix+1] = im;
      //Image[0][ipix][bin] = im;
    }
  }

  return 0;
}

int phcdetector::SaveData(string data_name, string file_name)
{
  int tmp_nx=NX, tmp_ny=NY;
  NX = ImageNx;
  NY = ImageNy;
  detectorarray::SaveData(data_name, file_name);
  NX = tmp_nx;
  NY = tmp_ny;

  return 0;
}

//////////////////////////////////////////////////////////////////////
// clear the detector pixel bin contents
//////////////////////////////////////////////////////////////////////
int phcdetector::Clear()
{
  for (int mode_idx=0; mode_idx<ModeNum; mode_idx++) {
    for (int iy=0; iy<ImageNy; iy++) {
      for (int ix=0; ix<ImageNx; ix++) {
	for (int ibin=0; ibin<NBins; ibin++) {
	  Image[mode_idx*NBins+ibin][iy][ix] = 0;
	}
      }
    }
  }

  return 0;
}

// begin event loop on clones
int phcdetector::CloneBegin()
{
  for (int i=0; i<THREAD_MAXNUM; i++) {
    SampleClones[i]->Begin();
  }

  return 0;
}

// next step on clone event loop
int phcdetector::CloneNext()
{
  for (int i=0; i<THREAD_MAXNUM; i++) {
    SampleClones[i]->Next();
  }

  return 0;
}

// end clone event loop
bool phcdetector::CloneEnd()
{
  return SampleClones[0]->End();
}

int phcdetector::CloneAcquisition(int ie, int thread_idx,
				  double **transm, double **propag,
				  randmt_t **rngs)
{
   const double KEV2ANGST=12.398;
   complex<double> I = complex<double>(0,1);
   int bin;
   vect3 Rp, DRp;
   //double Pgeom, signal;
   SampleClones[thread_idx]->Begin();
   for (int i=0; i<ie; i++) SampleClones[thread_idx]->Next();
     
   double E0 = SampleClones[thread_idx]->GetPhC_E0();
   //cout << "E0: " << E0 << endl;
   double Lambda0 = 1.e-8*KEV2ANGST/E0; // lambda in cm
   //printf("Lambda: %.20e\n", Lambda0);
   double K0 = 2.*PI/Lambda0;

   SampleClones[thread_idx]->Comp->Mu(E0);
   SampleClones[thread_idx]->Comp->Delta(E0);
   bin = 0;
   // check if energy binning is used
   if ((PixelType==2 || PixelType==3) && NBins > 1) {
     // evaluate the energy bin
     bin = ((int)trunc((E0-Emin)/(Emax-Emin)*NBins));
     if (bin < 0 || bin>= NBins) return 0;
   }
   
   for (int iy=0; iy<NY; iy++) { // loop on detector pixels
     for (int ix=0; ix<NX; ix++) {
       int ipix = NX*iy + ix;
       //int i1 = 2*ipix;
       if (ix>=NBx && ix<(NX-NBx) && iy>=NBy && iy<(NY-NBy)) { 
	 //cout << ipix << " " << iph << endl;
	 Rp = RandomPointOnPixel(ipix, rngs[thread_idx]);
	 
	 SampleClones[thread_idx]->Out_Phase_Photon_x1
	   (&PhotonArray[thread_idx], Rp);
	 double muL=0, deltaL=0; //, Pgeom=0;
	 if (PhotonArray[thread_idx].w != 0) { //check that weight is not 0
	   muL = SampleClones[thread_idx]->Path->MuL;
	   deltaL = SampleClones[thread_idx]->Path->DeltaL;
	   // if pixel is elliptical, correct the weight
	   if (Shape==1) PhotonArray[thread_idx].w *= PI/4;
	   
	   // IL PESO GEOMETRICO VA CORRETTO?
	   //DRp = Rp - PhotonArray[thread_idx].x;
	   //Pgeom = dOmega(DRp);// evaluates the factor dOmega
	 }
	 ImageWeight[thread_idx*NBins+bin][iy][ix] = PhotonArray[thread_idx].w;
	 complex<double> tr = exp(-I*K0*deltaL)*exp(-muL/2);
	 transm[iy][2*ix] = tr.real();
	 transm[iy][2*ix+1] = tr.imag();
	 //Image[0][ipix][bin] = transm.imag();
       }
       //else Image[0][ipix][bin] = 0;
       double x1, y1;
       if (ix < NX/2) x1 = dx1*ix;
       else x1 = dx1*(ix - NX);
       if (iy < NY/2) y1 = dy1*iy;
       else y1 = dy1*(iy - NY);
       double x1Csi = x1*Csix + y1*Csiy;
       complex<double> pr =
	 exp(I*K0/(2*Reff)*(x1*x1+y1*y1 - x1Csi*x1Csi));
       
       if (fabs(x1) > Lambda0*Lx1) {
	 double ratio = (fabs(x1) - Lambda0*Lx1)/Sigmax1/Lambda0;
	 pr *= exp(-ratio*ratio/2);
       }
       if (fabs(y1) > Lambda0*Ly1) {
	 double ratio = (fabs(y1) - Lambda0*Ly1)/Sigmay1/Lambda0;
	 pr *= exp(-ratio*ratio/2);
       }
       
       propag[iy][2*ix] = pr.real();
       propag[iy][2*ix+1] = pr.imag();
     }
   }
   FillBorders(transm);

   int Nmax = MAX(NY, NX);
   int *ip = alloc_1d_int(2 + (int) sqrt(Nmax + 0.5));
   Nmax = MAX(NY, 2*NX) * 3 / 2;
   double *w = alloc_1d_double(Nmax);
   ip[0] = 0;

   cdft2d(NY, 2*NX, 1, transm, NULL, ip, w);
   cdft2d(NY, 2*NX, 1, propag, NULL, ip, w);
   
   //fourn(transm[0]-1, NN1-1, 2, 1);
   //fourn(propag[0]-1, NN1-1, 2, 1);
   
   for (int iy=0; iy<NY; iy++) {
     for (int ix=0; ix<NX; ix++) {
       //int i1 = 2*(NX*iy + ix);
       complex<double> c1=complex<double>(transm[iy][2*ix],transm[iy][2*ix+1]);
       complex<double> c2=complex<double>(propag[iy][2*ix],propag[iy][2*ix+1]);
       c1 *= c2;
       transm[iy][2*ix] = dx1*dy1*c1.real()/NX/NY;
       transm[iy][2*ix+1] = dx1*dy1*c1.imag()/NX/NY;
     }
   }
   cdft2d(NY, 2*NX, -1, transm, NULL, ip, w);
   //fourn(transm[0]-1, NN1-1, 2, -1);
   
   for (int iy=0; iy<NY; iy++) {
     for (int ix=0; ix<NX; ix++) {
       //int ipix = NX*iy + ix;
       //int i1 = 2*ipix;
       // tr = complex<double>(propag[i1],propag[i1+1]);
       complex<double> tr=complex<double>(transm[iy][2*ix],transm[iy][2*ix+1])
	 * Csiz / (I*Lambda0*R1*R12);	  
       double signal = norm(tr)*ExpTime*PixelSizeX
	 *PixelSizeY/PhotonNum;
       int ix1 = ix - (NX - ImageNx)/2;
       int iy1 = iy - (NY - ImageNy)/2;
       
      if (ix1>=0 && ix1<ImageNx && iy1>=0 && iy1<ImageNy)
	#ifdef _OPENMP
	#pragma omp atomic
	#endif
	 Image[bin][iy1][ix1] += signal
	   *ImageWeight[thread_idx*NBins+bin][iy][ix];
       //Image[0][NX*iy+ix][bin] += signal;
     }
   }

   return 0;   
}
