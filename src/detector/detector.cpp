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
//     detector.cpp              //
//        31/01/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Methods of the class detectorarray
//
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "xrmc_math.h"
#include "xrmc_detector.h"
#include "xrmc_photon.h"
#include "xrmc_exception.h"
#include "xrmc_algo.h"
#include "xrmc_arrayNd.h"
#include "xrmc_sample.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace xrmc_algo;
using namespace arrayNd;

//////////////////////////////////////////////////////////////////////
// detectorarray initialization method
//////////////////////////////////////////////////////////////////////
int detectorarray::Init()
{
  double x, y;

  PixelSurf = PixelSizeX*PixelSizeY; // pixel surface
  N = NX*NY; // number of pixels
  PixelX = new vect3[N]; // array of pixel coordinates 
  ModeNum = Source->ModeNum(); // number of modes (scattering orders)
 
  Image = double_array3d(ModeNum, N, NBins); // allocate the image array

  double x0 = -(PixelSizeX*NX)/2 + PixelSizeX/2; // starting (local) coordinates
  double y0 = -(PixelSizeY*NY)/2 + PixelSizeY/2; // on the detector plane
  int i = 0;

  if (RunningFasterFlag==0) { // columns running faster than rows
    for (int iy=0; iy<NY; iy++) { // loop on detector pixels
      for (int ix=0; ix<NX; ix++) {
	x = x0 + PixelSizeX*ix; // local x coordinate of the pixel
	y = y0 + PixelSizeY*iy; // local y coordinate of the pixel
	PixelX[i] = X + ui*x + uj*y; // 3d absolute coordinates of the pixel
	i++; // increase pixel index
      }
    }
  }
  else { // rows running faster than columns
    for (int ix=0; ix<NX; ix++) { // loop on detector pixels
      for (int iy=0; iy<NY; iy++) {
	x = x0 + PixelSizeX*ix; // local x coordinate of the pixel
	y = y0 + PixelSizeY*iy; // local y coordinate of the pixel
	PixelX[i] = X + ui*x + uj*y; // 3d absolute coordinates of the pixel
	i++; // increase pixel index
      }
    }
  }
  
  return 0;
}

//////////////////////////////////////////////////////////////////////
// Generates a random point on the pixel surface
//////////////////////////////////////////////////////////////////////
vect3 detectorarray::RandomPointOnPixel(int i)
{
  double x, y, rx, ry;
  vect3 p;

  if (RandomPixelFlag == 1) { // check if random position on pixel is enabled
    do {
      rx=2.*Rnd()-1; // random real numbers between -1 and 1
      ry=2.*Rnd()-1;
    } while (Shape==1 && (rx*rx+ry*ry>1)); // if shape is elliptical check that
                                           // the point is inside the ellipse
    x = PixelSizeX*rx/2; // x,y coordinates of the point on pixel surface
    y = PixelSizeY*ry/2;
  }
  else {
    x = 0;
    y = 0;
  }
  p = PixelX[i] + ui*x + uj*y; // 3d absolute coordinates of the point

  return p;
}

vect3 detectorarray::RandomPointOnPixel(int i, randmt_t *rng)
{
  double x, y, rx, ry;
  vect3 p;

  if (RandomPixelFlag == 1) { // check if random position on pixel is enabled
    do {
      rx=2.*Rnd_r(rng)-1; // random real numbers between -1 and 1
      ry=2.*Rnd_r(rng)-1;
    } while (Shape==1 && (rx*rx+ry*ry>1)); // if shape is elliptical check that
                                           // the point is inside the ellipse
    x = PixelSizeX*rx/2; // x,y coordinates of the point on pixel surface
    y = PixelSizeY*ry/2;
  }
  else {
    x = 0;
    y = 0;
  }
  p = PixelX[i] + ui*x + uj*y; // 3d absolute coordinates of the point

  return p;
}

//////////////////////////////////////////////////////////////////////
// run the acquisition
//////////////////////////////////////////////////////////////////////
int detectorarray::Acquisition()
{
  int bin, mode_idx;
  vect3 Rp, DRp;
  double Pgeom, signal;
  photon Photon;
  const int ProgrUpdate=100000;

  Clear(); // clear the detector bin contents
  int event_idx=0, event_tot=EventMulti();
  if (event_tot>ProgrUpdate) {
    printf("\nProgress 0 %%\r");
    fflush(stdout);
  }

  Photon.rng = NULL;

#ifdef _OPENMP
  basesource **SourceClones = new basesource*[omp_get_max_threads()];
  for (int i = 0 ; i < omp_get_max_threads() ; i++) {
  	SourceClones[i] = Source->Clone(SourceName);
  }

  //initialize rngs
  randmt_t **rngs;
  photon *PhotonClones = new photon[omp_get_max_threads()];
  rngs = new randmt_t*[omp_get_max_threads()];
  for (int i = 0 ; i < omp_get_max_threads() ; i++) {
  	rngs[i] = new_randmt();
	if (rngs[i] == NULL)
    		throw xrmc_exception(string("Could not allocate new RNG\n"));
  	if (omp_get_max_threads() == 1)
		init_randmt_r(rngs[i], 54123UL);
	else
		init_randmt_auto_r(rngs[i]);
	dynamic_cast<sample*>(SourceClones[i])->rng = rngs[i];
	dynamic_cast<sample*>(SourceClones[i])->Path->rng = rngs[i];
	dynamic_cast<source*>(dynamic_cast<sample*>(SourceClones[i])->Source)->rng = rngs[i];
	dynamic_cast<source*>(dynamic_cast<sample*>(SourceClones[i])->Source)->Spectrum->rng = rngs[i];
	PhotonClones[i].rng = rngs[i];
  }
  //event_tot /= omp_get_max_threads();
  omp_lock_t print_lock;
  omp_init_lock(&print_lock);
#endif

#ifdef _OPENMP
#define Source SourceClones[omp_get_thread_num()]
#define Photon (PhotonClones[omp_get_thread_num()])
#pragma omp parallel for default(shared) private(DRp, Rp, mode_idx, Pgeom, signal, bin) schedule(guided) collapse(2)
#endif
  for (int ipix=0; ipix<N; ipix++) { // loop on detector pixels
    for (int iph=0; iph<PhotonNum; iph++) { // loop on event multiplicity
      // loop on events generated by the source device
      for (Source->Begin(); !Source->End(); Source->Next()) {
#ifdef _OPENMP
	Rp = RandomPointOnPixel(ipix, rngs[omp_get_thread_num()]); // extract random point on pixel
#else
	Rp = RandomPointOnPixel(ipix); // extract random point on pixel
#endif
	// get a photon from the source forcing it to terminate on the pixel
        // surface: 
	Source->Out_Photon_x1(&Photon, Rp, &mode_idx);
	if (Photon.w != 0) { // check that weight is not 0
	  if (Shape==1) Photon.w *= PI/4; // if pixel is elliptical, correct
	                                  // the weight
	  DRp = Rp - Photon.x;
	  Pgeom = dOmega(DRp);// evaluates the factor dOmega, related to
	  // the probability that the last photon trajectory crosses the pixel
	  signal = Photon.w*Pgeom*ExpTime/PhotonNum; // evaluate the signal
	  // associated to the single event
	  // Depending on pixel content type, multiply it by the energy
	  if (PixelType == 1 || PixelType == 3) signal *= Photon.E;
	  bin = 0;
	  // check if energy binning is used
	  if ((PixelType==2 || PixelType==3) && NBins > 1) {
	    // evaluate the energy bin
	    bin = ((int)trunc((Photon.E-Emin)/(Emax-Emin)*NBins));
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
#pragma omp atomic
	  Image[mode_idx][ipix][bin] += signal;
	}
#pragma omp atomic
	event_idx++;
#ifdef _OPENMP
	if (omp_get_thread_num() == 0)
#endif
	{
	  if (event_idx%ProgrUpdate==0 && event_tot>ProgrUpdate) {
	    printf("Progress %.3f %%\r", 100.*event_idx/event_tot);
	    fflush(stdout);
	  }
	}
      }
    }
  }
//end of big for loop
#undef Source


  // generate uncertainty on pixel count using Poisson statistic
  if (PoissonFlag == 1) {
#ifdef _OPENMP
	Poisson(rngs[0]);
#else
  	Poisson();
#endif
  }
  if (event_tot>ProgrUpdate) {
    printf("\rProgress 100.000 %%\n");
  }

#ifdef _OPENMP
  //free rngs
  for (int i = 0 ; i < omp_get_max_threads() ; i++) {
	free_randmt(rngs[i]);
  }
  delete [] rngs;

#endif
  return 0;
}
//////////////////////////////////////////////////////////////////////
// generate uncertainty on pixel count using Poisson statistic
//////////////////////////////////////////////////////////////////////
int detectorarray::Poisson()
{
  double signal, sigma;
  // loop on modes (scattering orders)
  for (int mode_idx=0; mode_idx<ModeNum; mode_idx++) {
    for (int ipix=0; ipix<N; ipix++) { // loop on detector pixels
      for (int ibin=0; ibin>NBins; ibin++) { // loop on energy bins
	signal = Image[mode_idx][ipix][ibin];
	sigma = sqrt(signal);          // approximate poisson deviation 
	signal += sigma*GaussRnd(); // using a random Gaussian distribution
	signal = signal>0 ? signal : 0; // threshold to 0
	if (RoundFlag == 1) signal = round(signal); // round to integer
	Image[mode_idx][ipix][ibin] = signal;
      }
    }
  }

  return 0;
}

int detectorarray::Poisson(randmt_t *rng)
{
  double signal, sigma;
  // loop on modes (scattering orders)
  for (int mode_idx=0; mode_idx<ModeNum; mode_idx++) {
    for (int ipix=0; ipix<N; ipix++) { // loop on detector pixels
      for (int ibin=0; ibin>NBins; ibin++) { // loop on energy bins
	signal = Image[mode_idx][ipix][ibin];
	sigma = sqrt(signal);          // approximate poisson deviation 
	signal += sigma*GaussRnd_r(rng); // using a random Gaussian distribution
	signal = signal>0 ? signal : 0; // threshold to 0
	if (RoundFlag == 1) signal = round(signal); // round to integer
	Image[mode_idx][ipix][ibin] = signal;
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
    for (int ipix=0; ipix<N; ipix++) {
      for (int ibin=0; ibin<NBins; ibin++) {
	Image[mode_idx][ipix][ibin] = 0;
      }
    }
  }

  return 0;
}

// event multiplicity
int detectorarray::EventMulti()
{
  return N*PhotonNum*Source->EventMulti();
}

//////////////////////////////////////////////////////////////////////
// method for importing a basesource device
/////////////////////////////////////////////////////////////////////
int detectorarray::ImportDevice(xrmc_device_map *dev_map)

{
  // check if the input device is defined
  xrmc_device_map::iterator it = dev_map->find(SourceName);

  // if not display error and exit
  if (it==dev_map->end())
    throw xrmc_exception(string("Device ") + SourceName
			 + " not found in device map\n");

  // get device pointer from the device map
  xrmc_device *dev_pt = (*it).second;
  // cast it to type basesource*
  Source = dynamic_cast<basesource*>(dev_pt);
  if (Source==0)
    throw xrmc_exception(string("Device ") + SourceName
			 + " cannot be casted to type basesource\n");

  Init();  // initialize the detectorarray 

  return 0;
}

//////////////////////////////////////////////////////////////////////
// evaluates the factor dOmega, related to
// the probability that the last photon trajectory crosses the pixel
//////////////////////////////////////////////////////////////////////
double detectorarray::dOmega(vect3 DRp)
{
  double r = DRp.Mod();        // distance from interaction point to the
                               // intersection of the trajectory with the pixel
  double num = -(DRp*uk)*PixelSurf; // numerator
  if (num<=0) return 0;
  double denom = r*r*r; // denominator
  if (num+denom == num) return dOmegaLim; // check for overflow in the ratio
  double dO = num/denom;
  if (dO>dOmegaLim) return dOmegaLim; // threshold to dOmegaLim

  return dO;
}
