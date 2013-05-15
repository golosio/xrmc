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
//     beamscreen.cpp            //
//        02/05/2013             //
//   author : Bruno Golosio      //
///////////////////////////////////
// Methods of the class detectorarray
//
#include <cmath>
#include <iostream>
#include "xrmc_math.h"
#include "beamscreen.h"
#include "xrmc_exception.h"
#include "xrmc_algo.h"
				    //#include "xrmc_arrayNd.h"
				    //#ifdef _OPENMP
				    //#include <omp.h>
				    //#endif

using namespace std;
using namespace xrmc_algo;
//using namespace arrayNd;

// destructor
beamscreen::~beamscreen() {
  if (Image!=NULL) delete[] Image;
  if (CumulImage!=NULL) delete[] CumulImage;
  if (CumulEnergy!=NULL) delete[] CumulEnergy;
  if (SumEnergyImage!=NULL) delete[] SumEnergyImage;
  if (CumulXY!=NULL) delete[] CumulXY;
  if (BinWeight!=NULL) delete[] BinWeight;
}

// constructor
beamscreen::beamscreen(string dev_name) {
  Runnable = false;
  NInputDevices = 0;

  PixelX=NULL;
  Image=NULL;
  CumulImage=NULL;
  CumulEnergy=NULL;
  SumEnergyImage=NULL;
  CumulXY=NULL;
  BinWeight=NULL;

  NX = NY = N = NBins = 0;
  Shape = 0;
  RandomPixelFlag = 1;
  RunningFasterFlag = 0;
  SetDevice(dev_name, "beamscreen");
}

//////////////////////////////////////////////////////////////////////
// beamscreen initialization before run method
//////////////////////////////////////////////////////////////////////
int beamscreen::RunInit()
{
  Init();

  if (Image==NULL || N*NBins<=0)
    throw xrmc_exception("Beamscreen image must be loaded before run.\n");

  if (CumulImage!=NULL) delete[] CumulImage;
  if (CumulEnergy!=NULL) delete[] CumulEnergy;
  if (SumEnergyImage!=NULL) delete[] SumEnergyImage;
  if (CumulXY!=NULL) delete[] CumulXY;
  if (BinWeight!=NULL) delete[] BinWeight;

  // allocate the cumulative image array
  CumulImage = new double [2*N*NBins];

  double sum=0;
  for(int i=0; i<2*N*NBins; i++) {
    CumulImage[i] = sum;
    double val = Image[i];
    if (val<0)
      throw xrmc_exception("Image bin contents must be nonnegative.\n");
    sum += val;
  }
  for(int i=0; i<2*N*NBins; i++) {
    Image[i] /= sum;
    CumulImage[i] /= sum;
  }

  if (LoopFlag==0) { // no loop mode
    // allocate the cumulative energy array
    CumulEnergy = new double [2*N*NBins];
    // allocate the probability distribution integrated over the energy
    SumEnergyImage = new double [N];
    
    for (int iy=0; iy<NY; iy++) {
      for (int ix=0; ix<NX; ix++) {
	int i0 = NBins*(NX*iy + ix);
	double sum=0;
	for(int iE=0; iE<NBins; iE++) {
	  CumulEnergy[2*i0+iE] = sum;
	  sum += Image[i0+iE];
	}
	for(int iE=0; iE<NBins; iE++) {
	  CumulEnergy[2*i0+NBins+iE] = sum;
	  sum += Image[N*NBins+i0+iE];
	}
	SumEnergyImage[NX*iy + ix] = sum;
	for(int iE=0; iE<2*NBins; iE++) {
	  CumulEnergy[2*i0+iE] /= sum;
	}
      }
    }
  }
  else { // loop mode
    BinWeight = new double [2*NBins]; // allocates energy bin weight
    CumulXY = new double [2*N*NBins]; // allocates cumulative function for x, y

    for(int ipol=0; ipol<2; ipol++) {
      for(int iE=0; iE<NBins; iE++) {
	int i0 = N*(NBins*ipol + iE);
	double sum=0; // initialize cumulative sum
	for (int iy=0; iy<NY; iy++) {
	  for (int ix=0; ix<NX; ix++) {
	    CumulXY[i0 + NX*iy + ix] = sum;  // cumulative function for x, y
	    int i = NBins*(N*ipol + NX*iy + ix) + iE;
	    sum += Image[i];
	  }
	}
	BinWeight[NBins*ipol+iE] = sum; // probab. weight associated to the bin
	for (int iy=0; iy<NY; iy++) {
	  for (int ix=0; ix<NX; ix++) {
	    CumulXY[i0 + NX*iy + ix] /= sum; // normalize cumulative sum to 1
	  }
	}
      }
    }
  }

  return 0;
}

// Generates a random pixel on the screen
int beamscreen::RandomPixel(int &iE, int &ix, int &iy, int &pol, randmt_t *rng)
{
  int i;
  double R = Rnd_r(rng);
  
  Locate(R, CumulImage, 2*N*NBins, &i);
  iE = i % NBins;
  i /= NBins;
  ix = i % NX;
  i /= NX;
  iy = i % NY;
  pol = i / NY;
  if (pol>1)
    throw xrmc_exception("Error in beam screen array dimensions.\n");

  return 0;
}

// Generates a random point and energy on the screen
vect3 beamscreen::RandomPoint(double &E, int &pol, double &w, randmt_t *rng)
{
  if (LoopFlag==1) { // loop mode
    pol = PolIdx;
    return RandomXY(E, w, Rng);
  }

  // no loop mode
  int iE, ix, iy;
  double rE, rx, ry;

  RandomPixel(iE, ix, iy, pol, rng);
  rE = Rnd_r(rng) - 0.5;
  rx = Rnd_r(rng) - 0.5;
  ry = Rnd_r(rng) - 0.5;

  // evaluates interpolation weight factor
  w = InterpolWeight(iE, ix, iy, pol, rE, rx, ry, 0);
  /*
  if (InterpolFlag==1) {
    if (rE>=0) {
      tE = rE;
      iE1 = iE+1;
      if (iE1>=NBins) iE1=NBins-1;
    }
    else {
      tE = -rE;
      iE1 = iE-1;
      if (iE1<0) iE1=0;
    }
    uE = 1. - tE;
    
    if (rx>=0) {
      tx = rx;
      ix1 = ix+1;
      if (ix1>=NX) ix1=NX-1;
    }
    else {
      tx = -rx;
      ix1 = ix-1;
      if (ix1<0) ix1=0;
    }
    ux = 1. - tx;
    
    if (ry>=0) {
      ty = ry;
      iy1 = iy+1;
      if (iy1>=NY) iy1=NY-1;
    }
    else {
      ty = -ry;
      iy1 = iy-1;
      if (iy1<0) iy1=0;
    }
    uy = 1. - ty;
    ipol = pol*NY*NX*NBins;
    i000 = ipol + NX*NBins*iy  + NBins*ix   + iE;
    i001 = ipol + NX*NBins*iy  + NBins*ix   + iE1;
    i010 = ipol + NX*NBins*iy  + NBins*ix1  + iE;
    i011 = ipol + NX*NBins*iy  + NBins*ix1  + iE1;
    i100 = ipol + NX*NBins*iy1 + NBins*ix   + iE;
    i101 = ipol + NX*NBins*iy1 + NBins*ix   + iE1;
    i110 = ipol + NX*NBins*iy1 + NBins*ix1  + iE;
    i111 = ipol + NX*NBins*iy1 + NBins*ix1  + iE1;
    
    c000 = Image[i000];
    c001 = Image[i001];
    if (c001>rmax*c000) c001 = rmax;
    else c001 /= c000;
    c010 = Image[i010];
    if (c010>rmax*c000) c010 = rmax;
    else c010 /= c000;
    c011 = Image[i011];
    if (c011>rmax*c000) c011 = rmax;
    else c011 /= c000;
    c100 = Image[i100];
    if (c100>rmax*c000) c100 = rmax;
    else c100 /= c000;
    c101 = Image[i101];
    if (c101>rmax*c000) c101 = rmax;
    else c101 /= c000;
    c110 = Image[i110];
    if (c110>rmax*c000) c110 = rmax;
    else c110 /= c000;
    c111 = Image[i111];
    if (c111>rmax*c000) c111 = rmax;
    else c111 /= c000;
    c000 = 1;
    
    w = c000*uy*ux*uE + c001*uy*ux*tE + c010*uy*tx*uE + c011*uy*tx*tE
      + c100*ty*ux*uE + c101*ty*ux*tE + c110*ty*tx*uE + c111*ty*tx*tE;
  }
  else w = 1;
  */
  // central bin energy
  E = Emin + (Emax - Emin)*(0.5 + iE + rE)/NBins;
 // local x, y coordinates of the pixel
  double x = PixelSizeX*(-0.5*NX + 0.5 + ix + rx);
  double y = PixelSizeY*(-0.5*NY + 0.5 + iy + ry);
  vect3 r = X + ui*x + uj*y; // 3d absolute coordinates of the pixel

  return r;
}

// Generates a random energy and polarization for a given trajectory
bool beamscreen::RandomEnergy(vect3 x0, vect3 u, double &E, int &pol,
			      double &w, randmt_t *rng)
{

  int iE, ix, iy;
  double tx, ty, t;
  if (!Intersect(x0, u, ix, iy, tx, ty, t)) {
    w = 0;
    return false;
  }
  double dO = dOmega(u*t);
  if (dO < 1e-10) {
    cout << "t: " << t << endl;
    cout << "u: " << u << endl;
    cout << PixelSurf << endl;
  }

  double w0;
  if (LoopFlag==1) { // loop mode
    iE = EnergyIdx;
    pol = PolIdx;
    w0 = Image[N*NBins*pol + NX*NBins*iy + NBins*ix + iE]/dOmega(u*t);
  }
  else {
    w0 = SumEnergyImage[NX*iy+ix]/dOmega(u*t);
    int i0 = NBins*(NX*iy + ix);
    double *cumul_func = &CumulEnergy[2*i0];
    double R = Rnd_r(rng);
    int i;
    Locate(R, cumul_func, 2*NBins, &i);
    iE = i % NBins;
    pol = i/NBins;
    if (pol>1)
      throw xrmc_exception("Error in cumulative energy array dimensions.\n");
  }
  double rE = Rnd_r(rng) - 0.5;
  double rx = tx - 0.5;
  double ry = ty - 0.5;
  w = w0*InterpolWeight(iE, ix, iy, pol, rE, rx, ry, 0);
  // central bin energy
  E = Emin + (Emax - Emin)*(0.5 + iE + rE)/NBins;

  return true;
}


double beamscreen::InterpolWeight(int iE, int ix, int iy, int pol,
        double rE, double rx, double ry, int ene_flag)
{
  const double rmax = 100;
  int iE1, ix1, iy1, ipol;
  int i000, i001, i010, i011, i100, i101, i110, i111;
  double c000, c001, c010, c011, c100, c101, c110, c111;
  double w, tE, tx, ty, uE, ux, uy;

  if (InterpolFlag==0) return 1;

  if (rE>=0) {
    tE = rE;
    iE1 = iE+1;
    if (iE1>=NBins) iE1=NBins-1;
  }
  else {
    tE = -rE;
    iE1 = iE-1;
    if (iE1<0) iE1=0;
  }
  uE = 1. - tE;
    
  if (rx>=0) {
    tx = rx;
    ix1 = ix+1;
    if (ix1>=NX) ix1=NX-1;
  }
  else {
    tx = -rx;
    ix1 = ix-1;
    if (ix1<0) ix1=0;
  }
  ux = 1. - tx;
    
  if (ry>=0) {
    ty = ry;
    iy1 = iy+1;
    if (iy1>=NY) iy1=NY-1;
  }
  else {
    ty = -ry;
    iy1 = iy-1;
    if (iy1<0) iy1=0;
  }
  uy = 1. - ty;
  ipol = pol*NY*NX*NBins;
  i000 = ipol + NX*NBins*iy  + NBins*ix   + iE;
  i001 = ipol + NX*NBins*iy  + NBins*ix   + iE1;
  i010 = ipol + NX*NBins*iy  + NBins*ix1  + iE;
  i011 = ipol + NX*NBins*iy  + NBins*ix1  + iE1;
  i100 = ipol + NX*NBins*iy1 + NBins*ix   + iE;
  i101 = ipol + NX*NBins*iy1 + NBins*ix   + iE1;
  i110 = ipol + NX*NBins*iy1 + NBins*ix1  + iE;
  i111 = ipol + NX*NBins*iy1 + NBins*ix1  + iE1;
    
  c000 = Image[i000];
  c001 = Image[i001];
  c010 = Image[i010];
  c011 = Image[i011];
  c100 = Image[i100];
  c101 = Image[i101];
  c110 = Image[i110];
  c111 = Image[i111];

  double norm;
  if (ene_flag==1) norm = c000*uy*ux + c010*uy*tx + c100*ty*ux + c110*ty*tx;
  else norm = c000; 
  if (c000>rmax*norm) c000 = rmax;
  else c000 /= norm;
  if (c001>rmax*norm) c001 = rmax;
  else c001 /= norm;
  if (c010>rmax*norm) c010 = rmax;
  else c010 /= norm;
  if (c011>rmax*norm) c011 = rmax;
  else c011 /= norm;
  if (c100>rmax*norm) c100 = rmax;
  else c100 /= norm;
  if (c101>rmax*norm) c101 = rmax;
  else c101 /= norm;
  if (c110>rmax*norm) c110 = rmax;
  else c110 /= norm;
  if (c111>rmax*norm) c111 = rmax;
  else c111 /= norm;
    
  w = c000*uy*ux*uE + c001*uy*ux*tE + c010*uy*tx*uE + c011*uy*tx*tE
    + c100*ty*ux*uE + c101*ty*ux*tE + c110*ty*tx*uE + c111*ty*tx*tE;
  //w = c000*ux + c010*tx;
  /*  if (iy==NY/2 && ix==2 && rx>0) {
    cout << "ix0 " << ix << endl;
    cout << "ix1 " << ix1 << endl;
    cout << "c0 " << c000 << endl;
    cout << "c1 " << c010 << endl;
    cout << "ux " << ux << endl;
    cout << "tx " << tx << endl;
    cout << "w " << w << endl;
    cout << "im0 " << Image[NBins*(NX*iy+ix)+iE] << endl;
    cout << "im1 " << Image[NBins*(NX*iy+ix1)+iE] << endl;

    }*/

  return w;
}

// initialize loop on events
int beamscreen::Begin()
{
  // check if flag for loop on all lines and on all intervals is enabled
  if (LoopFlag==0) LoopIdx = 0; // if not, set loop index to zero
  else { // if yes, initialize all nested loop indexes
    PolIdx = PhotonIdx = EnergyIdx = 0;
  }

  return 0;
}

// next step of the event loop
int beamscreen::Next()
{
  // check if flag for loop on all intervals is enabled
  if (LoopFlag==0) LoopIdx = 1;  // if not, set loop index to 1
  else { // if yes, update all nested loop indexes
    PhotonIdx++; // increase index of events within the interval 
    if (PhotonIdx == PhotonNum) { // check for last event
      PhotonIdx = 0; // if yes, reset index to zero ...
      EnergyIdx++; // go to the next energy bin
      if (EnergyIdx == NBins) { // check if the last bin is reached
	  // if yes reinitialize loop on intervals ...
	PhotonIdx = EnergyIdx = 0;
	PolIdx++; // with the other polarization type
      }
    }
  }

  return 0;
}

// check if the end of the event loop is reached
bool beamscreen::End()
{
  if (LoopFlag==0) return (LoopIdx==1);
  else return (PolIdx > 1); // end when all intervals have been
                            // completed with both polarization types 
}

// event multiplicity
long long beamscreen::EventMulti()
{
  if (LoopFlag==0) return 1;
  else {
    return 2*(NBins*PhotonNum);
  }
}

// Generates a random pixel x, y for a given polarization and enrgy bin
vect3 beamscreen::RandomXY(double &E, double &w, randmt_t *rng)
{
  int i0 = N*(NBins*PolIdx + EnergyIdx);
  double *cumul_func = &CumulXY[i0];
  double R = Rnd_r(rng);
  int i;
  Locate(R, cumul_func, N, &i);
  int ix = i % NX;
  int iy = i/NX;
  if (iy>NY)
    throw xrmc_exception("Error in cumulative x y array dimensions.\n");

  double w0 = BinWeight[NBins*PolIdx+EnergyIdx];

  double rE = Rnd_r(rng) - 0.5;
  double rx = Rnd_r(rng) - 0.5;
  double ry = Rnd_r(rng) - 0.5;

  // evaluates interpolation weight factor
  w = w0*InterpolWeight(EnergyIdx, ix, iy, PolIdx, rE, rx, ry, 0);
  // central bin energy
  E = Emin + (Emax - Emin)*(0.5 + EnergyIdx + rE)/NBins;

  // local x, y coordinates of the pixel
  double x = PixelSizeX*(-0.5*NX + 0.5 + ix + rx);
  double y = PixelSizeY*(-0.5*NY + 0.5 + iy + ry);
  vect3 r = X + ui*x + uj*y; // 3d absolute coordinates of the pixel

  return r;
}


beamscreen *beamscreen::Clone(string dev_name) {
  //cout << "Entering beamscreen::Clone\n";
  beamscreen *clone = new beamscreen(dev_name);

  /*
  clone->PhotonNum = PhotonNum;
  clone->PhotonIdx = PhotonIdx;
  clone->PolIdx = PolIdx;
  clone->EnergyIdx = EnergyIdx;
  clone->NX = NX;
  clone->NY = NY;
  clone->N = N;
  clone->NBins = NBins;
  */

  /* removed to avoid excessive ram usage
  if (Image!=NULL) {
    clone->Image = new double[2*N*NBins];
    memcpy(clone->Image, Image, sizeof(double)*2*N*NBins);
  }
  if (CumulImage!=NULL) {
    clone->CumulImage = new double[2*N*NBins];
    memcpy(clone->CumulImage, CumulImage, sizeof(double)*2*N*NBins);
  }
  if (CumulEnergy!=NULL) {
    clone->CumulEnergy = new double[2*N*NBins];
    memcpy(clone->CumulEnergy, CumulEnergy, sizeof(double)*2*N*NBins);
  }
  if (SumEnergyImage!=NULL) {
    clone->SumEnergyImage = new double[N];
    memcpy(clone->SumEnergyImage, SumEnergyImage, sizeof(double)*N);
  }
  if (BinWeight!=NULL) {
    clone->BinWeight = new double[2*NBins];
    memcpy(clone->BinWeight, BinWeight, sizeof(double)*2*NBins);
  }
  if (CumulXY!=NULL) {
    clone->CumulXY = new double[2*N*NBins];
    memcpy(clone->CumulXY, CumulXY, sizeof(double)*2*N*NBins);
  }
  if (PixelX!=NULL) {
    clone->PixelX = new vect3[N];
    memcpy(clone->PixelX, PixelX, sizeof(vect3)*N);
  }
  */

  /*
  clone->Image = Image;
  clone->CumulImage = CumulImage;
  clone->CumulEnergy = CumulEnergy;
  clone->SumEnergyImage = SumEnergyImage;
  clone->BinWeight = BinWeight;
  clone->CumulXY = CumulXY;
  clone->PixelX = PixelX;

  clone->TotalIntensity = TotalIntensity;
  clone->PixelSizeX = PixelSizeX;
  clone->PixelSizeY = PixelSizeY;
  clone->PixelSurf = PixelSurf 
  clone->dOmegaLim = dOmegaLim;
  clone->Shape = Shape;
  clone->RandomPixelFlag = RandomPixelFlag;
  clone->RunningFasterFlag = RunningFasterFlag;

  clone->PolarizedFlag = PolarizedFlag;
  clone->LoopFlag = LoopFlag;
  clone->EnergyBinFlag = EnergyBinFlag;
  clone->InterpolFlag = InterpolFlag;
  clone->Emin = Emin;
  clone->Emax = Emax;
  clone->X = X;
  clone->ui = ui;
  clone->uj = uj;
  clone->uk = uk;
  */

  *clone = *this;

  return clone;
}

