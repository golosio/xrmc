#include <stdio.h>
#include <math.h>
#include <complex>
#include "nr.h"
#include "xrmcglob.h"
#include "gridintegral3d.h"
#include "geom.h"
#include "mcxmath.h"
#include <time.h>
#include "sample.h"
#include "spectrum.h"
#include "photon.h"
#include "initphoton.h"
#include "signal.h"
#include "xraylib.h"
#include "phase.h"
using namespace std;

double Step[MAXNSTEPS], mu_step[MAXNSTEPS], sum_mu_s[MAXNSTEPS];
double CT[MAXNSTEPS];

void InitGrid()
{
  int i;

  XGrid = -(double)Nx*VoxelSize/2.0;
  YGrid = -(double)Ny*VoxelSize/2.0;
  ZGrid = -(double)Nz*VoxelSize/2.0;
  for (i=0; i<Nx*Ny*Nz; i++) {
    CTGrid[i] = AGrid*CTGrid[i] + BGrid;
  }
}

void CTMC()
{
  int PhotonIdx, PixelIdx, ScattOrderIdx, is, bin;
  double point_x[3], point_r[3]; //, point_x1[3], point_r1[3];
  double Pgeom, Psurv, integr;
  time_t time0, time1;
  double time_diff;
  clock_t cpu_time;
  double cpu_time_diff;
  double E0;
  int polar;
  photon Photon;
  int interaction_type, Z;

  int i1, ix1, iy1;
  double Lambda0, K0;
  double deltaL, muL;
  double x1, y1, x1Csi, ratio;
  complex<double> transm, propagator, c1, c2;

  printf("Start simulation...\n");
  time0 = time(NULL);
  cpu_time = clock();

  for (PhotonIdx=0; PhotonIdx<PhotonNum; PhotonIdx++) {

    // if ((PhotonIdx%10000)==0) printf("%d\n", PhotonIdx);
    for (PixelIdx=0; PixelIdx<PixelN; PixelIdx++) {
      Spectr1->RandomEnergy(&E0, &polar);
      Sample1->Comp->Mu(E0);
      RandomPointOnPixel(PixelIdx, point_x);
      // printf("pixel %d\t", PixelIdx);
      // printf("pt_x[] : %f\t%f\t%f\n", point_x[0], point_x[1], point_x[2]); 
      VectComb(1., point_x, -1., SourceX, point_r);
      Pgeom = dOmega(point_r);

      VectNormalize(point_r);
      // printf("pt_r[] : %f\t%f\t%f\n", point_r[0], point_r[1], point_r[2]); 
      
      integr = PathIntegral(SourceX, point_r);      
      Psurv = exp(-integr);
      // printf("POmega:%.3e\t Pgeom:%.3e\t Psurv:%.3e\n",  POmega(point_r),Pgeom,Psurv); 
      
      if (PixelEmax >=0 && PixelNbins > 1) 
	{                      
	  if (E0 >= PixelEmax )
	    {
	      bin = PixelNbins - 1;
	      printf("WARNING: energy greater than Emax, photon put in the last bin.\n");
	    }
	  else if (E0 < PixelEmin )
	    {
	      bin = 0;
	      printf("WARNING: energy smaller than Emin, photon put in the first bin.\n");
	    }
	  else
	    bin = ( (int) trunc((E0-PixelEmin)/(PixelEmax-PixelEmin)*PixelNbins)); //Pier
	}
      else                                                                         //Pier
	bin = 0;



      if (PixelType == 1 || PixelType == 3)                                     //Pier: va bene a Bruno?
	Image[PixelIdx][bin] += POmega(point_r)*Pgeom*Psurv*E0;                 //Pier
      else                                                                      //Pier
	Image[PixelIdx][bin] += POmega(point_r)*Pgeom*Psurv;                    //Pier
    }
    printf("%d / %d\n", PhotonIdx+1, PhotonNum);
  }

  for (ScattOrderIdx=1; ScattOrderIdx<=MaxScattOrder; ScattOrderIdx++) {
    for (PhotonIdx=0; PhotonIdx<ScattPhotonNum; PhotonIdx++) {
      
      for (PixelIdx=0; PixelIdx<PixelN; PixelIdx++) {
	Spectr1->RandomEnergy(&E0, &polar);
	Sample1->Comp->Mu(E0);
	InitPhoton(&Photon, E0, polar);
	// Sample1->RandomPoint(point_x);
	// printf("%f\t%f\t%f\n", Photon.uk[0], Photon.uk[1], Photon.uk[2]);
	for (is=1; is<=ScattOrderIdx; is++) {
	  Photon.MonteCarloStep(Sample1, &Z, &interaction_type);
	  if (Photon.w == 0) break;
	  if (is<ScattOrderIdx) {
	    Photon.Scatter(Z, interaction_type);
	    // TEST // if (is==1) printf("PhotonData %f\t%f\t%f\t%f\t%f\n",
	    // TEST // Photon.w, Photon.E, Photon.uk[0], Photon.uk[1],
            // TEST // Photon.uk[2]);
	    if (interaction_type==FLUORESCENCE || interaction_type==INCOHERENT)
	      Sample1->Comp->Mu(Photon.E);
	  }   
	}
	// printf("it1 %d\n", interaction_type);
	if (Photon.w != 0) ScattSignal(ScattOrderIdx, PixelIdx, &Photon, Z, 
				       interaction_type);
	// VectComb(1., point_x, -1., SourceX, point_r);
	
	// VectNormalize(point_r);
	// printf("pt_r[] : %f\t%f\t%f\n", point_r[0], point_r[1], point_r[2]); 
	
	// integr = PathIntegral(SourceX, point_r);
	// Psurv = exp(-integr);
	
	// RandomPointOnPixel(PixelIdx, point_x1);
	// printf("pixel %d\t", PixelIdx);
	// printf("pt_x[] : %f\t%f\t%f\n", point_x[0], point_x[1], point_x[2]); 
	// VectComb(1., point_x1, -1., point_x, point_r1);
	// Pgeom = dOmega(point_r1);
	
	// VectNormalize(point_r1);
	// printf("pt_r[] : %f\t%f\t%f\n", point_r[0], point_r[1], point_r[2]); 
	
	// integr = PathIntegral(point_x, point_r1);
	// Psurv = Psurv*exp(-integr);
	// ScattImage[PixelIdx] += POmega(point_r)*Pgeom*Psurv;
      }
    }
    printf("%d / %d\n", PhotonIdx+1, ScattPhotonNum);
  }


  for (PhotonIdx=0; PhotonIdx<PhasePhotonNum; PhotonIdx++) {    
    Spectr1->RandomEnergy(&E0, &polar);
    Lambda0 = 1.e-8*KEV2ANGST/E0; // lambda in cm
    //printf("Lambda: %.20e\n", Lambda0);
    K0 = 2.*PI/Lambda0;
    Sample1->Comp->Mu(E0);
    Sample1->Comp->Delta(E0);

    i1 = 0;
    for (ix1=0; ix1<Nx1; ix1++) {
      for (iy1=0; iy1<Ny1; iy1++) {
	if (ix1<PixelNX && iy1<PixelNY) {
	  PixelIdx = PixelNY*ix1 + iy1;
	  RandomPointOnPixel(PixelIdx, point_x);
	  //printf("\npoint_x:\t%f\t%f\t%f\n", point_x[0], point_x[1],
	  //point_x0[2]); 
	  VectComb(1., point_x, -1., SourceX, point_r);
	  VectNormalize(point_r);      
	  Sample1->LinearMuDelta(SourceX, point_r);
	  //printf("NSteps: %d\n", Sample1->Path->NSteps);
	  
	  muL = Sample1->Path->MuL;
	  deltaL = Sample1->Path->DeltaL;
	  transm = exp(-I*K0*deltaL)*exp(-muL/2);
	  Transm[i1] = transm.real();
	  Transm[i1+1] = transm.imag();
	}
	else {
	  Transm[i1] = 1;
	  Transm[i1+1] = 0;
	}
	if (ix1 < Nx1/2) x1 = dx1*ix1;
	else x1 = dx1*(ix1 - Nx1);
	if (iy1 < Ny1/2) y1 = dy1*iy1;
	else y1 = dy1*(iy1 - Ny1);
	x1Csi = x1*Csix + y1*Csiy;
	propagator = exp(I*K0/(2*Reff)*(x1*x1+y1*y1 - x1Csi*x1Csi));
	/* */
	if (fabs(x1) > Lambda0*Lx1) {
	  ratio = (fabs(x1) - Lambda0*Lx1)/Sigmax1/Lambda0;
	  propagator *= exp(-ratio*ratio/2);
	}
	if (fabs(y1) > Lambda0*Ly1) {
	  ratio = (fabs(y1) - Lambda0*Ly1)/Sigmay1/Lambda0;
	  propagator *= exp(-ratio*ratio/2);
	}
	/* */       
	Propagator[i1] = propagator.real();
	Propagator[i1+1] = propagator.imag();
	i1 += 2;
      }
    }
    
    fourn(Transm-1, NN1-1, 2, 1);
    fourn(Propagator-1, NN1-1, 2, 1);
    i1 = 0;
    for (ix1=0; ix1<Nx1; ix1++) {
      for (iy1=0; iy1<Ny1; iy1++) {
	c1 = complex<double>(Transm[i1],Transm[i1+1]);
	c2 = complex<double>(Propagator[i1],Propagator[i1+1]);
	c1 *= c2;
	// c1 = c2;
	Transm[i1] = dx1*dy1*c1.real()/Nx1/Ny1;
	Transm[i1+1] = dx1*dy1*c1.imag()/Nx1/Ny1;
	i1 += 2;
      }
    }
    fourn(Transm-1, NN1-1, 2, -1);

    PixelIdx = 0;
    for (ix1=0; ix1<PixelNX; ix1++) {
      for (iy1=0; iy1<PixelNY; iy1++) {
	i1 = 2*(Nx1*ix1 + iy1);
	// transm = complex<double>(Propagator[i1],Propagator[i1+1]);
	transm = complex<double>(Transm[i1],Transm[i1+1]);
	transm *= Csiz / (I*Lambda0*R1*R12);
      
	if (PixelEmax >=0 && PixelNbins > 1) {                      
	  if (E0 >= PixelEmax ) {
	    bin = PixelNbins - 1;
	    printf("WARNING: energy greater than Emax, photon put in the last bin.\n");
	  }
	  else if (E0 < PixelEmin ) {
	    bin = 0;
	    printf("WARNING: energy smaller than Emin, photon put in the first bin.\n");
	  }
	  else
	    bin = ( (int) trunc((E0-PixelEmin)/(PixelEmax-PixelEmin)*PixelNbins)); //Pier
	}
	else     //Pier
	  bin = 0;
	if (PixelType == 1 || PixelType == 3) //Pier: va bene a Bruno?
	  PhaseImage[PixelIdx][bin] += norm(transm)*E0;                //Pier
	else     //Pier
	  PhaseImage[PixelIdx][bin] += norm(transm);                   //Pier
	//printf("norm(F): %e\n", norm(transm));
	//printf("Pixel: %d / %d\n", PixelIdx+1, PixelN);
	PixelIdx++;
      }
    }
    printf("%d / %d\n", PhotonIdx+1, PhasePhotonNum);
  }

  //PhotonCounts();
  cpu_time_diff = (double)(clock() - cpu_time)/CLOCKS_PER_SEC;
  time1 = time(NULL);
  time_diff = difftime(time1, time0);
  printf("Elapsed time : %.0f sec\n", time_diff);
  printf("Elapsed CPU time : %.2f sec\n", cpu_time_diff);
}

double PathIntegral(double x[], double v_k[])
{

  Sample1->LinearAbsorption(x, v_k);

  return Sample1->Path->MuL;
}

double dOmega(double *point_r)
{
  double r;

  r = VectMod(point_r);
  return -ScalProd(point_r, Screen_k) * PixelSurf /r/r/r;
}

void PhotonCounts()
{
  double val; //, sigma;
  int PixelIdx, ScattOrderIdx;
  int PixelBin;

  for (PixelIdx=0; PixelIdx<PixelN; PixelIdx++) {
    for (PixelBin=0; PixelBin<PixelNbins; PixelBin++) {  //Pier
      val = Image[PixelIdx][PixelBin];                   //Pier
      val *= Spectr1->TotalIntensity*ExpTime/PhotonNum;
      // sigma = sqrt(val);
      // val += sigma*gasdev(&idum);
      // val = round(val);
      Image[PixelIdx][PixelBin] = val>0 ? val : 0;       //Pier
      val = PhaseImage[PixelIdx][PixelBin];                   //Pier
      val *= Spectr1->TotalIntensity*ExpTime/PhasePhotonNum;
      // sigma = sqrt(val);
      // val += sigma*gasdev(&idum);
      // val = round(val);
      PhaseImage[PixelIdx][PixelBin] = val>0 ? val : 0;       //Pier
    }
    for (ScattOrderIdx=1; ScattOrderIdx<=MaxScattOrder; ScattOrderIdx++) 
      for (PixelBin=0; PixelBin<PixelNbins; PixelBin++) {       //Pier
	val = ScattImage[ScattOrderIdx-1][PixelIdx][PixelBin];  //Pier
	val *= Spectr1->TotalIntensity*ExpTime/ScattPhotonNum;
	
	// sigma = sqrt(val);
	// val += sigma*gasdev(&idum);
	// val = round(val);
	ScattImage[ScattOrderIdx-1][PixelIdx][PixelBin] = val>0 ? val : 0;     //Pier
	//printf("PIER\t%d\t%d\t%d\t%f\n", ScattOrderIdx, PixelIdx, PixelBin,  ScattImage[ScattOrderIdx-1][PixelIdx][PixelBin]);
      }
  }
}

double POmega(double *point_r)
{
  double x, y, z, r;
  double cos_th, cos_th_l;
  
  r = sqrt(ScalProd(point_r, point_r));
  if (r==0) return 0;
  x = ScalProd(point_r, Source_i);
  y = ScalProd(point_r, Source_j);
  z = ScalProd(point_r, Source_k);
  cos_th_l = CosThLxy(x, y);
  cos_th = z / r;
  if (cos_th < cos_th_l) return 0;
  else return 1./Omega;
}








