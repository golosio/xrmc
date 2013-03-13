//---------------------------------------------------------------------------

#ifndef Detector_routinesH
#define Detector_routinesH
//---------------------------------------------------------------------------

int detector_convolution_SDD(double *spectrum_eff, char *spectrum_out,
			     double Emin, double Emax, int spectrum_size);

double *detector_efficiency(char *spectrum_in, int Z,double thickness,
			    int Z_window,double window_thick,
			    double Emin, double Emax, int spectrum_size);

#endif


