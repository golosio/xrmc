#include <stdio.h>
#include <stdlib.h>
#include "detector_routines.h"

int main(int argc, char *argv[])
{
  double *spectrum_eff, Emin=0, Emax=0;
  int spectrum_size=0;
  
  if (argc != 6) {
    printf("Usage: detector_response input-file output-file Emin Emax "
	   "spectrum-size\n");
    exit(0);
  }
  sscanf(argv[3], "%lf", &Emin);
  if (Emin<0) {
    printf("Error: Emin must be >= 0\n");
    exit(EXIT_FAILURE);
  }
  sscanf(argv[4], "%lf", &Emax);
  if (Emax<=Emin) {
    printf("Error: Emax must be > Emin\n");
    exit(EXIT_FAILURE);
  }
  sscanf(argv[5], "%d", &spectrum_size);
  if (spectrum_size<=0) {
    printf("Error: spectrum-size must be an integer > 0\n");
    exit(EXIT_FAILURE);
  }
  
  spectrum_eff=detector_efficiency(argv[1], 14, 0.05, 4, 0.00125, Emin, Emax,
				   spectrum_size);
  
  detector_convolution_SDD(spectrum_eff, argv[2], Emin, Emax, spectrum_size);

  return 0;
}
