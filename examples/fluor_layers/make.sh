IXRL=$(pkg-config --cflags libxrl)
LXRL=$(pkg-config --libs libxrl)
gcc -Wall -O3 -I ../include $IXRL $LXRL -o detector_response detector_response.c detector_routines.c -lm -lxrl
