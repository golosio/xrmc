IXRL=$(pkg-config --cflags libxrl)
LXRL=$(pkg-config --libs libxrl)
gcc -Wall $IXRL $LXRL -o theoretical theoretical.c -lm -lxrl
