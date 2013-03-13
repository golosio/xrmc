#include <stdio.h>

int main(int argc, char *argv[])
{
  FILE *fp;
  double x[3];

  fp = fopen(argv[1], "rb");
  fread(x, sizeof(double), 3, fp);
  printf("%le\n", x[2]);

  return 0;
}
