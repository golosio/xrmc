#ifndef FFTH
#define FFTH

#ifdef __cplusplus
extern "C" {
#endif

  void cdft(int, int, double *, int *, double *);
  void rdft(int, int, double *, int *, double *);
  void ddct(int, int, double *, int *, double *);
  void ddst(int, int, double *, int *, double *);
  void dfct(int, double *, double *, int *, double *);
  void dfst(int, double *, double *, int *, double *);
  void cdft2d(int, int, int, double **, double *, int *, double *);
  void rdft2d(int, int, int, double **, double *, int *, double *);
  void ddct2d(int, int, int, double **, double *, int *, double *);
  void ddst2d(int, int, int, double **, double *, int *, double *);

#ifdef __cplusplus
}
#endif

#endif
