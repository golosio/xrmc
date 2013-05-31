#ifndef FFTH
#define FFTH

#ifdef __cplusplus
extern "C" {
#endif

    void cdft2d(int, int, int, double **, double *, int *, double *);
    void rdft2d(int, int, int, double **, double *, int *, double *);
    void ddct2d(int, int, int, double **, double *, int *, double *);
    void ddst2d(int, int, int, double **, double *, int *, double *);

#ifdef __cplusplus
}
#endif

#endif
