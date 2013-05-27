xrmcincludedir = ${includedir}/xrmc
xrmcinclude_HEADERS = phcdetector.h phcdevice.h

noinst_LTLIBRARIES = libphcdetector.la

libphcdetector_la_SOURCES = phcdetector.cpp loadphcdetector.cpp phcdetector.h
libphcdetector_la_CPPFLAGS = -I${srcdir}/../xrmc_math -I${srcdir}/../xrmc_algo -I${srcdir}/../arrayNd -I${srcdir}/../device -I${srcdir}/../source -I${srcdir}/../photon -I${srcdir}/../spectrum -I${srcdir}/../exception -I${srcdir}/../main -I${srcdir}/../gettoken -I${srcdir}/../randmt -I${srcdir}/../sample -I${srcdir}/../composition -I${srcdir}/../geom3d  -I${srcdir}/../detector -I${srcdir}/../fft @OPENMP_CXXFLAGS@
