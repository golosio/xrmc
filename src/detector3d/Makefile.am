xrmcincludedir = ${includedir}/xrmc
xrmcinclude_HEADERS = xrmc_detector3d.h xrmc_loaddetector3d.h

noinst_LTLIBRARIES = libdetector3d.la

libdetector3d_la_SOURCES = detector3d.cpp loaddetector3d.cpp xrmc_detector3d.h xrmc_loaddetector3d.h
libdetector3d_la_CPPFLAGS = -I${srcdir}/../xrmc_math -I${srcdir}/../xrmc_algo -I${srcdir}/../arrayNd -I${srcdir}/../device -I${srcdir}/../source -I${srcdir}/../photon -I${srcdir}/../spectrum -I${srcdir}/../exception -I${srcdir}/../main -I${srcdir}/../gettoken -I${srcdir}/../randmt -I${srcdir}/../sample -I${srcdir}/../composition -I${srcdir}/../geom3d  -I${srcdir}/../phasecontrast @OPENMP_CXXFLAGS@
