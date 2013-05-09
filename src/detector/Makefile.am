xrmcincludedir = ${includedir}/xrmc
xrmcinclude_HEADERS = xrmc_detector.h xrmc_loaddetector.h xrmc_screen.h

noinst_LTLIBRARIES = libdetector.la

libdetector_la_SOURCES = detector.cpp loaddetector.cpp xrmc_screen.cpp xrmc_detector.h xrmc_loaddetector.h xrmc_screen.h
libdetector_la_CPPFLAGS = -I${srcdir}/../xrmc_math -I${srcdir}/../xrmc_algo -I${srcdir}/../arrayNd -I${srcdir}/../device -I${srcdir}/../source -I${srcdir}/../photon -I${srcdir}/../spectrum -I${srcdir}/../exception -I${srcdir}/../main -I${srcdir}/../gettoken -I${srcdir}/../randmt -I${srcdir}/../sample -I${srcdir}/../composition -I${srcdir}/../geom3d  @OPENMP_CXXFLAGS@
