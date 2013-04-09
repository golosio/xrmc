xrmcincludedir = ${includedir}/xrmc
xrmcinclude_HEADERS = xrmc_sample.h

noinst_LTLIBRARIES = libsample.la

libsample_la_SOURCES = loadsample.cpp sample.cpp xrmc_sample.h
libsample_la_CPPFLAGS = -I${srcdir}/../xrmc_math -I${srcdir}/../xrmc_algo -I${srcdir}/../arrayNd -I${srcdir}/../device -I${srcdir}/../source -I${srcdir}/../photon -I${srcdir}/../spectrum -I${srcdir}/../exception -I${srcdir}/../main -I${srcdir}/../gettoken @xraylib_CFLAGS@ -I${srcdir}/../sample -I${srcdir}/../composition -I${srcdir}/../geom3d -I${srcdir}/../randmt

