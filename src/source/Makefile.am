xrmcincludedir = ${includedir}/xrmc
xrmcinclude_HEADERS = xrmc_source.h

noinst_LTLIBRARIES = libsource.la

libsource_la_SOURCES = loadsource.cpp source.cpp xrmc_source.h
libsource_la_CPPFLAGS = -I${srcdir}/../xrmc_math -I${srcdir}/../xrmc_algo -I${srcdir}/../arrayNd -I${srcdir}/../device -I${srcdir}/../source -I${srcdir}/../photon -I${srcdir}/../spectrum -I${srcdir}/../exception -I${srcdir}/../main -I${srcdir}/../gettoken @xraylib_CFLAGS@ -I${srcdir}/../sample -I${srcdir}/../composition -I${srcdir}/../geom3d -I${srcdir}/../randmt -I${srcdir}/../phasecontrast

if ENABLE_XMIMSIM
libsource_la_CPPFLAGS += -I${srcdir}/../xmi-msim @xmi_msim_CFLAGS@
endif

