SUBDIRS = \
	arrayNd \
	composition \
	detector \
	device \
	exception \
	geom3d \
	gettoken \
	xrmc_algo \
	randmt \
	photon \
	sample \
	source \
	spectrum \
	xrmc_math \
	. \
	main

SUBLIBS = \
	arrayNd/libarrayNd.la \
	composition/libcomposition.la \
	detector/libdetector.la \
	device/libdevice.la \
	geom3d/libgeom3d.la \
	gettoken/libgettoken.la \
	xrmc_algo/libxrmcalgo.la \
	randmt/librandmt.la \
	photon/libphoton.la \
	sample/libsample.la \
	source/libsource.la \
	spectrum/libspectrum.la \
	xrmc_math/libxrmcmath.la

lib_LTLIBRARIES = libxrmc.la
libxrmc_la_SOURCES = 
libxrmc_la_LIBADD = $(SUBLIBS) @xraylib_LIBS@ -lm
libxrmc_la_LDFLAGS = -version-info @LIB_CURRENT@:@LIB_REVISION@:@LIB_AGE@ @OPENMP_CXXFLAGS@

