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
	beamsource

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
	xrmc_math/libxrmcmath.la \
	beamsource/libbeamsource.la

lib_LTLIBRARIES = libxrmc.la
libxrmc_la_SOURCES = 
libxrmc_la_LIBADD = $(SUBLIBS) @xraylib_LIBS@ -lm
libxrmc_la_LDFLAGS = -version-info @LIB_CURRENT@:@LIB_REVISION@:@LIB_AGE@ @OPENMP_CXXFLAGS@ @LDFLAGS_EXTRA@ 
libxrmc_la_LINK = $(LIBTOOL) --tag=CXX $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=link $(CCLD) $(AM_CFLAGS) $(CFLAGS) \
	$(libxrmc_la_LDFLAGS) $(LDFLAGS) -o $@

CCLD = $(CXX)

if ENABLE_XMIMSIM
SUBDIRS += xmi-msim
SUBLIBS += xmi-msim/libxrmc-xmimsim.la
libxrmc_la_LIBADD += @gmodule_LIBS@
endif

SUBDIRS += . main
