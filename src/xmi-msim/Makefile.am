noinst_LTLIBRARIES = libxmimsim.la


libxmimsim_la_SOURCES = xmimsim-detectorconvolute.cpp xrmc_detectorconvolute.h xmimsim-loaddetectorconvolute.cpp 
libxmimsim_la_CPPFLAGS = -I${srcdir}/../device \
		@xmi_msim_CFLAGS@ \
		@gmodule_CFLAGS@ \
		-I${srcdir}/../xrmc_math \
		-I${srcdir}/../xrmc_algo \
		-I${srcdir}/../arrayNd \
		-I${srcdir}/../device \
		-I${srcdir}/../source \
		-I${srcdir}/../photon \
		-I${srcdir}/../spectrum \
		-I${srcdir}/../exception \
		-I${srcdir}/../main \
		-I${srcdir}/../gettoken \
		-I${srcdir}/../randmt \
		-I${srcdir}/../detector \
		-I${srcdir}/../composition \
		-DXRMC_XMIMSIM_LIB=\"$(libdir)/xrmc\"


plugindir = $(libdir)/xrmc
plugin_LTLIBRARIES = xrmc-xmimsim.la
xrmc_xmimsim_la_CFLAGS = @xmi_msim_CFLAGS@ @gmodule_CFLAGS@ 
xrmc_xmimsim_la_LIBADD = @xmi_msim_LIBS@
xrmc_xmimsim_la_SOURCES = xrmc-xmimsim.c
xrmc_xmimsim_la_LDFLAGS = -module -avoid-version -shared
