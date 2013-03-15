bin_PROGRAMS = xrmc

xrmc_SOURCES = loaddevice.cpp loadparams.cpp xrmc.cpp       xrmc.h         xrmc_main.cpp
xrmc_LDADD = ../libxrmc.la
xrmc_CPPFLAGS = -I${srcdir}/../xrmc_math -I${srcdir}/../xrmc_algo -I${srcdir}/../arrayNd -I${srcdir}/../device -I${srcdir}/../source -I${srcdir}/../photon -I${srcdir}/../spectrum -I${srcdir}/../exception -I${srcdir}/../main -I${srcdir}/../gettoken @xraylib_CFLAGS@ -I${srcdir}/../sample -I${srcdir}/../composition -I${srcdir}/../geom3d -I${srcdir}/../detector -I${srcdir}/../randmt 
if ENABLE_XMIMSIM
xrmc_CPPFLAGS += -I${srcdir}/../xmi-msim @xmi_msim_CFLAGS@ @gmodule_CFLAGS@ \
		-DXRMC_XMIMSIM_LIB=\"$(libdir)/xrmc\"
endif
