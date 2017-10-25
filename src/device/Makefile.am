xrmcincludedir = ${includedir}/xrmc
xrmcinclude_HEADERS = xrmc_device.h

noinst_LTLIBRARIES = libdevice.la

libdevice_la_SOURCES = device.cpp loadnewdevice.cpp xrmc_device.h
libdevice_la_CPPFLAGS = -I${srcdir}/../main -I${srcdir}/../xrmc_math -I${srcdir}/../exception -I${srcdir}/../gettoken -I${srcdir}/../sample -I${srcdir}/../composition -I${srcdir}/../geom3d -I${srcdir}/../detector -I${srcdir}/../detector3d -I${srcdir}/../spectrum -I${srcdir}/../source -I${srcdir}/../randmt -I${srcdir}/../photon  -I${srcdir}/../arrayNd -I${srcdir}/../beamsource -I${srcdir}/../phasecontrast -I${srcdir}/../radionuclide @xraylib_CFLAGS@
if ENABLE_XMIMSIM
libdevice_la_CPPFLAGS += -I${srcdir}/../xmi-msim @xmi_msim_CFLAGS@ @gmodule_CFLAGS@ -DXRMC_XMIMSIM_LIB=\"$(libdir)/xrmc\"
endif

