xrmcincludedir = ${includedir}/xrmc
xrmcinclude_HEADERS = beamsource.h beamscreen.h anisotropicsource.h intensityscreen.h

noinst_LTLIBRARIES = libbeamsource.la

libbeamsource_la_SOURCES = loadbeamsource.cpp loadbeamscreen.cpp loadintensityscreen.cpp loadanisotropicsource.cpp beamsource.cpp beamscreen.cpp anisotropicsource.cpp intensityscreen.cpp beamsource.h beamscreen.h anisotropicsource.h intensityscreen.h    
libbeamsource_la_CPPFLAGS = -I${srcdir}/../xrmc_math -I${srcdir}/../xrmc_algo -I${srcdir}/../device -I${srcdir}/../source -I${srcdir}/../spectrum -I${srcdir}/../photon -I${srcdir}/../exception -I${srcdir}/../main -I${srcdir}/../gettoken -I${srcdir}/../detector @xraylib_CFLAGS@ -I${srcdir}/../randmt

