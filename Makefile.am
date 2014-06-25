pkgconfigdir=$(libdir)/pkgconfig
pkgconfig_DATA=libxrmc.pc

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src doc examples nsis

EXTRA_DIST = test libxrmc.pc.in xrmc.spec.in
#DISTCHECK_CONFIGURE_FLAGS = --disable-openmp --disable-xmi-msim
DISTCHECK_CONFIGURE_FLAGS = CXX=g++-mp-4.7 CC=gcc-mp-4.7
