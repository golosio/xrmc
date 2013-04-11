pkgconfigdir=$(libdir)/pkgconfig
pkgconfig_DATA=libxrmc.pc

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src doc examples

EXTRA_DIST = test libxrmc.pc.in xrmc.spec.in
DISTCHECK_CONFIGURE_FLAGS = --disable-openmp --disable-xmi-msim
