AC_INIT([xrmc],[6.3.3],[golosio@uniss.it])
AC_PREREQ([2.60])
AC_CONFIG_AUX_DIR([build-aux])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AS_INIT


LIB_CURRENT=0
LIB_REVISION=0
LIB_AGE=0
AC_SUBST(LIB_CURRENT)
AC_SUBST(LIB_REVISION)
AC_SUBST(LIB_AGE)

#m4_pattern_allow([AS_TR_SH])
AC_CONFIG_MACRO_DIR([m4])

AC_USE_SYSTEM_EXTENSIONS
AM_PROG_AR

#libtool
LT_PREREQ([2.0.0])

LT_INIT()
AC_PROG_LIBTOOL

AC_PROG_CC
AM_PROG_CC_C_O

#this next line may never be reached...
if test -z $CC ; then
        AC_MSG_ERROR([no C compiler was found on the system.])
fi

AC_PROG_CXX
AC_PROG_CXX_C_O

#this next line may never be reached...
if test -z $CXX ; then
        AC_MSG_ERROR([no C++ compiler was found on the system.])
fi

#look for xraylib
#initialize pkg-config
PKG_PROG_PKG_CONFIG

#search for xraylib and other modules
PKG_CHECK_MODULES([xraylib],libxrl >= 2.15.0)

AC_CONFIG_FILES([	Makefile 
			src/Makefile
			src/arrayNd/Makefile
			src/composition/Makefile
			src/detector/Makefile
			src/device/Makefile 
			src/geom3d/Makefile
			src/gettoken/Makefile
			src/xrmc_algo/Makefile
			src/photon/Makefile
			src/source/Makefile
			src/spectrum/Makefile 
			src/xrmc_math/Makefile
			src/exception/Makefile
			src/sample/Makefile
			src/main/Makefile
			src/randmt/Makefile
			doc/Makefile
			examples/Makefile
			examples/cylind_shell/Makefile
			examples/materials/Makefile
			examples/quadrics/Makefile
			examples/star1/Makefile
			examples/star2/Makefile
			examples/wheel/Makefile
			])
AC_CONFIG_HEADERS([config.h])


abs_top_builddir=`pwd -P`
AC_SUBST(abs_top_builddir)
abs_top_srcdir=`AS_DIRNAME([$0])`
cd $abs_top_srcdir
abs_top_srcdir=`pwd -P`
cd $abs_top_builddir
AC_SUBST(abs_top_srcdir)

AM_CONDITIONAL([ABS_SRC_BUILD_EQUAL],[test x$abs_top_builddir = x$abs_top_srcdir])


AC_OUTPUT

