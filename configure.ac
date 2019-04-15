AC_INIT([xrmc],[6.7.0],[golosio@uniss.it])
AC_PREREQ([2.60])
AC_CONFIG_SRCDIR([src/main/xrmc_main.cpp])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CANONICAL_HOST



LIB_CURRENT=0
LIB_REVISION=0
LIB_AGE=0
AC_SUBST(LIB_CURRENT)
AC_SUBST(LIB_REVISION)
AC_SUBST(LIB_AGE)

#m4_pattern_allow([AS_TR_SH])
AC_CONFIG_MACRO_DIR([m4])

AC_USE_SYSTEM_EXTENSIONS
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

#libtool
LT_PREREQ([2.0.0])

LT_INIT()
AC_PROG_LIBTOOL

AC_PROG_CC
if test `AS_BASENAME([$CC])` = $CC ; then
  AC_CHECK_PROG(CC_FULL, [$CC], $CC, [none])
  #this next line may never be reached...
  if test x$CC_FULL = "xnone" ; then
          AC_MSG_ERROR([no C compiler was found on the system.])
  fi
fi
AM_PROG_CC_C_O


AC_PROG_CXX
if test `AS_BASENAME([$CXX])` = $CXX ; then
  AC_CHECK_PROG(CXX_FULL, [$CXX], $CXX, [none])
  #this next line may never be reached...
  if test x$CXX_FULL = "xnone" ; then
          AC_MSG_ERROR([no C++ compiler was found on the system.])
  fi
fi
AC_PROG_CXX_C_O


LDFLAGS_EXTRA=""
case "$host" in
	*mingw*)
		LDFLAGS_EXTRA+="-no-undefined"
		;;
esac
AC_SUBST(LDFLAGS_EXTRA)

#look for openmp support in C++
AC_LANG_PUSH(C++)
AC_OPENMP
AC_LANG_POP(C++)
#look for openmp support in C
AC_LANG_PUSH(C)
AC_OPENMP
AC_LANG_POP(C)

AC_CHECK_LIB([m], [exp], [], [AC_MSG_ERROR([math library does not contain exp])])
AC_CHECK_FUNC(expm1,,[AC_MSG_ERROR([Required libc function not found])])
AC_CHECK_FUNC(log1p,,[AC_MSG_ERROR([Required libc function not found])])


#look for xraylib
#initialize pkg-config
PKG_PROG_PKG_CONFIG


#search for xraylib and other modules
PKG_CHECK_MODULES([xraylib],libxrl >= 3.99.0)


#determine if the xmi-msim should be built
AC_ARG_ENABLE([xmi-msim], AS_HELP_STRING([--disable-xmi-msim],[build without XMI-MSIM support]), [enable_xmi_msim=$enableval], [enable_xmi_msim=check])

if test x$enable_xmi_msim != "xno" ; then 
	PKG_CHECK_MODULES(xmi_msim, libxmimsim >= 7.2, [XMIMSIM_FOUND=yes], [XMIMSIM_FOUND=no])
	PKG_CHECK_MODULES(gmodule, gmodule-2.0, [GMODULE_FOUND=yes],[GMODULE_FOUND=no])

	if test "x$XMIMSIM_FOUND" = xyes && test "x$GMODULE_FOUND" = xyes ; then
		AC_MSG_NOTICE([Building with XMI-MSIM support])
		AC_DEFINE([HAVE_XMIMSIM], [], [XMI-MSIM support detected])
		AC_SUBST(HAVE_XMIMSIM, " -DHAVE_XMIMSIM")
	fi
else
	XMIMSIM_FOUND=no
	GMODULE_FOUND=no
fi


AM_CONDITIONAL(ENABLE_XMIMSIM,[test "x$XMIMSIM_FOUND" = xyes && test "x$GMODULE_FOUND" = xyes])

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
			src/beamsource/Makefile
			src/phasecontrast/Makefile
			src/spectrum/Makefile
			src/xrmc_math/Makefile
			src/exception/Makefile
			src/sample/Makefile
			src/main/Makefile
			src/randmt/Makefile
			src/fft/Makefile
			src/xmi-msim/Makefile
			src/radionuclide/Makefile
			doc/Makefile
			examples/Makefile
			examples/cylind_shell/Makefile
			examples/materials/Makefile
			examples/quadrics/Makefile
			examples/star1/Makefile
			examples/star2/Makefile
			examples/wheel/Makefile
			examples/xmi-msim/Makefile
			examples/xmi-msim2/Makefile
			examples/fluor_layers/Makefile
			examples/anisotropicsource/Makefile
			examples/beamsource/Makefile
			examples/phasecontrast/Makefile
			examples/radionuclide/Makefile
			libxrmc.pc
			xrmc.spec
			nsis/Makefile
			nsis/xrmc-win32.nsi
			nsis/xrmc-win64.nsi])
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

