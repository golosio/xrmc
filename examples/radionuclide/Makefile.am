xrmc_exampledir = $(datadir)/examples/xrmc/radionuclide

dist_xrmc_example_DATA = composition.dat detector.dat geom3d.dat input.dat quadric.dat sample.dat source.dat radionuclide.dat

if ENABLE_XMIMSIM
check_SCRIPTS = xrmc.sh
TESTS_ENVIRONMENT = XRMC_XMIMSIM_MODULE="${top_builddir}/src/xmi-msim/.libs"
TESTS = $(check_SCRIPTS)
endif

xrmc.sh:
if !ABS_SRC_BUILD_EQUAL
	cp ${srcdir}/*dat .
endif
	echo "${top_builddir}/src/main/xrmc input.dat" > xrmc.sh
	chmod a+x xrmc.sh

clean-local:
	rm -f xrmc.sh
	rm -f output.dat
	rm -f unconvoluted.dat
if !ABS_SRC_BUILD_EQUAL
	rm -f $(dist_xrmc_example_DATA)
endif
