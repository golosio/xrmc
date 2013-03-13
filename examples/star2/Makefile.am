xrmc_exampledir = $(datadir)/examples/xrmc/star2

dist_xrmc_example_DATA = composition.dat detector.dat geom3d.dat input.dat params.dat quadric.dat sample.dat source.dat spectrum.dat

check_SCRIPTS = xrmc.sh

TESTS = $(check_SCRIPTS)

xrmc.sh:
if !ABS_SRC_BUILD_EQUAL
	cp ${srcdir}/*dat .
endif
	echo "${top_builddir}/src/main/xrmc input.dat" > xrmc.sh
	chmod a+x xrmc.sh



clean-local:
	rm -f xrmc.sh
	rm -f image.dat
if !ABS_SRC_BUILD_EQUAL
	rm -f $(dist_xrmc_example_DATA)
endif
