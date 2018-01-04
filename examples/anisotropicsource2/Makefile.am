xrmc_exampledir = $(datadir)/examples/xrmc/anisotropicsource

dist_xrmc_example_DATA = anisotropicsource.dat composition.dat detector.dat geom3d.dat input.dat intensityscreen.dat quadric.dat sample.dat spectrum.dat intensityscreen_image.dat

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
