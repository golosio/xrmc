xrmc_exampledir = $(datadir)/xrmc/examples/phasecontrast

dist_xrmc_example_DATA = composition.dat phcdetector.dat geom3d.dat input.dat quadric.dat sample.dat source.dat spectrum.dat

#check_SCRIPTS = xrmc.sh

#TESTS = $(check_SCRIPTS)

xrmc.sh:
if !ABS_SRC_BUILD_EQUAL
	cp ${srcdir}/*dat .
endif
	echo "${top_builddir}/src/main/xrmc input.dat" > xrmc.sh
	chmod a+x xrmc.sh



clean-local:
	rm -f xrmc.sh
	rm -f image.dat
	rm -f convoluted_image.dat
if !ABS_SRC_BUILD_EQUAL
	rm -f $(dist_xrmc_example_DATA)
endif
