xrmc_exampledir = $(datadir)/examples/xrmc/fluor_layers

dist_xrmc_example_DATA = spectrumAg40kV.dat input.dat composition.dat quadric.dat continuous_spectrum.dat detector_0deg.dat geom3d.dat sample.dat source_20deg.dat detector_response.c detector_routines.c detector_routines.h experimental.dat make.sh oct_plot.m README

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
	rm -f output.dat
if !ABS_SRC_BUILD_EQUAL
	rm -f $(dist_xrmc_example_DATA)
endif
