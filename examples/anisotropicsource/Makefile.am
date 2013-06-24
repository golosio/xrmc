xrmc_exampledir = $(datadir)/examples/xrmc/anisotropicsource

dist_xrmc_example_DATA = anisotropicsource.dat composition.dat detector.dat geom3d.dat input.dat intensityscreen.dat intensityscreen_image.c quadric.dat sample.dat spectrum.dat make.sh README

clean-local:
	rm -f image.dat
if !ABS_SRC_BUILD_EQUAL
	rm -f $(dist_xrmc_example_DATA)
endif
