xrmc_exampledir = $(datadir)/examples/xrmc/beamsource

dist_xrmc_example_DATA = beamsource.dat composition.dat detector.dat geom3d.dat input.dat input_spectrum.dat beamscreen.dat beamscreen_image.c quadric.dat sample.dat make.sh README spectrum_detector.dat

clean-local:
	rm -f image.dat
if !ABS_SRC_BUILD_EQUAL
	rm -f $(dist_xrmc_example_DATA)
endif
