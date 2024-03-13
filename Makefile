.PHONY: analysis

propagation:
	gfortran src/propagate.f90 src/graphics.f90 src/dfft.f -o src/propagate

analysis:
	gfortran src/analysis.f90 src/dfft.f -o src/analysis

propagate_HA:
	@cd harmonic && ./../src/propagate && mv psi* data/

analyze_HA:
	@cd harmonic && ./../src/analysis

propagate_DW:
	@cd double-well && ./../src/propagate && mv psi* data/

analyze_DW:
	@cd double-well && ./../src/analysis

animation_HA:
	convert -set dispose previous -delay 20 harmonic/data/psi*.ps harmonic/harmonic.gif

animation_DW:
	convert -set dispose previous -delay 20 double-well/data/psi*.ps double-well/double-well.gif
