propagation:
	gfortran propagate.f90 graphics.f90 dfft.f -o propagate

analysis:
	gfortran analysis.f90 dfft.f -o analysis

propagate_HA:
	@cd harmonic && ./../propagate && mv psi* data/

analyze_HA:
	@cd harmonic && ./../analysis

propagate_DW:
	@cd double-well && ./../propagate && mv psi* data/

analyze_DW:
	@cd double-well && ./../analysis

animation_HA:
	convert -set dispose previous -delay 20 harmonic/data/psi*.ps harmonic/harmonic.gif

animation_DW:
	convert -set dispose previous -delay 20 double-well/data/psi*.ps double-well/double-well.gif
