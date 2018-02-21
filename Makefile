#####################################################
###                                               ###
###             Makefile for GRAPUS               ###
###                                               ###
###         Duncan H. Forgan (24/01/2018)         ###
###       				          ###
###                                               ###
#####################################################

# Compiler variable:
FC     = gfortran
VPATH = src/main/ src/disc/ src/embryo/ src/eos/ src/nbody/


# For serial runs use these flags
FFLAGS = -O3 -frecord-marker=4 -fdefault-real-8  -fbounds-check -Wunused

# For OpenMP runs
FFLAGS = -O3 -frecord-marker=4 -fdefault-real-8  -fbounds-check -fopenmp -Wunused

# Create object files:

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

SOURCESAF90 = star_module.f90 embryo_module.f90 eosmodule.f90  main.f90 \
	accrete_gas.f90 calc_core_formation.f90 \
	calc_core_radiative_feedback.f90 \
	calc_grain_growth.f90 calc_grain_sedimentation.f90 \
	calc_tidal_disruption.f90 check_embryo_thermal_state.f90 \
	check_for_dead_embryos.f90 eosread.f90 eos_cs.f90 \
	 evolve_disc_interpmodel.f90 evolve_radius.f90 \
	evolve_temperature.f90 evolve_embryos.f90 \
	generate_disc_interpmodel.f90 generate_embryos.f90 \
	generate_star.f90 initial.f90 \
	interpolate_1D.f90 interpolate_2D.f90 \
	move_embryos.f90 migration_timescales.f90 \
	nbody_integrate.f90 nbody_timestep.f90 nbody_acceleration.f90 \
	nbody_grav_acceleration.f90 nbody_drag_terms.f90 \
	nbody_deallocate_arrays.f90 nbody_orbits.f90 \
	nbody_output.f90 nbody_system_properties.f90  nbody_rk4.f90 \
	sample_gaussian.f90 timestep.f90 write_population_snapshot.f90 

OBJECTSA    = $(SOURCESAF90:.f90=.o)

# Create executable files:
build: grapus

grapus: $(OBJECTSA)
	$(FC) $(FFLAGS) -o $@ $(OBJECTSA)
 
# Clean statements:
clean: 
	\rm *.o *.mod grapus

# End Makefile
