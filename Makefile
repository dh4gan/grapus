#####################################################
###                                               ###
###     Makefile for the TD_synthesis code        ###
###                                               ###
###         Duncan H. Forgan (23/01/2012)         ###
###       				          ###
###                                               ###
#####################################################

# Compiler variable:
FC     = gfortran

# For serial runs use these flags
FFLAGS = -O3 -frecord-marker=4 -fdefault-real-8  -fbounds-check

# For OpenMP runs
FFLAGS = -O3 -frecord-marker=4 -fdefault-real-8  -fbounds-check -fopenmp -Wunused

# For files generated on stacpolly use these flags
#FC = ifort
#ZFFLAGS= -openmp -autodouble -O3
#ZZFFLAGS = ${ZFFLAGS} # -qflttrap=enable:invalid:zerodivide -g -C -qsigtrap -qf\
loat=nans
#FFLAGS = ${ZZFFLAGS} -fPIC -i-dynamic -xW

# For real files use these flags
#FFLAGS = -O3 -frecord-marker=4

# Create object files:
#%.o: %.f
#	$(FC) $(FFLAGS) -c $<

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

SOURCESAF90 = star_module.f90 embryo_module.f90 eosmodule.f90  main.f90 \
	 eosread.f90 eos_cs.f90 evolve_disc_interpmodel.f90 evolve_embryos.f90 \
	generate_disc_interpmodel.f90 generate_embryos.f90 \
	 generate_star.f90 interpolate_1D.f90 interpolate_2D.f90 \
	move_embryos.f90 migration_timescales.f90 \
	nbody_integrate.f90 nbody_timestep.f90 nbody_acceleration.f90 \
	nbody_grav_acceleration.f90 nbody_drag_terms.f90 \
	nbody_deallocate_arrays.f90 nbody_orbits.f90 \
	nbody_output.f90 nbody_system_properties.f90  nbody_rk4.f90 \
	sample_gaussian.f90 timestep.f90

OBJECTSA    = $(SOURCESAF90:.f90=.o)

# Create executable files:
build: TD_synth

TD_synth: $(OBJECTSA)
	$(FC) $(FFLAGS) -o $@ $(OBJECTSA)
 
#probedata: probedata.f90
#	 $(FC) $(FFLAGS) probedata.f90 -o probedata
#	rm -f probedata.o

# Clean statements:
clean: 
	\rm *.o *.mod TD_synth

# End Makefile
