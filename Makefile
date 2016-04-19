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
FFLAGS = -O3 -frecord-marker=4 -fdefault-real-8  -fbounds-check -fopenmp

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
	migration_timescales.f90 timestep.f90
	 

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
