SUBROUTINE evolve_embryos(t)
  !****************************************************************************
  ! Subroutine evolves the embryo population for a single timestep 
  ! This routine evolves embryos simultaneously on a constant timestep
  !****************************************************************************

  use stardata
  use embryodata
  use eosdata

  implicit none

  integer :: i,j
  real,intent(in) :: t

  !*************************************
  !1. Compute the motion of the embryos
  !*************************************

  ! Calculate migration timescales and gap opening criteria
  call migration_timescales

  ! Move embryos (either analytically or via N Body integration)
  call move_embryos

  !*************************************
  !2. Compute the internal structure of the embryos
  !*************************************

  !$OMP PARALLEL &
  !$OMP shared(embryo,nembryo,H_d,r_d) &
  !$OMP shared(t,dt,alpha_d,omega_d,mstar) &
  !$OMP private(i,j)
  !$OMP DO SCHEDULE(runtime)
  DO j=1,nembryo

     ! If embryo finished, skip to the next one
     IF(embryo(j)%finished==1) cycle

     ! Decide whether to accrete gas or not (depending on gap opening)
     ! Commented out - useful for future implementations
     ! In general, embryos open "cold gaps" as they are quite massive
     ! call accrete_gas

     ! Evolve the embryo radius and central temperature

     if(embryo(j)%idiss==0) then
        call evolve_radius(j,t)
        call evolve_temperature(j,t)
     endif

     ! Check for ice melting, grain vapourisation, H2 dissociation

     call check_embryo_thermal_state(j)

     ! If core not already formed by Jeans instability, AND
     ! If dust not vapourised, continue evolving the dust component

     IF(embryo(j)%ivap==0.and.embryo(j)%ijeans==0) THEN

        if(embryo(j)%igrown==0) then
           call calc_grain_growth(j,t)
        else if(embryo(j)%igrown==1) THEN
           call calc_grain_sedimentation(j,t)
        endif

        ! Check If grain cluster radius < jeans length -> core formation

        if(embryo(j)%iself==1) then
           call calc_core_formation(j,t)
        endif
     endif

     ! Check for tidal disruption of the embryo

     call calc_tidal_disruption(j)

     ! Check if embryo has been destroyed
     call check_for_dead_embryos(j)


  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL

  ! End of embryo loop       

END SUBROUTINE evolve_embryos
