  SUBROUTINE timestep
! Subroutine calculates minimum timestep based on evolutionary state of all active embryos
! Also checks to confirm whether all embryos are finished
use stardata
use embryodata
use eosdata

implicit none


integer :: j
real :: tcheck

dt = 1.0e30
finishcheck=1

! Migration timescales and N-Body timestep already accepted

DO j=1,nembryo  
   ! Use this checksum to see if all embryos finished
   finishcheck =finishcheck*embryo(j)%finished 
   tcheck = embryo(j)%t_cool0

   ! Only consider timescale if embryo has not evaporated grains
   IF(embryo(j)%finished==0.and.embryo(j)%ivap==0) THEN

      ! Find appropriate timescale given evolutionary state
      IF(embryo(j)%igrown==0) THEN
         tcheck = embryo(j)%t_grow
      else if(embryo(j)%iself==0) THEN
         tcheck=embryo(j)%t_sed/30.0 ! Factor of 30 ensures that self-gravitating collapse is resolved (critical value is 3)
      ENDIF
      IF(tcheck < dt) dt = tcheck   

      if(embryo(j)%tmig < dt) dt = embryo(j)%tmig

   ENDIF

ENDDO


if(nbody=='y') then
   ! Check against N Body timestep
   ! N Body units: 2pi units = 1 yr

   if(dt_nbody*yr/twopi<dt) dt=dt_nbody*yr/twopi
endif

! Prevent overly long timesteps
dt = min(dt,100.0*yr)

END SUBROUTINE timestep
