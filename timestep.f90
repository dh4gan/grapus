  SUBROUTINE timestep
! Subroutine calculates minimum timestep based on evolutionary state of all active embryos
! Also checks to confirm whether all embryos are finished
use stardata
use embryodata
use eosdata

implicit none


integer :: j,migtype
real :: tcheck, tmig,tgap,tcross

dt = 1.0e30
finishcheck=1

!print*, 'Calculating dt '

DO j=1,nembryo  
   ! Use this checksum to see if all embryos finished
   finishcheck =finishcheck*embryo(j)%finished 
   tcheck = embryo(j)%t_cool0
   ! Only consider timescale if embryo has not evaporated grains
   IF(embryo(j)%finished==0.and.embryo(j)%ivap==0) THEN
    !  print*,j, embryo(j)%igrown,embryo(j)%t_grow, embryo(j)%iself,embryo(j)%t_sed
      ! Find appropriate timescale given evolutionary state
      IF(embryo(j)%igrown==0) THEN
         tcheck = embryo(j)%t_grow
      else if(embryo(j)%iself==0) THEN
         tcheck=embryo(j)%t_sed/30.0 ! Factor of 30 ensures that self-gravitating collapse is resolved (critical value is 3)
      ENDIF
      IF(tcheck < dt) dt = tcheck   

      ! Calculate migration timescales - if too short, then reduce the timestep
      call migration_timescales(j,migtype, tmig, tgap,tcross)

      if(tmig < dt) dt = tmig
   ENDIF

ENDDO

! Prevent overly long timesteps
dt = min(dt,100.0*yr)

!print*, 'dt ', dt
END SUBROUTINE timestep
