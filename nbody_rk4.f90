subroutine nbody_rk4
! This subroutine drives the N-Body integration over a single timestep
! Integration is done in a separate set of arrays for speed
! Orbital elements are stored in GE_embryo type

! Do integration

logical :: withintolerance

withintolerance = .true.

do while(withintolerance .eqv. .false.)


call integrate(dt_nbody,pos,vel,newpos,newvel)
call nbody_timestep(newpos,newvel)

if(maxerror>tolerance) withintolerance=.false.

end do

pos = newpos
vel = newvel

end subroutine nbody_rk4
