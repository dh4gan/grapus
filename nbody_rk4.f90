subroutine nbody_rk4
! This subroutine drives the N-Body integration over a single timestep
! Integration is done in a separate set of arrays for speed
! Orbital elements are stored in GE_embryo type

! Do integration

use embryodata

logical :: withintolerance

withintolerance = .false.
print*, 'Attempting integration RK4'

do while(withintolerance .eqv. .false.)

newpos(:,:) = 0.0
newvel(:,:) = 0.0

call nbody_integrate(dt_nbody,pos,vel,newpos,newvel)
call nbody_timestep(newpos,newvel)

if(maxerror<tolerance) withintolerance=.true.

end do

pos = newpos
vel = newvel
print*, 'New positions:'

print*, pos

print*, 'New timestep: ',dt_nbody

call nbody_system_properties

end subroutine nbody_rk4
