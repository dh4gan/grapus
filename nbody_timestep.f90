subroutine nbody_timestep(position,velocity)
! Adjusts the timestep based on a standard step doubling algorithm
! Takes two half timesteps and compares to the current result

use stardata, only : istar ! Debug line (TODO-delete)
use eosdata, only : twopi ! Debug line also
use embryodata

implicit none

real,dimension(3,nbodies),intent(in) :: position,velocity

integer :: ibody,ix
real :: halfdt,error
real, dimension(3,nbodies) :: testpos1,testvel1,testpos2,testvel2

halfdt = dt_nbody/2

call nbody_integrate(halfdt,pos,vel,testpos1,testvel1)
call nbody_integrate(halfdt,testpos1,testvel1,testpos2,testvel2)

! Compute error between half timestep and full timestep for each particle

error = 0.0
maxerror = -1.0e30

if(istar==3) print*, 'HALF TIMESTEP', halfdt
do ibody=1,nbodies

   error = 0.0
   do ix=1,3

      error = error + (position(ix,ibody)-testpos2(ix,ibody))*(position(ix,ibody)-testpos2(ix,ibody))
      error = error + (velocity(ix,ibody)-testvel2(ix,ibody))*(velocity(ix,ibody)-testvel2(ix,ibody))

      if(istar==3) print*,ibody, ix, error, position(ix,ibody),&
           testpos1(ix,ibody),testpos2(ix,ibody)
   enddo

   error = sqrt(error)

   if(error>maxerror) maxerror = error
enddo

! Compute new deltat based on this error level
! If current error below user defined tolerance, then dt is increased
! If current error above tolerance, then dt is decreased

if(istar==3) then
   print*, 'TIMESTEP CHECK: dt, tol, maxerr: ',dt_nbody,tolerance,maxerror
  ! print*, position
  !print*, velocity
 !  print*, acc
endif

if(maxerror>small) then
   dt_nbody = dt_nbody*abs(tolerance/maxerror)**0.2
endif

if(istar==3) print*, 'NEW TIMESTEP: ',dt_nbody/twopi

if(maxerror>tolerance) then
   !print*,'error too large'
   !print*,tolerance,maxerror,dt_nbody
   dt_nbody = dt_nbody*0.1
endif

end subroutine nbody_timestep
