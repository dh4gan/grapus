subroutine move_embryos
! Handles the motion of each embryo in the disc
! This motion is either by simple drag calculations, or by N Body motion

use stardata
use embryodata
use eosdata

implicit none

integer :: j,ibody


if(nbody=='y') then

    ! Update the mass arrays for N Body calculation
    mass(1) = mstar/umass

    do ibody=2,nbodies
    mass(ibody) = embryo(ibody-1)%m/umass
    enddo

    ! Perform the nbody integration to find new positions
    call nbody_rk4

    ! Find the embryos grid index in the disc model
    call get_icurrent

else

   do j=1,nembryo

   ! Calculate number of timesteps required to traverse one grid
   embryo(j)%Nsteps = ((dr* embryo(j)%tmig/embryo(j)%a)-embryo(j)%t_spent)/dt

   ! Update time spent at this radius
   embryo(j)%t_spent = embryo(j)%t_spent + dt

   ! If core has spent long enough at this radius, then move it inwards one grid

   IF(embryo(j)%Nsteps<=0) THEN
      !print*, 'Moving embryo ', j, t/yr, embryo(j)%Nsteps, embryo(j)%t_spent/yr,embryo(j)%tmig/yr, dr/embryo(j)%a
      embryo(j)%a = embryo(j)%a-dr
      embryo(j)%icurrent = embryo(j)%icurrent-1
      embryo(j)%t_spent = 0.0
      IF(embryo(j)%icurrent==1) embryo(j)%finished=1

   ENDIF

   enddo

endif



end subroutine move_embryos


subroutine get_icurrent

use stardata
use embryodata

! For N Body runs, finds the grid index of the embryo
! Allows grid indices outside the model region

do j=1,nembryo
   embryo(j)%icurrent = int((embryo(j)%rmag-rin)/dr)+1

   ! embryos on the inner boundary are stopped
   IF(embryo(j)%icurrent<=1) embryo(j)%finished=1

   if(istar==3) print*, embryo(j)%semimaj, embryo(j)%ecc, embryo(j)%icurrent, embryo(j)%finished
enddo

end subroutine get_icurrent
