subroutine move_embryos
! Handles the motion of each embryo in the disc
! This motion is either by simple drag calculations, or by N Body motion

use stardata
use embryodata
use eosdata

implicit none

integer :: j


if(nbody=='y') then


else

   do j=1,nembryo

   ! Calculate number of timesteps required to traverse one grid
   embryo(j)%Nsteps = ((dr* embryo(j)%tmig/embryo(j)%a)-embryo(j)%t_spent)/dt

   ! Update time spent at this radius
   embryo(j)%t_spent = embryo(j)%t_spent + dt

   ! If core has spent long enough at this radius, then move it inwards one grid

   IF(embryo(j)%Nsteps<=0) THEN
      !print*, 'Moving embryo ', j, t/yr, embryo(j)%Nsteps, embryo(j)%t_spent/yr,tmig/yr, dr/embryo(j)%a
      embryo(j)%a = embryo(j)%a-dr
      embryo(j)%icurrent = embryo(j)%icurrent-1
      embryo(j)%t_spent = 0.0
      IF(embryo(j)%icurrent==1) embryo(j)%finished=1

   ENDIF

   enddo

endif

end subroutine move_embryos
