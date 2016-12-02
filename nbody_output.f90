subroutine nbody_output(t)
! Outputs data to file
! Currently writes each particle to separate file

use embryodata
use eosdata,only: yr,twopi

implicit none
real, intent(in) :: t
integer :: ibody

102 format (1P,23E15.5)
103 format (1P, 7E15.5)

call nbody_acceleration(pos,vel,acc)
call nbody_system_properties

! output individual bodies to separate files

   do ibody=1,nbodies
      write(ibody+ilog,102) t/yr, mass(ibody),pos(:,ibody), vel(:,ibody), acc(:,ibody),&
semimaj(ibody),ecc(ibody),inc(ibody),longascend(ibody),argper(ibody),trueanom(ibody), &
ekin(ibody),epot(ibody),etot(ibody),angmom(:,ibody)
      call flush(ibody)
   enddo

! Write log file containing global simulation data

write(ilog,103) t/yr, dt_nbody/twopi, 100.0*maxerror/tolerance, system_energy, dE, system_ang, dL


end subroutine nbody_output
