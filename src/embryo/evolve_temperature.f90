subroutine evolve_temperature(j,t)
!****************************************************************
! Evolve temperature of embryo j at time t
! Currently assumes a polytropic collapse (cf protostar formation)
!****************************************************************

use eosdata, only: Boltzmannk, mH
use embryodata
implicit none

integer,intent(in) :: j
real,intent(in) :: t


embryo(j)%T = embryo(j)%T0*(1.0+(1.0+p_kap)*t/embryo(j)%t_cool0)**(1.0/(1.0+p_kap))
embryo(j)%cs = sqrt(gamma*Boltzmannk*embryo(j)%T/(mu*mH))
  


end subroutine evolve_temperature
