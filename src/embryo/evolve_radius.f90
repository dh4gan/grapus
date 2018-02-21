subroutine evolve_radius(j,t)
!****************************************************************
! Evolve radius of embryo j at time t
! Currently assumes a polytropic collapse (cf protostar formation)
!****************************************************************

use eosdata,only: pi
use embryodata

implicit none

integer, intent(in) :: j
real, intent(in) :: t
real :: rchoose


IF(embryo(j)%itidal==0) rchoose = embryo(j)%R0
IF(embryo(j)%itidal==1) rchoose = embryo(j)%R 

embryo(j)%R = rchoose/(1.0 + 2.0*t/embryo(j)%t_cool0)**0.5
embryo(j)%rhoc = embryo(j)%M/(4.0*pi*theta_grad*embryo(j)%R**3)           


end subroutine evolve_radius
