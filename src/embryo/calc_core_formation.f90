subroutine calc_core_formation(j,t)

  !************************************************
  ! Computes collapse of grain cluster into a core
  !************************************************
  
  use embryodata
  use eosdata, only: pi
  implicit none
  
  integer, intent(in) :: j
  real, intent(in) :: t

  real :: l_jeans

  ! Compute jeans length of grain cluster

  l_jeans = embryo(j)%R*embryo(j)%fg**0.5

  IF(embryo(j)%rg < l_jeans .or. embryo(j)%rg < 0.0) THEN
     
     print*, "Embryo ", j,"forms a core ",embryo(j)%idiss
     embryo(j)%ijeans = 1
     
     ! Set up initial core parameters
     embryo(j)%mcore = embryo(j)%fg*embryo(j)%m
     embryo(j)%rcore = (3.0*embryo(j)%mcore/(4.0*pi*rho_s))**0.333
     embryo(j)%rg = embryo(j)%rcore
     
     ! If testing for radiative feedback, do it here
     if(core_feedback=='y') call calc_core_radiativefeedback(j)
     
  ELSE

     ! If rg > Jeans length, no core is formed
     embryo(j)%mcore = 0.0
     embryo(j)%rcore = 0.0
  ENDIF
 
end subroutine calc_core_formation
