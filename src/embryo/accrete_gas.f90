subroutine accrete_gas(j)
!****************************************
! Accrete gas from disc onto embryo j
! depends on if a gap is opened
! (Future development)
!****************************************

use embryodata
implicit none

integer, intent(in) :: j

if(embryo(j)%migtype==1) then


endif

end subroutine accrete_gas
