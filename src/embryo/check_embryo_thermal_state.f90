subroutine check_embryo_thermal_state(j)
!*********************************************************************
! Check if embryo has melted ices, vapourised grains or dissociated H2
!*********************************************************************

use embryodata

implicit none
integer, intent(in) :: j

! If the temperature exceeds Tmelt, ices are evaporated from the solids component
IF(embryo(j)%T > Tmelt.and.embryo(j)%imelt==0) THEN
   embryo(j)%imelt =1 
   embryo(j)%fg = embryo(j)%fg*0.333
   !print*, "Embryo ",j,"'s ices have melted", embryo(j)%fg*embryo(j)%m/mearth             
ENDIF

! If the temperature of the grains exceeds Tvap, dust vapourised

IF(embryo(j)%T > Tvap.and.embryo(j)%ivap==0) THEN
   embryo(j)%imelt = 1
   embryo(j)%ivap = 1
   embryo(j)%fg = 0.0
   IF(embryo(j)%ijeans==0) THEN
      embryo(j)%rcore = 0.0
      embryo(j)%rg = 0.0
      embryo(j)%mcore =0.0
   ENDIF
   !	print*, 'Embryo ', j, ' has grains vapourised: core growth ends'           
ENDIF

! If the temperature exceeds the dissociation temperature of H2, second core collapse begins - gas giant/brown dwarf forms

IF(embryo(j)%T>Tdiss .and.embryo(j)%idiss==0) THEN
   embryo(j)%imelt=1
   embryo(j)%ivap = 1
   embryo(j)%idiss = 1
   
   IF(embryo(j)%ijeans==0) THEN 
      embryo(j)%rcore = 0.0
      embryo(j)%rg = 0.0
      embryo(j)%mcore =0.0
   ENDIF

   ! Object is a brown dwarf - give it appropriate radius
   ! Radius derived from linear approximation to Burrows et al (1997)

   IF(embryo(j)%M/Mjup < 70.0 ) THEN
      embryo(j)%R = -2.8e-3*embryo(j)%M/Mjup + 1.0
   ELSE
      embryo(j)%R = 0.01* embryo(j)%M/Mjup + 0.2
   ENDIF
   embryo(j)%R = embryo(j)%R*Rjup
   !    print*, 'Embryo ', j, ' is a brown dwarf'

ENDIF
   
end subroutine check_embryo_thermal_state
