subroutine check_for_dead_embryos(j)
  !*************************************
  ! Check if embryo j has been destroyed
  !************************************

use stardata, only: r_d
use embryodata
implicit none

integer, intent(in) :: j

! If the embryo has disappeared, mark it as finished
IF(embryo(j)%r <1.0 .and. embryo(j)%rcore <1.0) THEN
   embryo(j)%finished = 1
ENDIF

! If the embryo has reached the inner disc radius, mark it as finished

IF(embryo(j)%a <r_d(2)) THEN
   embryo(j)%finished=1
ENDIF

! If the embryo is less than a thousandth of an earth mass, assume it is totally destroyed
! mark it as finished, and delete its properties appropriately

IF(embryo(j)%m/mearth < 1.0e-3) THEN
   embryo(j)%finished=1
   embryo(j)%r=0.0           
ENDIF

end subroutine check_for_dead_embryos
