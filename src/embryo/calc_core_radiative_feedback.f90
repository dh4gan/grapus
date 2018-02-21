subroutine calc_core_radiativefeedback(j)

  use embryodata
  implicit none

  integer, intent(in) :: j
  
  real :: core_energy,embryo_energy
  
  ! evaluate radiative feedback of core formation on envelope OpenmP FLAGS!
  core_energy = embryo(j)%mcore*embryo(j)%mcore/embryo(j)%rcore
  embryo_energy = embryo(j)%m*embryo(j)%m/embryo(j)%r
  
  ! If energy released in core formation greater than envelope binding energy, destroy the embryo
  if(core_energy>embryo_energy) then
     embryo(j)%r =0.0
     embryo(j)%M = 0.0                       
  endif
  
end subroutine calc_core_radiativefeedback
