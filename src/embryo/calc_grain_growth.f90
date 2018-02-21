subroutine calc_grain_growth(j,t)
!***********************************************
! Compute grain growth inside embryo j at time t
!***********************************************

use embryodata
implicit none

integer,intent(in) :: j
real,intent(in) :: t


! Radius of grain population = embryo radius if grains not fully grown
embryo(j)%rg = embryo(j)%R
   
! Once grains have fully grown, mark them as ready for sedimentation
   
if(t>embryo(j)%t_grow) then
   embryo(j)%igrown = 1
   embryo(j)%rg0 = embryo(j)%R                    
endif



end subroutine calc_grain_growth
