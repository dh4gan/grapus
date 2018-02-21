subroutine calc_grain_sedimentation(j,t)

!***********************************************
! Compute sedimentation of grown grains inside embryo j at time t
!***********************************************

use embryodata
use eosdata, only: pi
implicit none

integer,intent(in) :: j
real,intent(in) :: t
real :: factor

! Define initial sedimentation timescale
embryo(j)%t_sed = embryo(j)%t_sed0

! Define remaining mass of dust not stripped
embryo(j)%msolids = embryo(j)%fg*embryo(j)%M

! Check if grains travelling above critical velocity for destructive collisions

IF(embryo(j)%t_sed < embryo(j)%R/vfrag) THEN
   embryo(j)%t_sed = embryo(j)%R/vfrag
ENDIF

!
! If grains self-gravitate, sedimentation is boosted
!

IF(embryo(j)%iself==1) THEN
   factor = 1.0 - 3.0*(t-embryo(j)%tself)/embryo(j)%t_sed
   ! If growth rate too fast for this timestep, assume Jeans instability forms the entire core
   if(factor< 0.0) then
      ! Assume core fully formed at 1 g cm-3 constant density (this circumstance should be avoided by timestep correction)
      embryo(j)%rg = embryo(j)%msolids/(1.333*pi*1.689e6) 
   else
      ! Otherwise grow the core under self-gravity
      embryo(j)%rg = embryo(j)%rself*factor**0.333
   endif

else
   ! Case where no self gravity acts
   embryo(j)%rg = embryo(j)%rg*(1.0 - dt/embryo(j)%t_sed)
endif

if(embryo(j)%rg > embryo(j)%R) embryo(j)%rg = embryo(j)%R


! If Rg small enough that rho(grains) ~ rho(gas), mark core as self-gravitating, set rcore==Rg
if(embryo(j)%rg < (embryo(j)%fg**0.333)*embryo(j)%R .and. embryo(j)%iself==0.and.embryo(j)%ivap==0) then

   embryo(j)%iself=1                 
   embryo(j)%tself = t
   embryo(j)%rself = embryo(j)%rg
   embryo(j)%rcore = embryo(j)%rg

!   print*, 'Grains self-gravitating for embryo ',j,': ',embryo(j)%rself/udist,embryo(j)%tself/yr

endif



end subroutine calc_grain_sedimentation
