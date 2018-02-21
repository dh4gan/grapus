subroutine calc_tidal_disruption(j)
  !***************************************
  ! Computes tidal disruption of embryo j
  !***************************************

  use stardata,   only: mstar
  use embryodata
  use eosdata,    only: G,pi
  implicit none

  integer,intent(in) :: j

  real :: r_hill, hillcore, rstrip, orb

  ! Recalculate Hill Radius of embryo

  r_hill = embryo(j)%a*(embryo(j)%m/(3.0*mstar))**0.333
  
  embryo(j)%itidal = 0
  
  ! Check if tidal disruption is occurring
  
  if(embryo(j)%r > r_hill .and.embryo(j)%itidal==0) embryo(j)%itidal = 1
  
  
  if(embryo(j)%itidal==1.and.embryo(j)%idiss==0) then
     ! If we are in the tidal disruption regime, strip off the upper layers of gas and update properties

     ! Calculate the hill radius for the core only

     hillcore = embryo(j)%a*(embryo(j)%mcore/(3.0*mstar))**0.333

     ! Strip the outer layers on an orbital timescale
     
     if(embryo(j)%r> hillcore.and.r_hill < embryo(j)%r) then
        orb = sqrt(G*mstar/embryo(j)%a**3)
        orb = 2.0*pi/orb

        rstrip = max(r_hill,hillcore)
          
        embryo(j)%r = embryo(j)%r - (embryo(j)%r - rstrip)*(1.0-exp(-dt/(1.0*orb))) ! Gradual depletion of envelope over one orbital period

      
     endif
   
     embryo(j)%M = 4.0*pi*embryo(j)%rhoc*embryo(j)%r**3*theta_grad
      
     if(embryo(j)%m < embryo(j)%mcore) then
        embryo(j)%r = embryo(j)%rcore
        embryo(j)%m = embryo(j)%mcore
        
     endif
   

  endif



end subroutine calc_tidal_disruption
