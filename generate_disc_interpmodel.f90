SUBROUTINE generate_disc
! Subroutine selects a self-gravitating disc model from file

use stardata
use eosdata
use embryodata, only: fg

implicit none

integer :: i,j,nskips
real :: mtry, mdot_try,sigma_old,dT,fine
real :: rmax,sigma0,sigma_prof,Qfrag

mdisc = 0.0
i=0

imodel =imodel +1

IF(imodel==1) THEN

! Do pass of first model to check number of time dumps
   OPEN(12,file=discfile,status='unknown')
  ! read(12,*) nmodels,ntime_mod,nrad_mod

   ntime_mod = 900
   nrad_mod = 200
   nmodels = 100
! Allocate arrays

   allocate(timemod(ntime_mod))
   allocate(mstar_mod(ntime_mod))
   allocate(r_mod(ntime_mod,nrad_mod))
   allocate(sigma_mod(ntime_mod,nrad_mod))
   allocate(nu_mod(ntime_mod,nrad_mod))
   allocate(Tc_mod(ntime_mod,nrad_mod))
   allocate(alpha_mod(ntime_mod,nrad_mod))
   allocate(tau_mod(ntime_mod,nrad_mod))

ELSE IF(imodel == nmodels) THEN

   print*, 'All models tried'

   imodel=1

   rewind(12)

ENDIF

! Find initial conditions of next disc model

! Move back in the file a little to make sure next entry isn't missed
! (Only for imodel >1)

!if(imodel>1) then
!   nskips = (nrad_mod+1)*2
!   do i=1,nskips
!      backspace(12)
!   enddo

!endif

i = 1

timemod(i) = 1.0d0

Do While (timemod(i) .ne. 0.0d0)
   read(12,*) timemod(i), Lx
   do j = 1, nrad_mod
      read(12,*) r_mod(i,j), sigma_mod(i,j), nu_mod(i,j), &
                       Tc_mod(i,j),tau_mod(i,j),alpha_mod(i,j), Mstar_mod(i)
   enddo

EndDo

! Read entire model into memory



do i=2,ntime_mod
   read(12,*) timemod(i),Lx
   do j=1,nrad_mod
      read(12,*) r_mod(i,j), sigma_mod(i,j), nu_mod(i,j), &
                       Tc_mod(i,j),tau_mod(i,j),alpha_mod(i,j), Mstar_mod(i)
   enddo
enddo


! Print out initial conditions

!if (imodel==10) THEN
   write(75,*) imodel, Lx
   call flush(75)
!do i=1,ntime_mod
!   do j=1,nrad_mod

!      write(75,*) r_mod(1,j)/udist, sigma_mod(1,j), Tc_mod(1,j),nu_mod(1,j), tau!_mod(1,j),alpha_mod(1,j)
!   enddo
!enddo
!if(imodel==10) STOP
!ENDIF

! Set up "true" disc variables
IF(allocated(sigma_d)) deallocate(sigma_d)
IF(allocated(cs_d)) deallocate(cs_d)
IF(allocated(omega_d)) deallocate(omega_d)
IF(allocated(betac_d)) deallocate(betac_d)
IF(allocated(mjeans)) deallocate(mjeans)
IF(allocated(r_d)) deallocate(r_d)
IF(allocated(H_d)) deallocate(H_d)
IF(allocated(alpha_d)) deallocate(alpha_d)
IF(allocated(T_d)) deallocate(T_d)
IF(allocated(kappa_d)) deallocate(kappa_d)
IF(allocated(gamma_d)) deallocate(gamma_d)
IF(allocated(tau_d)) deallocate(tau_d)
IF(allocated(gamma_j)) deallocate(gamma_J)
IF(allocated(nu_d)) deallocate(nu_d)

allocate(sigma_d(nrannuli))
allocate(cs_d(nrannuli))
allocate(omega_d(nrannuli))
allocate(betac_d(nrannuli))
allocate(mjeans(nrannuli))
allocate(gamma_J(nrannuli))
allocate(r_d(nrannuli))
allocate(H_d(nrannuli))
allocate(alpha_d(nrannuli))
allocate(T_d(nrannuli))
allocate(kappa_d(nrannuli))
allocate(gamma_d(nrannuli))
allocate(tau_d(nrannuli))
allocate(nu_d(nrannuli))

! Set up disc radii

DO i=1,nrannuli
   r_d(i) = rin+ (i-1)*dr
ENDDO

! Interpolate t=0 variables in space to generate current ICs

CALL evolve_disc(0.0)

! Randomly sample the maximum disc radius between 50 and 100 AU
! Disc model ICs extend to 50 AU, assumes planets form outside this

if(truncate_disc=='y') then
   rmax = rtrunc*udist + (rtruncmax-rtrunc)*udist*ran2(iseed)
else
   rmax = r_d(nrannuli)
endif

rmax = rmax

! Find surface density profile at 50 AU for extrapolation out to rmax


!print*, 'Disc Model ',imodel
!print*, 'Mdisc is ', mdisc/umass
!print*, 'Mstar is ', mstar/umass
!print*, 'mass ratio is ', mdisc/mstar
!print*, 'Sigma profile is ', sigma_prof
!print*, 'Sigma_0 is ',sigma0
!print*, 'rmax is ',rmax/udist

print*, 'Generating Model ',imodel,Lx, mdisc/umass, mstar/umass, mdisc/mstar,rmax/udist, rout/udist, irout

sigma_prof = 1.0
sigma0 = sigma_d(irout-100)*(r_d(irout-100)/udist)**sigma_prof


! Calculate other disc variables from interpolated model
! and extend disc beyond model limit

fg = 0.01

DO i =1,nrannuli

   ! Exit if radii too large
   IF(r_d(i) > rmax) THEN
      rout = r_d(i-1)
      irout = i-1
      exit
   ENDIF

   ! Get omega from Keplerian assumption

   omega_d(i) = sqrt(G*mstar/(r_d(i)*r_d(i)*r_d(i)))

	!print*, omega_d(i), G, mstar, r_d(i)
   ! If we are outside the disc model radius, then extrapolate sigma
   ! Assume Q=1.5

   IF(i>= irout-100.and.r_d(i) < rmax) THEN
      sigma_d(i) = sigma0*(r_d(i)/udist)**(-sigma_prof)
      cs_d(i) = 1.5*pi*G*sigma_d(i)/omega_d(i)
   ELSE

      ! Otherwise calculate sound speed from nu

      IF(alpha_d(i)/=0.0) THEN
         cs_d(i) = sqrt(ABS(nu_d(i)*omega_d(i)/alpha_d(i)))
      ELSE
         cs_d(i) = 0.0
      ENDIF

   ENDIF

   ! H is now simply calculated

   H_d(i) = cs_d(i)/omega_d(i)

   ! as is midplane volume density

   IF(H_d(i)/=0.0) THEN
      rhomid = sigma_d(i)/(2.0*H_d(i))      
   ELSE
      rhomid = 0.0
   ENDIF

   ! Use EoS to calculate tau, gamma, betac

   IF(rhomid /=0.0) THEN
      CALL eos_cs(rhomid, cs_d(i))

      gamma_d(i) = gammamuT(1)
      kappa_d(i) = gammamuT(4)


      IF(i > irout-100) THEN
         tau_d(i) = sigma_d(i)*kappa_d(i)
         T_d(i) = gammamuT(3)
      ENDIF
      
      betac_d(i) = (tau_d(i)+1.0/tau_d(i))*cs_d(i)*cs_d(i)*omega_d(i)/&
           (sigma_SB*(T_d(i)**4.0)*gamma_d(i)*(gamma_d(i)-1.0))
	

      ! If outside the model disc radius, calculate alpha based on cooling time

      IF(i>irout-100) THEN
         alpha_d(i) = 4.0/(9.0*gamma_d(i)*(gamma_d(i)-1.0)*betac_d(i))
      ENDIF

      !	Calculate expected accretion rate
			
      mdotvisc = 3.0*pi*alpha_d(i)*cs_d(i)*cs_d(i)*sigma_d(i)/omega_d(i)
	
      ! Now calculate the local Jeans Mass

      mjeans(i) = 4.0*1.412*pi*pi*pi/(3.0*G)
      mjeans(i) = mjeans(i)*cs_d(i)*cs_d(i)*cs_d(i)/(omega_d(i)*(1.0+1.0/sqrt(betac_d(i))))
      
      IF(alpha_d(i)<0.1) THEN
         gamma_j(i) = 2.0*(9.0*alpha_d(i)*gamma_d(i)*(gamma_d(i)-1.0)/4.0 - 1/betac_d(i)) 
      ELSE
         gamma_j(i) = 2.0*(9.0*0.1*gamma_d(i)*(gamma_d(i)-1.0)/4.0 - 1/betac_d(i))
      ENDIF

      gamma_j(i) = 1.0/gamma_j(i)

      Qfrag = cs_d(i)*omega_d(i)/(pi*G*sigma_d(i))
!      IF(Qfrag > 1.5.or.i<irout-100) gamma_j(i)=10.0
      IF(Qfrag > 1.5) gamma_j(i) = 10.0      

      ! Calculate enclosed Mass and mass ratio

      !mdisc = mdisc + twopi*r_d(i)*sigma_d(i)*dr

   ELSE

      gamma_d(i) = 0.0
      kappa_d(i) = 0.0
      betac_d(i) = 0.0
      mjeans(i) = 0.0
      gamma_j(i) = 0.1
   ENDIF

!   write(76,*) r_d(i)/udist,sigma_d(i),T_d(i),tau_d(i),gamma_d(i),nu_d(i),betac_d(i),alpha_d(i),cs_d(i),gamma_j(i)

ENDDO

q_disc = mdisc/mstar

END SUBROUTINE generate_disc
