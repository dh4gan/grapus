SUBROUTINE generate_disc
! Subroutine creates a self-gravitating disc, constrained by M_d and dot(M_d)
! We fix Q=1, and follow a similar procedure to Clarke (2009)

use stardata
use eosdata
use embryodata, only: fg

implicit none

integer :: i
real :: mtry, mdot_try,sigma_old,dT,fine

mtry = 0.0
i=0

! TODO - randomly generate fg, q_disc in some way

fg = 0.01 + ran2(iseed)*0.01
q_disc = 0.5 +ran2(iseed)*0.5

mdisc = q_disc*mstar
mdotvisc = -8.5 + 3.5*ran2(iseed)
mdotvisc = 10.0**(mdotvisc)*umass/yr

print*, 'Generating disc for star ',istar,': ', fg,mstar, q_disc, mdotvisc*yr/umass

sigma_old = 1.0e6

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

!print*, i, nrannuli, mtry,mdisc

DO WHILE(mtry < mdisc .and. i<nrannuli)

i = i+1
r_d(i) = rin + (i-1)*dr

omega_d(i) = sqrt(G*mstar/(r_d(i)*r_d(i)*r_d(i)))

! Iterate to find sigma

	sigma_d(i) = sigma_old

	dT = 1.0e30
	ntries = 0
	fine = 0.01
	DO WHILE(ABS(dT)> tolerance)
	
!	Calculate sound speed assuming fixed Q

	cs_d(i) = pi*G*sigma_d(i)/(omega_d(i))

!	Calculate scale height

	H_d(i) = cs_d(i)/omega_d(i)

	rhomid = sigma_d(i)/(2.0*H_d(i))


!	Use EoS to calculate tau, gamma, betac

	CALL eos_cs(rhomid, cs_d(i))

	gamma_d(i) = gammamuT(1)
	T_d(i) = gammamuT(3)
	kappa_d(i) = gammamuT(4)

	tau_d(i) = sigma_d(i)*kappa_d(i)
	betac_d(i) = (tau_d(i)+1.0/tau_d(i))*cs_d(i)*cs_d(i)*omega_d(i)/&
			(sigma_SB*(T_d(i)**4.0)*gamma_d(i)*(gamma_d(i)-1.0))


!	Calculate alpha from this value --> accretion rate
				
	alpha_d(i) = 4.0/(9.0*gamma_d(i)*(gamma_d(i)-1.0)*betac_d(i))
	
!	print*, alpha_d(i), gamma_d(i),betac_d(i),cs_d(i),omega_d(i)		
!	Compare with imposed accretion rate
			
	mdot_try = 3.0*pi*alpha_d(i)*cs_d(i)*cs_d(i)*sigma_d(i)/omega_d(i)
	
	dT = (mdotvisc-mdot_try)/mdotvisc
		
	sigma_old = sigma_d(i)
!	sigma_d(i) = sigma_d(i)*(1.0 + dT*0.01)
	sigma_d(i) = sigma_d(i)*(1.0 +dT/(abs(dT))*fine)		
	
!	print*,i, sigma_d(i), sigma_old, dT, mdot_try, mdotvisc
	ntries = ntries + 1

	IF (ntries>500) THEN
	fine = fine/10.0
	ntries = 0
	ENDIF
	ENDDO	

!	Check for MRI activation
!	If so, then set alpha=0.01 and readjust sigma to maintain mdotvisc
			
	IF(T_d(i)>1000.0) THEN
!	print*, 'MRI active ',r/udist, T, alpha, sigma
	alpha_d(i) = 0.01
	sigma_d(i) = mdotvisc*omega_d(i)/(3.0*pi*alpha_d(i)*cs_d(i)*cs_d(i))
!	print*, 'MRI active ', T, alpha, sigma
	ENDIF			
			
 sigma_d(i+1) = sigma_old
 mjeans(i) = 4.0*1.412*pi*pi*pi/(3.0*G)
 
 mjeans(i) = mjeans(i)*cs_d(i)*cs_d(i)*cs_d(i)/(omega_d(i)*(1.0+1.0/sqrt(betac_d(i))))
 
 IF(alpha_d(i)<0.1) THEN
      gamma_j(i) = 2.0*(9.0*alpha_d(i)*gamma_d(i)*(gamma_d(i)-1.0)/4.0 - 1/betac_d(i)) 
   ELSE
      gamma_j(i) = 2.0*(9.0*0.1*gamma_d(i)*(gamma_d(i)-1.0)/4.0 - 1/betac_d(i))
   ENDIF

   gamma_j(i) = 1.0/gamma_j(i)

!	Calculate enclosed Mass and mass ratio

	mtry = mtry + twopi*r_d(i)*sigma_d(i)*dr

IF(i==nrannuli) THEN
print*, 'ERROR: nrannuli insufficient'
print*, 'Either recompile or increase dr in .params file'
STOP
ENDIF
ENDDO

rout = r_d(i)
irout =i

END SUBROUTINE generate_disc
