SUBROUTINE evolve_disc(t,dt)
  ! Subroutine updates disc properties according to simplified semi-analytic model
  ! Q and accretion rates are held constant, but mass is deducted from outer and inner regions

  use stardata
  use eosdata
  use embryodata
  
  implicit none
  integer :: i
  real ::  t,dt,dsigma, mlost,dtmin,counter
  real :: sigma_count

  ! Calculate mass lost in outer radius - if all mass lost, change irout

  dsigma = 1.1*sigma_d(irout)
  mlost = 0.0
  sigma_count = 0.0

  DO WHILE (dsigma > sigma_d(irout))

     dsigma = mdotvisc*dt/(2.0*pi*r_d(irout)*dr)

     
  IF(dsigma< sigma_d(irout)) THEN
     sigma_d(irout) = sigma_d(irout) - dsigma
     mlost= mlost + 2.0*pi*r_d(irout)*dsigma*dr

     exit
  ELSE IF(dsigma > sigma_d(irout)+sigma_count) THEN

! If mass loss greater than available mass in shell, pass it through shells until this is no longer true

     sigma_count = sigma_count + sigma_d(irout)
     mlost = mlost + 2.0*pi*r_d(irout)*sigma_d(irout)*dr

     sigma_d(irout) = 0.0
     irout = irout-1

  ELSE IF(dsigma < sigma_d(irout) + sigma_count) THEN
     exit
     ENDIF
     

  ENDDO    

! Update mass of star and disc
! If timestep too long, and disc disappears, set disc mass to zero and exit
  IF(mdotvisc*dt > mdisc) THEN
  mstar = mstar+ mdisc
  mdisc =0.0
  q_disc = 0.0
  return
  ENDIF

    mstar = mstar + mdotvisc*dt
    mdisc = mdisc - mdotvisc*dt

    q_disc = mdisc/mstar
    
    ! As Omega has changed, entire disc properties must be recalculated
    ! Sigma is known, so no iteration required

    DO i=1,irout
       r_d(i) = rin + (i-1)*dr

       omega_d(i) = sqrt(G*mstar/(r_d(i)*r_d(i)*r_d(i)))

       !	Calculate sound speed assuming fixed Q
       
       cs_d(i) = pi*G*sigma_d(i)/(omega_d(i))
       
       ! Calculate scale height
       
       H_d(i) = cs_d(i)/omega_d(i)
       rhomid = sigma_d(i)/(2.0*H_d(i))

       ! Use EoS to calculate tau, gamma, betac

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
			
!	mdotvisc = 3.0*pi*alpha_d(i)*cs_d(i)*cs_d(i)*sigma_d(i)/omega_d(i)			
!	Check for MRI activation
!	If so, then set alpha=0.01 and readjust sigma to maintain mdotvisc
			
	IF(T_d(i)>1000.0) THEN
!	print*, 'MRI active ',r/udist, T, alpha, sigma
	alpha_d(i) = 0.01
	sigma_d(i) = mdotvisc*omega_d(i)/(3.0*pi*alpha_d(i)*cs_d(i)*cs_d(i))
!	print*, 'MRI active ', T, alpha, sigma
	ENDIF			
			
 mjeans(i) = 4.0*1.412*pi*pi*pi/(3.0*G)
 mjeans(i) = mjeans(i)*cs_d(i)*cs_d(i)*cs_d(i)/(omega_d(i)*(1.0+1.0/sqrt(betac_d(i))))
 
 IF(alpha_d(i)<0.1) THEN
      gamma_j(i) = 2.0*(9.0*alpha_d(i)*gamma_d(i)*(gamma_d(i)-1.0)/4.0 - 1/betac_d(i)) 
   ELSE
      gamma_j(i) = 2.0*(9.0*0.1*gamma_d(i)*(gamma_d(i)-1.0)/4.0 - 1/betac_d(i))
   ENDIF

   gamma_j(i) = 1.0/gamma_j(i)


ENDDO


return
end subroutine evolve_disc

