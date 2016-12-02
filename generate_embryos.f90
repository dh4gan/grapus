SUBROUTINE generate_embryos
  ! Subroutine finds fragment radius, and then calculates fragment masses at equal log a spacing

  use stardata
  use embryodata
  use eosdata

  implicit none


  integer :: i,j
  real :: kappa,r_hill,rtest,exp1,exp2
  real, dimension(100) :: cspace

  i=0

  DO WHILE(i<nrannuli) 
     i=i+1

     IF(gamma_j(i)>-1.0 .and. gamma_j(i)<0.0) THEN
        rfrag = r_d(i)
	irfrag = i
        exit
     ENDIF
  ENDDO

! Move on rfrag by a random value, find new irfrag

!  rfrag = rfrag +10.0*ran2(iseed)
!  i = irfrag
!  DO WHILE(i<nrannuli)
!     i = i+1

!     IF(r_d(i) > rfrag) THEN
!        irfrag = i
!        exit
!     ENDIF
!  ENDDO

  IF(i==nrannuli) THEN

    ! print*, 'No fragments formed'
     IF(allocated(embryo)) deallocate(embryo)
     nembryo=0
     return
  ELSE
     print*, 'Fragmentation radius is ',rfrag/udist, ' AU ',irfrag
  ENDIF

  ! Generate number of possible embryos
  ! Embryo spacing is C Hill Radii, where C uniformly distributed
  ! between 1.5 and 3 (stored for later use)

  rtest = rfrag
  i=irfrag
  nembryo=0

  j = 0
  DO WHILE (i < irout)
     j =j+1
     cspace(j) = 1.5 + ran2(iseed)*1.5

     IF(r_d(i)>rmax) exit

     r_hill = rtest*(mjeans(i)/(3.0*mstar))**0.333

     DO WHILE (r_d(i) < rtest+cspace(j)*r_hill .and. i<irout)
        i=i+1
     END DO

     rtest = r_d(i)
     nembryo = nembryo+1
  ENDDO

  print*, 'There are ',nembryo, ' embryos'

  ! Allocate embryo type array

  IF(allocated(embryo)) deallocate(embryo)
  allocate(embryo(nembryo))

  ! Calculate initial embryo properties

  i=irfrag
  rtest =rfrag

  DO j=1,nembryo

     embryo(j)%m = mjeans(i)
     embryo(j)%a = rtest
     embryo(j)%iform = i


    ! If this is an n-body run, then set up positions,velocities
    embryo(j)%ecc = 0.0
    embryo(j)%inc = 0.0
    embryo(j)%longascend = 0.0
    embryo(j)%argper = 0.0
    embryo(j)%trueanom = twopi*ran2(iseed)



     exp1 = (1.0-n)/n
     exp2 = (3.0-n)/n

     m1 = embryo(j)%m/(0.01*umass)
     T1 = T_d(i)/10.0

     kappa_0 = kappa_d(i)/(T1)**p_kap
     kappa_star = kappa_0/0.01

     !T1 = 1.0
     !kappa_star = 1.0

     rho_ad = 5.0e-13*kappa_star**(-0.666)*(T1)**((4.0-2.0*p_kap)/3.0)

     embryo(j)%R0 = 17.5*udist*m1**(-0.333)*T1*(1.0e13*rho_ad)**(-0.666)
     embryo(j)%T0 = 146.0*m1**1.3333*T1**(-(1.0+4.0*p_kap)/9)*kappa_star**(-4.0/9.0)
     embryo(j)%cs0 = gamma_d(i)*Boltzmannk*embryo(j)%T0/(mu*mH)
     embryo(j)%cs0 = sqrt(embryo(j)%cs0)

     embryo(j)%R = embryo(j)%R0
     embryo(j)%rhoc = embryo(j)%M/(4.0*pi*theta_grad*embryo(j)%R**3)

     embryo(j)%RG = embryo(j)%R0
     embryo(j)%fg = fg
     embryo(j)%rcore = 0.0  ! No solid core formed yet
     embryo(j)%tself = 0.0
     embryo(j)%rself = 0.0

     !print*, T1, m1, 1.0e13*rho_ad, kappa_0, embryo(j)%R/udist

    ! Initialise state of embryo in terms of evolutionary phases

     embryo(j)%iself = 0
     embryo(j)%ivap = 0
     embryo(j)%imelt = 0
     embryo(j)%idiss = 0
     embryo(j)%igrown = 0
     embryo(j)%ijeans = 0
     embryo(j)%itidal = 0

     ! If embryo hot enough, ice can melt, grains will evaporate or H2 dissociates
     IF(embryo(j)%T0> Tmelt) embryo(j)%imelt = 1
     IF(embryo(j)%T0> Tvap) embryo(j)%ivap = 1
     IF(embryo(j)%T0> Tdiss) embryo(j)%idiss = 1

     kappa = kappa_d(i)

     IF(p_kap==1.0) THEN
        embryo(j)%t_cool0 = 380.0*yr*m1**0.666*T1**(-1.333)*kappa_star**(1.0/9.0)
     ELSE IF(p_kap==2.0) THEN
        embryo(j)%t_cool0 = 5700.0*yr*m1**2*T1**(-2.1111)*kappa_star**(-0.333)
     ENDIF

     ! Multiply by tunable parameter c_collapse
     embryo(j)%t_cool0 = c_collapse*embryo(j)%t_cool0

     mfp = 1.0e10*mfp0*embryo(j)%m/(embryo(j)%R0)**3

     ! Calculate critical grain size at which turbulence loses to sedimentation

     embryo(j)%scrit = -mfp + sqrt(mfp*mfp +6.0*mfp*alpha_dust*embryo(j)%M/(pi*embryo(j)%R0*embryo(j)%R0*rho_s))
     embryo(j)%scrit = embryo(j)%scrit/2.0

     embryo(j)%scrit = max(scrit_min, embryo(j)%scrit)
    ! embryo(j)%scrit = 10.0
     ! Calculate grain growth timescales for this critical size

     embryo(j)%t_grow0 = 3.0*embryo(j)%cs0*embryo(j)%R0*embryo(j)%R0/(embryo(j)%m*pi*fg*G)*log(embryo(j)%scrit/s0)

     ! Modify this growth timescale due to the core's contraction

     embryo(j)%t_grow = 1.0+embryo(j)%t_grow0*(2.5+p_kap)/embryo(j)%t_cool0
     embryo(j)%t_grow = embryo(j)%t_grow**(p_grow) -1
     embryo(j)%t_grow = embryo(j)%t_cool0*embryo(j)%t_grow/(1.0+p_kap)

     ! Calculate sedimentation timescale for this critical grain size

     embryo(j)%t_sed0 = 3.0*embryo(j)%cs0*mfp/(4.0*pi*G*rho_s*embryo(j)%scrit*(embryo(j)%scrit+mfp))
     embryo(j)%t_sed = embryo(j)%t_sed0

! Core mass, radius

	embryo(j)%mcore = 0.0
	embryo(j)%rcore = 0.0
	
     ! Now find location of next embryo

     r_hill = rtest*(mjeans(i)/(3.0*mstar))**0.333

     DO WHILE (r_d(i) < rtest+cspace(j)*r_hill .and. i<irout)
        i=i+1
     END DO

     rtest = r_d(i)

     !     print*, 'Embryo ', j
     !     print*, embryo(j)%m/umass, embryo(j)%a/udist,embryo(j)%R0/udist, embryo(j)%cs0, &
     !          embryo(j)%T0,embryo(j)%t_cool0/yr, embryo(j)%t_grow0/yr, embryo(j)%t_sed0/yr

     write(istart,'(8E18.10)') embryo(j)%a/udist, embryo(j)%m/mjup, embryo(j)%R0/rjup, & 
          embryo(j)%T0, embryo(j)%scrit, embryo(j)%t_cool0/yr, embryo(j)%t_grow/yr, embryo(j)%t_sed0/yr

     r_hill = embryo(j)%a*(embryo(j)%m/(3.0*mstar))**0.333
     WRITE(*, '(I2,1X,I4,6F18.10)') j, embryo(j)%iform, embryo(j)%m/mjup, embryo(j)%a/udist, embryo(j)%R/rjup,&
          embryo(j)%t_cool0/yr, embryo(j)%t_grow/yr, embryo(j)%t_sed0/yr

  ENDDO

! Write data referring to this star-disc-planets system to log file


  WRITE(ilog,'(5E18.10, I3)') mstar/umass, mdisc/umass, q_disc, rout/udist, rfrag/udist, nembryo
  RETURN

! If this is an N Body run, then create arrays for N body calculation (TODO)

if(nbody=='y') then

    allocate(pos(3,nembryo),vel(3,nembryo),acc(3,nembryo))
endif


END SUBROUTINE generate_embryos
