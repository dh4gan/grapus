SUBROUTINE evolve_disc(t_interp,timeup)
  ! Subroutine interpolates disc models in time (given input time value), and in space

  use stardata
  use eosdata
  use embryodata

  implicit none

  integer :: i,j,k,l,timeup

  real :: grad1, grad2, inter1,inter2, grad,inter, t_interp, rad
  real :: sdens1,sdens2,Tmid1,Tmid2, nu_pl1,nu_pl2,tau1,tau2
  real :: alphamid1, alphamid2,begintime,finishtime

  mdisc = 0.0

  ! Find time location first

  i=1
  DO WHILE(timemod(i) < t_interp)
     i = i+1
     IF(i> ntime_mod) THEN
        print*, 'Time up: t=',timemod(i-1)/yr
        timeup = 1
        return
     ENDIF
  ENDDO

  
  ! If at beginning of file , mstar equal to initial value
  IF(i==1) THEN

     mstar =mstar_mod(i)    
      
     ! Otherwise interpolate in time
     ELSE
        begintime = timemod(i-1)
        finishtime = timemod(i)
        
     CALL interpolate_1D(t_interp,mstar, begintime,finishtime,mstar_mod(i-1),mstar_mod(i))

  ENDIF

! Now do disc properties:

  IF(i==1) THEN         

     do k=1,nrannuli

        rad = r_d(k)
        
        j=1
        DO WHILE(r_mod(i,j) < rad.and. j<nrad_mod)
           j=j+1
        ENDDO

        ! If at the beginning spatial index, interpolate against zero

        IF(j==1) THEN

           ! Sigma
           CALL interpolate_1D(r_d(k), sigma_d(k),0.0,r_mod(i,j),0.0, sigma_mod(i,j))

           !print*, 'Sigma done for j=1'
           ! Nu
           CALL interpolate_1D(r_d(k), nu_d(k),0.0,r_mod(i,j),0.0, nu_mod(i,j))

           ! T
           CALL interpolate_1D(r_d(k), T_d(k),0.0,r_mod(i,j),0.0, Tc_mod(i,j))

           ! Alpha
           CALL interpolate_1D(r_d(k), alpha_d(k),0.0,r_mod(i,j),0.0, alpha_mod(i,j))

           ! Tau
           CALL interpolate_1D(r_d(k), tau_d(k),0.0,r_mod(i,j),0.0, tau_mod(i,j))


           ! Otherwise interpolate as normal in space
        ELSE

           
           ! Sigma
           CALL interpolate_1D(r_d(k), sigma_d(k),r_mod(i,j-1),r_mod(i,j),sigma_mod(i,j-1), sigma_mod(i,j))
          
           ! Nu
           CALL interpolate_1D(r_d(k), nu_d(k),r_mod(i,j-1),r_mod(i,j),nu_mod(i,j-1), nu_mod(i,j))

           ! T
           CALL interpolate_1D(r_d(k), T_d(k),r_mod(i,j-1),r_mod(i,j),Tc_mod(i,j-1), Tc_mod(i,j))

           ! Alpha
           CALL interpolate_1D(r_d(k), alpha_d(k),r_mod(i,j-1),r_mod(i,j),alpha_mod(i,j-1), alpha_mod(i,j))

           ! Tau
           CALL interpolate_1D(r_d(k), tau_d(k),r_mod(i,j-1),r_mod(i,j),tau_mod(i,j-1), tau_mod(i,j))


        ENDIF


        ! Calculate ancillary variables from interpolated data

        omega_d(k) = sqrt(G*mstar/r_d(k)*r_d(k)*r_d(k))
        H_d(k) = cs_d(k)/omega_d(k)

        mdisc = mdisc + twopi*r_d(k)*sigma_d(k)*dr
   !  If the maximum model disc radius is reached, record and exit the routine

        IF(j>nrad_mod) THEN
           !print*, 'Maximum model radius reached'
           DO l=k,nrannuli
              sigma_d(l) = 0.0
              nu_d(l) = 0.0
              T_d(l) = 0.0
              alpha_d(l) = 0.0
           ENDDO

           rout = r_d(k)
           irout = k

           return
        ENDIF

        ! If disc has disappeared, mark the outer edge and exit the routine
        IF(sigma_d(k) ==0.0) THEN
           !print*, 'Model disc empty!'
           DO l=k,nrannuli
              sigma_d(l) = 0.0
              nu_d(l) = 0.0
              T_d(l) = 0.0
              alpha_d(l) = 0.0
           ENDDO
           
           rout = r_d(k)
           irout = k
           return
        ENDIF

     ENDDO

     ! If not at beginning of file, interpolations must be 2D (generally)
  ELSE
     
     ! Now do 2D interpolations

    !print*, 'Interpolating in 2D'
     DO k=1,nrannuli

        rad = r_d(k)

        j=1
        DO WHILE(r_mod(i,j) < rad .and. j<nrad_mod)
           j=j+1
        ENDDO
     

        ! If at the very start of interpolation in space, interpolate against zero

        IF(j==1) THEN

           ! Sigma
           CALL interpolate_2D(rad,t_interp,sigma_d(k),0.0,0.0 &
                ,sigma_mod(i-1,j), sigma_mod(i,j),0.0,r_mod(i,j), begintime,finishtime)
           ! Nu
           CALL interpolate_2D(rad,t_interp,nu_d(k),0.0,0.0 &
                ,nu_mod(i-1,j), nu_mod(i,j),0.0,r_mod(i,j), begintime,finishtime)
           ! T
           CALL interpolate_2D(rad,t_interp,T_d(k),0.0,0.0 &
                ,Tc_mod(i-1,j), Tc_mod(i,j),0.0,r_mod(i,j), begintime,finishtime)
           ! Alpha
           CALL interpolate_2D(rad,t_interp,alpha_d(k),0.0,0.0 &
                ,alpha_mod(i-1,j), alpha_mod(i,j),0.0,r_mod(i,j), begintime,finishtime)
           ! Tau
           CALL interpolate_2D(rad,t_interp,tau_d(k),0.0,0.0 &
                ,tau_mod(i-1,j), tau_mod(i,j),0.0,r_mod(i,j), begintime,finishtime)

           ! Otherwise interpolate as normal
        ELSE

           ! Sigma
           CALL interpolate_2D(rad,t_interp,sigma_d(k),sigma_mod(i-1,j-1),sigma_mod(i,j-1) &
                ,sigma_mod(i-1,j), sigma_mod(i,j),r_mod(i,j-1),r_mod(i,j), begintime,finishtime)

           ! Nu
           CALL interpolate_2D(rad,t_interp,nu_d(k),nu_mod(i-1,j-1),nu_mod(i,j-1) &
                ,nu_mod(i-1,j), nu_mod(i,j),r_mod(i,j-1),r_mod(i,j), begintime,finishtime)

           ! T
           CALL interpolate_2D(rad,t_interp,T_d(k),Tc_mod(i-1,j-1),Tc_mod(i,j-1) &
                ,Tc_mod(i-1,j), Tc_mod(i,j),r_mod(i,j-1),r_mod(i,j), begintime,finishtime)

           ! Alpha
           CALL interpolate_2D(rad,t_interp,alpha_d(k),alpha_mod(i-1,j-1),alpha_mod(i,j-1) &
                ,alpha_mod(i-1,j), alpha_mod(i,j),r_mod(i,j-1),r_mod(i,j), begintime,finishtime)

           ! Tau
           CALL interpolate_2D(rad,t_interp,tau_d(k),tau_mod(i-1,j-1),tau_mod(i,j-1) &
                ,tau_mod(i-1,j), tau_mod(i,j),r_mod(i,j-1),r_mod(i,j), begintime,finishtime)
        ENDIF
              
        !write(74,*) r_d(k)/udist,sigma_d(k)

        mdisc = mdisc + twopi*r_d(k)*sigma_d(k)*dr
        !write(77,*), r_d(k)/udist,sigma_d(k), mdisc/umass

   !  If the maximum model disc radius is reached, record and exit the routine

        IF(j>nrad_mod) THEN
           !print*, 'Maximum model radius reached'
           DO l=k,nrannuli
              sigma_d(l) = 0.0
              nu_d(l) = 0.0
              T_d(l) = 0.0
              alpha_d(l) = 0.0
           ENDDO

           rout = r_d(k)
           irout = k

           return
        ENDIF

        ! If disc has disappeared, mark the outer edge and exit the routine
        IF(sigma_d(k) ==0.0) THEN
           !print*, 'Model disc empty!'
           DO l=k,nrannuli
              sigma_d(l) = 0.0
              nu_d(l) = 0.0
              T_d(l) = 0.0
              alpha_d(l) = 0.0
           ENDDO
           
           rout = r_d(k)
           irout = k
           return
        ENDIF

     ENDDO
  ENDIF

END SUBROUTINE evolve_disc
