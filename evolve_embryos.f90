SUBROUTINE evolve_embryos
  ! Routine takes population of nembryos embryos, and tracks their migration through the disc
  ! This routine evolves embryos simultaneously on a constant timestep
  use stardata
  use embryodata
  use eosdata

  implicit none

  integer :: i,j, jwrite,timeup
  real, parameter :: cmig = 1.0 ! Migration efficiency - low c, slow migration
  real :: M_t, r_hill, tmig,t,l_jeans,factor
  real :: vaptime, hillcore,rchoose,orb,rstrip

  ! Debug line - picks an embryo to write data to file on
  jwrite = 3

  ! Initialise all embryos


  DO j=1,nembryo

     embryo(j)%icurrent = embryo(j)%iform

     embryo(j)%R = embryo(j)%R0
     r_hill = embryo(j)%a*(embryo(j)%m/(3.0*mstar))**0.333

     embryo(j)%t_spent = 0.0

     embryo(j)%finished = 0
     vaptime = (embryo(j)%t_cool0)/(p_kap+1)*( (Tvap/embryo(j)%T0)**(p_kap+1) -1)
     !print*, 'Embryo ', j, ' has Vapourisation timescale is ', vaptime/yr

  ENDDO



  t = 0.0     

  ! Evolve the disc until it is no longer self-gravitating
  DO WHILE(q_disc > 0.05)

     ! Check to see if all embryos have finished, and update timestep


     CALL timestep

     !print*, 'dt ',dt/yr
     ! If all embryos have finished, exit the loop
     IF(finishcheck==1) exit

     ! Now begin main embryo loop

     !$OMP PARALLEL &
     !$OMP shared(embryo,nembryo,H_d,r_d) &
     !$OMP shared(t,dt,alpha_d,omega_d,mstar) &
     !$OMP private(i,j,M_t,orb,tmig,l_jeans,r_hill,hillcore)
     !$OMP DO SCHEDULE(runtime)
     DO j=1,nembryo

        i = embryo(j)%icurrent

        ! If embryo finished, skip to the next one
        IF(embryo(j)%finished==1) cycle

        ! Calculate current transition mass at this disc location

        M_t = 2.0*mstar*(H_d(i)/r_d(i))**3

        ! Calculate migration timescale (depending on regime)

        IF(embryo(j)%m<= M_t)THEN
           ! Type I
           IF(sigma_d(i)/=0.0) THEN
              tmig = H_d(i)*mstar/(omega_d(i)*embryo(j)%m*embryo(j)%a)
              tmig = tmig/cmig
           ELSE
              tmig = 1.0e10*yr
           ENDIF
           !  print*, 'Embryo ',j,' undergoes Type I', tmig/yr,sigma_d(i),mstar/umass,H_d(i)/udist,alpha_d(i),&
           !       omega_d(i), embryo(j)%m/mjup, embryo(j)%a/udist
        ELSE
           ! Type II
           IF(alpha_d(i)*omega_d(i)*H_d(i)/=0.0) THEN
              tmig = embryo(j)%a*embryo(j)%a/(alpha_d(i)*omega_d(i)*H_d(i)*H_d(i))
              tmig = tmig/cmig
           ELSE
              tmig = 1.0e10*yr
           ENDIF

           ! print*, 'Embryo ',j,' undergoes Type II', tmig/yr,sigma_d(i), mstar/umass,H_d(i)/udist,alpha_d(i),&
           !      omega_d(i),embryo(j)%m/mjup, embryo(j)%a/udist
        ENDIF

        ! Calculate number of timesteps required to traverse one grid

        embryo(j)%Nsteps = ((dr*tmig/embryo(j)%a)-embryo(j)%t_spent)/dt

        ! Update time spent at this radius
        embryo(j)%t_spent = embryo(j)%t_spent + dt

        ! If core has spent long enough at this radius, then move it inwards one grid

        IF(embryo(j)%Nsteps<=0) THEN
           !print*, 'Moving embryo ', j, t/yr, embryo(j)%Nsteps, embryo(j)%t_spent/yr,tmig/yr, dr/embryo(j)%a
           embryo(j)%a = embryo(j)%a-dr
           embryo(j)%icurrent = embryo(j)%icurrent-1
           embryo(j)%t_spent = 0.0
           IF(embryo(j)%icurrent==1) embryo(j)%finished=1

        ENDIF

        ! Decide whether to accrete gas or not (i.e. check for a gap opening)
        ! Commented out - useful for future implementations
        ! In general, embryos open "cold gaps" as they are quite massive

        ! r_hill = embryo(j)%a*(embryo(j)%m/(3.0*mstar))**0.333

        !IF(r_hill <= H_d(i)) THEN
        !   IF(embryo(j)%m/mstar > sqrt(3.0*pi*alpha_d(i)*(H_d(i)/r_d(i))**5)) THEN
        !      print*, 'Gap opened for embryo ', j
        !   ENDIF
        !ELSE IF(r_hill> H_d(i)) THEN
        !   IF(embryo(j)%m/mstar > pi*alpha_d(i)*(H_d(i)/r_d(i))**2) THEN
        !      print*, 'Cold Gap opened for embryo ', j
        !   ENDIF
        !ENDIF

        ! If not yet in tidal disruption regime, then continue as usual

        !IF(embryo(j)%itidal==0) THEN

        ! Evolve the radius and temperature

        IF(embryo(j)%idiss==0) THEN
           IF(embryo(j)%itidal==0) rchoose = embryo(j)%R0
           !IF(embryo(j)%itidal==1) rchoose = embryo(j)%R 

           embryo(j)%R = rchoose/(1.0 + 2.0*t/embryo(j)%t_cool0)**0.5
           embryo(j)%T = embryo(j)%T0*(1.0+(1.0+p_kap)*t/embryo(j)%t_cool0)**(1.0/(1.0+p_kap))
           embryo(j)%cs = sqrt(gamma*Boltzmannk*embryo(j)%T/(mu*mH))
           !IF(embryo(j)%itidal==0) embryo(j)%rhoc = embryo(j)%M/(4.0*pi*theta_grad*embryo(j)%R**3)
           embryo(j)%rhoc = embryo(j)%M/(4.0*pi*theta_grad*embryo(j)%R**3)           
        ENDIF

        ! If the temperature exceeds Tmelt, ices are evaporated from the solids component
        IF(embryo(j)%T > Tmelt.and.embryo(j)%imelt==0) THEN
           embryo(j)%imelt =1 
           embryo(j)%fg = embryo(j)%fg*0.333
           !print*, "Embryo ",j,"'s ices have melted", embryo(j)%fg*embryo(j)%m/mearth             
        ENDIF

        ! If the temperature of the grains exceeds Tvap, dust vapourised

        IF(embryo(j)%T > Tvap.and.embryo(j)%ivap==0) THEN
           embryo(j)%imelt = 1
           embryo(j)%ivap = 1
           embryo(j)%fg = 0.0
           IF(embryo(j)%ijeans==0) THEN
              embryo(j)%rcore = 0.0
              embryo(j)%rg = 0.0
              embryo(j)%mcore =0.0
           ENDIF
           !	print*, 'Embryo ', j, ' has grains vapourised: core growth ends'           
        ENDIF

        ! If the temperature exceeds the dissociation temperature of H2, second core collapse begins

        IF(embryo(j)%T>Tdiss .and.embryo(j)%idiss==0) THEN
           embryo(j)%imelt=1
           embryo(j)%ivap = 1
           embryo(j)%idiss = 1

           IF(embryo(j)%ijeans==0) THEN 
              embryo(j)%rcore = 0.0
              embryo(j)%rg = 0.0
              embryo(j)%mcore =0.0
           ENDIF

           ! Object is a brown dwarf - give it appropriate radius
           ! Radius derived from linear approximation to Burrows et al (1997)

           IF(embryo(j)%M/Mjup < 70.0 ) THEN
              embryo(j)%R = -2.8e-3*embryo(j)%M/Mjup + 1.0
           ELSE
              embryo(j)%R = 0.01* embryo(j)%M/Mjup + 0.2
           ENDIF
           embryo(j)%R = embryo(j)%R*Rjup
           !    print*, 'Embryo ', j, ' is a brown dwarf'
           ! embryo(j)%finished = 1

        ENDIF

        ! If core not already formed by Jeans instability, AND
        ! If dust not vapourised, continue evolving the dust component

        IF(embryo(j)%ivap==0.and.embryo(j)%ijeans==0) THEN

           ! If the dust has not yet grown to sedimentation, check on this timestep
           IF(embryo(j)%igrown==0) THEN             
              embryo(j)%rg = embryo(j)%R
              ! Once grains have fully grown, mark the core as ready for sedimentation
              IF(t>embryo(j)%t_grow) THEN
                 embryo(j)%igrown = 1
                 embryo(j)%rg0 = embryo(j)%R                    
                 !print*, 'Grains grown for embryo ',j,': ', embryo(j)%scrit, embryo(j)%t_sed0/yr
              ENDIF
           ENDIF

           IF(embryo(j)%igrown==1) THEN
              ! If grains already grown, calculate sedimentation timescales
              embryo(j)%t_sed = embryo(j)%t_sed0

              ! Define remaining mass of dust not stripped
              embryo(j)%mcore = fg*embryo(j)%M

              ! Check if grains travelling above critical velocity for destructive collisions

              IF(embryo(j)%t_sed < embryo(j)%R/vfrag) THEN
                 !   print*, 'Destruction limited sedimentation for embryo ',j
                 embryo(j)%t_sed = embryo(j)%R/vfrag
                 ! print*, t/yr, embryo(j)%t_sed/yr, embryo(j)%t_sed0/yr
              ENDIF





              ! If grains self-gravitate, sedimentation is boosted
              IF(embryo(j)%iself==1) THEN
                 factor = 1.0 - 3.0*(t-embryo(j)%tself)/embryo(j)%t_sed
                 ! If growth rate too fast for this timestep, assume Jeans instability forms the entire core
                 IF(factor< 0.0) THEN
                    embryo(j)%rg = embryo(j)%mcore/(1.333*pi*1.689e6) ! Assume core fully formed at 1 g cm-3 constant density (this circumstance should be avoided by timestep correction)
                 ELSE
                    ! Otherwise grow the core under self-gravity
                    embryo(j)%rg = embryo(j)%rself*factor**0.333
                 ENDIF

              ELSE
                 ! Case where no self gravity acts
                 embryo(j)%rg = embryo(j)%rg*(1.0 - dt/embryo(j)%t_sed)

              ENDIF

              IF(embryo(j)%rg > embryo(j)%R) embryo(j)%rg = embryo(j)%R
              ! print*,'Sedimenting ', t/yr,j, embryo(j)%rg/udist, (embryo(j)%R*fg**0.333)/udist

              ! If Rg small enough that rho(grains) ~ rho(gas), mark core as self-gravitating, set rcore==Rg
              IF(embryo(j)%rg < (fg**0.333)*embryo(j)%R .and. embryo(j)%iself==0.and.embryo(j)%ivap==0) THEN

                 embryo(j)%iself=1                 
                 embryo(j)%tself = t
                 embryo(j)%rself = embryo(j)%rg
                 embryo(j)%rcore = embryo(j)%rg

                 print*, 'Grains self-gravitating for embryo ',j,': ',embryo(j)%rself/udist,embryo(j)%tself/yr

              ENDIF

              ! If rg smaller than jeans length, core collapses

              IF(embryo(j)%iself==1.and.embryo(j)%ijeans==0.and.embryo(j)%ivap==0) THEN

                 !l_jeans = 3.0*Boltzmannk*embryo(j)%T*fg*embryo(j)%m/(4.0*pi*mu*mH*G*&
                 !     (embryo(j)%rg**3)*rho_s*rho_s)
                 !l_jeans = sqrt(l_jeans)

                 l_jeans = embryo(j)%R*fg**0.5
                 !print*, 'Testing for Jeans' ,j, embryo(j)%mcore/mearth, l_jeans/udist, embryo(j)%rg/udist

                 IF(embryo(j)%rg < l_jeans .or. embryo(j)%rg < 0.0) THEN

                    print*, "Embryo ", j,"forms a core ",embryo(j)%idiss
                    embryo(j)%ijeans = 1

                    ! Set up initial core parameters
                    embryo(j)%mcore = embryo(j)%fg*embryo(j)%m
                    embryo(j)%rcore = (3.0*embryo(j)%mcore/(4.0*pi*rho_s))**0.333
                    embryo(j)%rg = embryo(j)%rcore

                 ELSE
                    embryo(j)%mcore = 0.0
                    embryo(j)%rcore = 0.0
                 ENDIF

              ENDIF
           ENDIF
        ENDIF
        !        ENDIF

        ! Recalculate Hill Radius
        r_hill = embryo(j)%a*(embryo(j)%m/(3.0*mstar))**0.333

        IF(embryo(j)%r > r_hill .and.embryo(j)%itidal==0) THEN

           embryo(j)%itidal = 1

           ! print*, 'Tidal disruption begins for embryo ',j, embryo(j)%r/rjup, r_hill/rjup
        ENDIF
        
        IF(embryo(j)%itidal==1.and.embryo(j)%idiss==0) THEN
           ! If we are in the tidal disruption regime, strip off the upper layers of gas and update properties
           ! Calculate the hill radius for the core only

           hillcore = embryo(j)%a*(embryo(j)%mcore/(3.0*mstar))**0.333
           !IF(embryo(j)%ijeans==1) print*, 'CORE ',embryo(j)%mcore/mearth,hillcore/udist,r_hill/udist, &
           !     embryo(j)%a/udist,embryo(j)%mcore/(3.0*mstar)

           ! If no core has formed yet, then strip over a few orbits


           IF(embryo(j)%r> hillcore.and.r_hill < embryo(j)%r) THEN 
              orb = sqrt(G*mstar/embryo(j)%a**3)
              orb = 2.0*pi/orb

              rstrip = max(r_hill,hillcore)

              !print*, 'Stripping ', j,t/yr, embryo(j)%a/udist,embryo(j)%r/rjup, &
              !    r_hill/rjup, (1.0-exp(-dt/orb))

              embryo(j)%r = embryo(j)%r - (embryo(j)%r - rstrip)*(1.0-exp(-dt/(1.0*orb))) ! Gradual depletion of envelope over one orbital period
              !embryo(j)%r = r_hill
              !print*, 'New Radius ', embryo(j)%r/rjup

              ! If a core has formed, then immediately strip layers above hillcore

           ENDIF

           embryo(j)%M = 4.0*pi*embryo(j)%rhoc*embryo(j)%r**3*theta_grad

           ! If the grain radius is larger than the Hill Radius, need to deplete grains as well
           !IF(embryo(j)%rg > r_hill) THEN
           !   embryo(j)%rg = r_hill
           !   embryo(j)%mcore = 1.333*pi*embryo(j)%rg**3
           !   ENDIF


           IF(embryo(j)%m < embryo(j)%mcore) THEN
              embryo(j)%r = embryo(j)%rcore
              embryo(j)%m = embryo(j)%mcore

           ENDIF


        ! IF(t>embryo(j)%t_grow) THEN
        !    print*, j, 'GROWN'
        !    IF(t>embryo(j)%t_sed) print*, j, 'sedimented'
        ! ENDIF

     ENDIF

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

  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL

  ! End of embryo loop       

  t = t+dt
  !IF(t/yr > 100000)STOP

  ! Now evolve the disc for this timestep

  !print*, 'Evolving disc ', t/yr, dt/yr,q_disc
  IF(mstar<mdisc.or.t<dt.or.t/yr>1.0e6) exit 
  ! Evolve the disc

  ! Check for the end of the disc simulation using timeup
  timeup = 0
  CALL evolve_disc(t,dt,timeup)
  IF(timeup==1) exit

  i=0
  ! Recalculate ancillary variables
  q_disc = mdisc/mstar
  !print*, mstar/umass, mdisc/umass, q_disc
  DO WHILE(i < irout)
     i=i+1
     omega_d(i) = sqrt(G*mstar/(r_d(i)*r_d(i)*r_d(i)))

     IF(alpha_d(i)/=0.0.and.sigma_d(i)/=0.0) THEN
        cs_d(i) = sqrt(ABS(nu_d(i)*omega_d(i)/alpha_d(i)))
     ELSE
        cs_d(i) = 0.0
     ENDIF
     H_d(i) = cs_d(i)/omega_d(i)
  ENDDO

ENDDO

! Output data pertaining to all embryos

print*, 'Resulting Objects'
DO j=1,nembryo

  IF(embryo(j)%R > 1.0 .or.embryo(j)%rcore>1.0) THEN   
     WRITE(*,'("Embryo ", I2,": ",7I5,4F18.10)') j, embryo(j)%imelt, embryo(j)%ivap,embryo(j)%idiss, &
          embryo(j)%igrown, embryo(j)%iself, embryo(j)%ijeans, embryo(j)%itidal, &
          embryo(j)%M/Mjup, embryo(j)%R/Rjup, embryo(j)%mcore/mearth, embryo(j)%rcore/rearth

     WRITE(ifinal,'(7I5,5E18.10)') embryo(j)%imelt, embryo(j)%ivap, embryo(j)%idiss, &
          embryo(j)%igrown, embryo(j)%iself, embryo(j)%ijeans, embryo(j)%itidal,&
          embryo(j)%a/udist, embryo(j)%m/mjup, embryo(j)%r/rjup, & 
          embryo(j)%rcore/rearth, embryo(j)%mcore/mearth
  ENDIF

ENDDO

END SUBROUTINE evolve_embryos
