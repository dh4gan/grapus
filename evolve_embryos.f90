SUBROUTINE evolve_embryos
  ! Routine takes population of nembryos embryos, and tracks their migration through the disc
  ! This routine evolves embryos simultaneously on a constant timestep
  use stardata
  use embryodata
  use eosdata

  implicit none

  integer :: i,j, timeup,nsurvive,ifile
  real :: M_t, r_hill, t,l_jeans,factor
  real :: vaptime, hillcore,rchoose,orb,rstrip
  real :: core_energy,embryo_energy

  real :: tdump


  tsnap = 1.0e3*yr
  tdump = 0.0

  ! Initialise all embryos


  DO j=1,nembryo

     embryo(j)%icurrent = embryo(j)%iform

     embryo(j)%R = embryo(j)%R0
     !r_hill = embryo(j)%a*(embryo(j)%m/(3.0*mstar))**0.333

     embryo(j)%t_spent = 0.0

     embryo(j)%finished = 0
     vaptime = (embryo(j)%t_cool0)/(p_kap+1)*( (Tvap/embryo(j)%T0)**(p_kap+1) -1)
     !print*, 'Embryo ', j, ' has Vapourisation timescale is ', vaptime/yr

  ENDDO

  t = 0.0     

  ! Evolve the disc until it is no longer self-gravitating
  DO WHILE(q_disc > 0.05)

     !*************************************
     !1. Compute the motion of the embryos
     !*************************************

     ! Calculate migration timescales and gap opening criteria
     call migration_timescales

     ! Move embryos (either analytically or via N Body integration)
     call move_embryos

     !*************************************
     !2. Compute the internal structure of the embryos
     !*************************************

     !$OMP PARALLEL &
     !$OMP shared(embryo,nembryo,H_d,r_d) &
     !$OMP shared(t,dt,alpha_d,omega_d,mstar) &
     !$OMP private(i,j,M_t,orb,l_jeans,r_hill,hillcore) &
     !$OMP private(core_energy,embryo_energy) 
     !$OMP DO SCHEDULE(runtime)
     DO j=1,nembryo

        ! If embryo finished, skip to the next one
        IF(embryo(j)%finished==1) cycle

        ! Decide whether to accrete gas or not (i.e. check for a gap opening)
        ! Commented out - useful for future implementations
        ! In general, embryos open "cold gaps" as they are quite massive

!         r_hill = embryo(j)%a*(embryo(j)%m/(3.0*mstar))**0.333

        !IF(r_hill <= H_d(i)) THEN
        !   IF(embryo(j)%m/mstar > sqrt(3.0*pi*alpha_d(i)*(H_d(i)/r_d(i))**5)) THEN
        !      print*, 'Gap opened for embryo ', j
        !   ENDIF
        !ELSE IF(r_hill> H_d(i)) THEN
        !   IF(embryo(j)%m/mstar > pi*alpha_d(i)*(H_d(i)/r_d(i))**2) THEN
        !      print*, 'Cold Gap opened for embryo ', j
        !   ENDIF
        !ENDIF

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

                    if(core_feedback=='y') then
                       ! evaluate radiative feedback of core formation on envelope OpenmP FLAGS!
                       core_energy = embryo(j)%mcore*embryo(j)%mcore/embryo(j)%rcore
                       embryo_energy = embryo(j)%m*embryo(j)%m/embryo(j)%r

                       ! If energy released in core formation greater than envelope binding energy, destroy the embryo
                       if(core_energy>embryo_energy) then
                          embryo(j)%r =0.0
                          embryo(j)%M = 0.0                       
                       endif
                    endif

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

        embryo(j)%itidal = 0
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

  ! Check to see if all embryos have finished, and update global timestep

  CALL timestep

  ! If all embryos have finished, exit the loop
  IF(finishcheck==1) exit

  t = t+dt
  !IF(t/yr > 100000)STOP

  ! Now evolve the disc for this timestep

 
  IF(mstar<mdisc.or.t<dt.or.t/yr>tmax) exit 
  ! Evolve the disc

  ! Check for the end of the disc simulation using timeup
  timeup = 0
  CALL evolve_disc(t,dt,timeup)

  tdump = tdump + dt
  if(tdump>tsnap) then
     isnap = isnap +1
     if(debug=='y') write(*,'(A,1P,3e18.4)'), 't, dt, mdisc/mstar: ', t/yr, dt/yr,q_disc
         
     if(isnap<nsnaps) then
        ifile = isnapfile+isnap
        nsurvive = 0
               
        OPEN(ifile,file=snapshotfile(isnap),position='append')

        call write_population_snapshot(ifile,nsurvive)
     endif

     tdump = 0.0
  endif

  if(nbody=='y' .and.debug=='y') call nbody_output(t) ! Debug line - check nbody outputs

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

if(debug=='y') then
   print*, 'Resulting ', nembryo, ' Objects:'
   print*, '(istar, iproperties, a,e,i,m,r,rcore,mcore)'

   do j=1,nembryo
      WRITE(*,'("Run ",I5, ", Embryo ", I2,": ",7I2,7F12.4)') istar, j, &
           embryo(j)%imelt, &
           embryo(j)%ivap,embryo(j)%idiss, embryo(j)%igrown, &
           embryo(j)%iself, embryo(j)%ijeans, embryo(j)%itidal, &
           embryo(j)%a/udist, embryo(j)%ecc, embryo(j)%inc, &
           embryo(j)%M/Mjup, embryo(j)%R/Rjup, &
           embryo(j)%mcore/mearth, embryo(j)%rcore/rearth
   enddo
endif


nsurvive = 0

call write_population_snapshot(ifinal,nsurvive)

!!$DO j=1,nembryo
!!$   
!!$  IF(embryo(j)%R > 1.0 .or.embryo(j)%rcore>1.0 .or. embryo(j)%m/mearth >1.0e-3) THEN
!!$
!!$     nsurvive = nsurvive+1
!!$     if(debug=='y') then
!!$        WRITE(*,'("Run ",I5, ", Embryo ", I2,": ",7I2,7F12.4)') istar, j, embryo(j)%imelt, &
!!$             embryo(j)%ivap,embryo(j)%idiss, embryo(j)%igrown, &
!!$             embryo(j)%iself, embryo(j)%ijeans, embryo(j)%itidal, &
!!$             embryo(j)%a/udist, embryo(j)%ecc, embryo(j)%inc, &
!!$             embryo(j)%M/Mjup, embryo(j)%R/Rjup, &
!!$             embryo(j)%mcore/mearth, embryo(j)%rcore/rearth
!!$     endif
!!$
!!$     WRITE(ifinal,'(I6,7I2,1P,7E18.10)') istar, embryo(j)%imelt, embryo(j)%ivap, embryo(j)%idiss, &
!!$          embryo(j)%igrown, embryo(j)%iself, embryo(j)%ijeans, embryo(j)%itidal,&
!!$          embryo(j)%a/udist, embryo(j)%ecc, embryo(j)%inc, embryo(j)%m/mjup, embryo(j)%r/rjup, &
!!$          embryo(j)%rcore/rearth, embryo(j)%mcore/mearth
!!$  ENDIF
!!$
!!$ENDDO

write(*,'(A,I1)') 'Number of survivors: ',nsurvive
call flush(ifinal)

if(nbody=='y') call nbody_deallocate_arrays

END SUBROUTINE evolve_embryos
