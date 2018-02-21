SUBROUTINE evolve_embryos
  !****************************************************************************
  ! Subroutine takes population of nembryos embryos, and tracks their evolution
  ! both internally and in the disc
  ! This routine evolves embryos simultaneously on a constant timestep
  !****************************************************************************

  use stardata
  use embryodata
  use eosdata

  implicit none

  integer :: i,j, timeup,nsurvive,ifile
  real :: t, tdump

  tdump = 0.0
  isnap = 0

  ! Initialise all embryos

  do j=1,nembryo

     embryo(j)%icurrent = embryo(j)%iform
     embryo(j)%R = embryo(j)%R0
     embryo(j)%t_spent = 0.0
     embryo(j)%finished = 0
 
  enddo

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
     !$OMP private(i,j)
     !$OMP DO SCHEDULE(runtime)
     DO j=1,nembryo

        ! If embryo finished, skip to the next one
        IF(embryo(j)%finished==1) cycle

        ! Decide whether to accrete gas or not (depending on gap opening)
        ! Commented out - useful for future implementations
        ! In general, embryos open "cold gaps" as they are quite massive
        ! call accrete_gas

        ! Evolve the embryo radius and central temperature

        if(embryo(j)%idiss==0) then
           call evolve_radius(j,t)
           call evolve_temperature(j,t)
        endif

        ! Check for ice melting, grain vapourisation, H2 dissociation

        call check_embryo_thermal_state(j)


        ! If core not already formed by Jeans instability, AND
        ! If dust not vapourised, continue evolving the dust component

        IF(embryo(j)%ivap==0.and.embryo(j)%ijeans==0) THEN


           if(embryo(j)%igrown==0) then
              call calc_grain_growth(j,t)
           else if(embryo(j)%igrown==1) THEN
              call calc_grain_sedimentation(j,t)
           endif

              ! If grain cluster radius smaller than jeans length, collapses to form a core

              if(embryo(j)%iself==1) then
                 call calc_core_formation(j,t)
              endif
           endif        

        ! Check for tidal disruption of the embryo

         call calc_tidal_disruption(j)


         ! Check if embryo has been destroyed
         call check_for_dead_embryos(j)


      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL

  ! End of embryo loop       

  ! Check to see if all embryos have finished, and update global timestep

  CALL timestep

  ! If all embryos have finished, exit the loop
  IF(finishcheck==1) exit

  t = t+dt

  ! Now evolve the disc for this timestep
 
  IF(mstar<mdisc.or.t<dt.or.t/yr>tmax) exit 

  ! Evolve the disc

  ! Check for the end of the disc simulation using timeup
  timeup = 0
  CALL evolve_disc(t,dt,timeup)

  tdump = tdump + dt

  !********************************************************************
  ! If necessary, write a snapshot of the population at this timestep
  !********************************************************************

  if(tdump>tsnap) then
     isnap = isnap +1

     if(debug=='y') write(*,'(A,1P,3e18.4)'), 't, dt, mdisc/mstar: ', t/yr, dt/yr,q_disc
         
     if(isnap.le.nsnaps) then
        ifile = isnapfile+isnap
        nsurvive = 0
               
        OPEN(ifile,file=snapshotfile(isnap),form='formatted',status='old',position='append')

        call write_population_snapshot(ifile,nsurvive)
     endif

     tdump = 0.0
  endif

  ! If debugging, write nbody data to separate output files
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

!********************************
! Evolution of embryos complete
!******************************** 

!
! Print out data pertaining to all surviving embryos
!

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

! Write final population
call write_population_snapshot(ifinal,nsurvive)

write(*,'(A,I1)') 'Number of survivors: ',nsurvive
call flush(ifinal)

! Deallocate nbody arrays ready for next run
if(nbody=='y') call nbody_deallocate_arrays

END SUBROUTINE evolve_embryos
