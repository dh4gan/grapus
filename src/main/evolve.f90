SUBROUTINE evolve
  !****************************************************************************
  ! Subroutine evolves both the embryo and disc population simultaneously 
  ! This routine evolves embryos simultaneously on a shared timestep
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

     call evolve_embryos(t)

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
     !
     ! Recalculate ancillary variables
     ! (omega, cs, H)
     !
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
! Evolution of system complete
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

END SUBROUTINE evolve
