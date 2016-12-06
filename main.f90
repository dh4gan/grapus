PROGRAM TD_synthesis

  ! Program carries out MC exoplanet synthesis using the 
  ! tidal downsizing paradigm of Nayakshin et al 
  ! (see Forgan & Rice 2013, MNRAS, 432, pp 3168-3185 for v1.0)

  ! 6/12/16: Major upgrade - inclusion of RK4 N Body integrator for embryos
  
  use stardata
  use embryodata
  use eosdata

  IMPLICIT NONE

  INTERFACE

     SUBROUTINE generate_star
     END SUBROUTINE generate_star

     SUBROUTINE generate_disc
     END SUBROUTINE generate_disc

     SUBROUTINE generate_embryos
     END SUBROUTINE generate_embryos

     SUBROUTINE evolve_embryos
     END SUBROUTINE evolve_embryos

     SUBROUTINE evolve_disc
     END SUBROUTINE evolve_disc

  END INTERFACE

  !		Read in input parameters

  open(10, file='TD_synth.params', status='unknown')
  read(10,*) prefix                ! Output file prefix
  read(10,*) Nstar                 ! Number of star systems to simulate
  read(10,*) datafilepath          ! File path to location of disc file
  read(10,*) discfile              ! Disc file contains output semi-analytic disc models
  read(10,*) tmax                  ! Maximum runtime of the model (in years)
  read(10,*) iseed                 ! Random number seed
  read(10,*) rin                   ! Disc inner boundary
  read(10,*) dr                    ! Radial separation
  read(10,*) p_kap                 ! opacity index (kappa = kappa_0 T^p_kap)
  read(10,*) c_mig                 ! Migration Factor (tmig = c_mig*tmig)
  read(10,*) c_gap                 ! Gap opening factor (tgap = c_gap*tgap)
  read(10,*) c_collapse            ! Collapse timescale factor
  read(10,*) truncate_disc         ! Truncate disc or not (y/n)
  read(10,*) core_feedback         ! Radiative Feedback of Core Formation? (y/n)
  read(10,*) nbody                 ! Use N-Body integrator? (y/n)
  read(10,*) rtrunc                ! Truncation radius of disc
  read(10,*) rtruncmax             ! Maximum radius of disc models

  close(10)

  rin = rin*udist
  dr = dr*udist

  p_grow = (1.0+p_kap)/(2.5+p_kap)

  iseed = -abs(iseed)

  ! Disc model counter (only used when interpolating from file)
  imodel = 0

  ! Set up output files for initial and final parameters

  istart = 10
  ifinal = 20
  ilog = 30

  inbodylog= 40

  OPEN(istart, file=TRIM(prefix)//'.initial', status='unknown')
  OPEN(ifinal, file=TRIM(prefix)//'.final',status='unknown')
  OPEN(ilog, file=TRIM(prefix)//'.log', status='unknown')

  ! Read in Equation of State for SGD

  CALL eosread

  ! Loop over total number of stars

  DO istar=1,Nstar

     ! Create a T_Tauri star
     CALL generate_star

     ! Generate a self-consistent self-gravitating disc
     CALL generate_disc

     ! Where possible, generate embryos from fragmenting disc
     ! Subroutine outputs initial embryo M, R, a,  to file
     
     CALL generate_embryos 

     ! Evolve these embryos towards the inner disc
     ! Subroutine outputs final M,R,a to file
     IF(nembryo>0) CALL evolve_embryos

     STOP ! Debug line for N Body runs
  ENDDO

  close(istart)
  close(ifinal)
  close(ilog)

END PROGRAM TD_synthesis
