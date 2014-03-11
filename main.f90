PROGRAM TD_synthesis

  ! Program carries out MC exoplanet synthesis using the tidal downsizing paradigm of Nayakshin et al

  use stardata
  use embryodata
  use eosdata

  IMPLICIT NONE

  integer :: l,N_runs,start,k
  character(len=100) :: prefix

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
  read(10,*) prefix
  read(10,*) Nstar 
  read(10,*) discfile
  read(10,*) iseed
  read(10,*) rin
  read(10,*) dr
  read(10,*) p_kap
  read(10,*) cmig
  read(10,*) truncate_disc
  read(10,*) rtrunc
  read(10,*) rtruncmax

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

  ENDDO

  close(istart)
  close(ifinal)
  close(ilog)

END PROGRAM TD_synthesis
