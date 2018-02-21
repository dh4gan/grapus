PROGRAM grapus

  ! Program carries out population synthesis of gravitational instability 
  ! using the tidal downsizing paradigm of Nayakshin et al 
  ! see:
  ! Forgan & Rice 2013, MNRAS, 432, pp 3168-3185 (v1.0)
  ! Forgan et al, 2018, MNRAS, 474, pp 5036-5048  (v2.0: adding N body physics)
  !
  
  use stardata
  use embryodata
  use eosdata

  IMPLICIT NONE

  INTERFACE

     subroutine initial
     end subroutine initial

     SUBROUTINE generate_star
     END SUBROUTINE generate_star

     SUBROUTINE generate_disc
     END SUBROUTINE generate_disc

     SUBROUTINE generate_embryos
     END SUBROUTINE generate_embryos

     SUBROUTINE evolve
     END SUBROUTINE evolve

     SUBROUTINE evolve_disc
     END SUBROUTINE evolve_disc

  END INTERFACE


  ! Initialise run
  call initial

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
     IF(nembryo>0) CALL evolve

  ENDDO

  close(istart)
  close(ifinal)
  close(ilog)

END PROGRAM grapus
