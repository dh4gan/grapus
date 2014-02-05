SUBROUTINE generate_star
! 
use stardata
use eosdata
use embryodata

! For now, simply postulate star mass (in solar masses)
! (TODO - stellar mass function?)

mstar = (0.5 + ran2(iseed))*umass

END SUBROUTINE generate_star
