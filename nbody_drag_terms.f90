subroutine nbody_drag_terms(position,velocity,acceleration)
! Calculates approximation to 3 drag terms given a migration timescale
! Migration drag at timescale tmig
! Eccentricity damping on timescale tmig/10
! Inclincation damping on timescale tmig/10
! (cf Alibert et al 2013)

! Forces not calculated on body 1 (the star)

use embryodata
implicit none

integer :: ibody,ix
real,dimension(3,nbodies),intent(in) :: position,velocity
real,dimension(3,nbodies),intent(inout) :: acceleration
real,dimension(nbodies) :: vdotr,rmag

! Migration drag first

do ix=1,3
    do ibody=2,nbodies
    if(embryo(ibody-1)%tmig>small) then
    acceleration(ix,ibody) = acceleration(ix,ibody)-velocity(ix,ibody)/(2.0*embryo(ibody-1)%tmig)
    endif
    enddo
enddo

! Eccentricity damping
vdotr(:) =  0.0
do ix=1,3
   vdotr(:) = vdotr(:)+ position(ix,:)*velocity(ix,:)
enddo

rmag(:) = 0.0

rmag(:) = position(1,:)*position(1,:) + &
     position(2,:)*position(2,:) + &
     position(3,:)*position(3,:)

do ix=1,3

do ibody=2,nbodies
 if (embryo(ibody-1)%tmig*rmag(ibody)>small) then
acceleration(ix,:) = acceleration(ix,:) - 2.0*dampfac*vdotr(:)*position(ix,:)/(rmag(ibody)*embryo(ibody-1)%tmig)
endif
enddo
enddo

! Inclination damping

do ibody=2,nbodies
if(embryo(ibody-1)%tmig>small) then
    acceleration(3,:) = acceleration(3,:) - 2.0*dampfac*vel(3,:)/embryo(ibody-1)%tmig
endif
end do

end subroutine nbody_drag_terms
