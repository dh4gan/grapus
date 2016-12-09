subroutine nbody_grav_acceleration(position,acceleration)
! Calculates gravitational force between all bodies (brute force)
! given input positions
! IN THE HELIOCENTRIC FRAME (Keeps star at the centre)

! Calculated in units where G=1, M=msol, r=AU, t=2pi units/year
use stardata, only: istar ! Debug line (TODO)
use embryodata
implicit none

real, dimension(3,nbodies),intent(in) :: position
real, dimension(3,nbodies), intent(out) :: acceleration

integer :: ix,ibody, jbody
real :: relpos,magipos,magjpos
real, dimension(3) :: sep

! Skip body 1 as we fix the star at the origin
do ibody=2,nbodies

    ! Start by computing the force from the star
    magipos = sqrt(position(1,ibody)*position(1,ibody) + &
                position(2,ibody)*position(2,ibody)+ &
                position(3,ibody)*position(3,ibody))

    acceleration(:,ibody) = acceleration(:,ibody) - &
            (mass(1)+mass(ibody))*position(:,ibody)/(magipos*magipos*magipos)

    if(istar==3.and.ibody==2) print*, 'ACC STAR',acceleration(:,ibody)
    if(magipos<small) print*, 'AARGH Magipos'
    ! Again skip the central star
    do jbody=2,nbodies

       if(ibody==jbody) cycle ! Don't calculate force on itself

        do ix=1,3
            sep(ix) = position(ix,ibody) - position(ix,jbody)
        enddo

        relpos = sqrt(sep(1)*sep(1) + sep(2)*sep(2)+sep(3)*sep(3) +rsoft*rsoft)

           magjpos = sqrt(position(1,jbody)*position(1,jbody) + &
                position(2,jbody)*position(2,jbody)+ &
                position(3,jbody)*position(3,jbody))

           if(magjpos<small) print*, 'AARGH Magjpos'
           if(relpos<small) print*, 'AARGH relpos'
           do ix=1,3
              acceleration(ix,ibody) = acceleration(ix,ibody) - &
                   mass(jbody)*sep(ix)/(relpos*relpos*relpos) - &
                   mass(jbody)*position(ix,jbody)/(magjpos*magjpos*magjpos)
           enddo

           if(istar==3.and.ibody==2)then
              print*, 'ACC ',jbody,sep,position(1,jbody),acceleration(1,ibody)
           endif
        enddo

enddo

return
end subroutine nbody_grav_acceleration
