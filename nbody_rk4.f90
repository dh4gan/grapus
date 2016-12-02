subroutine nbody_rk4
! This subroutine drives the N-Body integration over a single timestep
! Integration is done in a separate set of arrays for speed
! Orbital elements are stored in GE_embryo type

! Do integration

call integrate(dt_nbody,pos,vel,newpos,newvel)

call nbody_timestep(newpos,newvel)

if(maxerror>tolerance) cycle

pos = newpos
vel = newvel

end subroutine nbody_rk4

! Pour data back into dummy arrays


end subroutine nbody_rk4

subroutine setup_nbody_arrays

pos(:,1) = 0.0
vel(:,1) = 0.0
acc(:,1) = 0.0

! Thought - is it even worth putting these arrays into GE_embryo?
pos(:,ibody) = embryo(ibody)%pos(:)
vel(:,ibody) = embryo(ibody)%vel(:)
acc(:,ibody) = embryo(ibody)%acc(:)


end subroutine setup_nbody_arrays
