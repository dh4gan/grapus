subroutine migration_timescales
! Subroutine calculates migration timescales 
! Returns the migration regime 'migtype' also

use embryodata
use stardata
use eosdata, only: yr,udist,umass

implicit none

integer :: i,j,migtype
real :: vmig, aspectratio,massratio, r_hill, pressure_crit


do j=1,nembryo

! Find embryo in disc

i = embryo(j)%icurrent

if(i>irout) then
    embryo(j)%tcross = 1.0e30*yr
    embryo(j)%tgap = 1.0e30*yr
    embryo(j)%tmig = 1.0e30*yr
    cycle

else

    aspectratio = H_d(i)/embryo(j)%a
    massratio = embryo(j)%m/mstar

    ! Calculate gap opening criteria

    ! First check Crida et al (2006) pressure criterion

    r_hill = embryo(j)%a*(massratio/3.0)**0.333
         
    pressure_crit = 0.75*H_d(i)/r_hill + &
     50.0*alpha_d(i)*cs_d(i)*aspectratio/(massratio*r_d(i)*omega_d(i))
  
    ! Check if gap opening time less than crossing time
    ! (assume type I migration timescale initially
         
    IF(sigma_d(i)>0.0) THEN
    embryo(j)%tmig = c_mig*aspectratio/(omega_d(i)*massratio)
    ELSE
    embryo(j)%tmig = 1.0e30*yr
    ENDIF

    vmig = r_d(i)/embryo(j)%tmig

    embryo(j)%tcross = 2.5*r_hill/vmig
    embryo(j)%tgap = c_gap*(aspectratio**5)/(omega_d(i)*massratio*massratio)

    ! Assume gap opens before testing - migration type II
    migtype = 2

    ! Failure to open a gap results in type I migration

    ! If crossing time too short to open gap, no gap opens
    if( embryo(j)%tcross <  embryo(j)%tgap)  migtype = 1
    ! If pressure criterion not met, gap doesn't open
    if(pressure_crit>1.0) migtype = 1
 
    ! Calculate migration timescale (depending on regime)

    IF(migtype==1)THEN
   ! Type I - Multiply migration timescale by tunable migration parameter
        IF(sigma_d(i)>0.0) THEN
        embryo(j)%tmig = c_mig*aspectratio/(omega_d(i)*massratio)
        ELSE
        embryo(j)%tmig = 1.0e30*yr
        ENDIF
   !  print*, 'Embryo ',j,' undergoes Type I', tmig/yr,sigma_d(i),mstar/umass,H_d(i)/udist,alpha_d(i),&
   !       omega_d(i), embryo(j)%m/mjup, embryo(j)%a/udist
    ELSE
        ! Type II - tunable migration parameter not effective here
        IF(alpha_d(i)*omega_d(i)*aspectratio >=0.0) THEN
        embryo(j)%tmig = 2.0/(3.0*alpha_d(i)*omega_d(i)*aspectratio*aspectratio)
        ELSE
        embryo(j)%tmig = 1.0e30*yr
        ENDIF
   
        ! print*, 'Embryo ',j,' undergoes Type II', tmig/yr,sigma_d(i), mstar/umass,H_d(i)/udist,alpha_d(i),&
   !      omega_d(i),embryo(j)%m/mjup, embryo(j)%a/udist
    ENDIF


    if(j==1) print*, 'Migration timescales: ', j, migtype,  embryo(j)%a/udist,embryo(j)%m/mjup, &
    pressure_crit,  embryo(j)%tmig/yr,  embryo(j)%tcross/yr,  embryo(j)%tgap/yr

endif

enddo

end subroutine migration_timescales
