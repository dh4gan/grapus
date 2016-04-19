subroutine migration_timescales(iembryo, migtype, tmig, tgap, tcross)
! Subroutine calculates migration timescales 
! Returns the migration regime 'migtype' also

use embryodata
use stardata
use eosdata, only: yr,udist,umass

implicit none

integer, intent(in) :: iembryo
integer, intent(inout) :: migtype
real, intent(inout) :: tmig,tgap,tcross

integer :: i
real :: vmig, aspectratio,massratio, r_hill, pressure_crit

! Find embryo in disc

i = embryo(iembryo)%icurrent

aspectratio = H_d(i)/embryo(iembryo)%a
massratio = embryo(iembryo)%m/mstar

! Calculate gap opening criteria

! First check Crida et al (2006) pressure criterion

r_hill = embryo(iembryo)%a*(massratio/3.0)**0.333
         
pressure_crit = 0.75*H_d(i)/r_hill + &
     50.0*alpha_d(i)*cs_d(i)*aspectratio/(massratio*r_d(i)*omega_d(i))
  
! Check if gap opening time less than crossing time
! (assume type I migration timescale initially
         
IF(sigma_d(i)>0.0) THEN
   tmig = aspectratio/(omega_d(i)*massratio)
ELSE
   tmig = 1.0e30*yr
ENDIF

tmig = c_mig*tmig

vmig = r_d(i)/tmig

tcross = 2.5*r_hill/vmig
tgap = c_gap*(aspectratio**5)/(omega_d(i)*massratio*massratio)

! Assume gap opens before testing - migration type II
migtype = 2

! Failure to open a gap results in type I migration

! If crossing time too short to open gap, no gap opens
if(tcross < tgap)  migtype = 1
! If pressure criterion not met, gap doesn't open
if(pressure_crit>1.0) migtype = 1

 
! Calculate migration timescale (depending on regime)

IF(migtype==1)THEN
   ! Type I
   IF(sigma_d(i)>0.0) THEN
      tmig = aspectratio/(omega_d(i)*massratio)
   ELSE
      tmig = 1.0e30*yr
   ENDIF
   !  print*, 'Embryo ',j,' undergoes Type I', tmig/yr,sigma_d(i),mstar/umass,H_d(i)/udist,alpha_d(i),&
   !       omega_d(i), embryo(j)%m/mjup, embryo(j)%a/udist
ELSE
   ! Type II
   IF(alpha_d(i)*omega_d(i)*aspectratio >=0.0) THEN
      tmig = 1.0/(alpha_d(i)*omega_d(i)*aspectratio*aspectratio)
   ELSE
      tmig = 1.0e30*yr
   ENDIF
   
   ! print*, 'Embryo ',j,' undergoes Type II', tmig/yr,sigma_d(i), mstar/umass,H_d(i)/udist,alpha_d(i),&
   !      omega_d(i),embryo(j)%m/mjup, embryo(j)%a/udist
ENDIF

! Multiply migration timescale by tunable migration parameter
tmig = tmig*c_mig

vmig = r_d(i)/tmig
tcross = 2.5*r_hill/vmig
if(iembryo==3) print*, 'Migration type: ', migtype, r_hill/udist, &
              dt/yr, tmig/yr, tcross/yr, tgap/yr, pressure_crit

end subroutine migration_timescales
