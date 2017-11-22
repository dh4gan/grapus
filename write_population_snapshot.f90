subroutine write_population_snapshot(ifile,nsurvive)

use embryodata
use eosdata, only: udist
use stardata,only: istar
implicit none

integer, intent(inout) :: ifile,nsurvive

integer :: j


DO j=1,nembryo

   
   IF(embryo(j)%R > 1.0 .or.embryo(j)%rcore>1.0 .or. embryo(j)%m/mearth >1.0e-3) THEN

      nsurvive = nsurvive+1     
      WRITE(ifile,'(I6,7I2,1P,7E18.10)') istar, embryo(j)%imelt, &
           embryo(j)%ivap, embryo(j)%idiss, embryo(j)%igrown, &
           embryo(j)%iself, embryo(j)%ijeans, embryo(j)%itidal, &
           embryo(j)%a/udist, embryo(j)%ecc, embryo(j)%inc, &
           embryo(j)%m/mjup, embryo(j)%r/rjup, &
           embryo(j)%rcore/rearth, embryo(j)%mcore/mearth
   ENDIF
   
ENDDO

call flush(ifile)

end subroutine write_population_snapshot
