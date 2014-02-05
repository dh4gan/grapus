      module eosdata

!-----------------------------------------------------------------------
! Data module for saving equation of state values
! PJC 22/05/2008
! DHF 09/09/2009
!-----------------------------------------------------------------------

      implicit none
      save

! Units
      integer :: nrhopoints,nUpoints
      real(kind=8),parameter :: Boltzmannk = 1.3807d-16
      real(kind=8),parameter :: mH = 1.6726d-24
      real(kind=8),parameter :: pi = 3.141592658285
	real(kind=8), parameter :: twopi = 2.0*pi
      real(kind=8),parameter :: G = 6.672041d-8
	real(kind=8), parameter :: sigma_SB = 5.675d-5
      real(kind=8),parameter :: udist = 1.50d13
      real(kind=8),parameter :: umass = 1.99d33
	real(kind=8),parameter :: yr = 3.15e7
      real(kind=8),parameter :: utime = sqrt((udist**3)/(G*umass))
      real(kind=8),parameter :: uergg = udist*udist/(utime*utime)

! Arrays
     real,dimension(5) :: gammamuT
      real,allocatable,dimension(:,:,:) :: eostable
	real,allocatable,dimension(:,:) :: cstab

      end module eosdata
