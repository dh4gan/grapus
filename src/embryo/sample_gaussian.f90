
!-------------------------------------------------------------------------------
!+
! Generates Gaussian samples N(mu,sigma^2) using the Box-Muller Transform (in polar form)
!+
!-------------------------------------------------------------------------------
subroutine sample_gaussian(x,mu,sigma)
use stardata, only: ran2, iseed
implicit none

real, intent(inout) :: x
real, intent(in) :: mu, sigma

real :: w,y1,y2

w = 10.0

do while(w>1.0)
   y1 = 2.0*ran2(iseed) - 1.0
   y2 = 2.0*ran2(iseed) - 1.0

   w = y1*y1 + y2*y2
enddo

! Gaussian variate of N(0,1)
x = y1*sqrt(-log(w)/w)

! Convert to N(mu, sigma^2)
x = x*sigma + mu

end subroutine sample_gaussian
