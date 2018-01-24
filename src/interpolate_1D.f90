SUBROUTINE interpolate_1D(x,y,x1,x2,y1,y2)
! Subroutine interpolates in 1 dimension, given
! Bounding independent coordinates x1,x2
! Corresponding dependent variables y1,y2
! Input independent coordinate x
! Desired interpolant y

implicit none

real :: x1,x2,y1,y2,x,y
real :: grad, inter

!print*, 'Calling interpolate 1D'

grad = (y2-y1)/(x2-x1)
inter = y2 - grad*x2

y = grad*x + inter


!print*, x,y,x1,x2,y1,y2, grad,inter

END SUBROUTINE interpolate_1D
