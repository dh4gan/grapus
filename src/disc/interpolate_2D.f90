SUBROUTINE interpolate_2D(x,y,z,z1_y1,z1_y2, z2_y1,z2_y2,x1,x2,y1,y2)
! Subroutine interpolates in two dimensions (x,y)
implicit none

real :: z1_y1, z2_y1 ! Values of x at y1
real :: z1_y2, z2_y2 ! Values of x at y2
real :: x1,x2,y1,y2,x,y,z, grad,inter
real :: grad1,grad2,inter1,inter2, z1,z2

!print*, 'Calling interpolate 2D'
!print*, x1,x2,y1,y2
!print*, z1_y1, z2_y1
!print*, z1_y2,z2_y2

grad1 = (z1_y2 - z1_y1)/(y2-y1)
inter1 = z1_y2-grad1*y2
z1 = grad1*y + inter1

grad2 = (z2_y2 - z2_y1)/(y2-y1)
inter2 = z2_y2-grad2*y2
z2 = grad2*y + inter2

grad = (z2-z1)/(x2-x1)
inter = z2 - grad*x2
z = grad*x + inter 

!write(74,*)x/1.496e13,y,z,z1,z1_y1,z1_y2, z2,z2_y1,z2_y2, grad,inter

END SUBROUTINE interpolate_2D
