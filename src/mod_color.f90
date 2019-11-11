!! This module contains color correction model
module mod_color
  !! This module contains color correction model
  use mod_constants
  use mod_math
  use mod_array

  implicit none
  
  private

  public :: poly_color, d_poly_color_dx, d_poly_color_dy

contains

  pure function poly_color(a, x, y, degree)
    implicit none

    real(xp), intent(in), dimension(:) :: a
    integer, intent(in) :: degree
    real(xp), intent(in) :: x
    real(xp), intent(in) :: y       

    real(xp) :: poly_color
    integer :: i,j,k
    
    poly_color = 0._xp
    
    k = 1
    do i=1, degree
       do j=1, degree+1-(i-1)
          poly_color = poly_color + (a(k) * x**(i-1) * y**(j-1))
          k = k + 1
       end do
    end do

  end function poly_color


  pure function d_poly_color_dx(a, x, y, degree)
    implicit none

    real(xp), intent(in), dimension(:) :: a
    integer, intent(in) :: degree
    real(xp), intent(in) :: x
    real(xp), intent(in) :: y       

    real(xp) :: d_poly_color_dx
    integer :: i,j,k
    
    d_poly_color_dx = 0._xp
    
    k = 1
    do i=1, degree
       do j=1, degree+1-(i-1)
          d_poly_color_dx = d_poly_color_dx + ((i-1) * a(k) * x**(i-2) * y**(j-1))
          k = k + 1
       end do
    end do

  end function d_poly_color_dx


  pure function d_poly_color_dy(a, x, y, degree)
    implicit none

    real(xp), intent(in), dimension(:) :: a
    integer, intent(in) :: degree
    real(xp), intent(in) :: x
    real(xp), intent(in) :: y       

    real(xp) :: d_poly_color_dy
    integer :: i,j,k
    
    d_poly_color_dy = 0._xp
    
    k = 1
    do i=1, degree
       do j=1, degree+1-(i-1)
          d_poly_color_dy = d_poly_color_dy + ((j-1) * a(k) * x**(i-1) * y**(j-2))
          k = k + 1
       end do
    end do

  end function d_poly_color_dy
  

end module mod_color
