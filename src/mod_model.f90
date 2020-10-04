!! This module contains optimization subroutine and parametric model
module mod_model
  !! This module contains optimization subroutine and parametric model
  use mod_constants

  implicit none
  
  private

  public :: gaussian, func_q, func_u

contains

  pure function gaussian(x, a, m, s)
    !! Gaussian function   
    implicit none
    
    integer, intent(in) :: x
    real(xp), intent(in) :: a, m, s
    real(xp) :: gaussian

    gaussian = a * exp(-( (real(x,xp) - m)**2 ) / (2._xp * s**2) )
  end function gaussian


  pure function func_q(x, a, mu, w, fn)
    !! Gaussian function   
    implicit none
    
    real(xp), intent(in) :: x
    real(xp), intent(in), dimension(:) :: w
    real(xp), intent(in) :: a, mu
    integer, intent(in) :: fn

    real(xp) :: func_q

    func_q = a / fn * sum(cos(-2_xp*(x-mu)*w**2_xp));
  end function func_q


  pure function func_u(x, a, mu, w, fn)
    !! Gaussian function   
    implicit none
    
    real(xp), intent(in) :: x
    real(xp), intent(in), dimension(:) :: w
    real(xp), intent(in) :: a, mu
    integer, intent(in) :: fn

    real(xp) :: func_u

    func_u = a / fn * sum(sin(-2_xp*(x-mu)*w**2_xp));
  end function func_u

end module mod_model
