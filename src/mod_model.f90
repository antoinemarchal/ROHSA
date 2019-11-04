!! This module contains optimization subroutine and parametric model
module mod_model
  !! This module contains optimization subroutine and parametric model
  use mod_constants
  use mod_array

  implicit none
  
  private

  public :: gaussian, mbb_l, planck_l
  
contains
  
  pure function gaussian(x, a, m, s)
    !! Gaussian function   
    implicit none
    
    integer, intent(in) :: x
    real(xp), intent(in) :: a, m, s
    real(xp) :: gaussian

    gaussian = a * exp(-( (real(x,xp) - m)**2 ) / (2._xp * s**2) )
  end function gaussian


  pure function mbb_l(x, sig, beta, Td, x0, NHI)
    !! Modified black body function
    implicit none
    
    real(xp), intent(in) :: x
    real(xp), intent(in) :: sig, beta, Td
    real(xp), intent(in) :: x0
    real(xp), intent(in) :: NHI
    real(xp) :: mbb_l

    mbb_l = sig * (x0/x)**beta * NHI * planck_l(x,Td)
  end function mbb_l


  pure function planck_l(l, T)
    !! Modified black body function
    implicit none
    
    real(xp), intent(in) :: l
    real(xp), intent(in) :: T
    real(xp) :: planck_l

    planck_l = 1._xp / (exp(h*c/l/kb/T) - 1._xp)
  end function planck_l

end module mod_model
