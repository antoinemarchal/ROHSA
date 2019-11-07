!! This module contains optimization subroutine and parametric model
module mod_model
  !! This module contains optimization subroutine and parametric model
  use mod_constants
  use mod_math
  use mod_array

  implicit none
  
  private

  public :: gaussian, mbb_l, planck_l, lumi_cst, d_lumi_cst_dsig, d_lumi_cst_dbeta, d_lumi_cst_dTd
  
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


  pure function lumi_cst(sig, beta, Td)
    !! LHI constant equation for regularization : sig * T^4 * Q_approx(beta)
    !! neglect the zeta function zeta(4+beta)/zeta(4)
    !! gamma(4) = 6.
    implicit none
    
    real(xp), intent(in) :: sig, beta, Td
    real(xp) :: a, q
    real(xp) :: lumi_cst

    a = 4._xp + beta
    q = 1._xp / a

    lumi_cst = Td * sig**q * gamma(a)**q
  end function lumi_cst


  pure function d_lumi_cst_dsig(sig, beta, Td)
    !! Derivative lumi_cst function with respect to sigma
    implicit none 

    real(xp), intent(in) :: sig, beta, Td
    real(xp) :: a, q
    real(xp) :: d_lumi_cst_dsig

    a = 4._xp + beta
    q = 1._xp / a

    d_lumi_cst_dsig = Td * q * sig**(q-1._xp) * gamma(a)**q
  end function d_lumi_cst_dsig


  function d_lumi_cst_dbeta(sig, beta, Td)
    !! Derivative lumi_cst function with respect to beta
    implicit none 

    real(xp), intent(in) :: sig, beta, Td
    real(xp) :: a, b, q
    real(xp) :: d_lumi_cst_dbeta
    integer :: ifault

    a = 4._xp + beta
    q = 1._xp / a

    ! d_lumi_cst_dbeta = - Td * sig**q * log(sig) / a**2._xp
    d_lumi_cst_dbeta = - Td * sig**q * gamma(a)**q * (log(sig) - a*digamma(a, ifault) + log(gamma(a))) / a**2._xp
  end function d_lumi_cst_dbeta


  pure function d_lumi_cst_dTd(sig, beta, Td)
    !! Derivative lumi_cst function with respect to Td
    implicit none 

    real(xp), intent(in) :: sig, beta, Td
    real(xp) :: a, q
    real(xp) :: d_lumi_cst_dTd

    a = 4._xp + beta
    q = 1._xp / a

    d_lumi_cst_dTd = sig**q * gamma(a)**q
  end function d_lumi_cst_dTd


end module mod_model
