!! This module contains optimization subroutine and parametric model
module mod_model
  !! This module contains optimization subroutine and parametric model
  use mod_constants
  use mod_math
  use mod_array
  use mod_color

  implicit none
  
  private

  public :: gaussian, mbb_l, d_mbb_l_dtau, d_mbb_l_db, d_mbb_l_dT, d_mbbcc_l_dtau, d_mbbcc_l_db, d_mbbcc_l_dT, &
       planck_l, lumi_cst, d_lumi_cst_dtau, d_lumi_cst_dbeta, d_lumi_cst_dTd, butterworth
  
contains

  subroutine butterworth(butter, k, H0, k0, n)
    !! Butterworth filter
    implicit none

    real(xp), intent(inout), dimension(:,:), allocatable :: butter
    real(xp), intent(in), dimension(:,:), allocatable :: k
    real(xp), intent(in) :: H0
    real(xp), intent(in) :: k0
    real(xp), intent(in) :: n

    butter = H0 / (1. + (k/k0)**(2*n))
    
  end subroutine butterworth
  

  pure function gaussian(x, a, m, s)
    !! Gaussian function   
    implicit none
    
    integer, intent(in) :: x
    real(xp), intent(in) :: a, m, s
    real(xp) :: gaussian

    gaussian = a * exp(-( (real(x,xp) - m)**2 ) / (2._xp * s**2) )
  end function gaussian


  pure function mbb_l(x, tau, beta, Td, x0)
    !! Modified black body function
    implicit none
    
    real(xp), intent(in) :: x
    real(xp), intent(in) :: tau, beta, Td
    real(xp), intent(in) :: x0
    real(xp) :: mbb_l

    mbb_l = tau * (x0/x)**beta * planck_l(x,Td)
  end function mbb_l


  pure function d_mbb_l_dtau(x, beta, Td, x0)
    !! Modified black body function derivative tauma
    implicit none
    
    real(xp), intent(in) :: x
    real(xp), intent(in) :: beta, Td
    real(xp), intent(in) :: x0
    real(xp) :: d_mbb_l_dtau

    d_mbb_l_dtau = (x0/x)**beta * planck_l(x,Td)
  end function d_mbb_l_dtau


  pure function d_mbb_l_db(x, tau, beta, Td, x0)
    !! Modified black body function derivative beta
    implicit none
    
    real(xp), intent(in) :: x
    real(xp), intent(in) :: tau, beta, Td
    real(xp), intent(in) :: x0
    real(xp) :: d_mbb_l_db

    d_mbb_l_db = tau * log(x0/x) * (x0/x)**beta * planck_l(x,Td)
  end function d_mbb_l_db


  pure function d_mbb_l_dT(x, tau, beta, Td, x0)
    !! Modified black body function derivative tauma
    implicit none
    
    real(xp), intent(in) :: x
    real(xp), intent(in) :: tau, beta, Td
    real(xp), intent(in) :: x0
    real(xp) :: d_mbb_l_dT

    d_mbb_l_dT = tau * (x0/x)**beta * (h*c/x/kb) * exp(h*c/x/kb/Td) / &
         (Td**2._xp * (exp(h*c/x/kb/(Td)) - 1._xp)**2._xp)                         
  end function d_mbb_l_dT


  ! function d_mbbcc_l_dtau(x, beta, Td, x0, NHI, color, degree)
  !   !! Modified black body function derivative tauam with colour correction
  !   implicit none
    
  !   real(xp), intent(in) :: x
  !   real(xp), intent(in) :: beta, Td
  !   real(xp), intent(in) :: x0
  !   real(xp), intent(in) :: NHI
  !   real(xp) :: d_mbbcc_l_dtau
  !   real(xp), intent(in), dimension(:) :: color
  !   integer, intent(in) :: degree

  !   d_mbbcc_l_dtau = d_mbb_l_dtau(x, beta, Td, x0, NHI) / poly_color(color, beta, Td, degree)
  ! end function d_mbbcc_l_dtau


  ! function d_mbbcc_l_db(x, tau, beta, Td, x0, NHI, color, degree)
  !   !! Modified black body function derivative beta with colour correction
  !   implicit none
    
  !   real(xp), intent(in) :: x
  !   real(xp), intent(in) :: tau, beta, Td
  !   real(xp), intent(in) :: x0
  !   real(xp), intent(in) :: NHI
  !   real(xp) :: d_mbbcc_l_db
  !   real(xp), intent(in), dimension(:) :: color
  !   integer, intent(in) :: degree

  !   d_mbbcc_l_db = (d_mbb_l_db(x, tau, beta, Td, x0, NHI) / poly_color(color, beta, Td, degree)) &
  !        - (mbb_l(x, tau, beta, Td, x0, NHI) * d_poly_color_dx(color, beta, Td, degree)) & 
  !        / poly_color(color, beta, Td, degree)**2._xp
  ! end function d_mbbcc_l_db


  ! function d_mbbcc_l_dT(x, tau, beta, Td, x0, NHI, color, degree)
  !   !! Modified black body function derivative temperature with colour correction
  !   implicit none
    
  !   real(xp), intent(in) :: x
  !   real(xp), intent(in) :: tau, beta, Td
  !   real(xp), intent(in) :: x0
  !   real(xp), intent(in) :: NHI
  !   real(xp) :: d_mbbcc_l_dT
  !   real(xp), intent(in), dimension(:) :: color
  !   integer, intent(in) :: degree

  !   d_mbbcc_l_dT = (d_mbb_l_dT(x, tau, beta, Td, x0, NHI) / poly_color(color, beta, Td, degree)) &
  !        - (mbb_l(x, tau, beta, Td, x0, NHI) * d_poly_color_dy(color, beta, Td, degree)) & 
  !        / poly_color(color, beta, Td, degree)**2._xp
  ! end function d_mbbcc_l_dT

  function d_mbbcc_l_dtau(x, beta, Td, x0, color, degree)
    !! Modified black body function derivative tauam with colour correction
    implicit none
    
    real(xp), intent(in) :: x
    real(xp), intent(in) :: beta, Td
    real(xp), intent(in) :: x0
    real(xp) :: d_mbbcc_l_dtau
    real(xp), intent(in), dimension(:) :: color
    integer, intent(in) :: degree

    d_mbbcc_l_dtau = d_mbb_l_dtau(x, beta, Td, x0) * poly_color(color, beta, Td, degree)
  end function d_mbbcc_l_dtau


  function d_mbbcc_l_db(x, tau, beta, Td, x0, color, degree)
    !! Modified black body function derivative beta with colour correction
    implicit none
    
    real(xp), intent(in) :: x
    real(xp), intent(in) :: tau, beta, Td
    real(xp), intent(in) :: x0
    real(xp) :: d_mbbcc_l_db
    real(xp), intent(in), dimension(:) :: color
    integer, intent(in) :: degree

    d_mbbcc_l_db = (d_mbb_l_db(x, tau, beta, Td, x0) * poly_color(color, beta, Td, degree)) &
         + (mbb_l(x, tau, beta, Td, x0) * d_poly_color_dx(color, beta, Td, degree))
  end function d_mbbcc_l_db


  function d_mbbcc_l_dT(x, tau, beta, Td, x0, color, degree)
    !! Modified black body function derivative temperature with colour correction
    implicit none
    
    real(xp), intent(in) :: x
    real(xp), intent(in) :: tau, beta, Td
    real(xp), intent(in) :: x0
    real(xp) :: d_mbbcc_l_dT
    real(xp), intent(in), dimension(:) :: color
    integer, intent(in) :: degree

    d_mbbcc_l_dT = (d_mbb_l_dT(x, tau, beta, Td, x0) * poly_color(color, beta, Td, degree)) &
         + (mbb_l(x, tau, beta, Td, x0) * d_poly_color_dy(color, beta, Td, degree)) 
  end function d_mbbcc_l_dT


  pure function planck_l(l, T)
    !! Modified black body function
    implicit none
    
    real(xp), intent(in) :: l
    real(xp), intent(in) :: T
    real(xp) :: planck_l

    planck_l = 1._xp / (exp(h*c/l/kb/T) - 1._xp)
  end function planck_l


  pure function lumi_cst(tau, beta, Td, l0)
    !! LHI constant equation for regularization :
    !! neglect the zeta function zeta(4+beta)/zeta(4)
    !! gamma(4) = 6.
    implicit none
    
    real(xp), intent(in) :: tau, beta, Td
    real(xp), intent(in) :: l0 !! reference wavelength
    real(xp) :: a, b
    real(xp) :: lumi_cst
    real(xp) :: norm

    a = 4._xp + beta
    b = kb * l0 / h / c
    norm = 1.d-5

    ! lumi_cst = Td * tau**q * gamma(a)**q
    lumi_cst = tau * Td**a * gamma(a) * b**beta * norm
    ! print*, lumi_cst
    ! stop
    ! lumi_cst = Td**a * gamma(a) * b**beta * norm
  end function lumi_cst


  pure function d_lumi_cst_dtau(beta, Td, l0)
    !! Derivative lumi_cst function with respect to tauma
    implicit none 

    real(xp), intent(in) :: beta, Td
    real(xp), intent(in) :: l0 !! reference wavelength
    real(xp) :: a, b
    real(xp) :: d_lumi_cst_dtau
    real(xp) :: norm

    a = 4._xp + beta
    b = kb * l0 / h / c
    norm = 1.d-5

    ! d_lumi_cst_dtau = Td * q * tau**(q-1._xp) * gamma(a)**q
    d_lumi_cst_dtau = Td**a * b**beta * gamma(a) * norm
  end function d_lumi_cst_dtau


  function d_lumi_cst_dbeta(tau, beta, Td, l0)
    !! Derivative lumi_cst function with respect to beta
    implicit none 

    real(xp), intent(in) :: tau, beta, Td
    real(xp), intent(in) :: l0 !! reference wavelength
    real(xp) :: a, b
    real(xp) :: d_lumi_cst_dbeta
    real(xp) :: norm
    integer :: ifault

    a = 4._xp + beta
    b = kb * l0 / h / c
    norm = 1.d-5

    ! d_lumi_cst_dbeta = - Td * tau**q * log(tau) / a**2._xp
    ! d_lumi_cst_dbeta = - Td * tau**q * gamma(a)**q * (log(tau) - a*digamma(a, ifault) + log(gamma(a))) / a**2._xp
    d_lumi_cst_dbeta = tau * b**beta * Td**a * gamma(a) * (log(b) + log(Td) + digamma(a, ifault)) * norm
    ! d_lumi_cst_dbeta = b**beta * Td**a * gamma(a) * (log(b) + log(Td) + digamma(a, ifault)) * norm
  end function d_lumi_cst_dbeta


  pure function d_lumi_cst_dTd(tau, beta, Td, l0)
    !! Derivative lumi_cst function with respect to Td
    implicit none 

    real(xp), intent(in) :: tau, beta, Td
    real(xp), intent(in) :: l0 !! reference wavelength
    real(xp) :: a, b
    real(xp) :: d_lumi_cst_dTd
    real(xp) :: norm

    a = 4._xp + beta
    b = kb * l0 / h / c
    norm = 1.d-5

    ! d_lumi_cst_dTd = tau**q * gamma(a)**q
    d_lumi_cst_dTd = tau * a * Td**(a-1._xp) * b**beta * gamma(a) * norm
    ! d_lumi_cst_dTd = a * Td**(a-1._xp) * b**beta * gamma(a) * norm
  end function d_lumi_cst_dTd


end module mod_model
