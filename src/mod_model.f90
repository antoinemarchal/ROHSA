!! This module contains optimization subroutine and parametric model
module mod_model
  !! This module contains optimization subroutine and parametric model
  use mod_constants
  use mod_math
  use mod_array
  use mod_color

  implicit none
  
  private

  public :: gaussian, mbb_l, d_mbb_l_dsig, d_mbb_l_db, d_mbb_l_dT, d_mbbcc_l_dsig, d_mbbcc_l_db, d_mbbcc_l_dT, &
       planck_l, lumi_cst, d_lumi_cst_dsig, d_lumi_cst_dbeta, d_lumi_cst_dTd
  
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


  pure function d_mbb_l_dsig(x, beta, Td, x0, NHI)
    !! Modified black body function derivative sigma
    implicit none
    
    real(xp), intent(in) :: x
    real(xp), intent(in) :: beta, Td
    real(xp), intent(in) :: x0
    real(xp), intent(in) :: NHI
    real(xp) :: d_mbb_l_dsig

    d_mbb_l_dsig = (x0/x)**beta * NHI * planck_l(x,Td)
  end function d_mbb_l_dsig


  pure function d_mbb_l_db(x, sig, beta, Td, x0, NHI)
    !! Modified black body function derivative beta
    implicit none
    
    real(xp), intent(in) :: x
    real(xp), intent(in) :: sig, beta, Td
    real(xp), intent(in) :: x0
    real(xp), intent(in) :: NHI
    real(xp) :: d_mbb_l_db

    d_mbb_l_db = sig * log(x0/x) * (x0/x)**beta * NHI * planck_l(x,Td)
  end function d_mbb_l_db


  pure function d_mbb_l_dT(x, sig, beta, Td, x0, NHI)
    !! Modified black body function derivative sigma
    implicit none
    
    real(xp), intent(in) :: x
    real(xp), intent(in) :: sig, beta, Td
    real(xp), intent(in) :: x0
    real(xp), intent(in) :: NHI
    real(xp) :: d_mbb_l_dT

    d_mbb_l_dT = sig * (x0/x)**beta * NHI * (h*c/x/kb) * exp(h*c/x/kb/Td) / &
         (Td**2._xp * (exp(h*c/x/kb/(Td)) - 1._xp)**2._xp)                         
  end function d_mbb_l_dT


  ! function d_mbbcc_l_dsig(x, beta, Td, x0, NHI, color, degree)
  !   !! Modified black body function derivative sigam with colour correction
  !   implicit none
    
  !   real(xp), intent(in) :: x
  !   real(xp), intent(in) :: beta, Td
  !   real(xp), intent(in) :: x0
  !   real(xp), intent(in) :: NHI
  !   real(xp) :: d_mbbcc_l_dsig
  !   real(xp), intent(in), dimension(:) :: color
  !   integer, intent(in) :: degree

  !   d_mbbcc_l_dsig = d_mbb_l_dsig(x, beta, Td, x0, NHI) / poly_color(color, beta, Td, degree)
  ! end function d_mbbcc_l_dsig


  ! function d_mbbcc_l_db(x, sig, beta, Td, x0, NHI, color, degree)
  !   !! Modified black body function derivative beta with colour correction
  !   implicit none
    
  !   real(xp), intent(in) :: x
  !   real(xp), intent(in) :: sig, beta, Td
  !   real(xp), intent(in) :: x0
  !   real(xp), intent(in) :: NHI
  !   real(xp) :: d_mbbcc_l_db
  !   real(xp), intent(in), dimension(:) :: color
  !   integer, intent(in) :: degree

  !   d_mbbcc_l_db = (d_mbb_l_db(x, sig, beta, Td, x0, NHI) / poly_color(color, beta, Td, degree)) &
  !        - (mbb_l(x, sig, beta, Td, x0, NHI) * d_poly_color_dx(color, beta, Td, degree)) & 
  !        / poly_color(color, beta, Td, degree)**2._xp
  ! end function d_mbbcc_l_db


  ! function d_mbbcc_l_dT(x, sig, beta, Td, x0, NHI, color, degree)
  !   !! Modified black body function derivative temperature with colour correction
  !   implicit none
    
  !   real(xp), intent(in) :: x
  !   real(xp), intent(in) :: sig, beta, Td
  !   real(xp), intent(in) :: x0
  !   real(xp), intent(in) :: NHI
  !   real(xp) :: d_mbbcc_l_dT
  !   real(xp), intent(in), dimension(:) :: color
  !   integer, intent(in) :: degree

  !   d_mbbcc_l_dT = (d_mbb_l_dT(x, sig, beta, Td, x0, NHI) / poly_color(color, beta, Td, degree)) &
  !        - (mbb_l(x, sig, beta, Td, x0, NHI) * d_poly_color_dy(color, beta, Td, degree)) & 
  !        / poly_color(color, beta, Td, degree)**2._xp
  ! end function d_mbbcc_l_dT

  function d_mbbcc_l_dsig(x, beta, Td, x0, NHI, color, degree)
    !! Modified black body function derivative sigam with colour correction
    implicit none
    
    real(xp), intent(in) :: x
    real(xp), intent(in) :: beta, Td
    real(xp), intent(in) :: x0
    real(xp), intent(in) :: NHI
    real(xp) :: d_mbbcc_l_dsig
    real(xp), intent(in), dimension(:) :: color
    integer, intent(in) :: degree

    d_mbbcc_l_dsig = d_mbb_l_dsig(x, beta, Td, x0, NHI) * poly_color(color, beta, Td, degree)
  end function d_mbbcc_l_dsig


  function d_mbbcc_l_db(x, sig, beta, Td, x0, NHI, color, degree)
    !! Modified black body function derivative beta with colour correction
    implicit none
    
    real(xp), intent(in) :: x
    real(xp), intent(in) :: sig, beta, Td
    real(xp), intent(in) :: x0
    real(xp), intent(in) :: NHI
    real(xp) :: d_mbbcc_l_db
    real(xp), intent(in), dimension(:) :: color
    integer, intent(in) :: degree

    d_mbbcc_l_db = (d_mbb_l_db(x, sig, beta, Td, x0, NHI) * poly_color(color, beta, Td, degree)) &
         + (mbb_l(x, sig, beta, Td, x0, NHI) * d_poly_color_dx(color, beta, Td, degree))
  end function d_mbbcc_l_db


  function d_mbbcc_l_dT(x, sig, beta, Td, x0, NHI, color, degree)
    !! Modified black body function derivative temperature with colour correction
    implicit none
    
    real(xp), intent(in) :: x
    real(xp), intent(in) :: sig, beta, Td
    real(xp), intent(in) :: x0
    real(xp), intent(in) :: NHI
    real(xp) :: d_mbbcc_l_dT
    real(xp), intent(in), dimension(:) :: color
    integer, intent(in) :: degree

    d_mbbcc_l_dT = (d_mbb_l_dT(x, sig, beta, Td, x0, NHI) * poly_color(color, beta, Td, degree)) &
         + (mbb_l(x, sig, beta, Td, x0, NHI) * d_poly_color_dy(color, beta, Td, degree)) 
  end function d_mbbcc_l_dT


  pure function planck_l(l, T)
    !! Modified black body function
    implicit none
    
    real(xp), intent(in) :: l
    real(xp), intent(in) :: T
    real(xp) :: planck_l

    planck_l = 1._xp / (exp(h*c/l/kb/T) - 1._xp)
  end function planck_l


  pure function lumi_cst(sig, beta, Td, l0)
    !! LHI constant equation for regularization :
    !! neglect the zeta function zeta(4+beta)/zeta(4)
    !! gamma(4) = 6.
    implicit none
    
    real(xp), intent(in) :: sig, beta, Td
    real(xp), intent(in) :: l0 !! reference wavelength
    real(xp) :: a, b
    real(xp) :: lumi_cst
    real(xp) :: norm

    a = 4._xp + beta
    b = kb * l0 / h / c
    norm = 1.d-5

    ! lumi_cst = Td * sig**q * gamma(a)**q
    lumi_cst = sig * Td**a * gamma(a) * b**beta * norm
    ! print*, lumi_cst
    ! stop
    ! lumi_cst = Td**a * gamma(a) * b**beta * norm
  end function lumi_cst


  pure function d_lumi_cst_dsig(sig, beta, Td, l0)
    !! Derivative lumi_cst function with respect to sigma
    implicit none 

    real(xp), intent(in) :: sig, beta, Td
    real(xp), intent(in) :: l0 !! reference wavelength
    real(xp) :: a, b
    real(xp) :: d_lumi_cst_dsig
    real(xp) :: norm

    a = 4._xp + beta
    b = kb * l0 / h / c
    norm = 1.d-5

    ! d_lumi_cst_dsig = Td * q * sig**(q-1._xp) * gamma(a)**q
    d_lumi_cst_dsig = Td**a * b**beta * gamma(a) * norm
  end function d_lumi_cst_dsig


  function d_lumi_cst_dbeta(sig, beta, Td, l0)
    !! Derivative lumi_cst function with respect to beta
    implicit none 

    real(xp), intent(in) :: sig, beta, Td
    real(xp), intent(in) :: l0 !! reference wavelength
    real(xp) :: a, b
    real(xp) :: d_lumi_cst_dbeta
    real(xp) :: norm
    integer :: ifault

    a = 4._xp + beta
    b = kb * l0 / h / c
    norm = 1.d-5

    ! d_lumi_cst_dbeta = - Td * sig**q * log(sig) / a**2._xp
    ! d_lumi_cst_dbeta = - Td * sig**q * gamma(a)**q * (log(sig) - a*digamma(a, ifault) + log(gamma(a))) / a**2._xp
    d_lumi_cst_dbeta = sig * b**beta * Td**a * gamma(a) * (log(b) + log(Td) + digamma(a, ifault)) * norm
    ! d_lumi_cst_dbeta = b**beta * Td**a * gamma(a) * (log(b) + log(Td) + digamma(a, ifault)) * norm
  end function d_lumi_cst_dbeta


  pure function d_lumi_cst_dTd(sig, beta, Td, l0)
    !! Derivative lumi_cst function with respect to Td
    implicit none 

    real(xp), intent(in) :: sig, beta, Td
    real(xp), intent(in) :: l0 !! reference wavelength
    real(xp) :: a, b
    real(xp) :: d_lumi_cst_dTd
    real(xp) :: norm

    a = 4._xp + beta
    b = kb * l0 / h / c
    norm = 1.d-5

    ! d_lumi_cst_dTd = sig**q * gamma(a)**q
    d_lumi_cst_dTd = sig * a * Td**(a-1._xp) * b**beta * gamma(a) * norm
    ! d_lumi_cst_dTd = a * Td**(a-1._xp) * b**beta * gamma(a) * norm
  end function d_lumi_cst_dTd


end module mod_model
