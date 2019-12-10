!! This module define physical and numerical constant used by ROHSA
module mod_constants
  !! This module define physical and numerical constant used by ROHSA
  use, intrinsic :: iso_fortran_env
  implicit none

  integer,  parameter, public :: xp = REAL64
  
  !physical constants FIXME CHECK VALUES
  real(xp), parameter, public :: G       = 6.67408e-8_xp             !gravitationnal cst in cgs
  real(xp), parameter, public :: c       = 2.99792458e10_xp          !speed of light in cgs
  real(xp), parameter, public :: pi      = 4.0_xp*atan(1.0_xp)
  real(xp), parameter, public :: M_sun   = 1.98855e33_xp             !mass of the sun in cgs
  real(xp), parameter, public :: cst_rad = 7.5657308531642009e-15_xp !radiation cst in cgs
  real(xp), parameter, public :: stefan  = (c * cst_rad) / 4.0_xp    !stefan cst in cgs
  real(xp), parameter, public :: kb      = 1.38064852e-16_xp         !boltzmann cst in cgs
  real(xp), parameter, public :: R       = 8.3144598e7_xp            !gas cst in csg = Boltzmann cst over proton mass
  real(xp), parameter, public :: gammag  = 5._xp / 3._xp             !adiabatic coefficient
  real(xp), parameter, public :: thomson = 6.6524587158e-25_xp       !thomson cross-section in cgs 
  real(xp), parameter, public :: h       = 6.62607004e-27_xp         !planck constant in cgs

  type parameters
     integer  :: n_mbb

     integer  :: m 
     integer  :: lstd, ustd
     integer  :: iprint, iprint_init
     integer  :: maxiter, maxiter_init

     real(xp) :: lambda_tau, lambda_beta, lambda_Td
     real(xp) :: lambda_var_tau, lambda_var_beta, lambda_var_Td
     real(xp) :: lambda_stefan

     real(xp) :: tau_init, beta_init, Td_init
     real(xp) :: tau_init_cib, beta_init_cib, Td_init_cib

     real(xp) :: lb_tau, ub_tau
     real(xp) :: lb_beta, ub_beta
     real(xp) :: lb_Td, ub_Td

     real(xp) :: lb_tau_cib, ub_tau_cib
     real(xp) :: lb_beta_cib, ub_beta_cib
     real(xp) :: lb_Td_cib, ub_Td_cib

     real(xp) :: l0
     integer  :: degree
     logical  :: noise, save_grid
     logical  :: cc
     logical  :: ciba

     character(len=512) :: filename, filename_NHI, filename_wavelength, filename_color, fileout, timeout, filename_noise
  end type parameters

  type indata
     real(xp), dimension(:,:,:), allocatable :: cube, std_cube, NHI  
     real(xp), dimension(:), allocatable     :: wavelength  !! wavelength Planck + IRAS
     real(xp), dimension(:,:), allocatable   :: color       !! polynomial coefficient for color correction
  end type indata
 
end module mod_constants
