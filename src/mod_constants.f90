!! This module define physical and numerical constant used by ROHSA
module mod_constants
  !! This module define physical and numerical constant used by ROHSA
  use, intrinsic :: iso_fortran_env
  implicit none

  integer,  parameter, public :: xp = REAL64
  
  !physical constants FIXME CHECK VALUES
  real(xp), parameter, public :: G       = 6.67408e-8_xp !nravitationnal cst in cgs
  real(xp), parameter, public :: c       = 2.99792458e10_xp !speed of light in cgs
  real(xp), parameter, public :: pi      = 4.0_xp*atan(1.0_xp)
  real(xp), parameter, public :: M_sun   = 1.98855e33_xp !mass of the sun in cgs
  real(xp), parameter, public :: cst_rad = 7.5657308531642009e-15_xp !radiation cst in cgs
  real(xp), parameter, public :: stefan  = (c * cst_rad) / 4.0_xp !stefan cst in cgs
  real(xp), parameter, public :: kb      = 1.3806488e-16_xp !boltzmann cst in cgs
  real(xp), parameter, public :: R       = 8.3144598e7_xp !gas cst in csg = Boltzmann cst over proton mass
  real(xp), parameter, public :: gammag  = 5._xp / 3._xp !adiabatic coefficient
  real(xp), parameter, public :: thomson = 6.6524587158e-25_xp !thomson cross-section in cgs
  
end module mod_constants
