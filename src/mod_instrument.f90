!! This module read the input user parameters (parameters.txt file / data / noise if true)
module mod_instrument
  !! This module read the input user parameters (parameters.txt file / data / noise if true)
  use mod_constants
  use mod_model
  use mod_read_parameters
  use mod_functions
  use mod_minimize

  ! real(xp), dimension(:), allocatable :: coeff

  private
  
  public :: compute_constant

contains

  subroutine compute_constant()
    implicit none

    real(xp), dimension(:), allocatable :: rmsf_norm
    
    allocate(rmsf_norm(params%rm_n))
    allocate(coeff(2*params%n_rmsf))

    ! coeff = 0.1_xp
    rmsf_norm = 1._xp * rmsf
    ! call fit_rmsf(coeff, rmsf_norm, params%rm_n)
    stop
    
  end subroutine compute_constant

end module mod_instrument
