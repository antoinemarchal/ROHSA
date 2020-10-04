!! This module read the input user parameters (parameters.txt file / data / noise if true)
module mod_instrument
  !! This module read the input user parameters (parameters.txt file / data / noise if true)
  use mod_constants
  use mod_model
  use mod_read_parameters

  real(xp), dimension(:), allocatable :: wl !wavelength
  real(xp), dimension(:), allocatable :: rm !rm axis
  real(xp), dimension(:), allocatable :: rmsf

  private
  
  public :: compute_constant, wl, rm, rmsf

contains

  subroutine compute_constant()
    implicit none
    
    real(xp), dimension(:), allocatable :: freq
    real(xp), dimension(:), allocatable :: rmsfQ, rmsfU
    
    integer :: i

    !Compute wavelength array
    params%freq_n = floor((params%freq_max - params%freq_min) / params%freq_step + 1)

    allocate(freq(params%freq_n))
    do i=1,params%freq_n
       freq(i) = (params%freq_min + ((i-1)*params%freq_step)) * 1e6_xp
    end do

    allocate(wl(params%freq_n))
    wl = c / 100_xp / freq !attention c divided by 100 because cgs to metric
    print*, 'largest Faraday scale = ', pi/minval(wl**2), 'rad/m2'
    print*, ""

    !Compute RM range and RMSF
    allocate(rm(params%rm_n))
    allocate(rmsfQ(params%rm_n))
    allocate(rmsfU(params%rm_n))
    allocate(rmsf(params%rm_n))

    do i=1,params%rm_n
       rm(i) = params%crval3 + params%cdelt3 * ((i-1) - params%crpix3)
       rmsfQ(i) = func_q(rm(i),1._xp,0._xp,wl,params%freq_n)
       rmsfU(i) = func_u(rm(i),1._xp,0._xp,wl,params%freq_n)
    end do

    rmsf = sqrt(rmsfQ**2._xp + rmsfU**2._xp)
    
  end subroutine compute_constant

end module mod_instrument
