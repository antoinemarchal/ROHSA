!! Regularized Optimization for Hypers-spectral Analysis (ROHSA)
program ROHSA

  use mod_constants
  use mod_start
  use mod_array
  use mod_fft
  use mod_inout
  use mod_rohsa
  use mod_optimize
  use mod_color
  use mod_model
  
  implicit none

  character(len=512) :: filename_parameters !! name of the parameters file (default parameters.txt)
  character(len=512) :: filename_fBm="fBm.dat"

  real(xp) :: start, finish

  real(xp), dimension(:,:,:), allocatable :: data        !! initial fits data
  real(xp), dimension(:), allocatable     :: wavelength  !! wavelength Planck + IRAS
  real(xp), dimension(:,:), allocatable   :: color       !! polynomial coefficient for color correction
  real(xp), dimension(:,:,:), allocatable :: std_cube    !! standard deviation cube
  real(xp), dimension(:,:,:), allocatable :: NHI         !! initial fits data NHI

  real(xp), dimension(:,:), allocatable    :: test_fft !! test fft
  real(xp), dimension(:,:), allocatable    :: tapper !! test fft
  real(xp), dimension(:,:), allocatable    :: butter !! test fft
  real(xp), dimension(:,:), allocatable    :: test_fft_shift !! test fft
  complex(xp), dimension(:,:), allocatable :: c_test_fft !! test fft
  complex(xp), dimension(:,:), allocatable :: c_test_fft2 !! test fft

  real(xp), dimension(:), allocatable :: x
  real(xp), dimension(:), allocatable :: y
  real(xp), dimension(:,:), allocatable :: xx
  real(xp), dimension(:,:), allocatable :: yy
  real(xp), dimension(:,:), allocatable :: kmat
  !
  !
  integer, dimension(3) :: dim_data !! number of frequencies
  integer, dimension(3) :: dim_NHI !! number of NHI maps

  call cpu_time(start)

  !Print header and get filename in argument
  call get_command_argument(1, filename_parameters)
  call get_parameters(filename_parameters)
 
  !Read parameters
  ! call read_parameters(filename_parameters, filename, filename_NHI, filename_wavelength, filename_color, fileout, &
  !      timeout, filename_noise, n_mbb, lambda_tau, lambda_beta, lambda_Td, lambda_var_tau, lambda_var_beta, &
  !      lambda_var_Td, lambda_stefan, tau_init, beta_init, Td_init, tau_init_cib, beta_init_cib, Td_init_cib, lb_tau, &
  !      ub_tau, lb_beta, ub_beta, lb_Td, ub_Td, lb_tau_cib, ub_tau_cib, lb_beta_cib, ub_beta_cib, lb_Td_cib, ub_Td_cib, &
  !      l0, maxiter_init, maxiter, m, noise, lstd, ustd, iprint, iprint_init, save_grid, degree, cc, ciba)

  !Call header
  call header()  

  print*, "filename = '",trim(params%filename),"'"

  ! !Load data
  call read_cube(params%filename, data)
  call read_array(params%filename_wavelength, wavelength)
  call read_map(params%filename_color, color)

  ! !Test fft
  ! call read_map(filename_fBm, test_fft)
  ! allocate(test_fft_shift(64,64))

  ! call shift(test_fft, test_fft_shift)
  ! allocate(c_test_fft(64,64), c_test_fft2(64,64))
  ! c_test_fft = cmplx(test_fft_shift,0._xp,xp)

  ! call cfft2d(64,64,c_test_fft,c_test_fft2)

  ! allocate(kmat(64,64))
  ! call kgrid(64,64,kmat)

  ! allocate(tapper(34,64))
  ! call apodize(tapper, 0.86_xp, 34,64)

  ! call butterworth(butter,kmat,1._xp,1._xp,2._xp)
  ! stop
  ! !

  if (params%noise .eqv. .false.) then
     print*, "no .false. option for rohsa-mbb, please provide a rms cube."
     stop
  end if
  call read_cube(params%filename_noise, std_cube)

  call read_cube(params%filename_NHI, NHI)

  ! !Check if n_mbb == number of NHI maps
  ! dim_NHI = shape(NHI)
  ! if (n_mbb .ne. dim_NHI(1)) then
  !    print*, "n_mbb_dust .ne. number of NHI maps / please correct your parameter file."
  !    stop
  ! end if
  
  !Call ROHSA subroutine
  call main_rohsa(data, wavelength, std_cube, NHI, params%fileout, params%timeout, params%n_mbb, params%lambda_tau, &
       params%lambda_beta, params%lambda_Td, params%lambda_var_tau, params%lambda_var_beta, params%lambda_var_Td, &
       params%lambda_stefan, params%tau_init, params%beta_init, params%Td_init, params%tau_init_cib, params%beta_init_cib, &
       params%Td_init_cib, params%lb_tau, params%ub_tau, params%lb_beta, params%ub_beta, params%lb_Td, params%ub_Td, &
       params%lb_tau_cib, params%ub_tau_cib, params%lb_beta_cib, params%ub_beta_cib, params%lb_Td_cib, params%ub_Td_cib, &
       params%l0, params%maxiter_init, params%maxiter, params%m, params%noise, params%lstd, params%ustd, params%iprint, &
       params%iprint_init, params%save_grid, color, params%degree, params%cc, params%ciba)  

  call ender()

  call cpu_time(finish)
  print '("Time = ",f6.3," seconds.")',finish-start
   
end program ROHSA
