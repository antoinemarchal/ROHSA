!! Regularized Optimization for Hypers-spectral Analysis (ROHSA)
program ROHSA

  use mod_constants
  use mod_start
  use mod_array
  use mod_fft
  use mod_read_parameters
  use mod_rw_data
  use mod_rohsa
  use mod_optimize
  use mod_color
  use mod_model
  
  implicit none

  character(len=512) :: filename_parameters !! name of the parameters file (default parameters.txt)
  real(xp)           :: start, finish
  integer, dimension(3) :: dim_data !! number of frequencies
  integer, dimension(3) :: dim_NHI !! number of NHI maps

  !
  character(len=512) :: filename_fBm="fBm.dat"

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

  call cpu_time(start)

  !Call header
  call header()  
  !Print header and get filename in argument
  call get_command_argument(1, filename_parameters)
  !Read parameters
  call get_parameters(filename_parameters) 

  !Load data
  call get_data()

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

  ! !Check if n_mbb == number of NHI maps
  ! dim_NHI = shape(NHI)
  ! if (n_mbb .ne. dim_NHI(1)) then
  !    print*, "n_mbb_dust .ne. number of NHI maps / please correct your parameter file."
  !    stop
  ! end if
  
  !Call ROHSA subroutine
  call main_rohsa(data, wavelength, std_cube, NHI, color)  

  call ender()

  call cpu_time(finish)
  print '("Time = ",f6.3," seconds.")',finish-start
   
end program ROHSA
