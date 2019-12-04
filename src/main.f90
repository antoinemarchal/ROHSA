!! Regularized Optimization for Hypers-spectral Analysis (ROHSA)
program ROHSA

  use mod_constants
  use mod_start
  use mod_fft
  use mod_inout
  use mod_rohsa
  use mod_optimize
  use mod_color
  
  implicit none

  logical :: noise           !! if false --> STD map computed by ROHSA with lstd and ustd (if true given by the user)
  logical :: save_grid       !! save grid of fitted parameters at each step of the multiresolution process
  logical :: cc              !! if true --> apply colour correction PLANCK+IRAS
  integer :: n_mbb           !! number of gaussian to fit
  integer :: m               !! number of corrections used in the limited memory matrix by LBFGS-B
  integer :: lstd            !! lower bound to compute the standard deviation map of the cube (if noise .eq. false)
  integer :: ustd            !! upper bound to compute the standrad deviation map of the cube (if noise .eq. false)
  integer :: iprint          !! print option 
  integer :: iprint_init     !! print option init
  integer :: maxiter         !! max iteration for L-BFGS-B alogorithm
  integer :: maxiter_init    !! max iteration for L-BFGS-B alogorithm (init mean spectrum)

  real(xp) :: lambda_sig     !! lambda for dust opacity parameter
  real(xp) :: lambda_beta    !! lamnda for spectral emissivity index parameter
  real(xp) :: lambda_Td      !! lambda for dust temperature parameter

  real(xp) :: lambda_var_sig  !! lambda for variance dust opacity parameter
  real(xp) :: lambda_var_beta !! lambda for variance spectral emissivity parameter
  real(xp) :: lambda_var_Td !! lambda for variance spectral emissivity parameter
  real(xp) :: lambda_stefan   !! lambda for variance dust temperature parameter

  real(xp) :: sig_fact_init !! times max siglitude of additional Gaussian

  real(xp) :: Td_init       !! dust opacity init
  real(xp) :: beta_init     !! 
  real(xp) :: sig_init      !!

  real(xp) :: lb_sig        !! lower bound 
  real(xp) :: ub_sig        !! upper bound
  real(xp) :: lb_beta       !! lower bound
  real(xp) :: ub_beta       !! upper bound
  real(xp) :: lb_Td         !! lower bound 
  real(xp) :: ub_Td         !! upper bound

  real(xp) :: l0 !! reference wavelength
  integer :: degree

  character(len=512) :: filename_parameters !! name of the parameters file (default parameters.txt)
  character(len=512) :: filename            !! name of the data file
  character(len=512) :: filename_NHI        !! name of the data file
  character(len=512) :: filename_wavelength !! name of the wavelength file
  character(len=512) :: filename_color      !! name of the wavelength file
  character(len=512) :: fileout             !! name of the output result
  character(len=512) :: timeout             !! name of the output result
  character(len=512) :: filename_noise      !! name of the file with STD map (if noise .eq. true)

  character(len=512) :: filename_fBm="fBm.dat"

  real(xp) :: start, finish

  real(xp), dimension(:,:,:), allocatable :: data        !! initial fits data
  real(xp), dimension(:), allocatable     :: wavelength  !! wavelength Planck + IRAS
  real(xp), dimension(:,:), allocatable   :: color       !! polynomial coefficient for color correction
  real(xp), dimension(:,:,:), allocatable :: std_cube    !! standard deviation cube
  real(xp), dimension(:,:,:), allocatable :: NHI         !! initial fits data NHI

  real(xp), dimension(:,:), allocatable    :: test_fft !! test fft
  complex(xp), dimension(:,:), allocatable :: c_test_fft !! test fft
  complex(xp), dimension(:,:), allocatable :: c_test_fft2 !! test fft
  !
  !
  integer, dimension(3) :: dim_data !! number of frequencies
  integer, dimension(3) :: dim_NHI !! number of NHI maps

  call cpu_time(start)

  !Print header and get filename in argument
  call get_command_argument(1, filename_parameters)

  !Default user parameters
  n_mbb = 2

  lambda_sig = 1._xp
  lambda_beta = 1._xp
  lambda_Td = 1._xp

  lambda_var_sig = 0._xp
  lambda_var_beta = 0._xp
  lambda_var_Td = 0._xp
  lambda_stefan = 1._xp

  sig_fact_init = 2._xp/3._xp

  sig_init = 1._xp
  beta_init = 1.7_xp
  Td_init = 17._xp

  !FIXME VALUE
  lb_sig = 0._xp     
  ub_sig = 100._xp
  lb_beta = 1._xp
  ub_beta = 2.5_xp
  lb_Td = 8.2_xp       		   
  ub_Td = 50._xp       

  maxiter_init = 15000
  maxiter = 800
  m = 10
  noise = .false.
  lstd = 1; ustd = 20
  iprint = -1
  iprint_init = -1
  save_grid = .true.
  cc = .false.
 
  !Read parameters
  call read_parameters(filename_parameters, filename, filename_NHI, filename_wavelength, filename_color, fileout, &
       timeout, filename_noise, n_mbb, lambda_sig, lambda_beta, lambda_Td, lambda_var_sig, lambda_var_beta, &
       lambda_var_Td, lambda_stefan, sig_fact_init, sig_init, beta_init, Td_init, lb_sig, ub_sig, lb_beta, ub_beta, &
       lb_Td, ub_Td, l0, maxiter_init, maxiter, m, noise, lstd, ustd, iprint, iprint_init, save_grid, degree, cc)

  !Call header
  call header()  

  print*, "filename = '",trim(filename),"'"

  !Load data
  call read_cube(filename, data)
  call read_array(filename_wavelength, wavelength)
  call read_map(filename_color, color)

  !Test fft
  ! call read_map(filename_fBm, test_fft)
  ! allocate(c_test_fft(64,64), c_test_fft2(64,64))
  ! c_test_fft = cmplx(test_fft,0._xp,xp)
  ! call cfft2d(64,64,c_test_fft,c_test_fft2)
  ! call icfft2d(64,64,c_test_fft2,c_test_fft)
  ! stop

  if (noise .eqv. .false.) then
     print*, "no .false. option for rohsa-mbb, please provide a rms cube."
     stop
  end if
  call read_cube(filename_noise, std_cube)

  call read_cube(filename_NHI, NHI)

  !Check if n_mbb == number of NHI maps
  dim_NHI = shape(NHI)
  if (n_mbb .ne. dim_NHI(1)) then
     print*, "n_mbb .ne. number of NHI maps / please correct your parameter file."
     stop
  end if
  
  !Call ROHSA subroutine
  call main_rohsa(data, wavelength, std_cube, NHI, fileout, timeout, n_mbb, lambda_sig, lambda_beta, lambda_Td, &
       lambda_var_sig, lambda_var_beta, lambda_var_Td, lambda_stefan, sig_fact_init, sig_init, beta_init, Td_init, &
       lb_sig, ub_sig, lb_beta, ub_beta, lb_Td, ub_Td, l0, maxiter_init, maxiter, m, noise, lstd, ustd, iprint, &
       iprint_init, save_grid, color, degree, cc)  

  call ender()

  call cpu_time(finish)
  print '("Time = ",f6.3," seconds.")',finish-start
   
end program ROHSA
