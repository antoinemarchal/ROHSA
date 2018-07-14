!! Regularized Optimization for Hypers-spectral Analysis (ROHSA)
program ROHSA

  use mod_constants
  use mod_read
  use mod_array
  use mod_functions
  use mod_start
  use mod_optimize
  
  implicit none

  real(xp), dimension(:,:,:), allocatable :: data        !! initial fits data
  real(xp), dimension(:,:,:), allocatable :: cube        !! reshape data with nside --> cube
  real(xp), dimension(:,:,:), allocatable :: cube_mean   !! mean cube over spatial axis
  real(xp), dimension(:,:,:), allocatable :: fit_params  !! parameters to optimize with cube mean at each iteration
  real(xp), dimension(:,:,:), allocatable :: grid_params !! parameters to optimize at final step (dim of initial cube)
  real(xp), dimension(:,:), allocatable :: std_map       !! standard deviation map fo the cube computed by ROHSA with lb and ub
  real(xp), dimension(:,:), allocatable :: std_cube      !! standard deviation map fo the cube is given by the user 

  integer, dimension(3) :: dim_data !! dimension of original data
  integer, dimension(3) :: dim_cube !! dimension of reshape cube

  real(xp), dimension(:,:), allocatable :: kernel !! convolution kernel 

  integer :: ios=0 !! ios integer
  integer :: i     !! loop index
  integer :: j     !! loop index
  integer :: k     !! loop index

  logical :: noise        !! if false --> STD map computed by ROHSA with lstd and ustd (if true given by the user)
  logical :: regul        !! if true --> activate regulation
  integer :: nside        !! size of the reshaped data \(2^{nside}\)
  integer :: n            !! loop index
  integer :: power        !! loop index
  integer :: n_gauss      !! number of gaussian to fit
  integer :: n_gauss_add  !! number of gaussian to add at each step
  integer :: m            !! number of corrections used in the limited memory matrix by LBFGS-B
  integer :: lstd         !! lower bound to compute the standard deviation map of the cube (if noise .eq. false)
  integer :: ustd         !! upper bound to compute the standrad deviation map of the cube (if noise .eq. false)
  integer :: iprint       !! print option 
  integer :: iprint_init  !! print option init
  integer :: maxiter      !! max iteration for L-BFGS-B alogorithm
  integer :: maxiter_init !! max iteration for L-BFGS-B alogorithm (init mean spectrum)
  real(xp) :: lambda_amp  !! lambda for amplitude parameter
  real(xp) :: lambda_mu   !! lamnda for mean position parameter
  real(xp) :: lambda_sig  !! lambda for dispersion parameter
  real(xp) :: lambda_var_sig  !! lambda for variance dispersion parameter

  character(len=512) :: filename_parameters !! name of the parameters file (default parameters.txt)
  character(len=512) :: filename            !! name of the data file
  character(len=512) :: fileout             !! name of the output result
  character(len=512) :: filename_noise      !! name of the file with STD map (if noise .eq. true)

  !Print header and get filename in argument
  call header()
  call get_command_argument(1, filename_parameters)
  print*, ""
  
!   call read_fits()

  !Default user parameters
  n_gauss = 1
  n_gauss_add = 1
  lambda_amp = 1._xp
  lambda_mu = 1._xp
  lambda_sig = 1._xp
  lambda_var_sig = 1._xp
  maxiter_init = 15000
  maxiter = 800
  m = 10
  noise = .false.
  regul = .true.
  lstd = 0; ustd = 20
  iprint = -1
  iprint_init = -1

  call read_parameters(filename_parameters, filename, fileout, filename_noise, n_gauss, n_gauss_add, &
       lambda_amp, lambda_mu, lambda_sig, lambda_var_sig, maxiter_init, maxiter, m, noise, regul, &
       lstd, ustd, iprint, iprint_init)
    
  print*, "filename = '",trim(filename),"'"
  print*, "fileout = '",trim(fileout),"'"

  print*,
  print*, "______Parameters_____"
  print*,
  print*, "n_gauss = ", n_gauss
  print*, "n_gauss_add = ", n_gauss_add
  print*, "lambda_amp = ", lambda_amp
  print*, "lambda_mu = ", lambda_mu
  print*, "lambda_sig = ", lambda_sig
  print*, "lambda_var_sig = ", lambda_var_sig
  print*, "maxiter_itit = ", maxiter_init
  print*, "maxiter = ", maxiter
  print*, "lstd = ", lstd
  print*, "ustd = ", ustd
  print*, "noise = ", noise
  print*, "regul = ", regul
  print*,
  
  allocate(kernel(3, 3))

  kernel(1,1) = 0._xp
  kernel(1,2) = -0.25_xp
  kernel(1,3) = 0._xp
  kernel(2,1) = -0.25_xp
  kernel(2,2) = 1._xp
  kernel(2,3) = -0.25_xp
  kernel(3,1) = 0._xp
  kernel(3,2) = -0.25_xp
  kernel(3,3) = 0._xp
  
  !Load data
  call read_cube(filename, data)
  
  if (noise .eqv. .true.) then
     if (filename_noise == " ") then
        print*, "--> noise = .true. (no input rms map)"
     end if
     call read_map(filename_noise, std_cube)
  end if
  
  write(*,*) ""
  write(*,*) "opening file and reading data"

  dim_data = shape(data)

  write(*,*) "dim_v, dim_y, dim_x = ", dim_data
  write(*,*) ""
  write(*,*) "number of los = ", dim_data(2)*dim_data(3)

  nside = dim2nside(dim_data)

  write(*,*) "nside = ", nside

  call dim_data2dim_cube(nside, dim_data, dim_cube)

  !Allocate moemory for cube
  allocate(cube(dim_cube(1), dim_cube(2), dim_cube(3)))

  !Reshape the data (new cube of size nside)
  print*,
  write(*,*) "Reshape cube, new dimensions :"
  write(*,*) "dim_v, dim_y, dim_x = ", dim_cube
  print*, 

  call reshape_up(data, cube, dim_data, dim_cube)

  !Allocate memory for fit_params array
  allocate(fit_params(3*n_gauss, 1, 1))

  print*, "                    Start iteration"
  print*,
  
  !Start iteration
  do n=0,nside
     power = 2**n
     
     allocate(cube_mean(dim_cube(1), power, power))
     
     call mean_array(power, cube, cube_mean)
     
     if (n == 0) then
        print*, "Init mean spectrum"        
        call init_spectrum(n_gauss, fit_params(:,1,1), dim_cube(1), cube_mean(:,1,1), maxiter_init, m, iprint_init)
     end if
     
     call go_up_level(fit_params)
     write(*,*) ""
     write(*,*) "Update parameters level ", n, ">", power

     if (regul .eqv. .false.) then
        call upgrade(cube_mean, fit_params, power, n_gauss, dim_cube(1), maxiter, m, iprint)
     end if

     if (regul .eqv. .true.) then
        if (n == 0) then
           call upgrade(cube_mean, fit_params, power, n_gauss, dim_cube(1), maxiter, m, iprint)
        end if
        
        if (n > 0 .and. n < nside) then
           allocate(std_map(power, power))
           
           if (noise .eqv. .true.) then
              call mean_map(power, std_cube, std_map)           
           else
              call set_stdmap(std_map, cube_mean, lstd, ustd)
           end if
           
           call update(cube_mean, fit_params, n_gauss, dim_cube(1), power, power, lambda_amp, lambda_mu, lambda_sig, &
                lambda_var_sig, maxiter, m, kernel, iprint, std_map)        
           deallocate(std_map)
        end if
     end if
     
     deallocate(cube_mean)
  enddo

  print*,
  write(*,*) "Reshape cube, restore initial dimensions :"
  write(*,*) "dim_v, dim_y, dim_x = ", dim_data

  allocate(grid_params(3*n_gauss, dim_data(2), dim_data(3)))

  call reshape_down(fit_params, grid_params,  (/ 3*n_gauss, dim_cube(2), dim_cube(3)/), (/ 3*n_gauss, dim_data(2), dim_data(3)/))

  print*,
  print*, "Update last level ..."
  print*,

  allocate(std_map(dim_data(2), dim_data(3)))
  
  if (noise .eqv. .true.) then
     std_map = std_cube
  else   
     call set_stdmap(std_map, data, lstd, ustd)
  end if

  if (regul .eqv. .true.) then
     call update(data, grid_params, n_gauss, dim_data(1), dim_data(2), dim_data(3), lambda_amp, lambda_mu, lambda_sig, &
          lambda_var_sig, maxiter, m, kernel, iprint, std_map)
  end if

  print*,
  print*, "_____ Write output file _____"
  print*, 

  ! Open file
  open(unit=12, file=fileout, action="write", iostat=ios)
  if (ios /= 0) stop "opening file error"
  
  ! Read cube dimension and compute the number of line
  write(12,fmt=*) "# "
  write(12,fmt=*) "# ______Parameters_____"
  write(12,fmt=*) "# "
  write(12,fmt=*) "# n_gauss = ", n_gauss
  write(12,fmt=*) "# lambda_amp = ", lambda_amp
  write(12,fmt=*) "# lambda_mu = ", lambda_mu
  write(12,fmt=*) "# lambda_sig = ", lambda_sig
  write(12,fmt=*) "# maxiter_itit = ", maxiter_init
  write(12,fmt=*) "# maxiter = ", maxiter
  write(12,fmt=*) "# lstd = ", lstd
  write(12,fmt=*) "# ustd = ", ustd
  write(12,fmt=*) "# noise = ", noise
  write(12,fmt=*) "# regul = ", regul
  write(12,fmt=*) "# "

  write(12,fmt=*) "# i, j, A, mean, sigma"

  do i=1, dim_data(2)
     do j=1, dim_data(3)
        do k=1, n_gauss
           write(12,fmt=*) i-1, j-1, grid_params(1+((k-1)*3),i,j), grid_params(2+((k-1)*3),i,j), grid_params(3+((k-1)*3),i,j)
        enddo
     enddo
  enddo
  close(12)
  
  call ender()
  
end program ROHSA
