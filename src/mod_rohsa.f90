!! This module contains ROHSA subrtoutine
module mod_rohsa
  !! This module contains ROHSA subrtoutine
  use mod_constants
  use mod_model
  use mod_array
  use mod_functions
  use mod_start
  use mod_optimize
  use mod_inout

  implicit none

  private
  
  public :: main_rohsa

contains

  subroutine main_rohsa(data, wavelength, std_cube, data_HI, fileout, timeout, n_mbb, lambda_sig, lambda_beta, &
       lambda_Td, lambda_var_sig, lambda_var_beta, lambda_var_Td, lambda_stefan, amp_fact_init, sig_init, beta_init, &
       Td_init, lb_sig, ub_sig, lb_beta, ub_beta, lb_Td, ub_Td, l0, maxiter_init, maxiter, m, noise, lstd, ustd, iprint, &
       iprint_init, save_grid, color, degree, cc)
    
    implicit none
    
    logical, intent(in) :: noise           !! if false --> STD map computed by ROHSA with lstd and ustd (if true given by the user)
    logical, intent(in) :: save_grid       !! save grid of fitted parameters at each step of the multiresolution process
    logical, intent(in) :: cc              !! if true --> apply colour correction PLANCK+IRAS
    integer, intent(in) :: m               !! number of corrections used in the limited memory matrix by LBFGS-B
    integer, intent(in) :: lstd            !! lower bound to compute the standard deviation map of the cube (if noise .eq. false)
    integer, intent(in) :: ustd            !! upper bound to compute the standrad deviation map of the cube (if noise .eq. false)
    integer, intent(in) :: iprint          !! print option 
    integer, intent(in) :: iprint_init     !! print option init
    integer, intent(in) :: maxiter         !! max iteration for L-BFGS-B alogorithm
    integer, intent(in) :: maxiter_init    !! max iteration for L-BFGS-B alogorithm (init mean spectrum)

    real(xp), intent(in) :: lambda_sig     !! lambda for amplitude parameter
    real(xp), intent(in) :: lambda_beta    !! lamnda for mean position parameter
    real(xp), intent(in) :: lambda_Td      !! lambda for dispersion parameter

    real(xp), intent(in) :: lambda_var_sig  !! lambda for amp dispersion parameter
    real(xp), intent(in) :: lambda_var_beta !! lambda for mean position dispersion parameter
    real(xp), intent(in) :: lambda_var_Td !! lambda for mean position dispersion parameter
    real(xp), intent(in) :: lambda_stefan   !! lambda for variance dispersion parameter

    real(xp), intent(in) :: amp_fact_init  !! times max amplitude of additional Gaussian

    real(xp), intent(in) :: sig_init      !! 
    real(xp), intent(in) :: beta_init     !! 
    real(xp), intent(in) :: Td_init       !! 

    real(xp), intent(in) :: lb_sig        !! lower bound
    real(xp), intent(in) :: ub_sig        !! upper bound
    real(xp), intent(in) :: lb_beta       !! lower bound
    real(xp), intent(in) :: ub_beta       !! upper bound
    real(xp), intent(in) :: lb_Td         !! lower bound
    real(xp), intent(in) :: ub_Td         !! upper bound

    real(xp), intent(in) :: l0 !! reference wavelength
    integer, intent(in) :: degree

    character(len=512), intent(in) :: fileout   !! name of the output result
    character(len=512), intent(in) :: timeout   !! name of the output result

    integer :: n_mbb      !! number of gaussian to fit
    integer :: nside        !! size of the reshaped data \(2^{nside}\)
    integer :: n            !! loop index
    integer :: power        !! loop index

    real(xp), intent(in), dimension(:,:,:), allocatable :: data          !! initial fits data
    real(xp), intent(in), dimension(:,:,:), allocatable :: data_HI       !! initial fits data data_HI
    real(xp), intent(in), dimension(:,:), allocatable   :: color         !! polynomial coefficient for color correction
    real(xp), intent(in), dimension(:), allocatable     :: wavelength    !! wavelength Planck + IRAS
    real(xp), intent(in), dimension(:,:,:), allocatable :: std_cube      !! standard deviation cube

    real(xp), dimension(:,:,:), allocatable :: cube            !! reshape data with nside --> cube
    real(xp), dimension(:,:,:), allocatable :: cube_HI         !! reshape data with nside --> cube
    real(xp), dimension(:,:,:), allocatable :: cube_mean       !! mean cube over spatial axis
    real(xp), dimension(:,:,:), allocatable :: cube_HI_mean    !! mean cube over spatial axis
    real(xp), dimension(:,:,:), allocatable :: fit_params      !! parameters to optimize with cube mean at each iteration
    real(xp), dimension(:,:,:), allocatable :: grid_params     !! parameters to optimize at final step (dim of initial cube)
    real(xp), dimension(:,:,:), allocatable :: std_cube_mean   !! standard deviation cube
    real(xp), dimension(:), allocatable :: b_params            !! unknow average sigma
    real(xp), dimension(:), allocatable :: c_params            !! unknow average beta
    real(xp), dimension(:), allocatable :: d_params            !! unknow average Td
    real(xp), dimension(:), allocatable :: stefan_params       !! unknow average propto Luminosity

    integer, dimension(3) :: dim_data !! dimension of original data
    integer, dimension(3) :: dim_cube !! dimension of reshape cube
    integer, dimension(3) :: dim_data_HI !! dimension of original data
    integer, dimension(3) :: dim_cube_HI !! dimension of reshape cube
    
    real(xp), dimension(:,:), allocatable :: kernel !! convolution kernel 

    real(xp) :: lctime, uctime
    
    integer :: ios=0 !! ios integer
    integer :: i     !! loop index
    integer :: j     !! loop index
    integer :: k     !! loop index
        
    print*, "fileout = '",trim(fileout),"'"
    print*, "timeout = '",trim(timeout),"'"
    
    print*, " "
    print*, "______Parameters_____"
    print*, "n_mbb = ", n_mbb

    print*, "lambda_sig = ", lambda_sig
    print*, "lambda_beta = ", lambda_beta
    print*, "lambda_Td = ", lambda_Td

    print*, "lambda_var_sig = ", lambda_var_sig
    print*, "lambda_var_beta = ", lambda_var_beta
    print*, "lambda_var_Td = ", lambda_var_Td
    print*, "lambda_stefan = ", lambda_stefan

    print*, "amp_fact_init = ", amp_fact_init

    print*, "sig_init = ", sig_init
    print*, "beta_init = ", beta_init
    print*, "Td_init = ", Td_init

    print*, "lb_sig = ", lb_sig
    print*, "ub_sig = ", ub_sig
    print*, "lb_beta = ", lb_beta
    print*, "ub_beta = ", ub_beta
    print*, "lb_Td = ", lb_Td
    print*, "ub_Td = ", ub_Td

    print*, "l0 = ", l0

    print*, "maxiter_init = ", maxiter_init
    print*, "maxiter = ", maxiter
    print*, "lstd = ", lstd
    print*, "ustd = ", ustd
    print*, "noise = ", noise
    print*, "save_grid = ", save_grid
    print*, "cc = ", cc

    print*, " "
    
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
        
    dim_data = shape(data)
    dim_data_HI = shape(data_HI)
    
    write(*,*) "dim_v, dim_y, dim_x = ", dim_data
    write(*,*) ""
    write(*,*) "number of los = ", dim_data(2)*dim_data(3)
    
    nside = dim2nside(dim_data)
    
    write(*,*) "nside = ", nside
    
    call dim_data2dim_cube(nside, dim_data, dim_cube)
    call dim_data2dim_cube(nside, dim_data_HI, dim_cube_HI)
    
    !Allocate moemory for cube
    allocate(cube(dim_cube(1), dim_cube(2), dim_cube(3)))
    allocate(cube_HI(dim_cube_HI(1), dim_cube_HI(2), dim_cube_HI(3)))
    
    !Reshape the data (new cube of size nside)
    print*, " "
    write(*,*) "Reshape cube, new dimensions :"
    write(*,*) "dim_v, dim_y, dim_x = ", dim_cube
    print*, " "
    
    print*, "Compute mean and std spectrum"
    allocate(b_params(n_mbb))
    allocate(c_params(n_mbb))
    allocate(d_params(n_mbb))
    allocate(stefan_params(n_mbb))

    call reshape_up(data, cube, dim_data, dim_cube)
    call reshape_up(data_HI, cube_HI, dim_data_HI, dim_cube_HI)
    
    !Allocate memory for parameters grids
    allocate(grid_params(3*n_mbb, dim_data(2), dim_data(3)))
    allocate(fit_params(3*n_mbb, 1, 1))
    
    print*, "                    Start iteration"
    print*, " "
    
    print*, "Start hierarchical descent"

    if (save_grid .eqv. .true.) then
       !Open file time step
       open(unit=11, file=timeout, status='replace', access='append', iostat=ios)
       write(11,fmt=*) "# size grid, Time (s)"
       close(11)
       call cpu_time(lctime)
    end if
    
    !Start iteration
    do n=0,nside-1
       power = 2**n
       
       allocate(cube_mean(dim_cube(1), power, power))
       allocate(cube_HI_mean(dim_cube_HI(1), power, power))
       allocate(std_cube_mean(dim_cube(1),power, power))          
       
       call mean_array(power, cube, cube_mean)
       call mean_array(power, cube_HI, cube_HI_mean)
       call mean_array(power, std_cube, std_cube_mean)           
       
       if (n == 0) then
          print*, "Init mean spectrum"        

          do i=1, n_mbb    
             fit_params(1+(3*(i-1)),1,1) = sig_init
             fit_params(2+(3*(i-1)),1,1) = beta_init
             fit_params(3+(3*(i-1)),1,1) = Td_init
          end do
                   
          call init_spectrum(n_mbb, fit_params(:,1,1), dim_cube(1), cube_mean(:,1,1), cube_HI_mean(:,1,1), &
               wavelength, lb_sig, ub_sig, lb_beta, ub_beta, lb_Td, ub_Td, l0, maxiter_init, m, iprint_init, &
               color, degree, std_cube_mean(:,1,1), cc)

          !Init b_params
          do i=1, n_mbb       
             b_params(i) = fit_params(1+(3*(i-1)),1,1)
             c_params(i) = fit_params(2+(3*(i-1)),1,1)
             d_params(i) = fit_params(3+(3*(i-1)),1,1)
             stefan_params(i) = lumi_cst(fit_params(1+(3*(i-1)),1,1), fit_params(2+(3*(i-1)),1,1), &
                  fit_params(3+(3*(i-1)),1,1), l0)
          end do
       end if
              
       if (n == 0) then                
          print*,  "Update level", n
          call upgrade(cube_mean, fit_params, cube_HI_mean, wavelength, power, n_mbb, dim_cube(1), lb_sig, ub_sig, &
               lb_beta, ub_beta, lb_Td, ub_Td, l0, maxiter, m, iprint, color, degree, std_cube_mean, cc)
       end if
              
       if (n > 0 .and. n < nside) then          
          ! Update parameters 
          print*,  "Update level", n, ">", power
          call update(cube_mean, cube_HI_mean, wavelength, fit_params, b_params, c_params, d_params, stefan_params, &
               n_mbb, dim_cube(1), power, power, lambda_sig, lambda_beta, lambda_Td, lambda_var_sig, lambda_var_beta, &
               lambda_var_Td, lambda_stefan, lb_sig, ub_sig, lb_beta, ub_beta, lb_Td, ub_Td, l0, maxiter, m, kernel, &
               iprint, std_cube_mean, color, degree, cc)
          
       end if
       
       deallocate(cube_mean)
       deallocate(cube_HI_mean)
       deallocate(std_cube_mean)
       
       ! Save grid in file
       if (save_grid .eqv. .true.) then
          print*, "Save grid parameters"
          call save_process(n, n_mbb, fit_params, power, fileout)
          !Save timestep
          if (n .ne. 0) then
             open(unit=11, file=timeout, status='unknown', access='append', iostat=ios)
             if (ios /= 0) stop "opening file error"
             call cpu_time(uctime)
             print*, dim_cube(1)
             ! write(11,fmt='(I6, f6.3)') power, uctime-lctime
             write(11,fmt=*) power, uctime-lctime
             close(11)
          end if
       end if
       
       ! Propagate solution on new grid (higher resolution)
       call go_up_level(fit_params)
       write(*,*) ""
       write(*,*) "Interpolate parameters level ", n!, ">", power
       
    enddo
    
    print*, " "
    write(*,*) "Reshape cube, restore initial dimensions :"
    write(*,*) "dim_v, dim_y, dim_x = ", dim_data
    
    call reshape_down(fit_params, grid_params,  (/ 3*n_mbb, dim_cube(2), dim_cube(3)/), &
         (/ 3*n_mbb, dim_data(2), dim_data(3)/))       
        
    !Update last level
    print*, " "
    print*, "Start updating last level."
    print*, " "
    
        
    call update(data, data_HI, wavelength, grid_params, b_params, c_params, d_params, stefan_params, n_mbb, &
         dim_data(1), dim_data(2), dim_data(3), lambda_sig, lambda_beta, lambda_Td, lambda_var_sig, &
         lambda_var_beta, lambda_var_Td, lambda_stefan, lb_sig, ub_sig, lb_beta, ub_beta, lb_Td, ub_Td, l0, maxiter, &
         m, kernel, iprint, std_cube, color, degree, cc)       
    
    print*, " "
    print*, "_____ Write output file _____"
    print*, " "
    
    ! Open file
    open(unit=12, file=fileout, action="write", iostat=ios)
    if (ios /= 0) stop "opening file error"
    
    write(12,fmt=*) "# "
    write(12,fmt=*) "# ______Parameters_____"
    write(12,fmt=*) "# "
    write(12,fmt=*) "# n_mbb = ", n_mbb
    write(12,fmt=*) "# lambda_sig = ", lambda_sig
    write(12,fmt=*) "# lambda_beta = ", lambda_beta
    write(12,fmt=*) "# lambda_Td = ", lambda_Td
    write(12,fmt=*) "# lambda_var_sig = ", lambda_var_sig
    write(12,fmt=*) "# lambda_var_beta = ", lambda_var_beta
    write(12,fmt=*) "# lambda_var_Td = ", lambda_var_Td
    write(12,fmt=*) "# lambda_stefan = ", lambda_stefan
    write(12,fmt=*) "# amp_fact_init = ", amp_fact_init
    write(12,fmt=*) "# sig_init = ", sig_init
    write(12,fmt=*) "# beta_init = ", beta_init
    write(12,fmt=*) "# Td_init = ", Td_init
    write(12,fmt=*) "# lb_sig = ", lb_sig
    write(12,fmt=*) "# ub_sig = ", ub_sig
    write(12,fmt=*) "# lb_beta = ", lb_beta
    write(12,fmt=*) "# ub_beta = ", ub_beta
    write(12,fmt=*) "# lb_Td = ", lb_Td
    write(12,fmt=*) "# ub_Td = ", ub_Td
    write(12,fmt=*) "# maxiter_itit = ", maxiter_init
    write(12,fmt=*) "# maxiter = ", maxiter
    write(12,fmt=*) "# lstd = ", lstd
    write(12,fmt=*) "# ustd = ", ustd
    write(12,fmt=*) "# noise = ", noise
    write(12,fmt=*) "# "
    
    write(12,fmt=*) "# i, j, Td, beta, Td"

    do i=1, dim_data(2)
       do j=1, dim_data(3)
          do k=1, n_mbb
             write(12,fmt=*) i-1, j-1, grid_params(1+((k-1)*3),i,j), grid_params(2+((k-1)*3),i,j), grid_params(3+((k-1)*3),i,j)
          enddo
       enddo
    enddo
    close(12)

    open(unit=11, file=timeout, status='unknown', access='append', iostat=ios)
    if (ios /= 0) stop "opening file error"
    call cpu_time(uctime)
    write(11,fmt=*) dim_cube(2), uctime-lctime
    close(11)
    
  end subroutine main_rohsa
  
end module mod_rohsa
