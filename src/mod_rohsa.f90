!! This module contains ROHSA subrtoutine
module mod_rohsa
  !! This module contains ROHSA subrtoutine
  use mod_constants
  use mod_array
  use mod_functions
  use mod_start
  use mod_optimize

  implicit none

  private
  
  public :: main_rohsa

contains

  subroutine main_rohsa(data, std_cube, fileout, n_gauss, n_gauss_add, lambda_amp, lambda_mu, lambda_sig, &
       lambda_var_amp, lambda_var_mu, lambda_var_sig, amp_fact_init, sig_init, maxiter_init, maxiter, m, &
       noise, regul, descent, lstd, ustd, init_option, iprint, iprint_init)
    
    implicit none
    
    logical, intent(in) :: noise           !! if false --> STD map computed by ROHSA with lstd and ustd (if true given by the user)
    logical, intent(in) :: regul           !! if true --> activate regulation
    logical, intent(in) :: descent         !! if true --> activate hierarchical descent to initiate the optimization
    integer, intent(in) :: n_gauss_add     !! number of gaussian to add at each step
    integer, intent(in) :: m               !! number of corrections used in the limited memory matrix by LBFGS-B
    integer, intent(in) :: lstd            !! lower bound to compute the standard deviation map of the cube (if noise .eq. false)
    integer, intent(in) :: ustd            !! upper bound to compute the standrad deviation map of the cube (if noise .eq. false)
    integer, intent(in) :: iprint          !! print option 
    integer, intent(in) :: iprint_init     !! print option init
    integer, intent(in) :: maxiter         !! max iteration for L-BFGS-B alogorithm
    integer, intent(in) :: maxiter_init    !! max iteration for L-BFGS-B alogorithm (init mean spectrum)
    real(xp), intent(in) :: lambda_amp     !! lambda for amplitude parameter
    real(xp), intent(in) :: lambda_mu      !! lamnda for mean position parameter
    real(xp), intent(in) :: lambda_sig     !! lambda for dispersion parameter
    real(xp), intent(in) :: lambda_var_amp !! lambda for amp dispersion parameter
    real(xp), intent(in) :: lambda_var_mu  !! lambda for mean position dispersion parameter
    real(xp), intent(in) :: lambda_var_sig !! lambda for variance dispersion parameter
    real(xp), intent(in) :: amp_fact_init  !! times max amplitude of additional Gaussian
    real(xp), intent(in) :: sig_init       !! dispersion of additional Gaussian

    character(len=8), intent(in) :: init_option !!Init ROHSA with the mean or the std spectrum    
    character(len=512), intent(in) :: fileout   !! name of the output result

    integer :: n_gauss         !! number of gaussian to fit
    integer :: nside        !! size of the reshaped data \(2^{nside}\)
    integer :: n            !! loop index
    integer :: power        !! loop index

    real(xp), intent(in), dimension(:,:,:), allocatable :: data        !! initial fits data
    real(xp), intent(in), dimension(:,:), allocatable   :: std_cube    !! standard deviation map fo the cube is given by the user 

    real(xp), dimension(:,:,:), allocatable :: cube        !! reshape data with nside --> cube
    real(xp), dimension(:,:,:), allocatable :: cube_mean   !! mean cube over spatial axis
    real(xp), dimension(:,:,:), allocatable :: fit_params  !! parameters to optimize with cube mean at each iteration
    real(xp), dimension(:,:,:), allocatable :: grid_params !! parameters to optimize at final step (dim of initial cube)
    real(xp), dimension(:,:), allocatable :: std_map       !! standard deviation map fo the cube computed by ROHSA with lb and ub
    real(xp), dimension(:), allocatable :: std_spect       !! std spectrum of the observation
    real(xp), dimension(:), allocatable :: max_spect       !! max spectrum of the observation
    real(xp), dimension(:), allocatable :: max_spect_norm  !! max spectrum of the observation normalized by the max of the mean spectrum
    real(xp), dimension(:), allocatable :: mean_spect      !! mean spectrum of the observation
    real(xp), dimension(:), allocatable :: guess_spect !! params obtain fi the optimization of the std spectrum of the observation
    
    integer, dimension(3) :: dim_data !! dimension of original data
    integer, dimension(3) :: dim_cube !! dimension of reshape cube
    
    real(xp), dimension(:,:), allocatable :: kernel !! convolution kernel 
    
    integer :: ios=0 !! ios integer
    integer :: i     !! loop index
    integer :: j     !! loop index
    integer :: k     !! loop index
    integer :: l     !! loop index

    call header()
    
    print*, "fileout = '",trim(fileout),"'"
    
    print*,
    print*, "______Parameters_____"
    print*,
    print*, "n_gauss = ", n_gauss
    print*, "n_gauss_add = ", n_gauss_add
    print*, "lambda_amp = ", lambda_amp
    print*, "lambda_mu = ", lambda_mu
    print*, "lambda_sig = ", lambda_sig
    print*, "lambda_var_amp = ", lambda_var_amp
    print*, "lambda_var_mu = ", lambda_var_mu
    print*, "lambda_var_sig = ", lambda_var_sig
    print*, "amp_fact_init = ", amp_fact_init
    print*, "sig_init = ", sig_init
    print*, "init_option = ", init_option
    print*, "maxiter_itit = ", maxiter_init
    print*, "init_option = ", init_option
    print*, "maxiter = ", maxiter
    print*, "lstd = ", lstd
    print*, "ustd = ", ustd
    print*, "noise = ", noise
    print*, "regul = ", regul
    print*, "descent = ", descent
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
    
    print*, "Compute mean and std spectrum"
    allocate(std_spect(dim_data(1)))
    allocate(max_spect(dim_data(1)), max_spect_norm(dim_data(1)))
    allocate(mean_spect(dim_data(1)))

    call std_spectrum(data, std_spect, dim_data(1), dim_data(2), dim_data(3))
    call mean_spectrum(data, mean_spect, dim_data(1), dim_data(2), dim_data(3))
    call max_spectrum(data, max_spect, dim_data(1), dim_data(2), dim_data(3))
    call max_spectrum(data, max_spect_norm, dim_data(1), dim_data(2), dim_data(3), maxval(mean_spect))
    
    call reshape_up(data, cube, dim_data, dim_cube)
    
    !Allocate memory for parameters grids
    if (descent .eqv. .true.) then
       allocate(grid_params(3*(n_gauss+(nside*n_gauss_add)), dim_data(2), dim_data(3)))
       allocate(fit_params(3*(n_gauss+(nside*n_gauss_add)), 1, 1))
    else 
       allocate(grid_params(3*(n_gauss+n_gauss_add), dim_data(2), dim_data(3)))
    end if
    
    print*, "                    Start iteration"
    print*,
    
    if (descent .eqv. .true.) then
       print*, "Start hierarchical descent"
       !Start iteration
       do n=0,nside
          power = 2**n
          
          allocate(cube_mean(dim_cube(1), power, power))
          
          call mean_array(power, cube, cube_mean)
          
          if (n == 0) then
             if (init_option .eq. "mean") then
                print*, "Init mean spectrum"        
                call init_spectrum(n_gauss, fit_params(:,1,1), dim_cube(1), cube_mean(:,1,1), amp_fact_init, sig_init, &
                     maxiter_init, m, iprint_init)
             elseif (init_option .eq. "std") then
                call init_spectrum(n_gauss, fit_params(:,1,1), dim_cube(1), std_spect, amp_fact_init, sig_init, &
                     maxiter_init, m, iprint_init)
             elseif (init_option .eq. "max") then
                call init_spectrum(n_gauss, fit_params(:,1,1), dim_cube(1), max_spect, amp_fact_init, sig_init, &
                     maxiter_init, m, iprint_init)
             elseif (init_option .eq. "maxnorm") then
                call init_spectrum(n_gauss, fit_params(:,1,1), dim_cube(1), max_spect_norm, amp_fact_init, sig_init, &
                     maxiter_init, m, iprint_init)
             else 
                print*, "init_option keyword should be 'mean' or 'std' or 'max' or 'maxnorm'"
                stop
             end if
          end if
          
          ! Propagate solution on new grid (higher resolution)
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

                ! Update parameters 
                call update(cube_mean, fit_params, n_gauss, dim_cube(1), power, power, lambda_amp, lambda_mu, &
                     lambda_sig, lambda_var_amp, lambda_var_mu, lambda_var_sig, maxiter, m, kernel, iprint, std_map)        

                if (n_gauss_add .ne. 0) then !FIXME
                   ! Add new Gaussian if one reduced chi square > 1 
                   call init_new_gauss(cube_mean, fit_params, std_map, n_gauss, dim_cube(1), power, power, amp_fact_init, &
                        sig_init)
                end if

                deallocate(std_map)
             end if
          end if
          
          deallocate(cube_mean)
       enddo
       
       print*,
       write(*,*) "Reshape cube, restore initial dimensions :"
       write(*,*) "dim_v, dim_y, dim_x = ", dim_data
              
       call reshape_down(fit_params, grid_params,  (/ 3*n_gauss, dim_cube(2), dim_cube(3)/), &
            (/ 3*n_gauss, dim_data(2), dim_data(3)/))       

       else
          allocate(guess_spect(3*(n_gauss+n_gauss_add)))
          if (init_option .eq. "mean") then
             print*, "Use of the mean spectrum to initialize each los"
             call init_spectrum(n_gauss, guess_spect, dim_cube(1), mean_spect, amp_fact_init, sig_init, &
                  maxiter_init, m, iprint_init)
          else if (init_option .eq. "std") then
             print*, "Use of the std spectrum to initialize each los"
             call init_spectrum(n_gauss, guess_spect, dim_cube(1), std_spect, amp_fact_init, sig_init, &
                  maxiter_init, m, iprint_init)
          else if (init_option .eq. "max") then
             print*, "Use of the max spectrum to initialize each los"
             call init_spectrum(n_gauss, guess_spect, dim_cube(1), max_spect, amp_fact_init, sig_init, &
                  maxiter_init, m, iprint_init)
          else if (init_option .eq. "maxnorm") then
             print*, "Use of the std spectrum to initialize each los"
             call init_spectrum(n_gauss, guess_spect, dim_cube(1), max_spect_norm, amp_fact_init, sig_init, &
                  maxiter_init, m, iprint_init)
          else
             print*, "init_option keyword should be 'mean' or 'std' or 'max'"
             stop
          end if
          call init_grid_params(grid_params, guess_spect, dim_data(2), dim_data(3))

          deallocate(guess_spect)      
    end if
    
    !Update last level
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
       call update(data, grid_params, n_gauss, dim_data(1), dim_data(2), dim_data(3), lambda_amp, lambda_mu, &
            lambda_sig, lambda_var_amp, lambda_var_mu, lambda_var_sig, maxiter, m, kernel, iprint, std_map)
       
       if (n_gauss_add .ne. 0) then !FIXME KEYWORD
          do l=1,n_gauss_add
             ! Add new Gaussian if one reduced chi square > 1 
             call init_new_gauss(data, grid_params, std_map, n_gauss, dim_data(1), dim_data(2), dim_data(3), amp_fact_init, &
                  sig_init)
             call update(data, grid_params, n_gauss, dim_data(1), dim_data(2), dim_data(3), lambda_amp, lambda_mu, &
                  lambda_sig, lambda_var_amp, lambda_var_mu, lambda_var_sig, maxiter, m, kernel, iprint, std_map)
          end do
       end if

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
    write(12,fmt=*) "# n_gauss_add = ", n_gauss_add
    write(12,fmt=*) "# lambda_amp = ", lambda_amp
    write(12,fmt=*) "# lambda_mu = ", lambda_mu
    write(12,fmt=*) "# lambda_sig = ", lambda_sig
    write(12,fmt=*) "# lambda_var_amp = ", lambda_var_amp
    write(12,fmt=*) "# lambda_var_mu = ", lambda_var_mu
    write(12,fmt=*) "# lambda_var_sig = ", lambda_var_sig
    write(12,fmt=*) "# amp_fact_init = ", amp_fact_init
    write(12,fmt=*) "# sig_init = ", sig_init
    write(12,fmt=*) "# init_option = ", init_option
    write(12,fmt=*) "# maxiter_itit = ", maxiter_init
    write(12,fmt=*) "# maxiter = ", maxiter
    write(12,fmt=*) "# lstd = ", lstd
    write(12,fmt=*) "# ustd = ", ustd
    write(12,fmt=*) "# noise = ", noise
    write(12,fmt=*) "# regul = ", regul
    write(12,fmt=*) "# descent = ", descent
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
    
  end subroutine main_rohsa
  
end module mod_rohsa
