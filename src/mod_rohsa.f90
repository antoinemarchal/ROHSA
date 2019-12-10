!! This module contains ROHSA subrtoutine
module mod_rohsa
  !! This module contains ROHSA subrtoutine
  use mod_constants
  use mod_model
  use mod_array
  use mod_functions
  use mod_start
  use mod_optimize
  use mod_read_parameters
  use mod_rw_data

  implicit none

  private
  
  public :: main_rohsa

contains

  subroutine main_rohsa()
    
    implicit none
    
    integer :: nside        !! size of the reshaped data \(2^{nside}\)
    integer :: n            !! loop index
    integer :: power        !! loop index

    real(xp), dimension(:,:,:), allocatable :: cube            !! reshape data with nside --> cube
    real(xp), dimension(:,:,:), allocatable :: cube_HI         !! reshape data with nside --> cube

    real(xp), dimension(:,:,:), allocatable :: cube_mean       !! mean cube over spatial axis
    real(xp), dimension(:,:,:), allocatable :: cube_HI_mean    !! mean cube over spatial axis

    real(xp), dimension(:,:,:), allocatable :: fit_params      !! parameters to optimize with cube mean at each iteration
    real(xp), dimension(:,:,:), allocatable :: grid_params     !! parameters to optimize at final step (dim of initial cube)

    real(xp), dimension(:,:,:), allocatable :: std_cube_mean   !! standard deviation cube

    real(xp), dimension(:), allocatable :: b_params            !! unknow average tauma
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
        
    dim_data = shape(data%cube)
    dim_data_HI = shape(data%NHI)
    
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
    allocate(b_params(params%n_mbb))
    allocate(c_params(params%n_mbb))
    allocate(d_params(params%n_mbb))
    allocate(stefan_params(params%n_mbb))

    call reshape_up(data%cube, cube, dim_data, dim_cube)
    call reshape_up(data%NHI, cube_HI, dim_data_HI, dim_cube_HI)
    
    !Allocate memory for parameters grids
    allocate(grid_params(3*params%n_mbb, dim_data(2), dim_data(3)))
    allocate(fit_params(3*params%n_mbb, 1, 1))
    
    print*, "                    Start iteration"
    print*, " "
    
    print*, "Start hierarchical descent"

    if (params%save_grid .eqv. .true.) then
       !Open file time step
       open(unit=11, file=params%timeout, status='replace', access='append', iostat=ios)
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
       
       call mean_array(power, data%cube, cube_mean)
       call mean_array(power, data%NHI, cube_HI_mean)
       call mean_array(power, data%std_cube, std_cube_mean)           
       
       if (n == 0) then
          print*, "Init mean spectrum"        
          
          do i=1, params%n_mbb    
             fit_params(1+(3*(i-1)),1,1) = params%tau_init
             fit_params(2+(3*(i-1)),1,1) = params%beta_init
             fit_params(3+(3*(i-1)),1,1) = params%Td_init                   
          end do
          !if ciba .true. __. replace lower and upper bounds
          if (params%ciba .eqv. .true.) then
             print*, "Your first MBB is attached to the CIBA fluctuations"
             fit_params(1,1,1) = params%tau_init_cib
             fit_params(2,1,1) = params%beta_init_cib
             fit_params(3,1,1) = params%Td_init_cib
          end if
          
          call init_spectrum(fit_params(:,1,1), dim_cube(1), cube_mean(:,1,1), data%wavelength, &
               data%color, std_cube_mean(:,1,1))

          !Init b_params
          do i=1, params%n_mbb       
             b_params(i) = fit_params(1+(3*(i-1)),1,1)
             c_params(i) = fit_params(2+(3*(i-1)),1,1)
             d_params(i) = fit_params(3+(3*(i-1)),1,1)
             stefan_params(i) = lumi_cst(fit_params(1+(3*(i-1)),1,1), fit_params(2+(3*(i-1)),1,1), &
                  fit_params(3+(3*(i-1)),1,1), params%l0)
          end do
       end if
              
       if (n == 0) then                
          print*,  "Update level", n
          call upgrade(cube_mean, fit_params, data%wavelength, power, dim_cube(1), data%color, &
               std_cube_mean)
       end if
              
       if (n > 0 .and. n < nside) then          
          ! Update parameters 
          print*,  "Update level", n, ">", power
          call update(cube_mean, cube_HI_mean, data%wavelength, fit_params, b_params, c_params, d_params, &
               stefan_params, dim_cube(1), power, power, kernel, std_cube_mean, data%color)
          
       end if
       
       deallocate(cube_mean)
       deallocate(cube_HI_mean)
       deallocate(std_cube_mean)
       
       ! Save grid in file
       if (params%save_grid .eqv. .true.) then
          print*, "Save grid parameters"
          call save_process(n, params%n_mbb, fit_params, power, params%fileout)
          !Save timestep
          if (n .ne. 0) then
             open(unit=11, file=params%timeout, status='unknown', access='append', iostat=ios)
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
    
    call reshape_down(fit_params, grid_params,  (/ 3*params%n_mbb, dim_cube(2), dim_cube(3)/), &
         (/ 3*params%n_mbb, dim_data(2), dim_data(3)/))       
        
    !Update last level
    print*, " "
    print*, "Start updating last level."
    print*, " "
    
        
    call update(data%cube, data%NHI, data%wavelength, grid_params, b_params, c_params, d_params, stefan_params, &
         dim_data(1), dim_data(2), dim_data(3), kernel, data%std_cube, data%color)       
    
    print*, " "
    print*, "_____ Write output file _____"
    print*, " "
    
    ! Open file
    open(unit=12, file=params%fileout, action="write", iostat=ios)
    if (ios /= 0) stop "opening file error"
    
    call print_parameters_unit(12)
    write(12,fmt=*) "# i, j, Td, beta, Td"

    do i=1, dim_data(2)
       do j=1, dim_data(3)
          do k=1, params%n_mbb
             write(12,fmt=*) i-1, j-1, grid_params(1+((k-1)*3),i,j), grid_params(2+((k-1)*3),i,j), grid_params(3+((k-1)*3),i,j)
          enddo
       enddo
    enddo
    close(12)

    open(unit=11, file=params%timeout, status='unknown', access='append', iostat=ios)
    if (ios /= 0) stop "opening file error"
    call cpu_time(uctime)
    write(11,fmt=*) dim_cube(2), uctime-lctime
    close(11)
    
  end subroutine main_rohsa
  
end module mod_rohsa
