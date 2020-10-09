!! This module contains ROHSA subrtoutine
module mod_rohsa
  !! This module contains ROHSA subrtoutine
  use mod_constants
  use mod_array
  use mod_functions
  use mod_start
  use mod_optimize
  use mod_rw_data
  use mod_read_parameters

  implicit none

  private
  
  public :: main_rohsa

contains

  subroutine main_rohsa()
    
    implicit none
    
    integer :: nside        !! size of the reshaped data \(2^{nside}\)
    integer :: n            !! loop index
    integer :: power        !! loop index

    type(indata) :: cube, cube_mean
    type(indata_s) :: mean

    real(xp), dimension(:,:,:), allocatable :: fit_params      !! parameters to optimize with cube mean at each iteration
    real(xp), dimension(:,:,:), allocatable :: grid_params     !! parameters to optimize at final step (dim of initial cube)

    integer, dimension(3) :: dim_data !! dimension of original data
    integer, dimension(3) :: dim_cube !! dimension of reshape cube
    
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
        
    dim_data = shape(data%q)
    
    write(*,*) "dim_v, dim_y, dim_x = ", dim_data
    write(*,*) ""
    write(*,*) "number of los = ", dim_data(2)*dim_data(3)
    
    nside = dim2nside(dim_data)
    
    write(*,*) "nside = ", nside
    
    call dim_data2dim_cube(nside, dim_data, dim_cube)
    
    !Allocate memory for cube
    allocate(cube%q(dim_cube(1), dim_cube(2), dim_cube(3)))
    allocate(cube%u(dim_cube(1), dim_cube(2), dim_cube(3)))
    allocate(cube%p(dim_cube(1), dim_cube(2), dim_cube(3)))
    allocate(cube%rms(dim_cube(2), dim_cube(3)))
    
    !Reshape the data (new cube of size nside)
    print*, " "
    write(*,*) "Reshape cube, new dimensions :"
    write(*,*) "dim_v, dim_y, dim_x = ", dim_cube
    print*, " "
    
    print*, "Compute mean and std spectrum"
    allocate(mean%q(dim_data(1)))
    allocate(mean%u(dim_data(1)))
    allocate(mean%p(dim_data(1)))

    call mean_spectrum(data%q, mean%q, dim_data(1), dim_data(2), dim_data(3))    
    call mean_spectrum(data%u, mean%u, dim_data(1), dim_data(2), dim_data(3))    
    call mean_spectrum(data%p, mean%p, dim_data(1), dim_data(2), dim_data(3))    

    call reshape_up(data%q, cube%q, dim_data, dim_cube)
    call reshape_up(data%u, cube%u, dim_data, dim_cube)
    call reshape_up(data%p, cube%p, dim_data, dim_cube)
    call reshape_noise_up(data%rms, cube%rms, dim_data, dim_cube)
    
    !Allocate memory for parameters grids
    allocate(grid_params(2*params%n, dim_data(2), dim_data(3)))
    allocate(fit_params(2*params%n, 1, 1))
    
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
       
       allocate(cube_mean%q(dim_cube(1), power, power))
       allocate(cube_mean%u(dim_cube(1), power, power))
       allocate(cube_mean%p(dim_cube(1), power, power))
       allocate(cube_mean%rms(power, power))
       
       call mean_array(power, cube%q, cube_mean%q)
       call mean_array(power, cube%u, cube_mean%u)
       call mean_array(power, cube%p, cube_mean%p)
       call mean_map(power, cube%rms, cube_mean%rms)                     
       
       if (n == 0) then
          print*, "Init mean spectrum"        

          do i=1, params%n
             fit_params(1+(2*(i-1)),1,1) = params%amp_init
             fit_params(2+(2*(i-1)),1,1) = params%mu_init
          end do

          call init_spectrum(fit_params(:,1,1), mean, dim_cube(1))
       end if
                     
       if (n > 0 .and. n < nside) then          
          ! Update parameters 
          print*,  "Update level", n, ">", power
          call update(cube_mean, fit_params, dim_cube(1), power, power, cube_mean%rms, kernel)                  
       end if
       
       deallocate(cube_mean%q)
       deallocate(cube_mean%u)
       deallocate(cube_mean%p)
       deallocate(cube_mean%rms)
       
       ! Save grid in file
       if (params%save_grid .eqv. .true.) then
          print*, "Save grid parameters"
          call save_process(n, params%n, fit_params, power, params%fileout)
          !Save timestep
          if (n .ne. 0) then
             open(unit=11, file=params%timeout, status='unknown', access='append', iostat=ios)
             if (ios /= 0) stop "opening file error"
             call cpu_time(uctime)
             print*, dim_cube(1)
             write(11,fmt=*) power, uctime-lctime
             close(11)
          end if
       end if
       
       ! Propagate solution on new grid (higher resolution)
       call go_up_level(fit_params)
       write(*,*) " "
       write(*,*) "Interpolate parameters level ", n!, ">", power
       
    enddo
       
    print*, " "
    write(*,*) "Reshape cube, restore initial dimensions :"
    write(*,*) "dim_v, dim_y, dim_x = ", dim_data
    
    call reshape_down(fit_params, grid_params,  (/ 2*params%n, dim_cube(2), dim_cube(3)/), &
         (/ 2*params%n, dim_data(2), dim_data(3)/))       
     
    !Update last level
    print*, " "
    print*, "Start updating last level."
    print*, " "
        
    call update(data, grid_params, dim_cube(1), dim_data(2), dim_data(3), data%rms, kernel)                  
    
    print*, " "
    print*, "_____ Write output file _____"
    print*, " "
    
    ! Open file
    open(unit=12, file=params%fileout, action="write", iostat=ios)
    if (ios /= 0) stop "opening file error"
    
    call print_parameters_unit(12)        
    write(12,fmt=*) "# i, j, a, mu"

    do i=1, dim_data(2)
       do j=1, dim_data(3)
          do k=1, params%n
             write(12,fmt=*) i-1, j-1, grid_params(1+((k-1)*2),i,j), grid_params(2+((k-1)*2),i,j)
          enddo
       enddo
    enddo
    close(12)
    
    open(unit=11, file=params%timeout, status='unknown', access='append', iostat=ios)
    if (ios /= 0) stop "opening file error"
    call cpu_time(uctime)
    print*, dim_cube(1)
    write(11,fmt=*) dim_cube(2), uctime-lctime
    close(11)
    
  end subroutine main_rohsa
  
end module mod_rohsa
