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
        
    nside = dim2nside(dim_data)
    
    write(*,*) "nside = ", nside
            
    !Allocate memory for parameters grids
    allocate(grid_params(3*params%n_mbb, dim_data(2), dim_data(3)))
    allocate(b_params(params%n_mbb))
    allocate(c_params(params%n_mbb))
    allocate(d_params(params%n_mbb))
    allocate(stefan_params(params%n_mbb))

    !Init
    grid_params(1,:,:) = params%tau_init
    grid_params(2,:,:) = params%beta_init
    grid_params(3,:,:) = params%Td_init
    grid_params(4,:,:) = params%tau_init
    grid_params(5,:,:) = params%beta_init
    grid_params(6,:,:) = params%Td_init    

    b_params = 0._xp
    c_params = 0._xp
    d_params = 0._xp
    stefan_params = 0._xp

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
