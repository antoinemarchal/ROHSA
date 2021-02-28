!! This module read the input user parameters (parameters.txt file / data / noise if true)
module mod_inout
  !! This module read the input user parameters (parameters.txt file / data / noise if true)
  use mod_constants
  use mod_convert
  
  implicit none
  
  private
  
  public :: read_cube, read_output, read_map, read_array, read_parameters, save_process

contains
  
  subroutine read_parameters(filename_parameters, filename, fileout, timeout, filename_noise, filename_init_spec, &
       n_gauss, lambda_amp, lambda_mu, lambda_sig, lambda_var_amp, lambda_var_mu, lambda_var_sig, lambda_lym_sig, &
       amp_fact_init, sig_init, lb_sig_init, ub_sig_init, lb_sig, ub_sig, init_option, maxiter_init, maxiter, m, noise, &
       regul, descent, lstd, ustd, iprint, iprint_init, save_grid, lym, init_grid, fileinit, init_spec)
    implicit none

    integer :: ios=0

    character(len=512), intent(in) :: filename_parameters

    integer, intent(inout) :: n_gauss
    integer, intent(inout) :: m 
    integer, intent(inout) :: lstd, ustd
    integer, intent(inout) :: iprint, iprint_init
    integer, intent(inout) :: maxiter, maxiter_init
    real(xp), intent(inout) :: lambda_amp, lambda_mu, lambda_sig
    real(xp), intent(inout) :: lambda_var_amp, lambda_var_mu, lambda_var_sig
    real(xp), intent(inout) :: lambda_lym_sig
    real(xp), intent(inout) :: amp_fact_init, sig_init
    real(xp), intent(inout) :: ub_sig_init, ub_sig
    real(xp), intent(inout) :: lb_sig_init, lb_sig
    logical, intent(inout) :: noise, regul, descent, save_grid
    logical, intent(inout) :: lym, init_spec
    logical, intent(inout) :: init_grid

    character(len=512), intent(inout) :: filename
    character(len=512), intent(inout) :: fileout
    character(len=512), intent(inout) :: timeout
    character(len=512), intent(inout) :: filename_noise
    character(len=512), intent(inout) :: filename_init_spec
    character(len=8), intent(inout) :: init_option
    character(len=512), intent(inout) :: fileinit

    namelist /user_parameters/ filename, fileout, timeout, filename_noise, filename_init_spec, n_gauss, lambda_amp, &
         lambda_mu, lambda_sig, lambda_var_amp, lambda_var_mu, lambda_var_sig, lambda_lym_sig, amp_fact_init, &
         sig_init, lb_sig_init, ub_sig_init, lb_sig, ub_sig, init_option, maxiter_init, maxiter, m, noise, regul, &
         descent, lstd, ustd, iprint, iprint_init, save_grid, lym, init_grid, fileinit, init_spec
    
    open(unit=11, file=filename_parameters, status="old", iostat=ios)
    if (ios /= 0) stop "opening file error"
    
    read(11, user_parameters)
    
    close(11)
  end subroutine read_parameters


  subroutine read_cube(filename, cube)
    implicit none
    integer           :: ios=0, i
    integer           :: v, y, x 
    real(xp)          :: val
    integer           :: nv, ny, nx !cube dimension
    integer           :: nl
    character(len=512), intent(in) :: filename
    real(xp), intent(inout), dimension(:,:,:), allocatable :: cube

    open(unit=11, file=filename, action="read", status="old", iostat=ios)
    if (ios /= 0) stop "opening file error"
    
    read(11,fmt=*) nv, ny, nx
    nl = nv*ny*nx

    allocate(cube(nv,ny,nx))

    do i=1,nl
       read(11,fmt=*) v, y, x, val
       cube(v+1,y+1,x+1) = val
    enddo
    
    close(11)
  end subroutine read_cube


  subroutine read_output(filename, cube, n_gauss, dim_data)
    implicit none
    integer, intent(in) :: n_gauss
    character(len=512), intent(in) :: filename
    integer, intent(in), dimension(3) :: dim_data !! dimension of original data
    real(xp), intent(inout), dimension(:,:,:), allocatable :: cube

    integer           :: ios=0, i, j, k
    integer           :: y, x 

    open(unit=11, file=filename, action="read", status="old", iostat=ios)
    if (ios /= 0) stop "opening file error"
    
    allocate(cube(3*n_gauss,dim_data(2),dim_data(3)))

    do i=1, dim_data(2)
       do j=1, dim_data(3)
          do k=1, n_gauss
             read(11,fmt=*) y, x, cube(1+((k-1)*3),i,j), cube(2+((k-1)*3),i,j), cube(3+((k-1)*3),i,j)
          enddo
       enddo
    enddo
    
    close(11)
  end subroutine read_output

  
  subroutine read_map(filename, map)
    implicit none

    integer           :: ios=0, i
    integer           :: y,x 
    real(xp)          :: val
    integer           :: ny,nx
    integer           :: nl
    character(len=512), intent(in) :: filename
    real(xp), intent(inout), dimension(:,:), allocatable :: map

    open(unit=11, file=filename, action="read", status="old", iostat=ios)
    if (ios /= 0) stop "opening file error"
    
    read(11,fmt=*) ny, nx
    nl = ny*nx

    allocate(map(ny,nx))

    do i=1,nl
       read(11,fmt=*) y, x, val
       map(y+1,x+1) = val
    enddo
    
    close(11)
  end subroutine read_map


  subroutine read_array(filename, array)
    implicit none
    integer           :: ios=0, i
    integer           :: nl
    real(xp)          :: val
    character(len=512), intent(in) :: filename
    real(xp), intent(inout), dimension(:), allocatable :: array

    open(unit=11, file=filename, action="read", status="old", iostat=ios)
    if (ios /= 0) stop "opening file error"
    
    read(11,fmt=*) nl

    allocate(array(nl))

    do i=1,nl
       read(11,fmt=*) val
       array(i) = val
    enddo
    
    close(11)
  end subroutine read_array
  

  subroutine save_process(nside, n_gauss, grid, dim_yx, fileout)
    implicit none

    integer, intent(in) :: n_gauss    !! number of gaussian to fit
    integer, intent(in) :: dim_yx     !! spatial dimension
    real(xp), intent(in), dimension(:,:,:), allocatable :: grid !! parameters to optimize with cube mean at each iteration
    character(len=512), intent(in) :: fileout   !! name of the output result

    integer :: nside                    !! size of the reshaped data \(2^{nside}\)
    character(len=512) :: fileout_nside   !! name of the output at level nside

    integer :: ios=0 !! ios integer
    integer :: i     !! loop index
    integer :: j     !! loop index
    integer :: k     !! loop index

    fileout_nside = trim(fileout(:len_trim(fileout)-4)) // "_nside_" // trim(str(nside)) // ".dat"

    ! Open file
    open(unit=12, file=fileout_nside, action="write", iostat=ios)
    if (ios /= 0) stop "opening file error"

    write(12,fmt=*) "# i, j, A, mean, sigma"

    do i=1, dim_yx
       do j=1, dim_yx
          do k=1, n_gauss
             write(12,fmt=*) i-1, j-1, grid(1+((k-1)*3),i,j), grid(2+((k-1)*3),i,j), grid(3+((k-1)*3),i,j)
          enddo
       end do
    end do

    close(12)
  end subroutine save_process

end Module mod_inout
