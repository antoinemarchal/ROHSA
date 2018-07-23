!! This module read the input user parameters (parameters.txt file / data / noise if true)
module mod_inout
  !! This module read the input user parameters (parameters.txt file / data / noise if true)
  use mod_constants
  
  implicit none
  
  private
  
  public :: read_cube, read_map, read_parameters!, read_fits_real_3d

contains
  
  subroutine read_parameters(filename_parameters, filename, fileout, filename_noise, n_gauss, n_gauss_add, &
       lambda_amp, lambda_mu, lambda_sig, lambda_var_sig, amp_fact_init, sig_init, maxiter_init, maxiter, m, &
       noise, regul, lstd, ustd, iprint, iprint_init)
    implicit none

    integer :: ios=0

    character(len=512), intent(in) :: filename_parameters

    integer, intent(inout) :: n_gauss, n_gauss_add
    integer, intent(inout) :: m 
    integer, intent(inout) :: lstd, ustd
    integer, intent(inout) :: iprint, iprint_init
    integer, intent(inout) :: maxiter, maxiter_init
    real(xp), intent(inout) :: lambda_amp, lambda_mu, lambda_sig, lambda_var_sig
    real(xp), intent(inout) :: amp_fact_init, sig_init
    logical, intent(inout) :: noise, regul

    character(len=512), intent(inout) :: filename
    character(len=512), intent(inout) :: fileout
    character(len=512), intent(inout) :: filename_noise

    namelist /user_parameters/ filename, fileout, filename_noise, n_gauss, n_gauss_add, lambda_amp, lambda_mu, &
         & lambda_sig, lambda_var_sig, amp_fact_init, sig_init, maxiter_init, maxiter, m, noise, regul, lstd, ustd, &
         iprint, iprint_init
    
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

end Module mod_inout
