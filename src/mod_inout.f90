!! This module read the input user parameters (parameters.txt file / data / noise if true)
module mod_inout
  !! This module read the input user parameters (parameters.txt file / data / noise if true)
  use mod_constants
  use mod_convert
  
  implicit none
  
  private
  
  public :: read_cube, read_map, read_array, read_parameters, save_process

contains
  
  subroutine read_parameters(filename_parameters, filename, filename_NHI, filename_wavelength, fileout, timeout, &
       filename_noise, n_mbb, lambda_sig, lambda_beta, lambda_Td, lambda_var_sig, lambda_var_beta, lambda_var_Td, &
       sig_fact_init, sig_init, beta_init, Td_init, lb_sig, ub_sig, lb_beta, ub_beta, lb_Td, ub_Td, &
       maxiter_init, maxiter, m, noise, lstd, ustd, iprint, iprint_init, save_grid)
    implicit none

    integer :: ios=0

    character(len=512), intent(in) :: filename_parameters

    integer, intent(inout)  :: n_mbb
    integer, intent(inout)  :: m 
    integer, intent(inout)  :: lstd, ustd
    integer, intent(inout)  :: iprint, iprint_init
    integer, intent(inout)  :: maxiter, maxiter_init
    real(xp), intent(inout) :: lambda_sig, lambda_beta, lambda_Td
    real(xp), intent(inout) :: lambda_var_sig, lambda_var_beta, lambda_var_Td
    real(xp), intent(inout) :: sig_fact_init
    real(xp), intent(inout) :: sig_init, beta_init, Td_init
    real(xp), intent(inout) :: lb_sig, ub_sig
    real(xp), intent(inout) :: lb_beta, ub_beta
    real(xp), intent(inout) :: lb_Td, ub_Td
    logical, intent(inout)  :: noise, save_grid

    character(len=512), intent(inout) :: filename
    character(len=512), intent(inout) :: filename_NHI
    character(len=512), intent(inout) :: filename_wavelength
    character(len=512), intent(inout) :: fileout
    character(len=512), intent(inout) :: timeout
    character(len=512), intent(inout) :: filename_noise

    namelist /user_parameters/ filename, filename_NHI, filename_wavelength, fileout, timeout, filename_noise, &
         n_mbb, lambda_sig, lambda_beta, lambda_Td, lambda_var_sig, lambda_var_beta, lambda_var_Td, &
         sig_fact_init, sig_init, beta_init, Td_init, lb_sig, ub_sig, lb_beta, ub_beta, lb_Td, ub_Td, &
         maxiter_init, maxiter, m, noise, lstd, ustd, iprint, iprint_init, save_grid
    
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

  
  subroutine save_process(nside, n_mbb, grid, dim_yx, fileout)
    implicit none

    integer, intent(in) :: n_mbb    !! number of gaussian to fit
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

    write(12,fmt=*) "# i, j, sig, beta, Td"

    do i=1, dim_yx
       do j=1, dim_yx
          do k=1, n_mbb
             write(12,fmt=*) i-1, j-1, grid(1+((k-1)*3),i,j), grid(2+((k-1)*3),i,j), grid(3+((k-1)*3),i,j)
          enddo
       end do
    end do

    close(12)
  end subroutine save_process

end Module mod_inout
