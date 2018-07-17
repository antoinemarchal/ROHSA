!! This module read the input user parameters (parameters.txt file / data / noise if true)
module mod_read
  !! This module read the input user parameters (parameters.txt file / data / noise if true)
  use mod_constants
  
  implicit none
  
  private
  
  public :: read_cube, read_map, read_parameters, read_fits_real_3d

contains
  
  subroutine read_parameters(filename_parameters, filename, fileout, filename_noise, n_gauss, n_gauss_add, &
       lambda_amp, lambda_mu, lambda_sig, lambda_var_sig, maxiter_init, maxiter, m, noise, regul, lstd, ustd, &
       iprint, iprint_init)
    implicit none

    integer :: ios=0

    character(len=512), intent(in) :: filename_parameters

    integer, intent(inout) :: n_gauss, n_gauss_add
    integer, intent(inout) :: m 
    integer, intent(inout) :: lstd, ustd
    integer, intent(inout) :: iprint, iprint_init
    integer, intent(inout) :: maxiter, maxiter_init
    real(xp), intent(inout) :: lambda_amp, lambda_mu, lambda_sig, lambda_var_sig
    logical, intent(inout) :: noise, regul

    character(len=512), intent(inout) :: filename
    character(len=512), intent(inout) :: fileout
    character(len=512), intent(inout) :: filename_noise

    namelist /user_parameters/ filename, fileout, filename_noise, n_gauss, n_gauss_add, lambda_amp, lambda_mu, &
         & lambda_sig, lambda_var_sig, maxiter_init, maxiter, m, noise, regul, lstd, ustd, iprint, iprint_init
    
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

  ! subroutine read_fits()
  !   implicit none
  !   integer status, unit, readwrite, blocksize, naxes, nfound, naxis 
  !   integer group, firstpix, nbuffer, npixels, i
  !   real datamin, datamax, nullval, buffer
  !   logical anynull
  !   character(len=52) :: filename   

  !   print*, "Start"
  !   status = 0
  !   filename = "GHIGLS_DFN_Tb.fits"
  !   readwrite = 0 
    
  !   call ftgiou(unit, status)
  !   call ftopen(unit, filename, readwrite, blocksize, status)
  !   print*, status
    
  ! end subroutine read_fits

  subroutine read_fits_real_3d (filename, arr, status, dim1, dim2, dim3, indfmt)
    !
    character(*),                intent(in)  :: filename
    real(xp), dimension(:,:,:), allocatable  :: arr
    real(xp), dimension(:,:), allocatable  :: foo
    integer,            intent(out) :: status
    integer,  optional, intent(out) :: dim1
    integer,  optional, intent(out) :: dim2
    integer,  optional, intent(out) :: dim3
    integer,  optional, intent(in)  :: indfmt
    !
    character(len=80)   :: comm
    integer :: unit, blk, i, d1, d2, d3, nele, ifmt
    logical          :: anyf
    ! ----------------------------------------------------------------------------
    ifmt = 0
    if (present(indfmt)) ifmt = indfmt
    !
    !			open the fits file.
    !
    status = 0
    i      = 0
    call ftgiou (unit, status)
    call ftopen (unit, filename, i, blk, status)
    !
    !			allocate space for the array in the primary hdu.
    !
    call ftgidm (unit, d1, status)
    if ((status .ne. 0) .or. (d1 .lt. 3)) then
       if (d1 .lt. 3) print 1, status, 'primary hdu has insufficient axes.'
       return
    end if
    !
    call ftgkyj (unit, 'naxis1', d1, comm, status)
    call ftgkyj (unit, 'naxis2', d2, comm, status)
    call ftgkyj (unit, 'naxis3', d3, comm, status)
    if (status .ne. 0) return
    if (present(dim1)) dim1 = d1
    if (present(dim2)) dim2 = d2
    if (present(dim3)) dim3 = d3
    !
    ! if (.not. associated(arr)) then
    !    select case (ifmt)
    !    case (7)
    !       allocate (arr(1:d1,     1:d2,     1:d3    ), stat=status)
    !    case (6)
    !       allocate (arr(0:(d1-1), 1:d2,     1:d3    ), stat=status)
    !    case (5)
    !       allocate (arr(1:d1,     0:(d2-1), 1:d3    ), stat=status)
    !    case (4)
    !       allocate (arr(0:(d1-1), 0:(d2-1), 1:d3    ), stat=status)
    !    case (3)
    !       allocate (arr(1:d1,     1:d2,     0:(d3-1)), stat=status)
    !    case (2)
    !       allocate (arr(0:(d1-1), 1:d2,     0:(d3-1)), stat=status)
    !    case (1)
    !       allocate (arr(1:d1,     0:(d2-1), 0:(d3-1)), stat=status)
    !    case default
    allocate (arr(0:(d1-1), 0:(d2-1), 0:(d3-1)), stat=status)
    !    end select
    !    if (status .ne. 0) return
    ! end if
    !
    !			fill the array.
    !
    allocate(foo(d1,d2))
    nele = d1 * d2
    call ftgpve (unit, 0, 1, nele, 0, foo, anyf, status)
    print*, foo
    !
    !			close the fits file.
    !
    call ftclos (unit, status)
    if (unit .ge. 50) call ftfiou (unit, status)
    !
    return
    ! ----------------------------------------------------------------------------
1   format ('read_fits_real_3d:  status = ',i10,3x,a)
    !
  end subroutine read_fits_real_3d
end Module mod_read
