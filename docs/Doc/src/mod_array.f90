module mod_array
  use mod_constants
  implicit none

  private

  public :: convolution_2D_mirror, ravel_3D, unravel_3D
  
contains

  subroutine convolution_2D_mirror(image, conv, dim_y, dim_x, kernel, dim_k)
    implicit none

    real(xp), intent(in), dimension(:,:), allocatable :: image
    real(xp), intent(inout), dimension(:,:), allocatable :: conv
    real(xp), intent(in), dimension(:,:), allocatable :: kernel
    integer, intent(in) :: dim_y, dim_x, dim_k

    integer :: i, j, ii, jj, m, n, mm, nn
    integer :: kCenterY, kCenterX
    real(xp), dimension(:,:), allocatable :: ext_conv, extended

    allocate(ext_conv(dim_y+4, dim_x+4))
    allocate(extended(dim_y+4, dim_x+4))

    ii = 0; jj = 0; mm = 0; nn = 0
    kCenterY = 0; kCenterX = 0
    ext_conv = 0._xp
    extended = 0._xp
    
    do j=1, dim_x
       do i=1, dim_y
          extended(2+i,2+j) = image(i,j)
       end do
    end do

    do j=1, 2
       do i=1, dim_y
          extended(2+i,j) = image(i,j);
       end do
    end do
    
    do i=1, 2
       do j=1, dim_x
          extended(i,2+j) = image(i,j);
       end do
    end do

    do j=dim_x+1, dim_x+2
       do i=1, dim_y
          extended(2+i,2+j) = image(i,j-2)
       end do
    end do

    do j=1, dim_x
       do i=dim_y+1, dim_y+2
          extended(2+i,2+j) = image(i-2,j)
       end do
    end do
    
    kCenterY = dim_k / 2 + 1
    kCenterX = kCenterY

    do j=1, dim_x+4
       do i=1, dim_y+4
          do m=1, dim_k
             mm = dim_k - m + 1
             
             do n=1, dim_k
                nn = dim_k - n + 1

                ii = i + (m - kCenterY)
                jj = j + (n - kCenterX)
                
                if( ii >= 1 .and. ii < dim_y+4 .and. jj >= 1 .and. jj < dim_x+4 ) then
                   ext_conv(i,j) = ext_conv(i,j) + extended(ii,jj) * kernel(mm,nn)
                end if
             end do
             
          end do
       end do
    end do

    do j=1, dim_x
       do i=1, dim_y
          conv(i,j) = ext_conv(2+i,2+j)
       end do
    end do    
    
  end subroutine convolution_2D_mirror


  ! Return a contiguous flattened 1D array from a 3D array
  subroutine ravel_3D(cube, vector, dim_v, dim_y, dim_x)
    implicit none

    integer, intent(in) :: dim_v, dim_y, dim_x
    real(xp), intent(in), dimension(:,:,:), allocatable :: cube
    real(xp), intent(inout), dimension(:), allocatable :: vector

    integer :: i, j, k, i__

    i__ = 1
    
    do k=1, dim_x
       do j=1, dim_y
          do i=1, dim_v
             vector(i__) = cube(i,j,k)
             i__ = i__ + 1
          end do
       end do
    end do
  end subroutine ravel_3D

  ! Return a 3D array from a contiguous flattened 1D array 
  subroutine unravel_3D(vector, cube, dim_v, dim_y, dim_x)
    implicit none

    integer, intent(in) :: dim_v, dim_y, dim_x
    real(xp), intent(in), dimension(:), allocatable :: vector
    real(xp), intent(inout), dimension(:,:,:), allocatable :: cube

    integer :: i, j, k, i__

    i__ = 1
    
    do k=1, dim_x
       do j=1, dim_y
          do i=1, dim_v
             cube(i,j,k) = vector(i__)
             i__ = i__ + 1
          end do
       end do
    end do
  end subroutine unravel_3D

end module mod_array
