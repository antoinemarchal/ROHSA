!! This module contains tools (convolution/ravel) to manipulate 2D/3Darray 
module mod_array
  !! This module contains tools (convolution/ravel) to manipulate 2D/3Darray 
  use mod_constants
  implicit none

  private

  public :: convolution_2D_mirror, ravel_2D, ravel_3D, unravel_3D, mean, std, std_2D, mean_2D, max_2D, &
       ravel_3D_abs, unravel_3D_abs, meshgrid, linspace, shift, apodize
  
  
contains

  subroutine apodize(tapper,radius,dimx,dimy)
    implicit none
    
    integer, intent(in) :: dimx, dimy
    real(xp), intent(in) :: radius
    real(xp), intent(inout), dimension(:,:), allocatable :: tapper

    integer :: nx=0, ny=0
    integer :: dni, dnj
    real(xp), dimension(:), allocatable :: tap1dx, tap1dy
    real(xp), dimension(:), allocatable :: xdni, ydnj

    real(xp) :: pi_over_2
    integer :: i

    nx = dimy
    ny = dimx

    allocate(tap1dx(nx),tap1dy(ny))

    tap1dx = 1._xp
    tap1dy = 1._xp

    if (radius .le. 0._xp .or. radius .ge. 1._xp) then
       print*, "Error: radius must be lower than 1 and greater than 0."
       stop
    end if
    
    dni = nx - int(radius*nx)
    dnj = ny - int(radius*ny)

    allocate(xdni(dni),ydnj(dnj))
    
    pi_over_2 = pi/2._xp

    call linspace(xdni,1._xp,real(dni,xp))
    call linspace(ydnj,1._xp,real(dnj,xp))

    xdni = xdni - 1._xp
    ydnj = ydnj - 1._xp

    tap1dx(1:dni) = (cos(3._xp * pi_over_2 + pi_over_2 * (1._xp * xdni / (dni-1))))
    tap1dx(nx-dni+1:) = (cos(0._xp + pi_over_2 * (1._xp * xdni / (dni-1))))
    tap1dy(1:dnj) = (cos(3._xp * pi_over_2 + pi_over_2 * (1._xp * ydnj / (dnj-1))))
    tap1dy(ny-dnj+1:) = (cos(0._xp + pi_over_2 * (1._xp * ydnj / (dnj-1))))

    do i=1,nx
       tapper(:,i) = tap1dy
    end do
    
    do i=1,ny
       tapper(i,:) = tapper(i,:) * tap1dx
    end do
  end subroutine apodize

  subroutine shift(x, xshift)
    implicit none

    real(xp), intent(in), dimension(:,:), allocatable :: x
    real(xp), intent(inout), dimension(:,:), allocatable :: xshift

    integer :: nx, ny
    integer :: i, j

    nx = size(x,2)
    ny = size(x,1)

    do i=1,nx
       do j=1,ny
          xshift(i,j) = x(i,j) * (-1._xp)**(real(i+j,xp))
       end do
    end do
    
  end subroutine shift


  subroutine meshgrid(x,y,xx,yy)
    implicit none

    real(xp), intent(in), dimension(:) :: x, y
    real(xp), intent(inout), dimension(:,:) :: xx, yy
    xx = spread(x, 1, size(y))
    yy = spread(y, 2, size(x))    
  end subroutine meshgrid

  
  subroutine linspace(vec,a,b)
    implicit none
    
    real(xp), intent(in) :: a,b
    real(xp), intent(inout), dimension(:) :: vec
    real(xp) :: h
    integer :: n
    integer :: i

    n=size(vec)

    h=(b-a)/(n-1)
    vec=a+h*(/(i,i=0,n-1)/)
  end subroutine linspace

  
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


  ! Return a contiguous flattened 1D array from a 2D array
  subroutine ravel_2D(map, vector, dim_y, dim_x)
    implicit none

    integer, intent(in) :: dim_y, dim_x
    real(xp), intent(in), dimension(:,:), allocatable :: map
    real(xp), intent(inout), dimension(:), allocatable :: vector

    integer :: j, k, i__

    i__ = 1
    
    do k=1, dim_x
       do j=1, dim_y
             vector(i__) = map(j,k)
             i__ = i__ + 1
       end do
    end do
  end subroutine ravel_2D

  
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


  ! Return a contiguous flattened 1D array from a 3D array
  subroutine ravel_3D_abs(cube, cube_abs, vector, dim_v, dim_y, dim_x)
    implicit none

    integer, intent(in) :: dim_v, dim_y, dim_x
    real(xp), intent(in), dimension(:,:,:), allocatable :: cube, cube_abs
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

    do k=1, dim_x
       do j=1, dim_y
          do i=1, dim_v
             vector(i__) = cube_abs(i,j,k)
             i__ = i__ + 1
          end do
       end do
    end do

  end subroutine ravel_3D_abs


  ! Return a 3D array from a contiguous flattened 1D array 
  subroutine unravel_3D_abs(vector, cube, cube_abs, dim_v, dim_y, dim_x)
    implicit none

    integer, intent(in) :: dim_v, dim_y, dim_x
    real(xp), intent(in), dimension(:), allocatable :: vector
    real(xp), intent(inout), dimension(:,:,:), allocatable :: cube, cube_abs

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

    do k=1, dim_x
       do j=1, dim_y
          do i=1, dim_v
             cube_abs(i,j,k) = vector(i__)
             i__ = i__ + 1
          end do
       end do
    end do

  end subroutine unravel_3D_abs


  pure function mean(array)
    !! Compute the mean of a 1D array
    implicit none

    real(xp), intent(in), dimension(:) :: array !! 1D array
    real(xp) :: mean

    mean = sum(array)/(max(1,size(array)))
    
    return
  end function mean    


  pure function std(array)
    !! Compute the STD of a 1D array
    implicit none

    real(xp), intent(in), dimension(:) :: array !! 1D array
    integer :: i
    integer :: n
    real(xp) :: std !! standard deviation 
    real(xp) :: mean
    real(xp) :: var

    mean = 0._xp; var = 0._xp
    std = 0._xp

    n = size(array)
    mean = sum(array) / n

    do i=1, n
       var = var + (array(i) - mean)**2._xp
    end do
    
    var = var / (n - 1)
    std = sqrt(var)
    
    return
  end function std


  function std_2D(map, dim_y, dim_x)
    !! Compute the STD of a 2D map
    implicit none

    integer, intent(in) :: dim_y !! dimension along spatial axis y
    integer, intent(in) :: dim_x !! dimension along spatial axis x
    real(xp), intent(in), dimension(:,:), allocatable :: map !! 2D array
    real(xp), dimension(:), allocatable :: vector !! 1D array 
    real(xp) :: std_2D

    allocate(vector(dim_y*dim_x))

    call ravel_2D(map, vector, dim_y, dim_x)
    std_2D = std(vector)

    deallocate(vector)

  end function std_2D


  function max_2D(map, dim_y, dim_x)
    !! Compute the MAX of a 2D map
    implicit none

    integer, intent(in) :: dim_y !! dimension along spatial axis y
    integer, intent(in) :: dim_x !! dimension along spatial axis x
    real(xp), intent(in), dimension(:,:), allocatable :: map !! 2D array
    real(xp), dimension(:), allocatable :: vector !! 1D array 
    real(xp) :: max_2D

    allocate(vector(dim_y*dim_x))

    call ravel_2D(map, vector, dim_y, dim_x)
    max_2D = maxval(vector)

    deallocate(vector)

  end function max_2D


  function mean_2D(map, dim_y, dim_x)
    !! Compute the MEAN of a 2D map
    implicit none

    integer, intent(in) :: dim_y !! dimension along spatial axis y
    integer, intent(in) :: dim_x !! dimension along spatial axis x
    real(xp), intent(in), dimension(:,:), allocatable :: map !! 2D array
    real(xp), dimension(:), allocatable :: vector !! 1D array 
    real(xp) :: mean_2D

    allocate(vector(dim_y*dim_x))

    call ravel_2D(map, vector, dim_y, dim_x)
    mean_2D = mean(vector)

    deallocate(vector)

  end function mean_2D


end module mod_array
