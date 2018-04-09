module mod_functions

  use mod_constants
  use mod_optimize
  use mod_array
  
  implicit none

  private
  
  public :: mean_array, mean_map, dim2nside, dim_data2dim_cube, reshape_up, reshape_down, go_up_level, init_spectrum, &
       upgrade, update, set_stdmap

contains
    
  ! Compute nside value depending on dim_y and dim_x 
  pure function dim2nside(dim_cube)
    implicit none

    integer :: dim2nside
    integer, intent(in), dimension(3) :: dim_cube
    dim2nside = max(0,int(ceiling(log(real(dim_cube(2))) / log(2.))), int(ceiling(log(real(dim_cube(3))) / log(2.))))
    return
  end function dim2nside

  
  subroutine dim_data2dim_cube(nside, dim_data, dim_cube)
    implicit none
    
    integer, intent(in), dimension(3) :: dim_data
    integer, intent(inout), dimension(3) :: dim_cube
    integer :: nside
    
    dim_cube(1) = dim_data(1)
    dim_cube(2) = 2**nside
    dim_cube(3) = dim_cube(2) 
  end subroutine dim_data2dim_cube


  ! Reshape the data in a grid of 2**nside 
  subroutine reshape_up(data, cube, dim_data, dim_cube)
    implicit none
    
    real(xp), intent(in), dimension(:,:,:), allocatable    :: data
    real(xp), intent(inout), dimension(:,:,:), allocatable :: cube
    integer, intent(in), dimension(3) :: dim_data
    integer, intent(in), dimension(3) :: dim_cube
 
    integer :: offset_w, offset_h
 
    offset_w = (dim_cube(2) - dim_data(2)) / 2
    offset_h = (dim_cube(3) - dim_data(3)) / 2
    cube(:, offset_w+1:offset_w+dim_data(2), offset_h+1:offset_h+dim_data(3)) = data
  end subroutine reshape_up

  ! Reshape the data in a grid with original dimension 
  subroutine reshape_down(cube, data, dim_cube, dim_data)
    implicit none

    real(xp), intent(in), dimension(:,:,:), allocatable :: cube
    real(xp), intent(inout), dimension(:,:,:), allocatable    :: data
    integer, intent(in), dimension(3) :: dim_data
    integer, intent(in), dimension(3) :: dim_cube
 
    integer :: offset_w, offset_h
 
    offset_w = (dim_cube(2) - dim_data(2)) / 2
    offset_h = (dim_cube(3) - dim_data(3)) / 2
    data = cube(:, offset_w+1:offset_w+dim_data(2), offset_h+1:offset_h+dim_data(3)) 
  end subroutine reshape_down


  ! Average cube depending on level n 
  subroutine mean_array(nside, cube, cube_mean)
    implicit none

    integer, intent(in) :: nside
    real(xp), intent(in), dimension(:,:,:), allocatable :: cube
    real(xp), intent(inout), dimension(:,:,:), allocatable :: cube_mean

    integer :: i, j, k, l, n
    real(xp), dimension(:), allocatable :: spectrum

    allocate(spectrum(size(cube,dim=1)))
    spectrum = 0.
    
    n = size(cube, dim=2) / nside

    do i=1,size(cube_mean,dim=2)
       do j=1,size(cube_mean,dim=3)
          do k=1,n
             do l=1,n
                spectrum = spectrum + cube(:,k+((i-1)*n),l+((j-1)*n))
             enddo
          enddo
          spectrum = spectrum / (n**2)
          cube_mean(:,i,j) = spectrum
          spectrum = 0.
       enddo
    enddo
    
  end subroutine mean_array

  
  ! Average map depending on level n 
  subroutine mean_map(nside, map, map_mean)
    implicit none

    integer, intent(in) :: nside
    real(xp), intent(in), dimension(:,:), allocatable :: map
    real(xp), intent(inout), dimension(:,:), allocatable :: map_mean

    integer :: i, j, k, l, n
    real(xp) :: val

    val = 0.
    
    n = size(map, dim=2) / nside

    do i=1,size(map_mean,dim=1)
       do j=1,size(map_mean,dim=2)
          do k=1,n
             do l=1,n
                val = val + map(k+((i-1)*n),l+((j-1)*n))
             enddo
          enddo
          val = val / (n**2)
          map_mean(i,j) = val
          val = 0.
       enddo
    enddo
    
  end subroutine mean_map


  ! Projection of the solution at next level 
  subroutine go_up_level(cube_params)
    implicit none

    real(xp), intent(inout), dimension(:,:,:), allocatable :: cube_params

    integer :: i, j, k, l
    real(xp), dimension(:,:,:), allocatable :: cube_params_down
    integer, dimension(3) :: dim

    dim = shape(cube_params)
    allocate(cube_params_down(dim(1), dim(2), dim(3)))
    cube_params_down = 0._xp
    
    cube_params_down = cube_params
    
    deallocate(cube_params)
    allocate(cube_params(dim(1),dim(2)*2, dim(3)*2))
    cube_params = 0._xp
    
    do i=1,size(cube_params_down,dim=2)
       do j=1,size(cube_params_down,dim=3)
          do k=1,2
             do l=1,2
                cube_params(:,k+((i-1)*2),l+((j-1)*2)) = cube_params_down(:,i,j)
             enddo
          enddo
       enddo
    enddo
    
  end subroutine go_up_level

  
  ! Init mean spreturm with N Gaussian
  subroutine init_spectrum(n_gauss, params, dim_v, line, maxiter, m, iprint)
    implicit none
    
    integer, intent(in) :: n_gauss, dim_v, maxiter, m, iprint
    real(xp), intent(in), dimension(dim_v) :: line
    real(xp), intent(inout), dimension(3*n_gauss)  :: params

    integer :: i, j, k, p
    real(xp), dimension(:), allocatable :: lb, ub
    real(xp), dimension(dim_v) :: model, residual
    real(xp), dimension(:), allocatable :: x
    
    do i=1, n_gauss
       allocate(lb(3*i), ub(3*i))
       model = 0._xp
       residual = 0._xp
       lb = 0._xp; ub=0._xp
       
       call init_bounds(line, i, dim_v, lb, ub)

       do j=1, i
          do k=1, dim_v
             model(k) = model(k) + gaussian(k, params(1+(3*(j-1))), params(2+(3*(j-1))), params(3+(3*(j-1))))
          end do
       enddo

       residual = model - line
       
       allocate(x(3*i))
       x = 0._xp
       
       do p=1, 3*(i-1)
          x(p) = params(p);
       end do
       
       x(2+(3*(i-1))) = minloc(residual, dim_v)
       x(1+(3*(i-1))) = line(int(x(2+(3*(i-1))))) * 2._xp/3._xp
       x(3+(3*(i-1))) = 5._xp;
       
       call minimize_spec(3*i, m, x, lb, ub, line, dim_v, i, maxiter, iprint)
       
       do p=1, 3*i   
          params(p) = x(p);
       end do
       
       deallocate(x)
       deallocate(lb, ub)
    enddo
  end subroutine init_spectrum
  

  ! Init bounds for optimization
  subroutine init_bounds(line, n_gauss, dim_v, lb, ub)
    implicit none
    
    integer, intent(in) :: n_gauss, dim_v
    real(xp), intent(in), dimension(dim_v) :: line    
    real(xp), intent(inout), dimension(3*n_gauss) :: lb, ub
    
    integer :: i
    real(xp) :: max_line

    max_line = 0._xp
    max_line = maxval(line)
    
    do i=1, n_gauss       
       ! amplitude bounds
       lb(1+(3*(i-1))) = 0._xp;
       ub(1+(3*(i-1))) = max_line;
       
       ! mean bounds 
       lb(2+(3*(i-1))) = 0._xp;
       ub(2+(3*(i-1))) = dim_v;
       
       ! sigma bounds 
       lb(3+(3*(i-1))) = 0.001_xp;
       ub(3+(3*(i-1))) = 100._xp;
    end do
  end subroutine init_bounds


  ! Upgrade parameters using minimize function (here based on L-BFGS-B optimization module)
  subroutine upgrade(cube, params, power, n_gauss, dim_v, maxiter, m, iprint)
    implicit none

    real(xp), intent(in), dimension(:,:,:), allocatable :: cube
    integer, intent(in) :: power, n_gauss, dim_v, maxiter, m, iprint
    real(xp), intent(inout), dimension(:,:,:), allocatable :: params

    integer :: i,j
    real(xp), dimension(:), allocatable :: line
    real(xp), dimension(:), allocatable :: x
    real(xp), dimension(:), allocatable :: lb, ub

    do i=1, power
       do j=1, power
          ! print*, (i-1)*power+j, " / ", power*power
          allocate(line(dim_v))
          allocate(x(3*n_gauss), lb(3*n_gauss), ub(3*n_gauss))

          line = cube(:,i,j)
          x = params(:,i,j)
          
          call init_bounds(line, n_gauss, dim_v, lb, ub)
          call minimize_spec(3*n_gauss, m, x, lb, ub, line, dim_v, n_gauss, maxiter, iprint)
          
          params(:,i,j) = x
          
          deallocate(line)
          deallocate(x, lb, ub)
       end do
    end do
  end subroutine upgrade

  ! Update parameters using minimize function (here based on L-BFGS-B optimization module)
  subroutine update(cube, params, n_gauss, dim_v, dim_y, dim_x, lambda_amp, lambda_mu, lambda_sig, maxiter, m, kernel, &
       iprint, std_map)
    implicit none

    real(xp), intent(in), dimension(:,:,:), allocatable :: cube
    real(xp), intent(in), dimension(:,:), allocatable :: std_map
    real(xp), intent(in), dimension(:,:), allocatable :: kernel
    integer, intent(in) :: dim_v, dim_y, dim_x
    integer, intent(in) :: n_gauss, maxiter, m, iprint
    real(xp), intent(in) :: lambda_amp, lambda_mu, lambda_sig
    real(xp), intent(inout), dimension(:,:,:), allocatable :: params
    
    integer :: i,j
    integer :: n_beta
    real(xp), dimension(:,:,:), allocatable :: lb_3D, ub_3D
    real(xp), dimension(:), allocatable :: lb, ub
    real(xp), dimension(:), allocatable :: beta

    n_beta = 3*n_gauss * dim_y * dim_x

    allocate(lb(n_beta), ub(n_beta), beta(n_beta))
    allocate(lb_3D(3*n_gauss,dim_y,dim_x), ub_3D(3*n_gauss,dim_y,dim_x))

    do j=1, dim_x
       do i=1, dim_y
          call init_bounds(cube(:,i,j), n_gauss, dim_v, lb_3D(:,i,j), ub_3D(:,i,j))
       end do
    end do

    call ravel_3D(lb_3D, lb, 3*n_gauss, dim_y, dim_x)
    call ravel_3D(ub_3D, ub, 3*n_gauss, dim_y, dim_x)
    call ravel_3D(params, beta, 3*n_gauss, dim_y, dim_x)

    call minimize(n_beta, m, beta, lb, ub, cube, n_gauss, dim_v, dim_y, dim_x, lambda_amp, lambda_mu, lambda_sig, maxiter, &
         kernel, iprint, std_map)

    call unravel_3D(beta, params, 3*n_gauss, dim_y, dim_x)
        
  end subroutine update


  ! Compute the STD map of a 3D array
  subroutine set_stdmap(std_map, cube, lb, ub)
    implicit none

    integer, intent(in) :: lb, ub
    real(xp), intent(in), dimension(:,:,:), allocatable :: cube
    real(xp), intent(inout), dimension(:,:), allocatable :: std_map
    real(xp), dimension(:), allocatable :: line
    integer, dimension(3) :: dim_cube
    integer :: i, j 

    dim_cube = shape(cube)

    do j=1, dim_cube(3)
       do i=1, dim_cube(2)
          line = cube(lb:ub,i,j)
          std_map(i,j) = std(line)
       end do
    end do
    
  end subroutine set_stdmap


  ! Compute the STD of a 1D array
  pure function std(array)
    implicit none

    real(xp) :: std
    real(xp), intent(in), dimension(:) :: array
    integer :: i
    integer :: n
    real(xp) :: mean, var

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
  
end module mod_functions
