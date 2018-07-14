!! This module contains the main routines of ROHSA
module mod_functions
  !! This module contains the main routines of ROHSA
  use mod_constants
  use mod_optimize
  use mod_array
  
  implicit none

  private
  
  public :: mean_array, mean_map, dim2nside, dim_data2dim_cube, reshape_up, reshape_down, go_up_level, init_spectrum, &
       upgrade, update, set_stdmap

contains
    
  pure function dim2nside(dim_cube)
  !! Compute nside value from \(dim_y\) and \(dim_x\) 
    implicit none

    integer :: dim2nside !! nside of the cube
    integer, intent(in), dimension(3) :: dim_cube !! cube dimension
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


  subroutine reshape_up(data, cube, dim_data, dim_cube)
  !! Reshape the data in a grid of \( 2^{nside} \)
    implicit none
    
    real(xp), intent(in), dimension(:,:,:), allocatable    :: data !! original cube
    real(xp), intent(inout), dimension(:,:,:), allocatable :: cube !! reshape cube
    integer, intent(in), dimension(3) :: dim_data !! original cube dimension
    integer, intent(in), dimension(3) :: dim_cube !! new cube dimension
 
    integer :: offset_w, offset_h
 
    offset_w = (dim_cube(2) - dim_data(2)) / 2
    offset_h = (dim_cube(3) - dim_data(3)) / 2
    cube(:, offset_w+1:offset_w+dim_data(2), offset_h+1:offset_h+dim_data(3)) = data
  end subroutine reshape_up

  subroutine reshape_down(cube, data, dim_cube, dim_data)
    !! Reshape the cube (\( 2^{nside} \)) into a grid with original dimension (opposite of reshape_up)
    implicit none

    real(xp), intent(in), dimension(:,:,:), allocatable :: cube !! original cube
    real(xp), intent(inout), dimension(:,:,:), allocatable :: data !! reshape cube
    integer, intent(in), dimension(3) :: dim_data !! original cube dimension
    integer, intent(in), dimension(3) :: dim_cube !! new cube dimension
 
    integer :: offset_w, offset_h
 
    offset_w = (dim_cube(2) - dim_data(2)) / 2
    offset_h = (dim_cube(3) - dim_data(3)) / 2
    data = cube(:, offset_w+1:offset_w+dim_data(2), offset_h+1:offset_h+dim_data(3)) 
  end subroutine reshape_down


  subroutine mean_array(nside, cube, cube_mean)
    !! Average cube along spatial axis depending on level n 
    implicit none

    integer, intent(in) :: nside !! nside of the cube
    real(xp), intent(in), dimension(:,:,:), allocatable :: cube !! cube
    real(xp), intent(inout), dimension(:,:,:), allocatable :: cube_mean !! average cube

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

  
  subroutine mean_map(nside, map, map_mean)
    !! Average map depending on level nside 
    implicit none

    integer, intent(in) :: nside !! nside
    real(xp), intent(in), dimension(:,:), allocatable :: map !! map
    real(xp), intent(inout), dimension(:,:), allocatable :: map_mean !! avarage map

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


  subroutine go_up_level(cube_params)
    !! Projection of the solution at next level (nside += 1)
    implicit none

    real(xp), intent(inout), dimension(:,:,:), allocatable :: cube_params !! cube of parameters

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

  
  subroutine init_spectrum(n_gauss, params, dim_v, line, maxiter, m, iprint)
    !! Initialization of the mean sprectrum with N Gaussian
    implicit none
    
    integer, intent(in) :: n_gauss !! number of Gaussian
    integer, intent(in) :: dim_v !! dimension along v axis
    integer, intent(in) :: maxiter !! Max number of iteration
    integer, intent(in) :: m !! number of corrections used in the limited memory matrix by LBFGS-B
    integer, intent(in) :: iprint !! print option

    real(xp), intent(in), dimension(dim_v) :: line !! spectrum
    real(xp), intent(inout), dimension(3*n_gauss)  :: params !! params to optimize

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
  

  subroutine init_bounds(line, n_gauss, dim_v, lb, ub)
    !! Initialize parameters bounds for optimization
    implicit none
    
    integer, intent(in) :: n_gauss !! number of Gaussian
    integer, intent(in) :: dim_v !! dimension along v axis
    real(xp), intent(in), dimension(dim_v) :: line !! spectrum   
    real(xp), intent(inout), dimension(3*n_gauss) :: lb !! lower bounds
    real(xp), intent(inout), dimension(3*n_gauss) :: ub !! upper bounds
    
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


  subroutine upgrade(cube, params, power, n_gauss, dim_v, maxiter, m, iprint)
    !! Upgrade parameters (spectra to spectra) using minimize function (here based on L-BFGS-B optimization module)
    implicit none

    real(xp), intent(in), dimension(:,:,:), allocatable :: cube !! cube
    integer, intent(in) :: power !! nside of the cube
    integer, intent(in) :: n_gauss !! number of Gaussian
    integer, intent(in) :: dim_v !! dimension along v axis
    integer, intent(in) :: maxiter !! max number of iteration
    integer, intent(in) :: m !! number of corrections used in the limited memory matrix by LBFGS-B
    integer, intent(in) :: iprint !! print option

    real(xp), intent(inout), dimension(:,:,:), allocatable :: params !! cube parameters to update

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


  subroutine update(cube, params, n_gauss, dim_v, dim_y, dim_x, lambda_amp, lambda_mu, lambda_sig, lambda_var_sig, maxiter, &
       m, kernel, iprint, std_map)
    !! Update parameters (entire cube) using minimize function (here based on L-BFGS-B optimization module)
    implicit none
    
    real(xp), intent(in), dimension(:,:,:), allocatable :: cube !! cube 
    real(xp), intent(in), dimension(:,:), allocatable :: std_map !! Standard deviation map 
    real(xp), intent(in), dimension(:,:), allocatable :: kernel !! convolution kernel
    integer, intent(in) :: dim_v !! dimension along v axis
    integer, intent(in) :: dim_y !! dimension along spatial axis y 
    integer, intent(in) :: dim_x !! dimension along spatial axis x
    integer, intent(in) :: n_gauss !! Number of Gaussian
    integer, intent(in) :: maxiter !! max number of iteration
    integer, intent(in) :: m !! number of corrections used in the limited memory matrix by LBFGS-B
    integer, intent(in) :: iprint !! print option
    real(xp), intent(in) :: lambda_amp !! lambda for amplitude parameter
    real(xp), intent(in) :: lambda_mu !! lambda for mean position parameter
    real(xp), intent(in) :: lambda_sig !! lambda for dispersion parameter
    real(xp), intent(in) :: lambda_var_sig !! lambda for variance dispersion parameter

    real(xp), intent(inout), dimension(:,:,:), allocatable :: params !! parameters cube to update
    
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

    call minimize(n_beta, m, beta, lb, ub, cube, n_gauss, dim_v, dim_y, dim_x, lambda_amp, lambda_mu, lambda_sig, &
         lambda_var_sig, maxiter, kernel, iprint, std_map)

    call unravel_3D(beta, params, 3*n_gauss, dim_y, dim_x)
        
  end subroutine update


  subroutine set_stdmap(std_map, cube, lb, ub)
    !! Compute the STD map of a 3D array
    implicit none

    integer, intent(in) :: lb !! lower bound 
    integer, intent(in) :: ub !! upper bound
    real(xp), intent(in), dimension(:,:,:), allocatable :: cube !! cube

    real(xp), intent(inout), dimension(:,:), allocatable :: std_map !! standard deviation map of the cube

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
  
end module mod_functions
