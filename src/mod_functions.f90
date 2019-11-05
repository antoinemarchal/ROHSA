!! This module contains the main routines of ROHSA
module mod_functions
  !! This module contains the main routines of ROHSA
  use mod_constants
  use mod_minimize
  use mod_optimize
  use mod_array
  use mod_model
  
  implicit none

  private
  
  public :: mean_array, mean_map, dim2nside, dim_data2dim_cube, reshape_up, reshape_down, go_up_level, init_spectrum, &
       upgrade, update, set_stdmap, std_spectrum, mean_spectrum, max_spectrum, init_grid_params

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

    deallocate(spectrum)
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

    deallocate(cube_params_down)
  end subroutine go_up_level

  
  subroutine init_spectrum(n_mbb, params, dim_v, line, NHI, wavelength, sig_fact_init, sig_init, beta_init, &
       Td_init, lb_sig, ub_sig, lb_beta, ub_beta, lb_Td, ub_Td, l0, maxiter, m, iprint)
    !! Initialization of the mean sprectrum with N Gaussian
    implicit none
    
    integer, intent(in) :: n_mbb !! number of Gaussian
    integer, intent(in) :: dim_v !! dimension along v axis
    integer, intent(in) :: maxiter !! Max number of iteration
    integer, intent(in) :: m !! number of corrections used in the limited memory matrix by LBFGS-B
    integer, intent(in) :: iprint !! print option

    real(xp), intent(in), dimension(dim_v) :: line !! spectrum
    real(xp), intent(in), dimension(dim_v)  :: wavelength  !! wavelength Planck + IRAS
    real(xp), intent(in), dimension(n_mbb) :: NHI !! spectrum
    real(xp), intent(in) :: sig_fact_init !! times max siglitude of additional Gaussian

    real(xp), intent(in) :: sig_init !! dispersion of additional Gaussian
    real(xp), intent(in) :: beta_init !! dispersion of additional Gaussian
    real(xp), intent(in) :: Td_init !! dispersion of additional Gaussian

    real(xp), intent(in) :: lb_sig !! lower bound sigma
    real(xp), intent(in) :: ub_sig !! upper bound sigma
    real(xp), intent(in) :: lb_beta !! lower bound beta
    real(xp), intent(in) :: ub_beta !! upper bound beta
    real(xp), intent(in) :: lb_Td !! lower bound Td
    real(xp), intent(in) :: ub_Td !! upper bound Td

    real(xp), intent(in) :: l0 !! reference wavelength

    real(xp), intent(inout), dimension(3*n_mbb)  :: params !! params to optimize

    integer :: i, j, k, p
    real(xp), dimension(:), allocatable :: lb, ub
    real(xp), dimension(dim_v) :: model, residual
    real(xp), dimension(:), allocatable :: x

    allocate(lb(3*n_mbb), ub(3*n_mbb))
    allocate(x(3*n_mbb))
    model = 0._xp
    residual = 0._xp
    lb = 0._xp; ub=0._xp
    x = 0._xp
            
    do p=1, 3*n_mbb
       x(p) = params(p)
    end do

    call init_bounds(line, n_mbb, dim_v, lb, ub, lb_sig, ub_sig, lb_beta, ub_beta, lb_Td, ub_Td)   
    call minimize_spec(3*n_mbb, m, x, lb, ub, line, NHI, wavelength, dim_v, n_mbb, l0, maxiter, iprint)

    do p=1, 3*n_mbb
       params(p) = x(p)
    end do
    
    deallocate(x)
    deallocate(lb, ub)
  end subroutine init_spectrum
  

  subroutine init_bounds(line, n_mbb, dim_v, lb, ub, lb_sig, ub_sig, lb_beta, ub_beta, lb_Td, ub_Td)
    !! Initialize parameters bounds for optimization
    implicit none
    
    integer, intent(in) :: n_mbb !! number of Gaussian
    integer, intent(in) :: dim_v !! dimension along v axis
    real(xp), intent(in) :: lb_sig !! lower bound sigma
    real(xp), intent(in) :: ub_sig !! upper bound sigma
    real(xp), intent(in) :: lb_beta !! lower bound beta
    real(xp), intent(in) :: ub_beta !! upper bound beta
    real(xp), intent(in) :: lb_Td !! lower bound Td
    real(xp), intent(in) :: ub_Td !! upper bound Td
    real(xp), intent(in), dimension(dim_v) :: line !! spectrum   
    real(xp), intent(inout), dimension(3*n_mbb) :: lb !! lower bounds
    real(xp), intent(inout), dimension(3*n_mbb) :: ub !! upper bounds
    
    integer :: i
    
    do i=1, n_mbb       
       ! sigma bounds
       lb(1+(3*(i-1))) = lb_sig
       ub(1+(3*(i-1))) = ub_sig
       
       ! beta bounds 
       lb(2+(3*(i-1))) = lb_beta
       ub(2+(3*(i-1))) = ub_beta
       
       ! Td bounds 
       lb(3+(3*(i-1))) = lb_Td
       ub(3+(3*(i-1))) = ub_Td
    end do
  end subroutine init_bounds


  subroutine upgrade(cube, params, NHI, wavelength, power, n_mbb, dim_v, lb_sig, ub_sig, lb_beta, ub_beta, &
       lb_Td, ub_Td, l0, maxiter, m, iprint)
    !! Upgrade parameters (spectra to spectra) using minimize function (here based on L-BFGS-B optimization module)
    implicit none

    real(xp), intent(in), dimension(:,:,:), allocatable :: cube !! cube
    real(xp), intent(in), dimension(dim_v)  :: wavelength  !! wavelength Planck + IRAS
    real(xp), intent(in), dimension(n_mbb) :: NHI !! spectrum

    real(xp), intent(in) :: lb_sig !! lower bound sigma
    real(xp), intent(in) :: ub_sig !! upper bound sigma
    real(xp), intent(in) :: lb_beta !! lower bound beta
    real(xp), intent(in) :: ub_beta !! upper bound beta
    real(xp), intent(in) :: lb_Td !! lower bound Td
    real(xp), intent(in) :: ub_Td !! upper bound Td
    integer, intent(in) :: power !! nside of the cube
    integer, intent(in) :: n_mbb !! number of Gaussian
    integer, intent(in) :: dim_v !! dimension along v axis
    integer, intent(in) :: maxiter !! max number of iteration
    integer, intent(in) :: m !! number of corrections used in the limited memory matrix by LBFGS-B
    integer, intent(in) :: iprint !! print option

    real(xp), intent(in) :: l0 !! reference wavelength

    real(xp), intent(inout), dimension(:,:,:), allocatable :: params !! cube parameters to update

    integer :: i,j
    real(xp), dimension(:), allocatable :: line
    real(xp), dimension(:), allocatable :: x
    real(xp), dimension(:), allocatable :: lb, ub

    do i=1, power
       do j=1, power
          ! print*, (i-1)*power+j, " / ", power*power
          allocate(line(dim_v))
          allocate(x(3*n_mbb), lb(3*n_mbb), ub(3*n_mbb))

          line = cube(:,i,j)
          x = params(:,i,j)
          
          call init_bounds(line, n_mbb, dim_v, lb, ub, lb_sig, ub_sig, lb_beta, ub_beta, lb_Td, ub_Td)
          call minimize_spec(3*n_mbb, m, x, lb, ub, line, NHI, wavelength, dim_v, n_mbb, l0, maxiter, iprint)
          
          params(:,i,j) = x
          
          deallocate(line)
          deallocate(x, lb, ub)
       end do
    end do
  end subroutine upgrade


  subroutine update(cube, cube_HI, wavelength, params, b_params, stefan_params, n_mbb, dim_v, dim_y, dim_x, &
       lambda_sig, lambda_beta, lambda_Td, lambda_var_sig, lambda_var_beta, lambda_stefan, lb_sig, ub_sig, &
       lb_beta, ub_beta, lb_Td, ub_Td, l0, maxiter, m, kernel, iprint, std_map)
    !! Update parameters (entire cube) using minimize function (here based on L-BFGS-B optimization module)
    implicit none
    
    real(xp), intent(in), dimension(:,:,:), allocatable :: cube !! cube 
    real(xp), intent(in), dimension(:,:,:), allocatable :: cube_HI !! cube HI
    real(xp), intent(in), dimension(:), allocatable     :: wavelength  !! wavelength Planck + IRAS
    real(xp), intent(in), dimension(:,:), allocatable :: std_map !! Standard deviation map 
    real(xp), intent(in), dimension(:,:), allocatable :: kernel !! convolution kernel
    integer, intent(in) :: dim_v !! dimension along v axis
    integer, intent(in) :: dim_y !! dimension along spatial axis y 
    integer, intent(in) :: dim_x !! dimension along spatial axis x
    integer, intent(in) :: n_mbb !! Number of Gaussian
    integer, intent(in) :: maxiter !! max number of iteration
    integer, intent(in) :: m !! number of corrections used in the limited memory matrix by LBFGS-B
    integer, intent(in) :: iprint !! print option

    real(xp), intent(in) :: lambda_sig !! lambda for siglitude parameter
    real(xp), intent(in) :: lambda_beta !! lambda for mean position parameter
    real(xp), intent(in) :: lambda_Td !! lambda for dispersion parameter

    real(xp), intent(in) :: lambda_var_sig !! lambda for sig dispersion parameter
    real(xp), intent(in) :: lambda_var_beta  !! lambda for mean position dispersion parameter
    real(xp), intent(in) :: lambda_stefan !! lambda for variance dispersion parameter

    real(xp), intent(in) :: lb_sig !! lower bound sigma
    real(xp), intent(in) :: ub_sig !! upper bound sigma
    real(xp), intent(in) :: lb_beta !! lower bound beta
    real(xp), intent(in) :: ub_beta !! upper bound beta
    real(xp), intent(in) :: lb_Td !! lower bound Td
    real(xp), intent(in) :: ub_Td !! upper bound Td

    real(xp), intent(in) :: l0 !! reference wavelength

    real(xp), intent(inout), dimension(:), allocatable :: b_params !! unknown average Tdma
    real(xp), intent(inout), dimension(:), allocatable :: stefan_params !! 
    real(xp), intent(inout), dimension(:,:,:), allocatable :: params !! parameters cube to update
    
    integer :: i,j
    integer :: n_beta
    integer :: n_cube

    real(xp), dimension(:,:,:), allocatable :: lb_3D, ub_3D
    real(xp), dimension(:), allocatable :: lb, ub
    real(xp), dimension(:), allocatable :: beta

    n_beta = (3*n_mbb * dim_y * dim_x) + (2*n_mbb)
    n_cube = (3*n_mbb * dim_y * dim_x)

    allocate(lb(n_beta), ub(n_beta), beta(n_beta))
    allocate(lb_3D(3*n_mbb,dim_y,dim_x), ub_3D(3*n_mbb,dim_y,dim_x))

    !Bounds
    do j=1, dim_x
       do i=1, dim_y
          call init_bounds(cube(:,i,j), n_mbb, dim_v, lb_3D(:,i,j), ub_3D(:,i,j), lb_sig, ub_sig, &
               lb_beta, ub_beta, lb_Td, ub_Td)
       end do
    end do

    call ravel_3D(lb_3D, lb, 3*n_mbb, dim_y, dim_x)
    call ravel_3D(ub_3D, ub, 3*n_mbb, dim_y, dim_x)
    call ravel_3D(params, beta, 3*n_mbb, dim_y, dim_x)

    do i=1,n_mbb
       lb(n_cube+i) = lb_sig
       ub(n_cube+i) = ub_sig
       beta(n_cube+i) = b_params(i)

       lb(n_cube+n_mbb+i) = lb_sig !FIXME MAYBE
       ub(n_cube+n_mbb+i) = ub_sig !FIXME MAYBE
       beta(n_cube+n_mbb+i) = stefan_params(i)
    end do

    call minimize(n_beta, m, beta, lb, ub, cube, cube_HI, n_mbb, dim_v, dim_y, dim_x, lambda_sig, lambda_beta, lambda_Td, &
         lambda_var_sig, lambda_var_beta, lambda_stefan, l0, maxiter, kernel, iprint, std_map, wavelength)

    !Unravel data
    call unravel_3D(beta, params, 3*n_mbb, dim_y, dim_x)
    do i=1,n_mbb
       b_params(i) = beta(n_cube+i)
       stefan_params(i) = beta(n_cube+n_mbb+i)
    end do        

    deallocate(lb, ub, beta)
    deallocate(lb_3D, ub_3D)
  end subroutine update


  subroutine set_stdmap(std_map, cube, lb, ub) !fixme allocate deallocate
    !! Compute the STD map of a 3D array
    implicit none

    integer, intent(in) :: lb !! lower bound 
    integer, intent(in) :: ub !! upper bound
    real(xp), intent(in), dimension(:,:,:), allocatable :: cube !! cube

    real(xp), intent(inout), dimension(:,:), allocatable :: std_map !! standard deviation map of the cube

    real(xp), dimension(:), allocatable :: line 
    integer, dimension(3) :: dim_cube
    integer :: i, j 

    allocate(line(ub-lb))

    dim_cube = shape(cube)

    do j=1, dim_cube(3)
       do i=1, dim_cube(2)
          line = cube(lb:ub,i,j)
          std_map(i,j) = std(line)
       end do
    end do

    deallocate(line)
  end subroutine set_stdmap

  
  subroutine std_spectrum(data, spectrum, dim_v, dim_y, dim_x)
    !! Compute the STD spectrum of a cube along the spatial axis
    implicit none
    
    real(xp), intent(in), dimension(:,:,:), allocatable :: data !! initial fits data
    integer, intent(in) :: dim_v !! dimension along v axis
    integer, intent(in) :: dim_y !! dimension along spatial axis y 
    integer, intent(in) :: dim_x !! dimension along spatial axis x

    real(xp), intent(inout), dimension(:), allocatable :: spectrum !! std_spectrum of the observation
    real(xp), dimension(:,:), allocatable :: map !! 2D array

    integer :: i !! loop index

    do i=1,dim_v
       allocate(map(dim_y,dim_x))
       map = data(i,:,:)
       spectrum(i) = std_2D(map, dim_y, dim_x)
       deallocate(map)
    end do
    
  end subroutine std_spectrum  


  subroutine max_spectrum(data, spectrum, dim_v, dim_y, dim_x, norm_value)
    !! Compute the MAX spectrum of a cube along the spatial axis
    implicit none

    real(xp), intent(in), dimension(:,:,:), allocatable :: data !! initial fits data
    real(xp), intent(in), optional :: norm_value !! max value of the mean spectrum to normalze the max spectrum if present
    integer, intent(in) :: dim_v !! dimension along v axis
    integer, intent(in) :: dim_y !! dimension along spatial axis y 
    integer, intent(in) :: dim_x !! dimension along spatial axis x

    real(xp), intent(inout), dimension(:), allocatable :: spectrum !! max_spectrum of the observation
    real(xp), dimension(:,:), allocatable :: map !! 2D array

    integer :: i !! loop index

    do i=1,dim_v
       allocate(map(dim_y,dim_x))
       map = data(i,:,:)
       spectrum(i) = max_2D(map, dim_y, dim_x) 
       deallocate(map)
    end do

    if (present(norm_value)) then
       spectrum = spectrum / (maxval(spectrum) / norm_value)
    end if    
  end subroutine max_spectrum  


  subroutine mean_spectrum(data, spectrum, dim_v, dim_y, dim_x)
    !! Compute the MEAN spectrum of a cube along the spatial axis
    implicit none
    
    real(xp), intent(in), dimension(:,:,:), allocatable :: data !! initial fits data
    integer, intent(in) :: dim_v !! dimension along v axis
    integer, intent(in) :: dim_y !! dimension along spatial axis y 
    integer, intent(in) :: dim_x !! dimension along spatial axis x

    real(xp), intent(inout), dimension(:), allocatable :: spectrum !! std_spectrum of the observation
    real(xp), dimension(:,:), allocatable :: map !! 2D array

    integer :: i !! loop index

    do i=1,dim_v
       allocate(map(dim_y,dim_x))
       map = data(i,:,:)
       spectrum(i) = mean_2D(map, dim_y, dim_x)
       deallocate(map)
    end do
    
  end subroutine mean_spectrum  


  subroutine init_grid_params(params, guess_spectrum, dim_y, dim_x)
    !! Set up a grid params array with std spectrum at each spatial position
    implicit none
 
    real(xp), intent(inout), dimension(:,:,:), allocatable :: params !! grid of paramters
    real(xp), intent(in), dimension(:), allocatable :: guess_spectrum !! std spectrum of the observation
    integer, intent(in) :: dim_y !! dimension along spatial axis y 
    integer, intent(in) :: dim_x !! dimension along spatial axis x

    integer :: i !! index loop
    integer :: j !! index loop

    do j=1, dim_x
       do i=1, dim_y
          params(:,i,j) = guess_spectrum
       end do
    end do

  end subroutine init_grid_params


end module mod_functions
