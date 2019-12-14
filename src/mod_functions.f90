!! This module contains the main routines of ROHSA
module mod_functions
  !! This module contains the main routines of ROHSA
  use mod_constants
  use mod_minimize
  use mod_optimize
  use mod_array
  use mod_model
  use mod_read_parameters
  use mod_fft
  
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

  
  subroutine init_spectrum(pars, dim_v, line, wavelength, color, std)
    !! Initialization of the mean sprectrum with N Gaussian
    implicit none
    
    integer, intent(in) :: dim_v !! dimension along v axis

    real(xp), intent(in), dimension(dim_v) :: line !! spectrum
    real(xp), intent(in), dimension(dim_v)  :: wavelength  !! wavelength Planck + IRAS
    real(xp), intent(in), dimension(:,:), allocatable :: color
    real(xp), intent(in), dimension(:) :: std

    real(xp), intent(inout), dimension(3*params%n_mbb)  :: pars !! params to optimize

    integer :: i
    real(xp), dimension(:), allocatable :: lb, ub
    real(xp), dimension(dim_v) :: model, residual
    real(xp), dimension(:), allocatable :: x

    allocate(lb(3*params%n_mbb), ub(3*params%n_mbb))
    allocate(x(3*params%n_mbb))
    model = 0._xp
    residual = 0._xp
    lb = 0._xp; ub=0._xp
    x = 0._xp
            
    do i=1, 3*params%n_mbb
       x(i) = pars(i)
    end do

    call init_bounds(lb, ub)   
    call minimize_spec(3*params%n_mbb, params%m, x, lb, ub, line, wavelength, dim_v, params%maxiter, &
         params%iprint, color, std)

    do i=1, 3*params%n_mbb
       pars(i) = x(i)
    end do
    
    deallocate(x)
    deallocate(lb, ub)
  end subroutine init_spectrum
  

  subroutine init_bounds(lb, ub)
    !! Initialize parameters bounds for optimization
    implicit none

    real(xp), intent(inout), dimension(3*params%n_mbb) :: lb !! lower bounds
    real(xp), intent(inout), dimension(3*params%n_mbb) :: ub !! upper bounds
    
    integer :: i
    
    do i=1, params%n_mbb       
       ! tau bounds
       lb(1+(3*(i-1))) = params%lb_tau
       ub(1+(3*(i-1))) = params%ub_tau

       ! beta bounds 
       lb(2+(3*(i-1))) = params%lb_beta
       ub(2+(3*(i-1))) = params%ub_beta

       ! Td bounds 
       lb(3+(3*(i-1))) = params%lb_Td
       ub(3+(3*(i-1))) = params%ub_Td
       ! end if
    end do

    if (params%ciba .eqv. .true.) then
          ! tau bounds cib
          lb(1) = params%lb_tau_cib
          ub(1) = params%ub_tau_cib

          ! beta bounds cib
          lb(2) = params%lb_beta_cib
          ub(2) = params%ub_beta_cib

          ! Td bounds cib
          lb(3) = params%lb_Td_cib
          ub(3) = params%ub_Td_cib
    end if

  end subroutine init_bounds


  subroutine upgrade(cube, pars, wavelength, power, dim_v, color, std_cube)
    !! Upgrade parameters (spectra to spectra) using minimize function (here based on L-BFGS-B optimization module)
    implicit none

    real(xp), intent(in), dimension(:,:,:), allocatable :: cube !! cube
    real(xp), intent(in), dimension(dim_v)  :: wavelength  !! wavelength Planck + IRAS
    real(xp), intent(in), dimension(:,:), allocatable :: color
    real(xp), intent(in), dimension(:,:,:), allocatable :: std_cube !! standard deviation cube

    integer, intent(in) :: power !! nside of the cube
    integer, intent(in) :: dim_v !! dimension along v axis

    real(xp), intent(inout), dimension(:,:,:), allocatable :: pars !! cube parameters to update

    integer :: i,j
    real(xp), dimension(:), allocatable :: line
    real(xp), dimension(:), allocatable :: x
    real(xp), dimension(:), allocatable :: lb, ub
    real(xp), dimension(:), allocatable :: std

    do i=1, power
       do j=1, power
          ! print*, (i-1)*power+j, " / ", power*power
          allocate(line(dim_v))
          allocate(x(3*params%n_mbb), lb(3*params%n_mbb), ub(3*params%n_mbb))

          line = cube(:,i,j)
          x = pars(:,i,j)
          std = std_cube(:,i,j)
          
          call init_bounds(lb, ub)
          call minimize_spec(3*params%n_mbb, params%m, x, lb, ub, line, wavelength, dim_v, params%maxiter_init, &
               params%iprint_init, color, std)
          
          pars(:,i,j) = x
          
          deallocate(line)
          deallocate(x, lb, ub)
       end do
    end do
  end subroutine upgrade


  subroutine update(cube, cube_HI, wavelength, pars, b_pars, c_pars, d_pars, stefan_pars, &
       dim_v, dim_y, dim_x, kernel, std_cube, color)
    !! Update parameters (entire cube) using minimize function (here based on L-BFGS-B optimization module)
    implicit none
    
    real(xp), intent(in), dimension(:,:,:), allocatable :: cube        !! cube 
    real(xp), intent(in), dimension(:,:,:), allocatable :: cube_HI     !! cube HI
    real(xp), intent(in), dimension(:), allocatable     :: wavelength  !! wavelength Planck + IRAS
    real(xp), intent(in), dimension(:,:), allocatable   :: color       !! polynomial coefficient for color correction
    real(xp), intent(in), dimension(:,:,:), allocatable :: std_cube    !! Standard deviation cube
    real(xp), intent(in), dimension(:,:), allocatable :: kernel !! convolution kernel

    integer, intent(in) :: dim_v !! dimension along v axis
    integer, intent(in) :: dim_y !! dimension along spatial axis y 
    integer, intent(in) :: dim_x !! dimension along spatial axis x

    real(xp), intent(inout), dimension(:), allocatable :: b_pars !! unknown average Tdma
    real(xp), intent(inout), dimension(:), allocatable :: c_pars !! unknown average Tdma
    real(xp), intent(inout), dimension(:), allocatable :: d_pars !! unknown average Tdma
    real(xp), intent(inout), dimension(:), allocatable :: stefan_pars !! 
    real(xp), intent(inout), dimension(:,:,:), allocatable :: pars !! parameters cube to update
    
    integer :: i,j
    integer :: n_beta
    integer :: n_cube

    real(xp), dimension(:,:,:), allocatable :: lb_3D, ub_3D
    real(xp), dimension(:), allocatable :: lb, ub
    real(xp), dimension(:), allocatable :: beta

    real(xp), dimension(:,:), allocatable :: tapper 
    real(xp), dimension(:,:), allocatable :: kmat
    real(xp), dimension(:,:), allocatable :: butter, filter

    n_beta = (3*params%n_mbb * dim_y * dim_x) + (4*params%n_mbb)
    n_cube = (3*params%n_mbb * dim_y * dim_x)

    allocate(lb(n_beta), ub(n_beta), beta(n_beta))
    allocate(lb_3D(3*params%n_mbb,dim_y,dim_x), ub_3D(3*params%n_mbb,dim_y,dim_x))

    !Bounds
    do j=1, dim_x
       do i=1, dim_y
          call init_bounds(lb_3D(:,i,j), ub_3D(:,i,j))
       end do
    end do
    
    call ravel_3D(lb_3D, lb, 3*params%n_mbb, dim_y, dim_x)
    call ravel_3D(ub_3D, ub, 3*params%n_mbb, dim_y, dim_x)
    call ravel_3D(pars, beta, 3*params%n_mbb, dim_y, dim_x)

    do i=1,params%n_mbb
       lb(n_cube+(0*params%n_mbb)+i) = params%lb_tau
       ub(n_cube+(0*params%n_mbb)+i) = params%ub_tau

       lb(n_cube+(1*params%n_mbb)+i) = params%lb_tau !FIXME MAYBE
       ub(n_cube+(1*params%n_mbb)+i) = params%ub_tau !FIXME MAYBE

       lb(n_cube+(2*params%n_mbb)+i) = params%lb_beta
       ub(n_cube+(2*params%n_mbb)+i) = params%ub_beta

       lb(n_cube+(3*params%n_mbb)+i) = params%lb_Td
       ub(n_cube+(3*params%n_mbb)+i) = params%ub_Td

       beta(n_cube+(0*params%n_mbb)+i) = b_pars(i)
       beta(n_cube+(1*params%n_mbb)+i) = stefan_pars(i)
       beta(n_cube+(2*params%n_mbb)+i) = c_pars(i)
       beta(n_cube+(3*params%n_mbb)+i) = d_pars(i)
    end do

    if (params%ciba .eqv. .true.) then
       lb(n_cube+(0*params%n_mbb)+1) = params%lb_tau_cib
       ub(n_cube+(0*params%n_mbb)+1) = params%ub_tau_cib

       lb(n_cube+(1*params%n_mbb)+1) = params%lb_tau_cib
       ub(n_cube+(1*params%n_mbb)+1) = params%ub_tau_cib

       lb(n_cube+(2*params%n_mbb)+1) = params%lb_beta_cib
       ub(n_cube+(2*params%n_mbb)+1) = params%ub_beta_cib

       lb(n_cube+(3*params%n_mbb)+1) = params%lb_Td_cib
       ub(n_cube+(3*params%n_mbb)+1) = params%ub_Td_cib
    end if

    !Compute tapper and kmat
    allocate(tapper(dim_y,dim_x))
    allocate(kmat(dim_y,dim_x))
    
    !Compute tapper and kmat
    call kgrid(dim_y,dim_x,kmat)
    call apodize(tapper, params%radius_tapper,dim_y,dim_x)

    !compute butterworth map and filter
    call butterworth(butter,kmat,1._xp,1._xp,0.5_xp)
    filter = 1._xp / butter

    !run minimization
    call minimize(n_beta, params%m, beta, lb, ub, cube, cube_HI, dim_v, dim_y, dim_x, params%maxiter, &
         kernel, params%iprint, std_cube, wavelength, color, filter, tapper)

    !Unravel data
    call unravel_3D(beta, pars, 3*params%n_mbb, dim_y, dim_x)
    do i=1,params%n_mbb
       b_pars(i) = beta(n_cube+(0*params%n_mbb)+i)
       stefan_pars(i) = beta(n_cube+(1*params%n_mbb)+i)
       c_pars(i) = beta(n_cube+(2*params%n_mbb)+i)       
       d_pars(i) = beta(n_cube+(3*params%n_mbb)+i)       
    end do        

    ! print*, b_pars
    ! print*, stefan_pars
    ! print*, c_pars
    ! print*, d_pars

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


  subroutine init_grid_params(pars, guess_spectrum, dim_y, dim_x)
    !! Set up a grid params array with std spectrum at each spatial position
    implicit none
 
    real(xp), intent(inout), dimension(:,:,:), allocatable :: pars !! grid of paramters
    real(xp), intent(in), dimension(:), allocatable :: guess_spectrum !! std spectrum of the observation
    integer, intent(in) :: dim_y !! dimension along spatial axis y 
    integer, intent(in) :: dim_x !! dimension along spatial axis x

    integer :: i !! index loop
    integer :: j !! index loop

    do j=1, dim_x
       do i=1, dim_y
          pars(:,i,j) = guess_spectrum
       end do
    end do

  end subroutine init_grid_params


end module mod_functions
