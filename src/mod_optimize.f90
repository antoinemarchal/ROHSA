!! This module contains optimization subroutine and parametric model
module mod_optimize
  !! This module contains optimization subroutine and parametric model
  use mod_constants
  use mod_array
  use mod_read_parameters
  use mod_model

  implicit none
  
  private

  public :: myfunc_spec, f_g_cube_fast, myresidual, mygrad_spec
  
contains
    
  ! Compute the residual between model and data
  subroutine myresidual(pars, line, residual, dim_v, n)
    implicit none

    integer, intent(in) :: n
    integer, intent(in) :: dim_v
    real(xp), intent(in), dimension(dim_v) :: line
    real(xp), intent(in), dimension(2*params%n) :: pars

    real(xp), dimension(:), allocatable :: residual
    real(xp), dimension(:), allocatable :: model

    integer :: i, k

    allocate(model(dim_v))

    model = 0._xp;
    
    do i=1, n
       do k=1, dim_v
          model(k) = model(k) + f_rmsf(rm(k), pars(1+(2*(i-1))), pars(2+(2*(i-1))), coeff)
       enddo
    enddo

    residual = model - line
  end subroutine myresidual


  ! Objective function to minimize for a spectrum
  pure function  myfunc_spec(residual)
    implicit none
    
    real(xp), intent(in), dimension(:), allocatable :: residual
    real(xp) :: myfunc_spec
    
    myfunc_spec = 0._xp
    
    myfunc_spec = 0.5_xp * sum(residual**2._xp)
  end function  myfunc_spec


  ! Gradient of the objective function to minimize for a spectrum
  subroutine mygrad_spec(n, gradient, residual, pars, dim_v)
    implicit none

    integer, intent(in) :: n, dim_v
    real(xp), intent(in), dimension(2*n) :: pars
    real(xp), intent(in), dimension(:), allocatable :: residual
    real(xp), intent(inout), dimension(2*n) :: gradient

    integer :: i, k
    real(xp) :: g

    real(xp), dimension(:,:), allocatable :: dL

    allocate(dL(2*n, dim_v))

    g = 0._xp
    dL = 0._xp
    gradient = 0._xp

    do i=1, n
       do k=1, dim_v          
          dL(1+(2*(i-1)),k) = dL(1+(2*(i-1)),k) +&
               df_rmsf_da(rm(real(k,xp)), pars(2+(2*(i-1))), coeff)

          dL(2+(2*(i-1)),k) = dL(2+(2*(i-1)),k) +&
               df_rmsf_dm(rm(real(k,xp)), pars(1+(2*(i-1))), pars(2+(2*(i-1))), coeff)
       enddo
    enddo
    
    do i=1, dim_v
       do k=1, 2*n
          gradient(k) = gradient(k) + dL(k,i) * residual(i)
       end do
    end do

    deallocate(dL)
  end subroutine mygrad_spec

  
  ! Compute the objective function for a cube and the gradient of the obkective function
  subroutine f_g_cube_fast(f, g, cube, beta, dim_v, dim_y, dim_x, kernel, std_map)
    implicit none

    integer, intent(in) :: dim_v, dim_y, dim_x
    real(xp), intent(in), dimension(:), allocatable :: beta
    real(xp), intent(in), dimension(:,:), allocatable :: kernel
    real(xp), intent(in), dimension(:,:), allocatable :: std_map
    real(xp), intent(in), dimension(:,:,:), allocatable :: cube
    real(xp), intent(inout) :: f
    real(xp), intent(inout), dimension(:), allocatable :: g

    integer :: i, j, k, l
    integer :: n_beta
    real(xp), dimension(:,:,:), allocatable :: residual
    real(xp), dimension(:), allocatable :: residual_1D
    real(xp), dimension(:,:,:), allocatable :: pars
    real(xp), dimension(:,:), allocatable :: conv_amp, conv_mu
    real(xp), dimension(:,:), allocatable :: conv_conv_amp, conv_conv_mu
    real(xp), dimension(:,:), allocatable :: image_amp, image_mu
    real(xp), dimension(:,:,:), allocatable :: deriv
    real(xp), dimension(:), allocatable :: model

    allocate(deriv(2*params%n, dim_y, dim_x))
    allocate(residual(dim_v, dim_y, dim_x))
    allocate(pars(2*params%n, dim_y, dim_x))
    allocate(conv_amp(dim_y, dim_x), conv_mu(dim_y, dim_x))
    allocate(conv_conv_amp(dim_y, dim_x), conv_conv_mu(dim_y, dim_x))
    allocate(image_amp(dim_y, dim_x), image_mu(dim_y, dim_x))
    allocate(model(dim_v))

    deriv = 0._xp
    f = 0._xp
    g = 0._xp
    residual = 0._xp    
    pars = 0._xp
    model = 0._xp
    
    n_beta = (2*params%n * dim_y * dim_x)

    call unravel_3D(beta, pars, 2*params%n, dim_y, dim_x)    

    ! Compute the objective function and the gradient
    do j=1, dim_x
       do i=1, dim_y
          allocate(residual_1D(dim_v))
          residual_1D = 0._xp
          call myresidual(pars(:,i,j), cube(:,i,j), residual_1D, dim_v, params%n) !FIXMEEEEEEEE
          residual(:,i,j) = residual_1D
          if (std_map(i,j) > 0._xp) then
             f = f + (myfunc_spec(residual_1D)/std_map(i,j)**2._xp) !FIXMEEEEEEEEEEE
          end if
          deallocate(residual_1D)
       end do
    end do

    ! Compute the objective function and the gradient
    do i=1, params%n
       !
       conv_amp = 0._xp; conv_mu = 0._xp
       conv_conv_amp = 0._xp; conv_conv_mu = 0._xp
       image_amp = 0._xp; image_mu = 0._xp
       
       image_amp = pars(1+(2*(i-1)),:,:)
       image_mu = pars(2+(2*(i-1)),:,:)
       
       call convolution_2D_mirror(image_amp, conv_amp, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(image_mu, conv_mu, dim_y, dim_x, kernel, 3)
       
       call convolution_2D_mirror(conv_amp, conv_conv_amp, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(conv_mu, conv_conv_mu, dim_y, dim_x, kernel, 3)
       
       do l=1, dim_x
          do j=1, dim_y
             !Regularization
             f = f + (0.5_xp * params%lambda_amp * conv_amp(j,l)**2)
             f = f + (0.5_xp * params%lambda_mu * conv_mu(j,l)**2)                          
             !
             do k=1, dim_v                          
                if (std_map(j,l) > 0._xp) then
                   deriv(1+(2*(i-1)),j,l) = deriv(1+(2*(i-1)),j,l) +&
                        df_rmsf_da(rm(real(k,xp)), pars(2+(2*(i-1)),j,l), coeff) &
                   * (residual(k,j,l)/std_map(j,l)**2._xp) 

                   deriv(2+(2*(i-1)),j,l) = deriv(2+(2*(i-1)),j,l) +&
                        df_rmsf_dm(rm(real(k,xp)), pars(1+(2*(i-1)),j,l), pars(2+(2*(i-1)),j,l), coeff) &
                   * (residual(k,j,l)/std_map(j,l)**2._xp) 
                end if
             end do

             deriv(1+(2*(i-1)),j,l) = deriv(1+(2*(i-1)),j,l) + (params%lambda_amp * conv_conv_amp(j,l))
             deriv(2+(2*(i-1)),j,l) = deriv(2+(2*(i-1)),j,l) + (params%lambda_mu * conv_conv_mu(j,l))
          end do
          !
       end do
    end do        
    
    call ravel_3D(deriv, g, 2*params%n, dim_y, dim_x)

    deallocate(deriv)
    deallocate(residual)
    deallocate(pars)
    deallocate(conv_amp, conv_mu)
    deallocate(conv_conv_amp, conv_conv_mu)
    deallocate(image_amp, image_mu)

  end subroutine f_g_cube_fast


end module mod_optimize

  ! ! Compute the residual between model and data
  ! subroutine myresidual_rmsf(pars, line, residual, dim_v)
  !   implicit none

  !   integer, intent(in) :: dim_v
  !   real(xp), intent(in), dimension(:), allocatable :: line
  !   real(xp), intent(in), dimension(2*params%n_rmsf) :: pars
  !   real(xp), intent(inout), dimension(:), allocatable :: residual

  !   real(xp), dimension(:), allocatable :: model

  !   integer :: i, k

  !   allocate(model(dim_v))

  !   model = 0._xp
    
  !   do i=1, params%n_rmsf
  !      do k=1, dim_v
  !         model(k) = model(k) + gaussian(rm(k), pars(1+(2*(i-1))), pars(2+(2*(i-1))))
  !      enddo
  !   enddo

  !   residual = model - line
  ! end subroutine myresidual_rmsf

  ! ! Gradient of the objective function to minimize for a spectrum
  ! subroutine mygrad_rmsf(gradient, residual, pars, dim_v)
  !   implicit none

  !   integer, intent(in) :: dim_v
  !   real(xp), intent(in), dimension(2*params%n_rmsf) :: pars
  !   real(xp), intent(in), dimension(:), allocatable :: residual
  !   real(xp), intent(inout), dimension(2*params%n_rmsf) :: gradient

  !   integer :: i, k
  !   real(xp) :: g

  !   real(xp), dimension(:,:), allocatable :: dL

  !   allocate(dL(2*params%n_rmsf, dim_v))

  !   g = 0._xp
  !   dL = 0._xp
  !   gradient = 0._xp
    
  !   do i=1, params%n_rmsf
  !      do k=1, dim_v          
  !         dL(1+(2*(i-1)),k) = dL(1+(2*(i-1)),k) +&
  !              exp( ( - rm(k)**2._xp) / (2._xp * pars(2+(2*(i-1)))**2._xp))

  !         dL(2+(2*(i-1)),k) = dL(2+(2*(i-1)),k) +&
  !              pars(1+(2*(i-1))) * rm(k)**2._xp / pars(2+(2*(i-1)))**3._xp *&
  !              exp( ( - rm(k)**2._xp) / (2._xp * pars(2+(2*(i-1)))**2._xp))
  !      enddo
  !   enddo
    
  !   do i=1, dim_v
  !      do k=1, 2*params%n_rmsf
  !         gradient(k) = gradient(k) + dL(k,i) * residual(i)
  !      end do
  !   end do

  !   deallocate(dL)
  ! end subroutine mygrad_rmsf
