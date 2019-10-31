!! This module contains optimization subroutine and parametric model
module mod_optimize
  !! This module contains optimization subroutine and parametric model
  use mod_constants
  use mod_array

  implicit none
  
  private

  public :: gaussian, myfunc_spec, myresidual, mygrad_spec, f_g_cube_fast
  
contains
  
  pure function gaussian(x, a, m, s)
    !! Gaussian function   
    implicit none
    
    integer, intent(in) :: x
    real(xp), intent(in) :: a, m, s
    real(xp) :: gaussian

    gaussian = a * exp(-( (real(x,xp) - m)**2 ) / (2._xp * s**2) );
  end function gaussian
  
  ! Compute the residual between model and data
  subroutine myresidual(params, line, residual, n_mbb, dim_v)
    implicit none

    integer, intent(in) :: dim_v, n_mbb
    real(xp), intent(in), dimension(dim_v) :: line
    real(xp), intent(in), dimension(3*n_mbb) :: params
    real(xp), intent(inout), dimension(:), allocatable :: residual

    integer :: i, k
    real(xp) :: g    
    real(xp), dimension(dim_v) :: model

    g = 0._xp
    model = 0._xp
    
    do i=1, n_mbb
       do k=1, dim_v
          g = gaussian(k, params(1+(3*(i-1))), params(2+(3*(i-1))), params(3+(3*(i-1))))
          model(k) = model(k) + g
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

  
  ! Griadient of the objective function to minimize for a spectrum
  subroutine mygrad_spec(n_mbb, gradient, residual, params, dim_v)
    implicit none

    integer, intent(in) :: n_mbb, dim_v
    real(xp), intent(in), dimension(3*n_mbb) :: params
    real(xp), intent(in), dimension(:), allocatable :: residual
    real(xp), intent(inout), dimension(3*n_mbb) :: gradient

    integer :: i, k
    real(xp) :: g

    real(xp), dimension(:,:), allocatable :: dF_over_dB

    allocate(dF_over_dB(3*n_mbb, dim_v))

    g = 0._xp
    dF_over_dB = 0._xp
    gradient = 0._xp

    do i=1, n_mbb
       do k=1, dim_v          
          dF_over_dB(1+(3*(i-1)),k) = dF_over_dB(1+(3*(i-1)),k) +&
               exp( ( -(real(k,xp) - params(2+(3*(i-1))))**2._xp) / (2._xp * params(3+(3*(i-1)))**2._xp))

          dF_over_dB(2+(3*(i-1)),k) = dF_over_dB(2+(3*(i-1)),k) +&
               params(1+(3*(i-1))) * ( real(k,xp) - params(2+(3*(i-1))) ) / (params(3+(3*(i-1)))**2._xp) *&
               exp( ( -(real(k,xp) - params(2+(3*(i-1))))**2._xp) / (2._xp * params(3+(3*(i-1)))**2._xp))
          
          dF_over_dB(3+(3*(i-1)),k) = dF_over_dB(3+(3*(i-1)),k) +&
               params(1+(3*(i-1))) * ( real(k,xp) - params(2+(3*(i-1))) )**2._xp / (params(3+(3*(i-1)))**3._xp) *&
               exp( ( -(real(k,xp) - params(2+(3*(i-1))))**2._xp) / (2._xp * params(3+(3*(i-1)))**2._xp))
       enddo
    enddo
    
    do i=1, dim_v
       do k=1, 3*n_mbb
          gradient(k) = gradient(k) + dF_over_dB(k,i) * residual(i)
       end do
    end do

    deallocate(dF_over_dB)
  end subroutine mygrad_spec

  
  ! Compute the objective function for a cube and the gradient of the obkective function
  subroutine f_g_cube_fast(f, g, cube, beta, dim_v, dim_y, dim_x, n_mbb, kernel, lambda_sig, lambda_beta, lambda_Td, &
       lambda_var_sig, lambda_var_beta, lambda_var_Td, std_map)
    implicit none

    integer, intent(in) :: n_mbb
    integer, intent(in) :: dim_v, dim_y, dim_x
    real(xp), intent(in) :: lambda_sig, lambda_beta, lambda_Td
    real(xp), intent(in) :: lambda_var_sig, lambda_var_beta, lambda_var_Td
    real(xp), intent(in), dimension(:), allocatable :: beta
    real(xp), intent(in), dimension(:,:,:), allocatable :: cube
    real(xp), intent(in), dimension(:,:), allocatable :: kernel
    real(xp), intent(in), dimension(:,:), allocatable :: std_map
    real(xp), intent(inout) :: f
    real(xp), intent(inout), dimension(:), allocatable :: g

    integer :: i, j, k, l
    integer :: n_beta
    real(xp), dimension(:,:,:), allocatable :: residual
    real(xp), dimension(:), allocatable :: residual_1D
    real(xp), dimension(:,:,:), allocatable :: params
    real(xp), dimension(:), allocatable :: b_params
    real(xp), dimension(:,:), allocatable :: conv_sig, conv_beta, conv_Td
    real(xp), dimension(:,:), allocatable :: conv_conv_sig, conv_conv_beta, conv_conv_Td
    real(xp), dimension(:,:), allocatable :: image_sig, image_beta, image_Td
    real(xp), dimension(:,:,:), allocatable :: deriv
    real(xp), dimension(:), allocatable :: model
    real(xp) :: gauss

    allocate(deriv(3*n_mbb, dim_y, dim_x))
    allocate(residual(dim_v, dim_y, dim_x))
    allocate(b_params(n_mbb))
    allocate(params(3*n_mbb, dim_y, dim_x))
    allocate(conv_sig(dim_y, dim_x), conv_beta(dim_y, dim_x), conv_Td(dim_y, dim_x))
    allocate(conv_conv_sig(dim_y, dim_x), conv_conv_beta(dim_y, dim_x), conv_conv_Td(dim_y, dim_x))
    allocate(image_sig(dim_y, dim_x), image_beta(dim_y, dim_x), image_Td(dim_y, dim_x))
    allocate(model(dim_v))

    deriv = 0._xp
    f = 0._xp
    g = 0._xp
    residual = 0._xp    
    params = 0._xp
    model = 0._xp
    gauss = 0._xp
    
    n_beta = (3*n_mbb * dim_y * dim_x) + n_mbb

    call unravel_3D(beta, params, 3*n_mbb, dim_y, dim_x)    
    do i=1,n_mbb
       b_params(i) = beta((n_beta-n_mbb)+i)
    end do

    ! Compute the objective function and the gradient
    do j=1, dim_x
       do i=1, dim_y
          allocate(residual_1D(dim_v))
          residual_1D = 0._xp
          call myresidual(params(:,i,j), cube(:,i,j), residual_1D, n_mbb, dim_v)
          residual(:,i,j) = residual_1D
          if (std_map(i,j) > 0._xp) then
             f = f + (myfunc_spec(residual_1D)/std_map(i,j)**2._xp)
          end if
          deallocate(residual_1D)
       end do
    end do

    ! Compute the objective function and the gradient
    do i=1, n_mbb
       !
       conv_sig = 0._xp; conv_beta = 0._xp; conv_Td = 0._xp
       conv_conv_sig = 0._xp; conv_conv_beta = 0._xp; conv_conv_Td = 0._xp
       image_sig = 0._xp; image_beta = 0._xp; image_Td = 0._xp
       
       image_sig = params(1+(3*(i-1)),:,:)
       image_beta = params(2+(3*(i-1)),:,:)
       image_Td = params(3+(3*(i-1)),:,:)
       
       call convolution_2D_mirror(image_sig, conv_sig, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(image_beta, conv_beta, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(image_Td, conv_Td, dim_y, dim_x, kernel, 3)
       
       call convolution_2D_mirror(conv_sig, conv_conv_sig, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(conv_beta, conv_conv_beta, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(conv_Td, conv_conv_Td, dim_y, dim_x, kernel, 3)
       
       do l=1, dim_x
          do j=1, dim_y
             !Regularization
             f = f + (0.5_xp * lambda_sig * conv_sig(j,l)**2)
             f = f + (0.5_xp * lambda_beta * conv_beta(j,l)**2)
             f = f + (0.5_xp * lambda_Td * conv_Td(j,l)**2) + (0.5_xp * lambda_var_Td * (image_Td(j,l) - b_params(i))**2._xp)
             
             g((n_beta-n_mbb)+i) = g((n_beta-n_mbb)+i) - (lambda_var_Td * (image_Td(j,l) - b_params(i)))        
             
             !
             do k=1, dim_v                          
                if (std_map(j,l) > 0._xp) then
                   deriv(1+(3*(i-1)),j,l) = deriv(1+(3*(i-1)),j,l) + (exp( ( -(real(k,xp) - params(2+(3*(i-1)),j,l))**2._xp) &
                        / (2._xp * params(3+(3*(i-1)),j,l)**2._xp))) &
                        * (residual(k,j,l)/std_map(j,l)**2._xp) 

                   deriv(2+(3*(i-1)),j,l) = deriv(2+(3*(i-1)),j,l) + (params(1+(3*(i-1)),j,l) * &
                        ( real(k,xp) - params(2+(3*(i-1)),j,l) ) / (params(3+(3*(i-1)),j,l)**2._xp) * &
                        exp( ( -(real(k,xp) - params(2+(3*(i-1)),j,l))**2._xp) &
                        / (2._xp * params(3+(3*(i-1)),j,l)**2._xp))) &
                        * (residual(k,j,l)/std_map(j,l)**2._xp) 

                   deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + (params(1+(3*(i-1)),j,l) * &
                        ( real(k,xp) - params(2+(3*(i-1)),j,l) )**2._xp / (params(3+(3*(i-1)),j,l)**3._xp) * &
                        exp( ( -(real(k,xp) - params(2+(3*(i-1)),j,l))**2._xp) / (2._xp * params(3+(3*(i-1)),j,l)**2._xp))) &
                        * (residual(k,j,l)/std_map(j,l)**2._xp)
                end if
             end do

             deriv(1+(3*(i-1)),j,l) = deriv(1+(3*(i-1)),j,l) + (lambda_sig * conv_conv_sig(j,l))
             deriv(2+(3*(i-1)),j,l) = deriv(2+(3*(i-1)),j,l) + (lambda_beta * conv_conv_beta(j,l))
             deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + (lambda_Td * conv_conv_Td(j,l) + &
                  (lambda_var_Td * (image_Td(j,l) - b_params(i))))
          end do
          !
       end do
    end do        
    
    call ravel_3D(deriv, g, 3*n_mbb, dim_y, dim_x)

    deallocate(deriv)
    deallocate(residual)
    deallocate(b_params)
    deallocate(params)
    deallocate(conv_sig, conv_beta, conv_Td)
    deallocate(conv_conv_sig, conv_conv_beta, conv_conv_Td)
    deallocate(image_sig, image_beta, image_Td)

  end subroutine f_g_cube_fast


end module mod_optimize
