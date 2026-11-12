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
  subroutine myresidual(params, line, residual, n_gauss, dim_v)
    implicit none

    integer, intent(in) :: dim_v, n_gauss
    real(xp), intent(in), dimension(dim_v) :: line
    real(xp), intent(in), dimension(3*n_gauss) :: params
    real(xp), intent(inout), dimension(:), allocatable :: residual

    integer :: i, k
    real(xp) :: g    
    real(xp), dimension(dim_v) :: model

    g = 0._xp
    model = 0._xp
    
    do i=1, n_gauss
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
    
    myfunc_spec = 0.5_xp * sum(residual**2)    
  end function  myfunc_spec

  
  ! Griadient of the objective function to minimize for a spectrum
  subroutine mygrad_spec(n_gauss, gradient, residual, params, dim_v)
    implicit none

    integer, intent(in) :: n_gauss, dim_v
    real(xp), intent(in), dimension(3*n_gauss) :: params
    real(xp), intent(in), dimension(:), allocatable :: residual
    real(xp), intent(inout), dimension(3*n_gauss) :: gradient

    integer :: i, k
    real(xp) :: g

    real(xp), dimension(:,:), allocatable :: dF_over_dB

    allocate(dF_over_dB(3*n_gauss, dim_v))

    g = 0._xp
    dF_over_dB = 0._xp
    gradient = 0._xp

    do i=1, n_gauss
       do k=1, dim_v          
          dF_over_dB(1+(3*(i-1)),k) = dF_over_dB(1+(3*(i-1)),k) +&
               exp( ( -(real(k,xp) - params(2+(3*(i-1))))**2) / (2._xp * params(3+(3*(i-1)))**2))

          dF_over_dB(2+(3*(i-1)),k) = dF_over_dB(2+(3*(i-1)),k) +&
               params(1+(3*(i-1))) * ( real(k,xp) - params(2+(3*(i-1))) ) / (params(3+(3*(i-1)))**2) *&
               exp( ( -(real(k,xp) - params(2+(3*(i-1))))**2) / (2._xp * params(3+(3*(i-1)))**2))
          
          dF_over_dB(3+(3*(i-1)),k) = dF_over_dB(3+(3*(i-1)),k) +&
               params(1+(3*(i-1))) * ( real(k,xp) - params(2+(3*(i-1))) )**2 / (params(3+(3*(i-1)))**3) *&
               exp( ( -(real(k,xp) - params(2+(3*(i-1))))**2) / (2._xp * params(3+(3*(i-1)))**2))
       enddo
    enddo
    
    do i=1, dim_v
       do k=1, 3*n_gauss
          gradient(k) = gradient(k) + dF_over_dB(k,i) * residual(i)
       end do
    end do

    deallocate(dF_over_dB)
  end subroutine mygrad_spec

  
  !! Data term of one spectrum: its share of f and the data part of the gradient of its 3*n_gauss
  !! parameters (d). Each Gaussian is evaluated once per channel (ee) and used for the model and for
  !! the three derivatives, so no residual cube is stored. The model, the residual and every
  !! derivative are summed in the same order as in the separate loops this replaces. Explicit-shape
  !! arguments let gfortran optimise the loops inside the OpenMP region.
  subroutine data_term(dim_v, n_gauss, p, line, std, fp, d)
    implicit none

    integer, intent(in) :: dim_v, n_gauss
    real(xp), intent(in) :: p(3*n_gauss), line(dim_v), std
    real(xp), intent(out) :: fp, d(3*n_gauss)

    integer :: k, l
    real(xp) :: rr
    real(xp) :: ee(dim_v, n_gauss), res(dim_v)

    res = 0._xp
    do l=1, n_gauss
       do k=1, dim_v
          ee(k,l) = exp( ( -(real(k,xp) - p(2+(3*(l-1))))**2) / (2._xp * p(3+(3*(l-1)))**2))
          res(k) = res(k) + p(1+(3*(l-1))) * ee(k,l)
       end do
    end do
    res = res - line

    d = 0._xp
    fp = 0._xp
    if (std > 0._xp) then
       fp = (0.5_xp * sum(res**2))/std**2
       do l=1, n_gauss
          do k=1, dim_v
             rr = res(k)/std**2
             d(1+(3*(l-1))) = d(1+(3*(l-1))) + ee(k,l) * rr
             d(2+(3*(l-1))) = d(2+(3*(l-1))) + (p(1+(3*(l-1))) * &
                  ( real(k,xp) - p(2+(3*(l-1))) ) / (p(3+(3*(l-1)))**2) * ee(k,l)) * rr
             d(3+(3*(l-1))) = d(3+(3*(l-1))) + (p(1+(3*(l-1))) * &
                  ( real(k,xp) - p(2+(3*(l-1))) )**2 / (p(3+(3*(l-1)))**3) * ee(k,l)) * rr
          end do
       end do
    end if
  end subroutine data_term


  ! Compute the objective function for a cube and the gradient of the obkective function
  !! The parameters are read from beta and the gradient is written to g in place: params(p,i,j)
  !! of the original code is beta(p + np*((i-1) + dim_y*(j-1))) with np = 3*n_gauss, and deriv
  !! has the same layout in g, so unravel_3D/ravel_3D were plain copies.
  !! OpenMP: f and the b_params gradient are summed in an order that does not depend on the number
  !! of threads (pixel terms of the data part in the original serial order, regularisation terms
  !! column by column and then over columns), so the result is the same for any OMP_NUM_THREADS.
  subroutine f_g_cube_fast(f, g, cube, beta, dim_v, dim_y, dim_x, n_gauss, kernel, lambda_amp, lambda_mu, lambda_sig, &
       lambda_var_amp, lambda_var_mu, lambda_var_sig, std_map)
    implicit none

    integer, intent(in) :: n_gauss
    integer, intent(in) :: dim_v, dim_y, dim_x
    real(xp), intent(in) :: lambda_amp, lambda_mu, lambda_sig
    real(xp), intent(in) :: lambda_var_amp, lambda_var_mu, lambda_var_sig
    real(xp), intent(in), dimension(:), allocatable :: beta
    real(xp), intent(in), dimension(:,:,:), allocatable :: cube
    real(xp), intent(in), dimension(:,:), allocatable :: kernel
    real(xp), intent(in), dimension(:,:), allocatable :: std_map
    real(xp), intent(inout) :: f
    real(xp), intent(inout), dimension(:), allocatable :: g

    integer :: i, j, l, np
    integer(8) :: o
    integer :: n_beta
    real(xp), dimension(:), allocatable :: b_params
    real(xp), dimension(:,:), allocatable :: conv_amp, conv_mu, conv_sig
    real(xp), dimension(:,:), allocatable :: conv_conv_amp, conv_conv_mu, conv_conv_sig
    real(xp), dimension(:,:), allocatable :: image_amp, image_mu, image_sig
    real(xp), dimension(:,:), allocatable :: fpix
    real(xp), dimension(:), allocatable :: frcol, gbcol
    real(xp) :: fr, gb

    np = 3*n_gauss
    n_beta = (np * dim_y * dim_x) + n_gauss

    allocate(b_params(n_gauss))
    allocate(conv_amp(dim_y, dim_x), conv_mu(dim_y, dim_x), conv_sig(dim_y, dim_x))
    allocate(conv_conv_amp(dim_y, dim_x), conv_conv_mu(dim_y, dim_x), conv_conv_sig(dim_y, dim_x))
    allocate(image_amp(dim_y, dim_x), image_mu(dim_y, dim_x), image_sig(dim_y, dim_x))
    allocate(fpix(dim_y, dim_x), frcol(dim_x), gbcol(dim_x))

    do i=1,n_gauss
       b_params(i) = beta((n_beta-n_gauss)+i)
    end do

    ! Data term, spectra in parallel; each one sets its own np entries of g.
    !$omp parallel do private(i, o) schedule(static)
    do j=1, dim_x
       do i=1, dim_y
          o = int(np,8) * (int(i-1,8) + int(dim_y,8)*int(j-1,8))
          call data_term(dim_v, n_gauss, beta(o+1:o+np), cube(:,i,j), std_map(i,j), fpix(i,j), g(o+1:o+np))
       end do
    end do
    !$omp end parallel do

    f = 0._xp
    do j=1, dim_x
       do i=1, dim_y
          if (std_map(i,j) > 0._xp) f = f + fpix(i,j)
       end do
    end do

    ! Regularisation, one Gaussian at a time
    do i=1, n_gauss
       !$omp parallel do private(j, o) schedule(static)
       do l=1, dim_x
          do j=1, dim_y
             o = int(np,8) * (int(j-1,8) + int(dim_y,8)*int(l-1,8)) + 3*(i-1)
             image_amp(j,l) = beta(o+1)
             image_mu(j,l) = beta(o+2)
             image_sig(j,l) = beta(o+3)
          end do
       end do
       !$omp end parallel do

       call convolution_2D_mirror(image_amp, conv_amp, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(image_mu, conv_mu, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(image_sig, conv_sig, dim_y, dim_x, kernel, 3)

       call convolution_2D_mirror(conv_amp, conv_conv_amp, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(conv_mu, conv_conv_mu, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(conv_sig, conv_conv_sig, dim_y, dim_x, kernel, 3)

       ! The regularisation terms of f and of the b_params gradient are summed per column, then
       ! over the columns in order.
       !$omp parallel do private(j, o, fr, gb) schedule(static)
       do l=1, dim_x
          fr = 0._xp
          gb = 0._xp
          do j=1, dim_y
             fr = fr + (0.5_xp * lambda_amp * conv_amp(j,l)**2)
             fr = fr + (0.5_xp * lambda_mu * conv_mu(j,l)**2)
             fr = fr + (0.5_xp * lambda_sig * conv_sig(j,l)**2) + (0.5_xp * lambda_var_sig * (image_sig(j,l) - b_params(i))**2)

             gb = gb - (lambda_var_sig * (image_sig(j,l) - b_params(i)))

             o = int(np,8) * (int(j-1,8) + int(dim_y,8)*int(l-1,8)) + 3*(i-1)
             g(o+1) = g(o+1) + (lambda_amp * conv_conv_amp(j,l))
             g(o+2) = g(o+2) + (lambda_mu * conv_conv_mu(j,l))
             g(o+3) = g(o+3) + (lambda_sig * conv_conv_sig(j,l) + &
                  (lambda_var_sig * (image_sig(j,l) - b_params(i))))
          end do
          frcol(l) = fr
          gbcol(l) = gb
       end do
       !$omp end parallel do

       fr = 0._xp
       gb = 0._xp
       do l=1, dim_x
          fr = fr + frcol(l)
          gb = gb + gbcol(l)
       end do
       f = f + fr
       g((n_beta-n_gauss)+i) = gb
    end do

    deallocate(b_params)
    deallocate(conv_amp, conv_mu, conv_sig)
    deallocate(conv_conv_amp, conv_conv_mu, conv_conv_sig)
    deallocate(image_amp, image_mu, image_sig)
    deallocate(fpix, frcol, gbcol)

  end subroutine f_g_cube_fast


 ! Compute the objective function for a cube and the gradient of the obkective function
  subroutine f_g_cube(f, g, cube, beta, dim_v, dim_y, dim_x, n_gauss, kernel, lambda_amp, lambda_mu, lambda_sig, &
       lambda_var_amp, lambda_var_mu, lambda_var_sig, std_map)
    implicit none

    integer, intent(in) :: n_gauss
    integer, intent(in) :: dim_v, dim_y, dim_x
    real(xp), intent(in) :: lambda_amp, lambda_mu, lambda_sig
    real(xp), intent(in) :: lambda_var_amp, lambda_var_mu, lambda_var_sig
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
    real(xp), dimension(:,:), allocatable :: conv_amp, conv_mu, conv_sig
    real(xp), dimension(:,:), allocatable :: conv_conv_amp, conv_conv_mu, conv_conv_sig
    real(xp), dimension(:,:), allocatable :: image_amp, image_mu, image_sig
    real(xp), dimension(:,:,:), allocatable :: g_3D
    real(xp), dimension(:,:,:,:), allocatable :: dF_over_dB
    real(xp), dimension(:,:,:), allocatable :: dR_over_dB
    real(xp), dimension(:,:,:), allocatable :: deriv

    allocate(dR_over_dB(3*n_gauss, dim_y, dim_x))
    allocate(dF_over_dB(3*n_gauss, dim_v, dim_y, dim_x))
    allocate(deriv(3*n_gauss, dim_y, dim_x))
    allocate(g_3D(3*n_gauss, dim_y, dim_x))
    allocate(residual(dim_v, dim_y, dim_x))
    allocate(b_params(n_gauss))
    allocate(params(3*n_gauss, dim_y, dim_x))
    allocate(conv_amp(dim_y, dim_x), conv_mu(dim_y, dim_x), conv_sig(dim_y, dim_x))
    allocate(conv_conv_amp(dim_y, dim_x), conv_conv_mu(dim_y, dim_x), conv_conv_sig(dim_y, dim_x))
    allocate(image_amp(dim_y, dim_x), image_mu(dim_y, dim_x), image_sig(dim_y, dim_x))
    
    dR_over_dB = 0._xp
    dF_over_dB = 0._xp
    deriv = 0._xp
    f = 0._xp
    g = 0._xp
    g_3D = 0._xp
    residual = 0._xp    
    params = 0._xp
    
    n_beta = (3*n_gauss * dim_y * dim_x) + n_gauss

    call unravel_3D(beta, params, 3*n_gauss, dim_y, dim_x)    
    do i=1,n_gauss
       b_params(i) = beta((n_beta-n_gauss)+i)
    end do

    ! Compute the objective function and the gradient
    do j=1, dim_x
       do i=1, dim_y
          allocate(residual_1D(dim_v))
          residual_1D = 0._xp
          call myresidual(params(:,i,j), cube(:,i,j), residual_1D, n_gauss, dim_v)
          residual(:,i,j) = residual_1D
          if (std_map(i,j) > 0._xp) then
             f = f + (myfunc_spec(residual_1D)/std_map(i,j)**2)
          end if
          deallocate(residual_1D)
       end do
    end do

    do l=1, dim_x
       do j=1, dim_y
          do i=1, n_gauss
             do k=1, dim_v          
                dF_over_dB(1+(3*(i-1)),k,j,l) = dF_over_dB(1+(3*(i-1)),k,j,l) +&
                     exp( ( -(real(k,xp) - params(2+(3*(i-1)),j,l))**2) / (2._xp * params(3+(3*(i-1)),j,l)**2))
                
                dF_over_dB(2+(3*(i-1)),k,j,l) = dF_over_dB(2+(3*(i-1)),k,j,l) +&
                     params(1+(3*(i-1)),j,l) * ( real(k,xp) - params(2+(3*(i-1)),j,l) ) / (params(3+(3*(i-1)),j,l)**2) *&
                     exp( ( -(real(k,xp) - params(2+(3*(i-1)),j,l))**2) / (2._xp * params(3+(3*(i-1)),j,l)**2))
                
                dF_over_dB(3+(3*(i-1)),k,j,l) = dF_over_dB(3+(3*(i-1)),k,j,l) +&
                     params(1+(3*(i-1)),j,l) * ( real(k,xp) - params(2+(3*(i-1)),j,l) )**2 / (params(3+(3*(i-1)),j,l)**3) *&
                     exp( ( -(real(k,xp) - params(2+(3*(i-1)),j,l))**2) / (2._xp * params(3+(3*(i-1)),j,l)**2))
             enddo
          enddo
       end do
    end do
    
    do k=1, dim_v
       do j=1, dim_x
          do i=1, dim_y
             do l=1, 3*n_gauss
                if (std_map(i,j) > 0._xp) then
                   deriv(l,i,j) = deriv(l,i,j) + dF_over_dB(l,k,i,j) * (residual(k,i,j)/std_map(i,j)**2)
                end if
             end do
          end do
       end do
    end do
    
    do k=1, n_gauss
       conv_amp = 0._xp; conv_mu = 0._xp; conv_sig = 0._xp
       conv_conv_amp = 0._xp; conv_conv_mu = 0._xp; conv_conv_sig = 0._xp
       image_amp = 0._xp; image_mu = 0._xp; image_sig = 0._xp

       image_amp = params(1+(3*(k-1)),:,:)
       image_mu = params(2+(3*(k-1)),:,:)
       image_sig = params(3+(3*(k-1)),:,:)
       
       call convolution_2D_mirror(image_amp, conv_amp, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(image_mu, conv_mu, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(image_sig, conv_sig, dim_y, dim_x, kernel, 3)

       call convolution_2D_mirror(conv_amp, conv_conv_amp, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(conv_mu, conv_conv_mu, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(conv_sig, conv_conv_sig, dim_y, dim_x, kernel, 3)

       do j=1, dim_x
          do i=1, dim_y
             f = f + (0.5_xp * lambda_amp * conv_amp(i,j)**2)
             f = f + (0.5_xp * lambda_mu * conv_mu(i,j)**2)
             f = f + (0.5_xp * lambda_sig * conv_sig(i,j)**2) + (0.5_xp * lambda_var_sig * (image_sig(i,j) - b_params(k))**2)
                          
             dR_over_dB(1+(3*(k-1)),i,j) = lambda_amp * conv_conv_amp(i,j)
             dR_over_dB(2+(3*(k-1)),i,j) = lambda_mu * conv_conv_mu(i,j)
             dR_over_dB(3+(3*(k-1)),i,j) = lambda_sig * conv_conv_sig(i,j) + (lambda_var_sig * (image_sig(i,j) - b_params(k)))

             g((n_beta-n_gauss)+k) = g((n_beta-n_gauss)+k) - (lambda_var_sig * (image_sig(i,j) - b_params(k)))        
             ! g((n_beta-n_gauss)+k) = - lambda_var_sig * (image_sig(i,j) - b_params(k))
          end do
       end do       
    end do
    
    g_3D = deriv + dR_over_dB
    call ravel_3D(g_3D, g, 3*n_gauss, dim_y, dim_x)

  end subroutine f_g_cube

end module mod_optimize
