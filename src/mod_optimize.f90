!! This module contains optimization subroutine and parametric model
module mod_optimize
  !! This module contains optimization subroutine and parametric model
  use mod_constants
  use mod_array

  implicit none
  
  private

  public :: gaussian, minimize_spec, minimize, myfunc_spec, myresidual, minimize_abs
  
contains
  
  pure function gaussian(x, a, m, s)
    !! Gaussian function   
    implicit none
    
    integer, intent(in) :: x
    real(xp), intent(in) :: a, m, s
    real(xp) :: gaussian

    gaussian = a * exp(-( (real(x,xp) - m)**2 ) / (2._xp * s**2) );
  end function gaussian


  subroutine minimize_spec(n, m, x, lb, ub, line, dim_v, n_gauss, maxiter, iprint)
    !! Minimize algorithn for a specturm
    implicit none      

    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: dim_v
    integer, intent(in) :: n_gauss, maxiter
    integer, intent(in) :: iprint

    real(xp), intent(in), dimension(:), allocatable :: lb, ub
    real(xp), intent(in), dimension(dim_v) :: line

    real(xp), intent(in), dimension(:), allocatable :: x
    
    real(xp), parameter    :: factr  = 1.0d+7, pgtol  = 1.0d-5
    
    character(len=60)      :: task, csave
    logical                :: lsave(4)
    integer                :: isave(44)
    real(xp)               :: f
    real(xp)               :: dsave(29)
    integer,  dimension(:), allocatable  :: nbd, iwa
    real(xp), dimension(:), allocatable  :: g, wa

    real(xp), dimension(:), allocatable  :: residual
    
    !     Allocate dynamic arrays
    allocate(nbd(n), g(n))
    allocate(iwa(3*n))
    allocate(wa(2*m*n + 5*n + 11*m*m + 8*m))

    allocate(residual(dim_v))
    residual = 0._xp

    ! Init nbd
    nbd = 2
    
    !     We now define the starting point.
    !     We start the iteration by initializing task.
    task = 'START'
    
    !     The beginning of the loop
    do while(task(1:2).eq.'FG'.or. task.eq.'NEW_X' .or. task.eq.'START') 
       
       !     This is the call to the L-BFGS-B code.
       call setulb (n, m, x, lb, ub, nbd, f, g, factr, pgtol, wa, iwa, task, iprint, csave, lsave, isave, dsave)
       
       if (task(1:2) .eq. 'FG') then          
          !     Compute function f and gradient g for the sample problem.
          call myresidual(x, line, residual, n_gauss, dim_v)
          f = myfunc_spec(residual)          
          call mygrad_spec(n_gauss, g, residual, x, dim_v)
          
       elseif (task(1:5) .eq. 'NEW_X') then
          !        1) Terminate if the total number of f and g evaluations
          !             exceeds maxiter.
          if (isave(34) .ge. maxiter) &
               task='STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT'
          
          !        2) Terminate if  |proj g|/(1+|f|) < 1.0d-10.            
          if (dsave(13) .le. 1.d-10*(1.0d0 + abs(f))) &
               task='STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL'
       endif
    !     end of loop do while       
    end do
  end subroutine minimize_spec

  
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
    
    myfunc_spec = 0.5_xp * sum(residual**2._xp)    
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
       do k=1, 3*n_gauss
          gradient(k) = gradient(k) + dF_over_dB(k,i) * residual(i)
       end do
    end do
  end subroutine mygrad_spec

  
  ! Minimize algorithn for a cube with regularization
  subroutine minimize(n, m, x, lb, ub, cube, n_gauss, dim_v, dim_y, dim_x, lambda_amp, lambda_mu, lambda_sig, &
       lambda_var_amp, lambda_var_mu, lambda_var_sig, maxiter, kernel, iprint, std_map, mean_amp, mean_mu, mean_sig)
    implicit none      

    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: dim_v, dim_y, dim_x
    integer, intent(in) :: n_gauss, maxiter
    integer, intent(in) :: iprint
    
    real(xp), intent(in) :: lambda_amp, lambda_mu, lambda_sig
    real(xp), intent(in) :: lambda_var_amp, lambda_var_mu, lambda_var_sig
    real(xp), intent(in), dimension(:), allocatable :: lb, ub
    real(xp), intent(in), dimension(:,:,:), allocatable :: cube
    real(xp), intent(in), dimension(:,:), allocatable :: kernel
    real(xp), intent(in), dimension(:,:), allocatable :: std_map
    real(xp), intent(in), dimension(:), allocatable :: mean_amp, mean_mu, mean_sig    

    real(xp), intent(in), dimension(:), allocatable :: x
    
    real(xp), parameter    :: factr  = 1.0d+7, pgtol  = 1.0d-5
    
    character(len=60)      :: task, csave
    logical                :: lsave(4)
    integer                :: isave(44)
    real(xp)               :: f
    real(xp)               :: dsave(29)
    integer,  dimension(:), allocatable  :: nbd, iwa
    real(xp), dimension(:), allocatable  :: g, wa

    real(xp), dimension(:,:,:), allocatable  :: residual
    
    !     Allocate dynamic arrays
    allocate(nbd(n), g(n))
    allocate(iwa(3*n))
    allocate(wa(2*m*n + 5*n + 11*m*m + 8*m))

    allocate(residual(dim_v, dim_y, dim_x))

    residual = 0._xp
    f = 0._xp
    g = 0._xp

    ! Init nbd
    nbd = 2
    
    !     We now define the starting point.
    !     We start the iteration by initializing task.
    task = 'START'
    
    !     The beginning of the loop
    do while(task(1:2).eq.'FG'.or. task.eq.'NEW_X' .or. task.eq.'START') 
       
       !     This is the call to the L-BFGS-B code.
       call setulb (n, m, x, lb, ub, nbd, f, g, factr, pgtol, wa, iwa, task, iprint, csave, lsave, isave, dsave)
       
       if (task(1:2) .eq. 'FG') then          
          !     Compute function f and gradient g for the sample problem.
          call f_g_cube(f, g, cube, x, dim_v, dim_y, dim_x, n_gauss, kernel, lambda_amp, lambda_mu, lambda_sig, &
               lambda_var_amp, lambda_var_mu, lambda_var_sig, std_map, mean_amp, mean_mu, mean_sig)
          
       elseif (task(1:5) .eq. 'NEW_X') then
          !        1) Terminate if the total number of f and g evaluations
          !             exceeds maxiter.
          if (isave(34) .ge. maxiter) &
               task='STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT'
          
          !        2) Terminate if  |proj g|/(1+|f|) < 1.0d-10.            
          if (dsave(13) .le. 1.d-10*(1.0d0 + abs(f))) &
               task='STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL'
       endif
    !     end of loop do while       
    end do
  end subroutine minimize


  ! Compute the objective function for a cube and the gradient of the obkective function
  subroutine f_g_cube(f, g, cube, beta, dim_v, dim_y, dim_x, n_gauss, kernel, lambda_amp, lambda_mu, lambda_sig, &
       lambda_var_amp, lambda_var_mu, lambda_var_sig, std_map, mean_amp, mean_mu, mean_sig)
    implicit none

    integer, intent(in) :: n_gauss
    integer, intent(in) :: dim_v, dim_y, dim_x
    real(xp), intent(in) :: lambda_amp, lambda_mu, lambda_sig
    real(xp), intent(in) :: lambda_var_amp, lambda_var_mu, lambda_var_sig
    real(xp), intent(in), dimension(:), allocatable :: beta
    real(xp), intent(in), dimension(:,:,:), allocatable :: cube
    real(xp), intent(in), dimension(:,:), allocatable :: kernel
    real(xp), intent(in), dimension(:,:), allocatable :: std_map
    real(xp), intent(in), dimension(:), allocatable :: mean_amp, mean_mu, mean_sig    
    real(xp), intent(inout) :: f
    real(xp), intent(inout), dimension(:), allocatable :: g

    integer :: i, j, k, l
    real(xp), dimension(:,:,:), allocatable :: residual
    real(xp), dimension(:), allocatable :: residual_1D
    real(xp), dimension(:,:,:), allocatable :: params
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
    
    call unravel_3D(beta, params, 3*n_gauss, dim_y, dim_x)    

    ! Compute the objective function and the gradient
    do j=1, dim_x
       do i=1, dim_y
          allocate(residual_1D(dim_v))
          residual_1D = 0._xp
          call myresidual(params(:,i,j), cube(:,i,j), residual_1D, n_gauss, dim_v)
          residual(:,i,j) = residual_1D
          if (std_map(i,j) > 0._xp) then
             f = f + (myfunc_spec(residual_1D)/std_map(i,j)**2._xp)
          end if
          deallocate(residual_1D)
       end do
    end do

    do l=1, dim_x
       do j=1, dim_y
          do i=1, n_gauss
             do k=1, dim_v          
                dF_over_dB(1+(3*(i-1)),k,j,l) = dF_over_dB(1+(3*(i-1)),k,j,l) +&
                     exp( ( -(real(k,xp) - params(2+(3*(i-1)),j,l))**2._xp) / (2._xp * params(3+(3*(i-1)),j,l)**2._xp))
                
                dF_over_dB(2+(3*(i-1)),k,j,l) = dF_over_dB(2+(3*(i-1)),k,j,l) +&
                     params(1+(3*(i-1)),j,l) * ( real(k,xp) - params(2+(3*(i-1)),j,l) ) / (params(3+(3*(i-1)),j,l)**2._xp) *&
                     exp( ( -(real(k,xp) - params(2+(3*(i-1)),j,l))**2._xp) / (2._xp * params(3+(3*(i-1)),j,l)**2._xp))
                
                dF_over_dB(3+(3*(i-1)),k,j,l) = dF_over_dB(3+(3*(i-1)),k,j,l) +&
                     params(1+(3*(i-1)),j,l) * ( real(k,xp) - params(2+(3*(i-1)),j,l) )**2._xp / (params(3+(3*(i-1)),j,l)**3._xp) *&
                     exp( ( -(real(k,xp) - params(2+(3*(i-1)),j,l))**2._xp) / (2._xp * params(3+(3*(i-1)),j,l)**2._xp))
             enddo
          enddo
       end do
    end do
    
    do k=1, dim_v
       do j=1, dim_x
          do i=1, dim_y
             do l=1, 3*n_gauss
                if (std_map(i,j) > 0._xp) then
                   deriv(l,i,j) = deriv(l,i,j) + dF_over_dB(l,k,i,j) * (residual(k,i,j)/std_map(i,j)**2._xp)
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

       !New term on sig
       do j=1, dim_x
          do i=1, dim_y
             f = f + (0.5_xp * lambda_amp * conv_amp(i,j)**2) + (0.5_xp * lambda_var_amp * (image_amp(i,j) - mean_amp(k))**2._xp)
             f = f + (0.5_xp * lambda_mu * conv_mu(i,j)**2) + (0.5_xp * lambda_var_mu * (image_mu(i,j) - mean_mu(k))**2._xp)
             f = f + (0.5_xp * lambda_sig * conv_sig(i,j)**2) + (0.5_xp * lambda_var_sig * (image_sig(i,j) - mean_sig(k))**2._xp)
                          
             dR_over_dB(1+(3*(k-1)),i,j) = lambda_amp * conv_conv_amp(i,j) + (lambda_var_amp * (image_amp(i,j) - mean_amp(k)))
             dR_over_dB(2+(3*(k-1)),i,j) = lambda_mu * conv_conv_mu(i,j) + (lambda_var_mu * (image_mu(i,j) - mean_mu(k)))
             dR_over_dB(3+(3*(k-1)),i,j) = lambda_sig * conv_conv_sig(i,j) + (lambda_var_sig * (image_sig(i,j) - mean_sig(k)))
          end do
       end do       
    end do
    
    g_3D = deriv + dR_over_dB

    call ravel_3D(g_3D, g, 3*n_gauss, dim_y, dim_x)
  end subroutine f_g_cube


  ! Minimize algorithn for a cube with regularization
  subroutine minimize_abs(n, m, x, lb, ub, cube, cube_abs, n_gauss, dim_v, dim_y, dim_x, lambda_amp, lambda_mu, lambda_sig, &
       lambda_abs_amp, lambda_abs_mu, lambda_abs_sig, lambda_abs_tot, lambda_abs_amp_joint, lambda_abs_mu_joint, &
       lambda_abs_sig_joint, lambda_var_amp, lambda_var_mu, lambda_var_sig, maxiter, kernel, iprint, std_map, std_map_abs, &
       mean_amp, mean_mu, mean_sig)
    implicit none      

    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: dim_v, dim_y, dim_x
    integer, intent(in) :: n_gauss, maxiter
    integer, intent(in) :: iprint
    
    real(xp), intent(in) :: lambda_amp, lambda_mu, lambda_sig
    real(xp), intent(in) :: lambda_abs_amp, lambda_abs_mu, lambda_abs_sig
    real(xp), intent(in) :: lambda_abs_tot, lambda_abs_amp_joint, lambda_abs_mu_joint, lambda_abs_sig_joint
    real(xp), intent(in) :: lambda_var_amp, lambda_var_mu, lambda_var_sig
    real(xp), intent(in), dimension(:), allocatable :: lb, ub
    real(xp), intent(in), dimension(:,:,:), allocatable :: cube, cube_abs
    real(xp), intent(in), dimension(:,:), allocatable :: kernel
    real(xp), intent(in), dimension(:,:), allocatable :: std_map, std_map_abs
    real(xp), intent(in), dimension(:), allocatable :: mean_amp, mean_mu, mean_sig    

    real(xp), intent(in), dimension(:), allocatable :: x
    
    real(xp), parameter    :: factr  = 1.0d+7, pgtol  = 1.0d-5
    
    character(len=60)      :: task, csave
    logical                :: lsave(4)
    integer                :: isave(44)
    real(xp)               :: f
    real(xp)               :: dsave(29)
    integer,  dimension(:), allocatable  :: nbd, iwa
    real(xp), dimension(:), allocatable  :: g, wa

    real(xp), dimension(:,:,:), allocatable  :: residual
    
    !     Allocate dynamic arrays
    allocate(nbd(n), g(n))
    allocate(iwa(3*n))
    allocate(wa(2*m*n + 5*n + 11*m*m + 8*m))

    allocate(residual(dim_v, dim_y, dim_x))

    residual = 0._xp
    f = 0._xp
    g = 0._xp

    ! Init nbd
    nbd = 2
    
    !     We now define the starting point.
    !     We start the iteration by initializing task.
    task = 'START'
    
    !     The beginning of the loop
    do while(task(1:2).eq.'FG'.or. task.eq.'NEW_X' .or. task.eq.'START') 
       
       !     This is the call to the L-BFGS-B code.
       call setulb (n, m, x, lb, ub, nbd, f, g, factr, pgtol, wa, iwa, task, iprint, csave, lsave, isave, dsave)
       
       if (task(1:2) .eq. 'FG') then          
          !     Compute function f and gradient g for the sample problem.
          call f_g_cube_abs(f, g, cube, cube_abs, x, dim_v, dim_y, dim_x, n_gauss, kernel, lambda_amp, lambda_mu, lambda_sig, &
               lambda_abs_amp, lambda_abs_mu, lambda_abs_sig, lambda_abs_tot, lambda_abs_amp_joint, lambda_abs_mu_joint, &
               lambda_abs_sig_joint, lambda_var_amp, lambda_var_mu, lambda_var_sig, std_map, std_map_abs, mean_amp, mean_mu, &
               mean_sig)
          
       elseif (task(1:5) .eq. 'NEW_X') then
          !        1) Terminate if the total number of f and g evaluations
          !             exceeds maxiter.
          if (isave(34) .ge. maxiter) &
               task='STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT'
          
          !        2) Terminate if  |proj g|/(1+|f|) < 1.0d-10.            
          if (dsave(13) .le. 1.d-10*(1.0d0 + abs(f))) &
               task='STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL'
       endif
    !     end of loop do while       
    end do
  end subroutine minimize_abs

  ! Compute the objective function for a cube and the gradient of the obkective function
  subroutine f_g_cube_abs(f, g, cube, cube_abs, beta, dim_v, dim_y, dim_x, n_gauss, kernel, lambda_amp, lambda_mu, lambda_sig, &
        lambda_abs_amp, lambda_abs_mu, lambda_abs_sig, lambda_abs_tot, lambda_abs_amp_joint, lambda_abs_mu_joint, &
        lambda_abs_sig_joint, lambda_var_amp, lambda_var_mu, lambda_var_sig, std_map, std_map_abs, mean_amp, mean_mu, mean_sig)
    implicit none

    integer, intent(in) :: n_gauss
    integer, intent(in) :: dim_v, dim_y, dim_x
    real(xp), intent(in) :: lambda_amp, lambda_mu, lambda_sig
    real(xp), intent(in) :: lambda_abs_amp, lambda_abs_mu, lambda_abs_sig
    real(xp), intent(in) :: lambda_abs_tot, lambda_abs_amp_joint, lambda_abs_mu_joint, lambda_abs_sig_joint
    real(xp), intent(in) :: lambda_var_amp, lambda_var_mu, lambda_var_sig
    real(xp), intent(in), dimension(:), allocatable :: beta
    real(xp), intent(in), dimension(:,:,:), allocatable :: cube, cube_abs
    real(xp), intent(in), dimension(:,:), allocatable :: kernel
    real(xp), intent(in), dimension(:,:), allocatable :: std_map, std_map_abs
    real(xp), intent(in), dimension(:), allocatable :: mean_amp, mean_mu, mean_sig    
    real(xp), intent(inout) :: f
    real(xp), intent(inout), dimension(:), allocatable :: g

    integer :: i, j, k, l
    real(xp), dimension(:,:,:), allocatable :: residual, residual_abs
    real(xp), dimension(:), allocatable :: residual_1D
    real(xp), dimension(:,:,:), allocatable :: params, params_abs
    real(xp), dimension(:,:), allocatable :: conv_amp, conv_mu, conv_sig
    real(xp), dimension(:,:), allocatable :: conv_abs_amp, conv_abs_mu, conv_abs_sig
    real(xp), dimension(:,:), allocatable :: conv_conv_amp, conv_conv_mu, conv_conv_sig
    real(xp), dimension(:,:), allocatable :: conv_conv_abs_amp, conv_conv_abs_mu, conv_conv_abs_sig
    real(xp), dimension(:,:), allocatable :: image_amp, image_mu, image_sig
    real(xp), dimension(:,:), allocatable :: image_amp_abs, image_mu_abs, image_sig_abs
    real(xp), dimension(:,:,:), allocatable :: g_3D, g_3D_abs
    real(xp), dimension(:,:,:,:), allocatable :: dF_over_dB, dF_over_dB_abs
    real(xp), dimension(:,:,:), allocatable :: dR_over_dB, dR_over_dB_abs
    real(xp), dimension(:,:,:), allocatable :: deriv, deriv_abs

    allocate(dR_over_dB(3*n_gauss, dim_y, dim_x))
    allocate(dF_over_dB(3*n_gauss, dim_v, dim_y, dim_x))
    allocate(dR_over_dB_abs(3*n_gauss, dim_y, dim_x))
    allocate(dF_over_dB_abs(3*n_gauss, dim_v, dim_y, dim_x))
    allocate(deriv(3*n_gauss, dim_y, dim_x))
    allocate(deriv_abs(3*n_gauss, dim_y, dim_x))
    allocate(g_3D(3*n_gauss, dim_y, dim_x))
    allocate(g_3D_abs(3*n_gauss, dim_y, dim_x))
    allocate(residual(dim_v, dim_y, dim_x))
    allocate(residual_abs(dim_v, dim_y, dim_x))
    allocate(params(3*n_gauss, dim_y, dim_x))
    allocate(params_abs(3*n_gauss, dim_y, dim_x))
    allocate(conv_amp(dim_y, dim_x), conv_mu(dim_y, dim_x), conv_sig(dim_y, dim_x))
    allocate(conv_abs_amp(dim_y, dim_x), conv_abs_mu(dim_y, dim_x), conv_abs_sig(dim_y, dim_x))
    allocate(conv_conv_amp(dim_y, dim_x), conv_conv_mu(dim_y, dim_x), conv_conv_sig(dim_y, dim_x))
    allocate(conv_conv_abs_amp(dim_y, dim_x), conv_conv_abs_mu(dim_y, dim_x), conv_conv_abs_sig(dim_y, dim_x))
    allocate(image_amp(dim_y, dim_x), image_mu(dim_y, dim_x), image_sig(dim_y, dim_x))
    allocate(image_amp_abs(dim_y, dim_x), image_mu_abs(dim_y, dim_x), image_sig_abs(dim_y, dim_x))
    
    dR_over_dB = 0._xp; dR_over_dB_abs = 0._xp
    dF_over_dB = 0._xp; dF_over_dB_abs = 0._xp
    deriv = 0._xp; deriv_abs = 0._xp
    f = 0._xp
    g = 0._xp
    g_3D = 0._xp; g_3D_abs = 0._xp
    residual = 0._xp    
    params = 0._xp
    params_abs = 0._xp
    
    call unravel_3D_abs(beta, params, params_abs, 3*n_gauss, dim_y, dim_x)    

    ! Compute the objective function and the gradient
    do j=1, dim_x
       do i=1, dim_y
          allocate(residual_1D(dim_v))
          residual_1D = 0._xp
          call myresidual(params(:,i,j), cube(:,i,j), residual_1D, n_gauss, dim_v)
          residual(:,i,j) = residual_1D
          if (std_map(i,j) > 0._xp) then
             f = f + (myfunc_spec(residual_1D)/std_map(i,j)**2._xp)
          end if
          residual_1D = 0._xp          
          call myresidual(params_abs(:,i,j), cube_abs(:,i,j), residual_1D, n_gauss, dim_v)
          residual_abs(:,i,j) = residual_1D
          if (std_map_abs(i,j) > 0._xp) then
             f = f + (lambda_abs_tot * (myfunc_spec(residual_1D)/std_map_abs(i,j)**2._xp))
          end if          
          deallocate(residual_1D)
       end do
    end do

    do l=1, dim_x
       do j=1, dim_y
          do i=1, n_gauss
             do k=1, dim_v          
                dF_over_dB(1+(3*(i-1)),k,j,l) = dF_over_dB(1+(3*(i-1)),k,j,l) +&
                     exp( ( -(real(k,xp) - params(2+(3*(i-1)),j,l))**2._xp) / (2._xp * params(3+(3*(i-1)),j,l)**2._xp))
                
                dF_over_dB(2+(3*(i-1)),k,j,l) = dF_over_dB(2+(3*(i-1)),k,j,l) +&
                     params(1+(3*(i-1)),j,l) * ( real(k,xp) - params(2+(3*(i-1)),j,l) ) / (params(3+(3*(i-1)),j,l)**2._xp) *&
                     exp( ( -(real(k,xp) - params(2+(3*(i-1)),j,l))**2._xp) / (2._xp * params(3+(3*(i-1)),j,l)**2._xp))
                
                dF_over_dB(3+(3*(i-1)),k,j,l) = dF_over_dB(3+(3*(i-1)),k,j,l) +&
                     params(1+(3*(i-1)),j,l) * ( real(k,xp) - params(2+(3*(i-1)),j,l) )**2._xp / (params(3+(3*(i-1)),j,l)**3._xp) *&
                     exp( ( -(real(k,xp) - params(2+(3*(i-1)),j,l))**2._xp) / (2._xp * params(3+(3*(i-1)),j,l)**2._xp))

                dF_over_dB_abs(1+(3*(i-1)),k,j,l) = dF_over_dB_abs(1+(3*(i-1)),k,j,l) +&
                     exp( ( -(real(k,xp) - params_abs(2+(3*(i-1)),j,l))**2._xp) / (2._xp * params_abs(3+(3*(i-1)),j,l)**2._xp))
                
                dF_over_dB_abs(2+(3*(i-1)),k,j,l) = dF_over_dB_abs(2+(3*(i-1)),k,j,l) +&
                     params_abs(1+(3*(i-1)),j,l) * ( real(k,xp) - params_abs(2+(3*(i-1)),j,l) ) / &
                     (params_abs(3+(3*(i-1)),j,l)**2._xp) * exp( ( -(real(k,xp) - params_abs(2+(3*(i-1)),j,l))**2._xp) &
                     / (2._xp * params_abs(3+(3*(i-1)),j,l)**2._xp))
                
                dF_over_dB_abs(3+(3*(i-1)),k,j,l) = dF_over_dB_abs(3+(3*(i-1)),k,j,l) +&
                     params_abs(1+(3*(i-1)),j,l) * ( real(k,xp) - params_abs(2+(3*(i-1)),j,l) )**2._xp &
                     / (params_abs(3+(3*(i-1)),j,l)**3._xp) * exp( ( -(real(k,xp) - params_abs(2+(3*(i-1)),j,l))**2._xp) &
                     / (2._xp * params_abs(3+(3*(i-1)),j,l)**2._xp))

             enddo
          enddo
       end do
    end do
    
    do k=1, dim_v
       do j=1, dim_x
          do i=1, dim_y
             do l=1, 3*n_gauss
                if (std_map(i,j) > 0._xp) then
                   deriv(l,i,j) = deriv(l,i,j) + dF_over_dB(l,k,i,j) * (residual(k,i,j)/std_map(i,j)**2._xp)
                end if
                if (std_map_abs(i,j) > 0._xp) then
                   deriv_abs(l,i,j) = deriv_abs(l,i,j) + (lambda_abs_tot * dF_over_dB_abs(l,k,i,j) * &
                        (residual_abs(k,i,j)/std_map_abs(i,j)**2._xp))
                end if
             end do
          end do
       end do
    end do

    do k=1, n_gauss
       conv_amp = 0._xp; conv_mu = 0._xp; conv_sig = 0._xp
       conv_conv_amp = 0._xp; conv_conv_mu = 0._xp; conv_conv_sig = 0._xp
       conv_abs_amp = 0._xp; conv_abs_mu = 0._xp; conv_abs_sig = 0._xp
       conv_conv_abs_amp = 0._xp; conv_conv_abs_mu = 0._xp; conv_conv_abs_sig = 0._xp
       image_amp = 0._xp; image_mu = 0._xp; image_sig = 0._xp
       image_amp_abs = 0._xp; image_mu_abs = 0._xp; image_sig_abs = 0._xp

       image_amp = params(1+(3*(k-1)),:,:)
       image_mu = params(2+(3*(k-1)),:,:)
       image_sig = params(3+(3*(k-1)),:,:)

       image_amp_abs = params_abs(1+(3*(k-1)),:,:)
       image_mu_abs = params_abs(2+(3*(k-1)),:,:)
       image_sig_abs = params_abs(3+(3*(k-1)),:,:)
       
       call convolution_2D_mirror(image_amp, conv_amp, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(image_mu, conv_mu, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(image_sig, conv_sig, dim_y, dim_x, kernel, 3)

       call convolution_2D_mirror(conv_amp, conv_conv_amp, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(conv_mu, conv_conv_mu, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(conv_sig, conv_conv_sig, dim_y, dim_x, kernel, 3)

       call convolution_2D_mirror(image_amp_abs, conv_abs_amp, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(image_mu_abs, conv_abs_mu, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(image_sig_abs, conv_abs_sig, dim_y, dim_x, kernel, 3)

       call convolution_2D_mirror(conv_abs_amp, conv_conv_abs_amp, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(conv_abs_mu, conv_conv_abs_mu, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(conv_abs_sig, conv_conv_abs_sig, dim_y, dim_x, kernel, 3)

       !New term on sig
       do j=1, dim_x
          do i=1, dim_y
             f = f + (0.5_xp * lambda_amp * conv_amp(i,j)**2) + (0.5_xp * lambda_var_amp * (image_amp(i,j) - mean_amp(k))**2._xp)
             f = f + (0.5_xp * lambda_mu * conv_mu(i,j)**2) + (0.5_xp * lambda_var_mu * (image_mu(i,j) - mean_mu(k))**2._xp) 
             f = f + (0.5_xp * lambda_sig * conv_sig(i,j)**2) + (0.5_xp * lambda_var_sig * (image_sig(i,j) - mean_sig(k))**2._xp)

             f = f + (0.5_xp * lambda_abs_amp * conv_abs_amp(i,j)**2)
             f = f + (0.5_xp * lambda_abs_mu * conv_abs_mu(i,j)**2)
             f = f + (0.5_xp * lambda_abs_sig * conv_abs_sig(i,j)**2)

             f = f + (0.5_xp * lambda_abs_mu_joint * ((image_mu_abs(i,j) / image_mu(i,j)) - 1._xp)**2._xp)
             f = f + (0.5_xp * lambda_abs_sig_joint * ((image_sig_abs(i,j) / image_sig(i,j)) - 1._xp)**2._xp)
             
                          
             dR_over_dB(1+(3*(k-1)),i,j) = lambda_amp * conv_conv_amp(i,j) + (lambda_var_amp * (image_amp(i,j) - mean_amp(k))) 
             dR_over_dB(2+(3*(k-1)),i,j) = lambda_mu * conv_conv_mu(i,j) + (lambda_var_mu * (image_mu(i,j) - mean_mu(k)))
             dR_over_dB(3+(3*(k-1)),i,j) = lambda_sig * conv_conv_sig(i,j) + (lambda_var_sig * (image_sig(i,j) - mean_sig(k)))

             dR_over_dB_abs(1+(3*(k-1)),i,j) = lambda_abs_amp * conv_conv_abs_amp(i,j)
             dR_over_dB_abs(2+(3*(k-1)),i,j) = lambda_abs_mu * conv_conv_abs_mu(i,j)
             dR_over_dB_abs(3+(3*(k-1)),i,j) = lambda_abs_sig * conv_conv_abs_sig(i,j)

             !Gradient correlation velocity fields
             dR_over_dB(2+(3*(k-1)),i,j) = dR_over_dB(2+(3*(k-1)),i,j) - (lambda_abs_mu_joint * &
                  image_mu_abs(i,j) / image_mu(i,j)**2._xp & 
                  * (image_mu_abs(i,j) / image_mu(i,j) - 1._xp))

             dR_over_dB_abs(2+(3*(k-1)),i,j) = dR_over_dB_abs(2+(3*(k-1)),i,j) + (lambda_abs_mu_joint / image_mu(i,j)&
                  * (image_mu_abs(i,j) / image_mu(i,j) - 1._xp))

             !Gradient correlation velocity dispersion fields
             dR_over_dB(3+(3*(k-1)),i,j) = dR_over_dB(3+(3*(k-1)),i,j) - (lambda_abs_sig_joint * &
                  image_sig_abs(i,j) / image_sig(i,j)**2._xp & 
                  * (image_sig_abs(i,j) / image_sig(i,j) - 1._xp))

             dR_over_dB_abs(3+(3*(k-1)),i,j) = dR_over_dB_abs(3+(3*(k-1)),i,j) + (lambda_abs_sig_joint / image_sig(i,j)&
                  * (image_sig_abs(i,j) / image_sig(i,j) - 1._xp))

          end do
       end do       
    end do
    
    g_3D = deriv + dR_over_dB
    g_3D_abs = deriv_abs + dR_over_dB_abs

    call ravel_3D_abs(g_3D, g_3D_abs, g, 3*n_gauss, dim_y, dim_x)
  end subroutine f_g_cube_abs


  ! ! Minimize algorithn for a cube with regularization
  ! subroutine minimize_new_abs(n, m, x, lb, ub, params_mirror, cube, n_gauss, dim_v, dim_y, dim_x, lambda_amp, lambda_mu, &
  !      lambda_sig, lambda_mu_mirror, lambda_var_amp, lambda_var_mu, lambda_var_sig, maxiter, kernel, iprint, std_map, mean_amp, &
  !      mean_mu, mean_sig)
  !   implicit none      

  !   integer, intent(in) :: n
  !   integer, intent(in) :: m
  !   integer, intent(in) :: dim_v, dim_y, dim_x
  !   integer, intent(in) :: n_gauss, maxiter
  !   integer, intent(in) :: iprint
    
  !   real(xp), intent(in) :: lambda_amp, lambda_mu, lambda_sig
  !   real(xp), intent(in) :: lambda_mu_mirror    
  !   real(xp), intent(in) :: lambda_var_amp, lambda_var_mu, lambda_var_sig
  !   real(xp), intent(in), dimension(:), allocatable :: lb, ub
  !   real(xp), intent(in), dimension(:,:,:), allocatable :: cube
  !   real(xp), intent(in), dimension(:,:), allocatable :: kernel
  !   real(xp), intent(in), dimension(:,:), allocatable :: std_map
  !   real(xp), intent(in), dimension(:), allocatable :: mean_amp, mean_mu, mean_sig    
  !   real(xp), intent(in), dimension(:,:,:), allocatable :: params_mirror

  !   real(xp), intent(in), dimension(:), allocatable :: x
    
  !   real(xp), parameter    :: factr  = 1.0d+7, pgtol  = 1.0d-5
    
  !   character(len=60)      :: task, csave
  !   logical                :: lsave(4)
  !   integer                :: isave(44)
  !   real(xp)               :: f
  !   real(xp)               :: dsave(29)
  !   integer,  dimension(:), allocatable  :: nbd, iwa
  !   real(xp), dimension(:), allocatable  :: g, wa

  !   real(xp), dimension(:,:,:), allocatable  :: residual
    
  !   !     Allocate dynamic arrays
  !   allocate(nbd(n), g(n))
  !   allocate(iwa(3*n))
  !   allocate(wa(2*m*n + 5*n + 11*m*m + 8*m))

  !   allocate(residual(dim_v, dim_y, dim_x))

  !   residual = 0._xp
  !   f = 0._xp
  !   g = 0._xp

  !   ! Init nbd
  !   nbd = 2
    
  !   !     We now define the starting point.
  !   !     We start the iteration by initializing task.
  !   task = 'START'
    
  !   !     The beginning of the loop
  !   do while(task(1:2).eq.'FG'.or. task.eq.'NEW_X' .or. task.eq.'START') 
       
  !      !     This is the call to the L-BFGS-B code.
  !      call setulb (n, m, x, lb, ub, nbd, f, g, factr, pgtol, wa, iwa, task, iprint, csave, lsave, isave, dsave)
       
  !      if (task(1:2) .eq. 'FG') then          
  !         !     Compute function f and gradient g for the sample problem.
  !         call f_g_cube_new_abs(f, g, cube, x, dim_v, dim_y, dim_x, n_gauss, kernel, lambda_amp, lambda_mu, lambda_sig, &
  !              lambda_mu_mirror, lambda_var_amp, lambda_var_mu, lambda_var_sig, std_map, mean_amp, mean_mu, mean_sig, &
  !              params_mirror)
          
  !      elseif (task(1:5) .eq. 'NEW_X') then
  !         !        1) Terminate if the total number of f and g evaluations
  !         !             exceeds maxiter.
  !         if (isave(34) .ge. maxiter) &
  !              task='STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT'
          
  !         !        2) Terminate if  |proj g|/(1+|f|) < 1.0d-10.            
  !         if (dsave(13) .le. 1.d-10*(1.0d0 + abs(f))) &
  !              task='STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL'
  !      endif
  !   !     end of loop do while       
  !   end do
  ! end subroutine minimize_new_abs


  ! ! Compute the objective function for a cube and the gradient of the obkective function
  ! subroutine f_g_cube_new_abs(f, g, cube, beta, dim_v, dim_y, dim_x, n_gauss, kernel, lambda_amp, lambda_mu, lambda_sig, &
  !      lambda_mu_mirror, lambda_var_amp, lambda_var_mu, lambda_var_sig, std_map, mean_amp, mean_mu, mean_sig, params_mirror)
  !   implicit none

  !   integer, intent(in) :: n_gauss
  !   integer, intent(in) :: dim_v, dim_y, dim_x
  !   real(xp), intent(in) :: lambda_amp, lambda_mu, lambda_sig
  !   real(xp), intent(in) :: lambda_mu_mirror
  !   real(xp), intent(in) :: lambda_var_amp, lambda_var_mu, lambda_var_sig
  !   real(xp), intent(in), dimension(:), allocatable :: beta
  !   real(xp), intent(in), dimension(:,:,:), allocatable :: cube
  !   real(xp), intent(in), dimension(:,:), allocatable :: kernel
  !   real(xp), intent(in), dimension(:,:), allocatable :: std_map
  !   real(xp), intent(in), dimension(:), allocatable :: mean_amp, mean_mu, mean_sig    
  !   real(xp), intent(in), dimension(:,:,:), allocatable :: params_mirror
  !   real(xp), intent(inout) :: f
  !   real(xp), intent(inout), dimension(:), allocatable :: g

  !   integer :: i, j, k, l
  !   real(xp), dimension(:,:,:), allocatable :: residual
  !   real(xp), dimension(:), allocatable :: residual_1D
  !   real(xp), dimension(:,:,:), allocatable :: params
  !   real(xp), dimension(:,:), allocatable :: conv_amp, conv_mu, conv_sig
  !   real(xp), dimension(:,:), allocatable :: conv_conv_amp, conv_conv_mu, conv_conv_sig
  !   real(xp), dimension(:,:), allocatable :: image_amp, image_mu, image_sig
  !   real(xp), dimension(:,:), allocatable :: image_mu_mirror
  !   real(xp), dimension(:,:,:), allocatable :: g_3D
  !   real(xp), dimension(:,:,:,:), allocatable :: dF_over_dB
  !   real(xp), dimension(:,:,:), allocatable :: dR_over_dB
  !   real(xp), dimension(:,:,:), allocatable :: deriv

  !   allocate(dR_over_dB(3*n_gauss, dim_y, dim_x))
  !   allocate(dF_over_dB(3*n_gauss, dim_v, dim_y, dim_x))
  !   allocate(deriv(3*n_gauss, dim_y, dim_x))
  !   allocate(g_3D(3*n_gauss, dim_y, dim_x))
  !   allocate(residual(dim_v, dim_y, dim_x))
  !   allocate(params(3*n_gauss, dim_y, dim_x))
  !   allocate(conv_amp(dim_y, dim_x), conv_mu(dim_y, dim_x), conv_sig(dim_y, dim_x))
  !   allocate(conv_conv_amp(dim_y, dim_x), conv_conv_mu(dim_y, dim_x), conv_conv_sig(dim_y, dim_x))
  !   allocate(image_amp(dim_y, dim_x), image_mu(dim_y, dim_x), image_sig(dim_y, dim_x))
    
  !   dR_over_dB = 0._xp
  !   dF_over_dB = 0._xp
  !   deriv = 0._xp
  !   f = 0._xp
  !   g = 0._xp
  !   g_3D = 0._xp
  !   residual = 0._xp    
  !   params = 0._xp
    
  !   call unravel_3D(beta, params, 3*n_gauss, dim_y, dim_x)    

  !   ! Compute the objective function and the gradient
  !   do j=1, dim_x
  !      do i=1, dim_y
  !         allocate(residual_1D(dim_v))
  !         residual_1D = 0._xp
  !         call myresidual(params(:,i,j), cube(:,i,j), residual_1D, n_gauss, dim_v)
  !         residual(:,i,j) = residual_1D
  !         if (std_map(i,j) > 0._xp) then
  !            f = f + (myfunc_spec(residual_1D)/std_map(i,j)**2._xp)
  !         end if
  !         deallocate(residual_1D)
  !      end do
  !   end do

  !   do l=1, dim_x
  !      do j=1, dim_y
  !         do i=1, n_gauss
  !            do k=1, dim_v          
  !               dF_over_dB(1+(3*(i-1)),k,j,l) = dF_over_dB(1+(3*(i-1)),k,j,l) +&
  !                    exp( ( -(real(k,xp) - params(2+(3*(i-1)),j,l))**2._xp) / (2._xp * params(3+(3*(i-1)),j,l)**2._xp))
                
  !               dF_over_dB(2+(3*(i-1)),k,j,l) = dF_over_dB(2+(3*(i-1)),k,j,l) +&
  !                    params(1+(3*(i-1)),j,l) * ( real(k,xp) - params(2+(3*(i-1)),j,l) ) / (params(3+(3*(i-1)),j,l)**2._xp) *&
  !                    exp( ( -(real(k,xp) - params(2+(3*(i-1)),j,l))**2._xp) / (2._xp * params(3+(3*(i-1)),j,l)**2._xp))
                
  !               dF_over_dB(3+(3*(i-1)),k,j,l) = dF_over_dB(3+(3*(i-1)),k,j,l) +&
  !                    params(1+(3*(i-1)),j,l) * ( real(k,xp) - params(2+(3*(i-1)),j,l) )**2._xp / (params(3+(3*(i-1)),j,l)**3._xp) *&
  !                    exp( ( -(real(k,xp) - params(2+(3*(i-1)),j,l))**2._xp) / (2._xp * params(3+(3*(i-1)),j,l)**2._xp))
  !            enddo
  !         enddo
  !      end do
  !   end do
    
  !   do k=1, dim_v
  !      do j=1, dim_x
  !         do i=1, dim_y
  !            do l=1, 3*n_gauss
  !               if (std_map(i,j) > 0._xp) then
  !                  deriv(l,i,j) = deriv(l,i,j) + dF_over_dB(l,k,i,j) * (residual(k,i,j)/std_map(i,j)**2._xp)
  !               end if
  !            end do
  !         end do
  !      end do
  !   end do
    
  !   do k=1, n_gauss
  !      conv_amp = 0._xp; conv_mu = 0._xp; conv_sig = 0._xp
  !      conv_conv_amp = 0._xp; conv_conv_mu = 0._xp; conv_conv_sig = 0._xp
  !      image_amp = 0._xp; image_mu = 0._xp; image_sig = 0._xp

  !      image_amp = params(1+(3*(k-1)),:,:)
  !      image_mu = params(2+(3*(k-1)),:,:)
  !      image_sig = params(3+(3*(k-1)),:,:)

  !      image_mu_mirror = params_mirror(2+(3*(k-1)),:,:)
       
  !      call convolution_2D_mirror(image_amp, conv_amp, dim_y, dim_x, kernel, 3)
  !      call convolution_2D_mirror(image_mu, conv_mu, dim_y, dim_x, kernel, 3)
  !      call convolution_2D_mirror(image_sig, conv_sig, dim_y, dim_x, kernel, 3)

  !      call convolution_2D_mirror(conv_amp, conv_conv_amp, dim_y, dim_x, kernel, 3)
  !      call convolution_2D_mirror(conv_mu, conv_conv_mu, dim_y, dim_x, kernel, 3)
  !      call convolution_2D_mirror(conv_sig, conv_conv_sig, dim_y, dim_x, kernel, 3)

  !      !New term on sig
  !      do j=1, dim_x
  !         do i=1, dim_y
  !            f = f + (0.5_xp * lambda_amp * conv_amp(i,j)**2) + (0.5_xp * lambda_var_amp * (image_amp(i,j) - mean_amp(k))**2._xp)
  !            f = f + (0.5_xp * lambda_mu * conv_mu(i,j)**2) + (0.5_xp * lambda_var_mu * (image_mu(i,j) - mean_mu(k))**2._xp) &
  !                 + (0.5_xp * lambda_mu_mirror * (image_mu(i,j) - image_mu_mirror(i,j))**2._xp)
  !            f = f + (0.5_xp * lambda_sig * conv_sig(i,j)**2) + (0.5_xp * lambda_var_sig * (image_sig(i,j) - mean_sig(k))**2._xp)
                          
  !            dR_over_dB(1+(3*(k-1)),i,j) = lambda_amp * conv_conv_amp(i,j) + (lambda_var_amp * (image_amp(i,j) - mean_amp(k)))
  !            dR_over_dB(2+(3*(k-1)),i,j) = lambda_mu * conv_conv_mu(i,j) + (lambda_var_mu * (image_mu(i,j) - mean_mu(k))) & 
  !                 + (lambda_mu_mirror * (image_mu(i,j) - image_mu_mirror(i,j)))
  !            dR_over_dB(3+(3*(k-1)),i,j) = lambda_sig * conv_conv_sig(i,j) + (lambda_var_sig * (image_sig(i,j) - mean_sig(k)))
  !         end do
  !      end do       
  !   end do
    
  !   g_3D = deriv + dR_over_dB

  !   call ravel_3D(g_3D, g, 3*n_gauss, dim_y, dim_x)
  ! end subroutine f_g_cube_new_abs

end module mod_optimize
