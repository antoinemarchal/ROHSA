!! This module contains optimization subroutine and parametric model
module mod_optimize
  !! This module contains optimization subroutine and parametric model
  use mod_constants
  use mod_array
  use mod_model
  use mod_color
  use mod_read_parameters
  use mod_fft

  implicit none
  
  private

  public :: myfunc_spec, myresidual, mygrad_spec, f_g_cube_fast
  
contains
  
  ! Compute the residual between model and data
  subroutine myresidual(pars, line, residual, dim_v, wavelength, color, std)
    implicit none

    integer, intent(in) :: dim_v
    real(xp), intent(in), dimension(dim_v) :: line
    real(xp), intent(in), dimension(dim_v) :: wavelength
    real(xp), intent(in), dimension(3*params%n_mbb) :: pars
    real(xp), intent(in), dimension(:,:) :: color
    real(xp), intent(in), dimension(:) :: std
    real(xp), intent(inout), dimension(:), allocatable :: residual

    integer :: i, k
    real(xp) :: g    
    real(xp), dimension(dim_v) :: model

    g = 0._xp
    model = 0._xp

    do i=1, params%n_mbb
       do k=1, dim_v
          if (params%cc .eqv. .true.) then
             g = mbb_l(wavelength(k), pars(1+(3*(i-1))), pars(2+(3*(i-1))), pars(3+(3*(i-1))), params%l0) &
                  * poly_color(color(k,:), pars(2+(3*(i-1))), pars(3+(3*(i-1))), params%degree) 
          else
             g = mbb_l(wavelength(k), pars(1+(3*(i-1))), pars(2+(3*(i-1))), pars(3+(3*(i-1))), params%l0)
          end if
          model(k) = model(k) + g
       enddo
    enddo

    residual = (model - line)/std
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
  subroutine mygrad_spec(gradient, residual, pars, dim_v, wavelength, color, std)
    implicit none

    integer, intent(in) :: dim_v
    real(xp), intent(in), dimension(3*params%n_mbb) :: pars
    real(xp), intent(in), dimension(:), allocatable :: residual
    real(xp), intent(in), dimension(dim_v) :: wavelength
    real(xp), intent(in), dimension(:,:) :: color
    real(xp), intent(inout), dimension(3*params%n_mbb) :: gradient
    real(xp), intent(in), dimension(:) :: std

    integer :: i, k
    real(xp) :: g

    real(xp), dimension(:,:), allocatable :: dF_over_dB

    allocate(dF_over_dB(3*params%n_mbb, dim_v))

    g = 0._xp
    dF_over_dB = 0._xp
    gradient = 0._xp

    if (params%cc .eqv. .true.) then
       do i=1, params%n_mbb
          do k=1, dim_v          
             dF_over_dB(1+(3*(i-1)),k) = dF_over_dB(1+(3*(i-1)),k) + ( &
                  d_mbbcc_l_dtau(wavelength(k), pars(2+(3*(i-1))), pars(3+(3*(i-1))), params%l0, &
                  color(k,:), params%degree) &
                  )
             
             dF_over_dB(2+(3*(i-1)),k) = dF_over_dB(2+(3*(i-1)),k) + ( &
                  d_mbbcc_l_db(wavelength(k), pars(1+(3*(i-1))), pars(2+(3*(i-1))), pars(3+(3*(i-1))), &
                  params%l0, color(k,:), params%degree) &               
                  )
             
             dF_over_dB(3+(3*(i-1)),k) = dF_over_dB(3+(3*(i-1)),k) + ( &
                  d_mbbcc_l_dT(wavelength(k), pars(1+(3*(i-1))), pars(2+(3*(i-1))), pars(3+(3*(i-1))), &
                  params%l0, color(k,:), params%degree) &                              
                  )
          enddo
       enddo
    else       
       do i=1, params%n_mbb
          do k=1, dim_v          
             dF_over_dB(1+(3*(i-1)),k) = dF_over_dB(1+(3*(i-1)),k) + ( &
                  d_mbb_l_dtau(wavelength(k), pars(2+(3*(i-1))), pars(3+(3*(i-1))), params%l0) &               
                  )
             
             dF_over_dB(2+(3*(i-1)),k) = dF_over_dB(2+(3*(i-1)),k) + ( &
                  d_mbb_l_db(wavelength(k), pars(1+(3*(i-1))), pars(2+(3*(i-1))), pars(3+(3*(i-1))), &
                  params%l0) &               
                  )
             
             dF_over_dB(3+(3*(i-1)),k) = dF_over_dB(3+(3*(i-1)),k) + ( &
                  d_mbb_l_dT(wavelength(k), pars(1+(3*(i-1))), pars(2+(3*(i-1))), pars(3+(3*(i-1))), &
                  params%l0) &                              
                  )
          enddo
       enddo
    end if
    
    do i=1, dim_v
       do k=1, 3*params%n_mbb
          gradient(k) = gradient(k) + dF_over_dB(k,i) * residual(i)/std(i)
       end do
    end do

    deallocate(dF_over_dB)
  end subroutine mygrad_spec

  
  ! Compute the objective function for a cube and the gradient of the obkective function
  subroutine f_g_cube_fast(f, g, cube, cube_HI, beta, dim_v, dim_y, dim_x, kernel, std_cube, &
       wavelength, color, filter, tapper)
    implicit none

    integer, intent(in) :: dim_v, dim_y, dim_x
    real(xp), intent(in), dimension(:), allocatable :: beta
    real(xp), intent(in), dimension(:,:,:), allocatable :: cube
    real(xp), intent(in), dimension(:,:,:), allocatable :: cube_HI
    real(xp), intent(in), dimension(:), allocatable :: wavelength
    real(xp), intent(in), dimension(:,:), allocatable :: color
    real(xp), intent(in), dimension(:,:), allocatable :: kernel
    real(xp), intent(in), dimension(:,:,:), allocatable :: std_cube

    real(xp), intent(in), dimension(:,:), allocatable :: filter
    real(xp), intent(in), dimension(:,:), allocatable :: tapper

    real(xp), intent(inout) :: f
    real(xp), intent(inout), dimension(:), allocatable :: g

    integer :: i, j, k, l
    integer :: n_beta
    integer :: n_cube
    real(xp), dimension(:,:,:), allocatable :: residual
    real(xp), dimension(:), allocatable :: residual_1D
    real(xp), dimension(:,:,:), allocatable :: pars
    real(xp), dimension(:), allocatable :: b_pars
    real(xp), dimension(:), allocatable :: c_pars
    real(xp), dimension(:), allocatable :: d_pars
    real(xp), dimension(:), allocatable :: stefan_pars
    real(xp), dimension(:,:), allocatable :: conv_tau, conv_beta, conv_Td
    real(xp), dimension(:,:), allocatable :: conv_conv_tau, conv_conv_beta, conv_conv_Td
    real(xp), dimension(:,:), allocatable :: image_tau, image_beta, image_Td
    real(xp), dimension(:,:,:), allocatable :: deriv
    real(xp), dimension(:), allocatable :: model
    real(xp) :: gauss

    real(xp), dimension(:,:), allocatable :: tau_ciba
    complex(xp), dimension(:,:), allocatable :: c_tau_ciba
    complex(xp), dimension(:,:), allocatable :: tf_tau_ciba    
    complex(xp), dimension(:,:), allocatable :: filtered_tf
    complex(xp), dimension(:,:), allocatable :: filtered_itf
    complex(xp), dimension(:,:), allocatable :: d_filtered_tf
    complex(xp), dimension(:,:), allocatable :: d_filtered_itf

    !Cross correlation
    real(xp), dimension(:,:), allocatable :: corr, ra, rb, corr_grad

    allocate(deriv(3*params%n_mbb, dim_y, dim_x))
    allocate(residual(dim_v, dim_y, dim_x))
    allocate(b_pars(params%n_mbb))
    allocate(c_pars(params%n_mbb))
    allocate(d_pars(params%n_mbb))
    allocate(stefan_pars(params%n_mbb))
    allocate(pars(3*params%n_mbb, dim_y, dim_x))
    allocate(conv_tau(dim_y, dim_x), conv_beta(dim_y, dim_x), conv_Td(dim_y, dim_x))
    allocate(conv_conv_tau(dim_y, dim_x), conv_conv_beta(dim_y, dim_x), conv_conv_Td(dim_y, dim_x))
    allocate(image_tau(dim_y, dim_x), image_beta(dim_y, dim_x), image_Td(dim_y, dim_x))
    allocate(model(dim_v))    
    allocate(tau_ciba(dim_y,dim_x), c_tau_ciba(dim_y,dim_x), tf_tau_ciba(dim_y,dim_x))
    allocate(filtered_tf(dim_y,dim_x), filtered_itf(dim_y,dim_x))
    allocate(d_filtered_tf(dim_y,dim_x), d_filtered_itf(dim_y,dim_x))
    allocate(corr(dim_y,dim_x), corr_grad(dim_y,dim_x), ra(dim_y,dim_x), rb(dim_y,dim_x))
    
    deriv = 0._xp
    f = 0._xp
    g = 0._xp
    residual = 0._xp    
    pars = 0._xp
    model = 0._xp
    gauss = 0._xp
    
    n_beta = (3*params%n_mbb * dim_y * dim_x) + (4*params%n_mbb)
    n_cube = (3*params%n_mbb * dim_y * dim_x)

    call unravel_3D(beta, pars, 3*params%n_mbb, dim_y, dim_x)    
    do i=1,params%n_mbb
       b_pars(i) = beta(n_cube+(0*params%n_mbb)+i)
       stefan_pars(i) = beta(n_cube+(1*params%n_mbb)+i)
       c_pars(i) = beta(n_cube+(2*params%n_mbb)+i)
       d_pars(i) = beta(n_cube+(3*params%n_mbb)+i)
    end do

    ! print*, c_pars
    ! print*, d_pars
    ! print*, sum(pars(1,:,:))/size(pars(1,:,:))
    
    ! Compute the objective function and the gradient
    do j=1, dim_x
       do i=1, dim_y
          allocate(residual_1D(dim_v))
          residual_1D = 0._xp
          call myresidual(pars(:,i,j), cube(:,i,j), residual_1D, dim_v, wavelength, &
               color, std_cube(:,i,j))
          residual(:,i,j) = residual_1D
          !Attache aux donnees
          f = f + (myfunc_spec(residual_1D))
          deallocate(residual_1D)
       end do
    end do

    ! Compute the objective function and the gradient
    do i=1, params%n_mbb
       !
       conv_tau = 0._xp; conv_beta = 0._xp; conv_Td = 0._xp
       conv_conv_tau = 0._xp; conv_conv_beta = 0._xp; conv_conv_Td = 0._xp
       image_tau = 0._xp; image_beta = 0._xp; image_Td = 0._xp
       tau_ciba = 0._xp
       
       image_tau = pars(1+(3*(i-1)),:,:)
       image_beta = pars(2+(3*(i-1)),:,:)
       image_Td = pars(3+(3*(i-1)),:,:)
       
       !CIBA FIXME
       ! if (params%ciba .eqv. .true.) then 
       !    if (dim_y .gt. 4 .and. i .eq. 1) then
       !       !Shift real image
       !       call shift(image_tau, tau_ciba)
       !       !Prepare complex array
       !       ! c_tau_ciba = cmplx(tapper*tau_ciba,0._xp,xp)
       !       c_tau_ciba = cmplx(tau_ciba,0._xp,xp)
       !       !Compute centered FFT with unitary transform using fftpack
       !       call cfft2d(dim_y,dim_x,c_tau_ciba,tf_tau_ciba)
       !       filtered_tf = filter * tf_tau_ciba
       !       d_filtered_tf = filter**2._xp * tf_tau_ciba
       !       call icfft2d(dim_y,dim_x,filtered_tf,filtered_itf)
       !       call icfft2d(dim_y,dim_x,d_filtered_tf,d_filtered_itf)
       !    end if
       ! end if

       call convolution_2D_mirror(image_tau, conv_tau, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(image_beta, conv_beta, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(image_Td, conv_Td, dim_y, dim_x, kernel, 3)
       
       call convolution_2D_mirror(conv_tau, conv_conv_tau, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(conv_beta, conv_conv_beta, dim_y, dim_x, kernel, 3)
       call convolution_2D_mirror(conv_Td, conv_conv_Td, dim_y, dim_x, kernel, 3)       
       
       do l=1, dim_x
          do j=1, dim_y
             !Attache aux donnees
             if (params%cc .eqv. .true.) then
                do k=1, dim_v                          
                   deriv(1+(3*(i-1)),j,l) = deriv(1+(3*(i-1)),j,l) + ( &
                        d_mbbcc_l_dtau(wavelength(k), pars(2+(3*(i-1)),j,l), pars(3+(3*(i-1)),j,l), params%l0, &
                        color(k,:), params%degree) &
                        
                        * residual(k,j,l) / std_cube(k,j,l) &
                        )
                   
                   deriv(2+(3*(i-1)),j,l) = deriv(2+(3*(i-1)),j,l) + ( &
                        d_mbbcc_l_db(wavelength(k), pars(1+(3*(i-1)),j,l), pars(2+(3*(i-1)),j,l), pars(3+(3*(i-1)),j,l), &
                        params%l0, color(k,:), params%degree) &
                        
                        * residual(k,j,l) / std_cube(k,j,l) &
                        )
                   
                   deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + ( &
                        d_mbbcc_l_dT(wavelength(k), pars(1+(3*(i-1)),j,l), pars(2+(3*(i-1)),j,l), pars(3+(3*(i-1)),j,l), &
                        params%l0, color(k,:), params%degree) &
                        
                        * residual(k,j,l) / std_cube(k,j,l)&
                        )
                end do
             else
                do k=1, dim_v                          
                   deriv(1+(3*(i-1)),j,l) = deriv(1+(3*(i-1)),j,l) + ( &
                        d_mbb_l_dtau(wavelength(k), pars(2+(3*(i-1)),j,l), pars(3+(3*(i-1)),j,l), params%l0 &
                        ) &                        
                        * residual(k,j,l) / std_cube(k,j,l) &
                        )
                   
                   deriv(2+(3*(i-1)),j,l) = deriv(2+(3*(i-1)),j,l) + ( &
                        d_mbb_l_db(wavelength(k), pars(1+(3*(i-1)),j,l), pars(2+(3*(i-1)),j,l), pars(3+(3*(i-1)),j,l), &
                        params%l0) &                        
                        * residual(k,j,l) / std_cube(k,j,l) &
                        )
                   
                   deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + ( &
                        d_mbb_l_dT(wavelength(k), pars(1+(3*(i-1)),j,l), pars(2+(3*(i-1)),j,l), pars(3+(3*(i-1)),j,l), &
                        params%l0) &
                        * residual(k,j,l) / std_cube(k,j,l)&
                        )
                end do
             end if

             !Regularization function
             if (i .eq. 1) then 
                !Smoothmess 
                f = f + (0.5_xp * params%lambda_tau_cib * conv_tau(j,l)**2._xp) 
                f = f + (0.5_xp * params%lambda_beta_cib * conv_beta(j,l)**2._xp)
                f = f + (0.5_xp * params%lambda_Td_cib * conv_Td(j,l)**2._xp) 

                !Variance
                f = f + (0.5_xp * params%lambda_var_beta * (image_beta(j,l) - 1._xp)**2._xp)
                f = f + (0.5_xp * params%lambda_var_tau_cib * ((image_tau(j,l)))**2._xp)
                ! f = f + (0.5_xp * params%lambda_var_beta_cib * (image_beta(j,l) - c_pars(i))**2._xp)
                f = f + (0.5_xp * params%lambda_var_Td_cib * (image_Td(j,l) - d_pars(i))**2._xp)

             else
                !Smoothmess DUST
                f = f + (0.5_xp * params%lambda_tau * conv_tau(j,l)**2._xp) 
                f = f + (0.5_xp * params%lambda_beta * conv_beta(j,l)**2._xp)
                f = f + (0.5_xp * params%lambda_Td * conv_Td(j,l)**2._xp) 

                !Correlation HI + variance beta and T
                ! f = f + (0.5_xp * params%lambda_var_tau_cib * ((image_tau(j,l)) - b_pars(i))**2._xp)
                f = f + (0.5_xp * params%lambda_var_tau * ((image_tau(j,l)/cube_HI(i,j,l)) - b_pars(i))**2._xp)                
                f = f + (0.5_xp * params%lambda_var_beta * (image_beta(j,l) - c_pars(i))**2._xp)
                f = f + (0.5_xp * params%lambda_var_Td * (image_Td(j,l) - d_pars(i))**2._xp)
             end if
                
             !Regularization gradient
             !Variance
             if (i .eq. 1) then
                deriv(1+(3*(i-1)),j,l) = deriv(1+(3*(i-1)),j,l) + (params%lambda_tau_cib * conv_conv_tau(j,l))
                deriv(2+(3*(i-1)),j,l) = deriv(2+(3*(i-1)),j,l) + (params%lambda_beta_cib * conv_conv_beta(j,l))
                deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + (params%lambda_Td_cib * conv_conv_Td(j,l)) 

                !Variance
                ! deriv(1+(3*(i-1)),j,l) = deriv(1+(3*(i-1)),j,l) + (params%lambda_var_tau_cib * (image_tau(j,l) - b_pars(i)))
                deriv(1+(3*(i-1)),j,l) = deriv(1+(3*(i-1)),j,l) + (params%lambda_var_tau_cib * (image_tau(j,l)))
                deriv(2+(3*(i-1)),j,l) = deriv(2+(3*(i-1)),j,l) + (params%lambda_var_beta_cib * (image_beta(j,l) - 1._xp))
                ! deriv(2+(3*(i-1)),j,l) = deriv(2+(3*(i-1)),j,l) + (params%lambda_var_beta_cib * (image_beta(j,l) - c_pars(i)))
                deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + (params%lambda_var_Td_cib * (image_Td(j,l) - d_pars(i)))

                g(n_cube+(0*params%n_mbb)+i) = 0.
                g(n_cube+(2*params%n_mbb)+i) = g(n_cube+(2*params%n_mbb)+i) - (params%lambda_var_beta_cib &
                     * (image_beta(j,l) - 1._xp))
                ! g(n_cube+(2*params%n_mbb)+i) = g(n_cube+(2*params%n_mbb)+i) - (params%lambda_var_beta_cib &
                !      * (image_beta(j,l) - c_pars(i)))
                g(n_cube+(3*params%n_mbb)+i) = g(n_cube+(3*params%n_mbb)+i) - (params%lambda_var_Td_cib &
                * (image_Td(j,l) - d_pars(i)))                   
                ! g(n_cube+(0*params%n_mbb)+i) = g(n_cube+(0*params%n_mbb)+i) - &
                !(params%lambda_var_tau_cib * (image_tau(j,l) - b_pars(i)))  
             else
                deriv(1+(3*(i-1)),j,l) = deriv(1+(3*(i-1)),j,l) + (params%lambda_tau * conv_conv_tau(j,l))
                deriv(2+(3*(i-1)),j,l) = deriv(2+(3*(i-1)),j,l) + (params%lambda_beta * conv_conv_beta(j,l))
                deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + (params%lambda_Td * conv_conv_Td(j,l)) 

                !Correlation HI + variance
                deriv(1+(3*(i-1)),j,l) = deriv(1+(3*(i-1)),j,l) + (params%lambda_var_tau * &
                     (image_tau(j,l)/cube_HI(i,j,l) - b_pars(i)) / cube_HI(i,j,l))                
                deriv(2+(3*(i-1)),j,l) = deriv(2+(3*(i-1)),j,l) + (params%lambda_var_beta * (image_beta(j,l) - c_pars(i)))
                deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + (params%lambda_var_Td * (image_Td(j,l) - d_pars(i)))

                g(n_cube+(0*params%n_mbb)+i) = g(n_cube+(0*params%n_mbb)+i) - (params%lambda_var_tau &
                     * (image_tau(j,l)/cube_HI(i,j,l) - b_pars(i))) 
                g(n_cube+(2*params%n_mbb)+i) = g(n_cube+(2*params%n_mbb)+i) - (params%lambda_var_beta &
                     * (image_beta(j,l) - c_pars(i)))
                g(n_cube+(3*params%n_mbb)+i) = g(n_cube+(3*params%n_mbb)+i) - (params%lambda_var_Td &
                * (image_Td(j,l) - d_pars(i)))                   
             end if

             ! !Stefan
             ! if (params%lambda_stefan .ne. 0._xp) then                          
             !    f = f + (0.5_xp * params%lambda_stefan * &
             !         (lumi_cst(image_tau(j,l),image_beta(j,l),image_Td(j,l),params%l0) - stefan_pars(i))**2._xp)
             ! end if

             !CIBA
             ! if (params%ciba .eqv. .true.) then 
             !    if (dim_y .gt. 4 .and. i .eq. 1) then
             !       f = f + (0.5_xp * params%lambda_tau_ciba * abs(filtered_itf(j,l))**2._xp)   
             !       ! print*, abs(tf_tau_ciba(j,l))
             !       ! stop
             !    end if
             ! end if

             ! !Stefan
             ! if (params%lambda_stefan .ne. 0._xp) then
             !    g(n_cube+(1*params%n_mbb)+i) = g(n_cube+(1*params%n_mbb)+i) - params%lambda_stefan * ( &
             !         lumi_cst(image_tau(j,l),image_beta(j,l),image_Td(j,l),params%l0) - stefan_pars(i))
             ! end if
             
             ! !Stefan
             ! if (params%lambda_stefan .ne. 0._xp) then
             !    deriv(1+(3*(i-1)),j,l) = deriv(1+(3*(i-1)),j,l) + params%lambda_stefan * ( &
             !         d_lumi_cst_dtau(image_beta(j,l),image_Td(j,l),params%l0) * &
             !         (lumi_cst(image_tau(j,l),image_beta(j,l),image_Td(j,l),params%l0) - stefan_pars(i)) &
             !         )
                
             !    deriv(2+(3*(i-1)),j,l) = deriv(2+(3*(i-1)),j,l) + params%lambda_stefan * ( &
             !         d_lumi_cst_dbeta(image_tau(j,l),image_beta(j,l),image_Td(j,l),params%l0) * &
             !         (lumi_cst(image_tau(j,l),image_beta(j,l),image_Td(j,l),params%l0) - stefan_pars(i)) &
             !         )
             !    deriv(3+(3*(i-1)),j,l) = deriv(3+(3*(i-1)),j,l) + params%lambda_stefan * ( &
             !         d_lumi_cst_dTd(image_tau(j,l),image_beta(j,l),image_Td(j,l),params%l0) * &
             !         (lumi_cst(image_tau(j,l),image_beta(j,l),image_Td(j,l),params%l0) - stefan_pars(i)) &
             !         )
             ! end if

             !CIBA
             ! if (params%ciba .eqv. .true.) then 
             !    if (dim_y .gt. 4 .and. i .eq. 1) then
             !       deriv(1+(3*(i-1)),j,l) = deriv(1+(3*(i-1)),j,l) + (params%lambda_tau_ciba * &
             !            abs(d_filtered_itf(j,l)))             
             !    end if
             ! end if
             
          end do
          !
       end do
    end do        

    ! !CROSS CORRELATION
    ra = cube_HI(2,:,:)-(sum(cube_HI(2,:,:))/size(cube_HI(2,:,:)))
    do i=1,3 !FIXME FOR EACH CIBA PARAMETER MAP BTW 1 AND 3
       rb = pars(i,:,:)
       call crosscorrel(ra,rb,corr,corr_grad)
       do l=1, dim_x
          do j=1, dim_y
             f = f + (0.5_xp * params%lambda_cross * corr(j,l)**2._xp) 
             deriv(i,j,l) = deriv(i,j,l) + (params%lambda_cross * corr_grad(j,l))
          end do
       end do
    end do
    !    
    
    call ravel_3D(deriv, g, 3*params%n_mbb, dim_y, dim_x)

    deallocate(deriv)
    deallocate(residual)
    deallocate(b_pars)
    deallocate(c_pars)
    deallocate(stefan_pars)
    deallocate(pars)
    deallocate(conv_tau, conv_beta, conv_Td)
    deallocate(conv_conv_tau, conv_conv_beta, conv_conv_Td)
    deallocate(image_tau, image_beta, image_Td)

  end subroutine f_g_cube_fast


end module mod_optimize
