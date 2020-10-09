!! This module contains optimization subroutine and parametric model
module mod_model
  !! This module contains optimization subroutine and parametric model
  use mod_constants

  implicit none
  
  private

  public :: gaussian, func_q, func_u, f_rmsf, df_rmsf_da, df_rmsf_dm

contains

  pure function f_rmsf(x, a, m, cc)
    !! Approximation rmsf
    implicit none
    
    real(xp), intent(in) :: x
    real(xp), intent(in) :: a, m
    real(xp), intent(in), dimension(8) :: cc
    real(xp) :: f_rmsf

    integer :: i

    f_rmsf = 0._xp

    do i=1, 4
       f_rmsf = f_rmsf + a * (cc(1+(2*(i-1))) * exp(-(x-m)**2._xp / (2._xp * cc(2+(2*(i-1)))**2)))
    end do
    
  end function f_rmsf


  pure function df_rmsf_da(x, m, cc)
    !! Derivative with respect to a of f_rmsf
    implicit none
    
    real(xp), intent(in) :: x
    real(xp), intent(in) :: m
    real(xp), intent(in), dimension(8) :: cc
    real(xp) :: df_rmsf_da

    integer :: i

    df_rmsf_da = 0._xp

    do i=1, 4
       df_rmsf_da = df_rmsf_da + (cc(1+(2*(i-1))) * exp(-(x-m)**2._xp / (2._xp * cc(2+(2*(i-1)))**2)))
    end do
    
  end function df_rmsf_da


  pure function df_rmsf_dm(x, a, m, cc)
    !! Derivative with respect to m of f_rmsf
    implicit none
    
    real(xp), intent(in) :: x
    real(xp), intent(in) :: a, m
    real(xp), intent(in), dimension(8) :: cc
    real(xp) :: df_rmsf_dm

    integer :: i

    df_rmsf_dm = 0._xp

    do i=1, 4
       df_rmsf_dm = df_rmsf_dm + a*(x-m)/cc(2+(2*(i-1)))**2 * (cc(1+(2*(i-1))) * exp(-(x-m)**2._xp / (2._xp * cc(2+(2*(i-1)))**2)))
    end do
    
  end function df_rmsf_dm


  pure function gaussian(x, a, s)
    !! Gaussian function   
    implicit none
    
    real(xp), intent(in) :: x
    real(xp), intent(in) :: a, s
    real(xp) :: gaussian

    gaussian = a * exp(- x**2._xp / (2._xp * s**2))
  end function gaussian


  pure function func_q(x, a, mu, w, fn)
    !! Gaussian function   
    implicit none
    
    real(xp), intent(in) :: x
    real(xp), intent(in), dimension(:) :: w
    real(xp), intent(in) :: a, mu
    integer, intent(in) :: fn

    real(xp) :: func_q

    func_q = a / fn * sum(cos(-2_xp*(x-mu)*w**2_xp));
  end function func_q


  pure function func_u(x, a, mu, w, fn)
    !! Gaussian function   
    implicit none
    
    real(xp), intent(in) :: x
    real(xp), intent(in), dimension(:) :: w
    real(xp), intent(in) :: a, mu
    integer, intent(in) :: fn

    real(xp) :: func_u

    func_u = a / fn * sum(sin(-2_xp*(x-mu)*w**2_xp));
  end function func_u

end module mod_model
