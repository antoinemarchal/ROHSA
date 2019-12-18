!! This module contains optimization subroutine and parametric model
module mod_fft
  !! This module contains optimization subroutine and parametric model
  use mod_constants
  use mod_array
  use mod_random

  implicit none
  
  private

  public :: cfft2d, icfft2d, kgrid, crosscorrel

contains

  subroutine crosscorrel(a, b, corr, corr_grad)
    implicit none

    real(xp), intent(in), dimension(:,:), allocatable :: a, b 
    real(xp), intent(inout), dimension(:,:), allocatable :: corr 
    real(xp), intent(inout), dimension(:,:), allocatable :: corr_grad 

    complex(xp), dimension(:,:), allocatable :: ca, cb
    complex(xp), dimension(:,:), allocatable :: fftca, fftcb
    complex(xp), dimension(:,:), allocatable :: cross, ifftcross
    complex(xp), dimension(:,:), allocatable :: cross_grad, ifftcross_grad

    integer :: dimy, dimx

    dimy = size(a,1)
    dimx = size(a,2)

    allocate(ca(dimy,dimx),cb(dimy,dimx),fftca(dimy,dimx),fftcb(dimy,dimx))
    allocate(cross(dimy,dimx), ifftcross(dimy,dimx))

    ca = cmplx(a,0._xp,xp)
    cb = cmplx(b,0._xp,xp)

    call cfft2d(dimy,dimx,ca,fftca)
    call cfft2d(dimy,dimx,cb,fftcb)

    cross = conjg(fftca) * fftcb
    cross_grad = fftca * conjg(fftca) * fftcb

    call icfft2d(dimy,dimx,cross,ifftcross)
    call icfft2d(dimy,dimx,cross_grad,ifftcross_grad)
    
    corr = real(ifftcross,xp)
    corr_grad = real(ifftcross_grad,xp)
    
  end subroutine crosscorrel

  subroutine kgrid(nx,ny,kmat)
    implicit none
    
    integer, intent(in) :: nx, ny
    real(xp), intent(inout), dimension(:,:), allocatable :: kmat

    real(xp), dimension(:), allocatable :: x, y
    real(xp), dimension(:,:), allocatable :: xx, yy

    real(xp), dimension(:,:), allocatable :: kx
    real(xp), dimension(:,:), allocatable :: ky

    integer :: na, nb

    na = ny
    nb = nx

    allocate(x(na), y(nb))
    allocate(xx(nb,na), yy(nb,na))
    allocate(kx(nb,na), ky(nb,na))

    call linspace(x,1._xp,real(na,xp))
    call linspace(y,1._xp,real(nb,xp))
    call meshgrid(x,y,xx,yy)

    xx = xx - 1._xp
    yy = yy - 1._xp

    if (mod(na,2) .eq. 0) then
       kx = (1._xp*xx - ((na)/2._xp) ) / na
    else
       kx = (1._xp*xx - (na-1)/2._xp)/ na
    end if

    if (mod(nb,2) .eq. 0) then
       ky = (1._xp*yy - ((nb)/2._xp) ) / nb
    else
       ky = (1._xp*yy - (nb-1)/2._xp)/ nb
    end if
    
    kmat = sqrt(kx**2._xp + ky**2._xp)

  end subroutine kgrid

  subroutine cfft2d(l,m,data,cfft)
    implicit none
    
    complex(xp), intent(in), allocatable, dimension(:,:) :: data

    integer, intent(in) :: l
    integer, intent(in) :: m
    
    integer ( kind = 4 ) ier
    integer ( kind = 4 ) ldim
    integer ( kind = 4 ) lensav
    integer ( kind = 4 ) lenwrk
    real ( kind = 8 ), allocatable, dimension ( : ) :: work
    real ( kind = 8 ), allocatable, dimension ( : ) :: wsave

    complex(xp), allocatable, dimension(:,:) :: c
    complex(xp), intent(inout), allocatable, dimension(:,:) :: cfft

    allocate(c(l,m))
    
    !  Allocate work arrays.
    lenwrk = 2 * l * m

    lensav = 2 * l + int ( log ( real ( l, kind = 8 ) ) / log ( 2.0D+00 ) ) &
         + 2 * m + int ( log ( real ( m, kind = 8 ) ) / log ( 2.0D+00 ) ) &
         + 8 

    allocate(work(1:lenwrk))
    allocate(wsave(1:lensav))

    call cfft2i ( l, m, wsave, lensav, ier )

    !data input in working complex array c
    c = data

    !  Compute the FFT coefficients.
    ldim = l

    call cfft2f(ldim, l, m, c, wsave, lensav, work, lenwrk, ier)

    !normalization unitary transforms
    c = c * sqrt(real((l*m),xp))

    !working complex array c in output cfft
    cfft = c

    deallocate(c)
    deallocate(work)
    deallocate(wsave)

    return
  end subroutine cfft2d

  subroutine icfft2d(l,m,data,icfft)
    implicit none
    
    complex(xp), intent(in), allocatable, dimension(:,:) :: data

    integer, intent(in) :: l
    integer, intent(in) :: m
    
    integer ( kind = 4 ) ier
    integer ( kind = 4 ) ldim
    integer ( kind = 4 ) lensav
    integer ( kind = 4 ) lenwrk
    real ( kind = 8 ), allocatable, dimension ( : ) :: work
    real ( kind = 8 ), allocatable, dimension ( : ) :: wsave

    complex(xp), allocatable, dimension(:,:) :: c
    complex(xp), intent(inout), allocatable, dimension(:,:) :: icfft

    allocate(c(l,m))
    
    !  Allocate work arrays.
    lenwrk = 2 * l * m

    lensav = 2 * l + int ( log ( real ( l, kind = 8 ) ) / log ( 2.0D+00 ) ) &
         + 2 * m + int ( log ( real ( m, kind = 8 ) ) / log ( 2.0D+00 ) ) &
         + 8 

    allocate ( work(1:lenwrk) )
    allocate ( wsave(1:lensav) )

    call cfft2i ( l, m, wsave, lensav, ier )

    !data input in working complex array c
    c = data

    !  Compute the inverse FFT coefficients.
    ldim = l

    call cfft2b(ldim, l, m, c, wsave, lensav, work, lenwrk, ier)

    !normalization unitary transforms
    c = c / sqrt(real((l*m),xp))

    !working complex array c in output cfft
    icfft = c

    deallocate(c)
    deallocate(work)
    deallocate(wsave)

    return
  end subroutine icfft2d

  ! subroutine cfft2d_test(l,m,data,cfft)
  !   implicit none
    
  !   real(xp), intent(in), allocatable, dimension(:,:) :: data

  !   integer, intent(in) :: l
  !   integer, intent(in) :: m
    
  !   integer ( kind = 4 ) ier
  !   integer ( kind = 4 ) ldim
  !   integer ( kind = 4 ) lensav
  !   integer ( kind = 4 ) lenwrk
  !   ! integer ( kind = 4 ) seed
  !   real ( kind = 8 ), allocatable, dimension ( : ) :: work
  !   real ( kind = 8 ), allocatable, dimension ( : ) :: wsave

  !   complex(xp), allocatable, dimension(:,:) :: c
  !   complex(xp), intent(inout), allocatable, dimension(:,:) :: cfft

  !   allocate(c(l,m))

  !   ! call c8mat_uniform_01(l, m, seed, c)
    
  !   write ( *, '(a)' ) ' '
  !   write ( *, '(a)' ) 'TEST02'
  !   write ( *, '(a)' ) '  For complex double precision fast Fourier transforms, 2D,'
  !   write ( *, '(a)' ) '  CFFT2I initializes the transform,'
  !   write ( *, '(a)' ) '  CFFT2F does a forward transform;'
  !   write ( *, '(a)' ) '  CFFT2B does a backward transform.'
  !   write ( *, '(a)' ) ' '
  !   write ( *, '(a)' ) '  The data is stored in an L by M array, with'
  !   write ( *, '(a,i8)' ) '  L = ', l
  !   write ( *, '(a,i8)' ) '  M = ', m
  !   !
  !   !  Allocate work arrays.
  !   !
  !   lenwrk = 2 * l * m

  !   lensav = 2 * l + int ( log ( real ( l, kind = 8 ) ) / log ( 2.0D+00 ) ) &
  !        + 2 * m + int ( log ( real ( m, kind = 8 ) ) / log ( 2.0D+00 ) ) &
  !        + 8 

  !   write ( *, '(a,i8)' ) '  LENSAV = ', lensav
  !   write ( *, '(a,i8)' ) '  LENWRK = ', lenwrk

  !   allocate ( work(1:lenwrk) )
  !   allocate ( wsave(1:lensav) )

  !   call cfft2i ( l, m, wsave, lensav, ier )
  !   !
  !   !  Set the data values.
  !   !
  !   ! seed = 1973

  !   !data input in working complex array c
  !   c = cmplx(data,0._xp,xp)

  !   call c8mat_print_some ( l, m, c, 1, 1, 5, 5, &
  !        '  Part of the original data:' )
  !   !
  !   !  Compute the FFT coefficients.
  !   !
  !   ldim = l

  !   call cfft2f(ldim, l, m, c, wsave, lensav, work, lenwrk, ier)

  !   call c8mat_print_some(l, m, c, 1, 1, 5, 5, &
  !        '  Part of the FFT coefficients:')
  !   !
  !   !  Compute inverse FFT of coefficients.  Should get back the
  !   !  original data.
  !   !
  !   call cfft2b( ldim, l, m, c, wsave, lensav, work, lenwrk, ier)

  !   call c8mat_print_some( l, m, c, 1, 1, 5, 5, '  Part of the retrieved data:')

  !   !working complex array c in output cfft
  !   cfft = c

  !   deallocate(c)
  !   deallocate(work)
  !   deallocate(wsave)

  !   return
  ! end subroutine cfft2d_test
  
end module mod_fft


