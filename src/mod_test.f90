!! 
module mod_test
  !! 
  use mod_constants
  use mod_convert
  use mod_read_parameters
  use mod_array
  use mod_fft
  use mod_model
  use mod_rw_data
  
  implicit none

  private

  public :: test
  
contains
  
  subroutine test()
    implicit none
    
    character(len=512) :: filename_fBm="fBm.dat"
    
    real(xp), dimension(:,:), allocatable    :: test_fft
    real(xp), dimension(:,:), allocatable    :: tapper 
    real(xp), dimension(:,:), allocatable    :: butter 
    real(xp), dimension(:,:), allocatable    :: test_fft_shift 
    complex(xp), dimension(:,:), allocatable :: c_test_fft 
    complex(xp), dimension(:,:), allocatable :: c_test_fft2 
    complex(xp), dimension(:,:), allocatable :: c_test_fft3
    real(xp), dimension(:,:), allocatable :: kmat

    call read_map(filename_fBm, test_fft)
    allocate(test_fft_shift(64,64))

    call shift(test_fft, test_fft_shift)
    allocate(c_test_fft(64,64), c_test_fft2(64,64))
    c_test_fft = cmplx(test_fft_shift,0._xp,xp)

    call cfft2d(64,64,c_test_fft,c_test_fft2)
    call icfft2d(64,64,c_test_fft2,c_test_fft3)    
    ! print*, c_test_fft(1,:)
    ! print*, real(c_test_fft3(1,:),xp)

    allocate(kmat(64,64))
    call kgrid(64,64,kmat)

    allocate(tapper(34,64))
    call apodize(tapper, 0.86_xp, 34,64)

    !test normalized FFT unity transform
    ! print*, c_test_fft2(1,:)
    ! print*, tapper(1,:)

    call butterworth(butter,kmat,1._xp,1._xp,2._xp)

  end subroutine test


end module mod_test
