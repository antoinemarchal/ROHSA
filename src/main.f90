!! Regularized Optimization for Hypers-spectral Analysis (ROHSA)
program ROHSA

  use mod_constants
  use mod_start
  use mod_read_parameters
  use mod_instrument
  use mod_rw_data
  use mod_rohsa
  
  implicit none

  character(len=512) :: filename_parameters
  real(xp) :: start, finish

  call cpu_time(start)
  call header()  
  call get_command_argument(1, filename_parameters)
  call get_parameters(filename_parameters) 
  ! call compute_constant()
  call get_data()
  call main_rohsa() 
  call ender()

  call cpu_time(finish)
  print '("Time = ",f6.3," seconds.")',finish-start
      
end program ROHSA
