!! This module read the input user parameters (parameters.txt file / data / noise if true)
module mod_read_parameters
  !! This module read the input user parameters (parameters.txt file / data / noise if true)
  use mod_constants
  use mod_convert
  
  implicit none

  type(parameters) :: params
  
  private
  
  public :: get_parameters, print_parameters, print_parameters_unit, params

contains

  subroutine get_parameters(filename_parameters)
    implicit none 
    
    character(len=512), intent(in) :: filename_parameters
    integer :: ios=0

    namelist /user_parameters/ params 

    open(unit=11, file=filename_parameters, status="old", iostat=ios)
    if (ios /= 0) stop "opening file error"
    read(11, user_parameters)    
    close(11)

    ! Display parameters
    call print_parameters()
        
  end subroutine get_parameters


  subroutine print_parameters()
    implicit none 

    print*, "filename_r = '",trim(params%filename_q),"'"
    print*, "filename_i = '",trim(params%filename_u),"'"

    print*, "fileout = '",trim(params%fileout),"'"
    print*, "timeout = '",trim(params%timeout),"'"
    
    print*, " "
    print*, "______Parameters_____"
    print*, "n = ", params%n

    print*, "lambda_amp = ", params%lambda_amp
    print*, "lambda_mu = ", params%lambda_mu

    print*, " "
    print*, "amp_init = ", params%amp_init
    print*, "mu_init = ", params%mu_init

    print*, " "
    print*, "lb_amp = ", params%lb_amp
    print*, "ub_amp = ", params%ub_amp
    print*, "lb_mu = ", params%lb_mu
    print*, "ub_mu = ", params%ub_mu

    print*, " "
    print*, "maxiter_init = ", params%maxiter_init
    print*, "maxiter = ", params%maxiter
    print*, "lstd = ", params%lstd
    print*, "ustd = ", params%ustd
    print*, "noise = ", params%noise
    print*, "save_grid = ", params%save_grid

    print*, " "
    
  end subroutine print_parameters

  subroutine print_parameters_unit(unit)
    implicit none 

    integer, intent(in) :: unit

    write(unit,fmt=*) "#filename_r = '",trim(params%filename_q),"'"
    write(unit,fmt=*) "#filename_i = '",trim(params%filename_u),"'"

    write(unit,fmt=*) "#fileout = '",trim(params%fileout),"'"
    write(unit,fmt=*) "#timeout = '",trim(params%timeout),"'"
    
    write(unit,fmt=*) "# "
    write(unit,fmt=*) "#______Parameters_____"
    write(unit,fmt=*) "#n = ", params%n

    write(unit,fmt=*) "#lambda_amp = ", params%lambda_amp
    write(unit,fmt=*) "#lambda_mu = ", params%lambda_mu

    write(unit,fmt=*) "# "
    write(unit,fmt=*) "#amp_init = ", params%amp_init
    write(unit,fmt=*) "#mu_init = ", params%mu_init

    write(unit,fmt=*) "# "
    write(unit,fmt=*) "#lb_amp = ", params%lb_amp
    write(unit,fmt=*) "#ub_amp = ", params%ub_amp
    write(unit,fmt=*) "#lb_mu = ", params%lb_mu
    write(unit,fmt=*) "#ub_mu = ", params%ub_mu

    write(unit,fmt=*) "# "
    write(unit,fmt=*) "#maxiter_init = ", params%maxiter_init
    write(unit,fmt=*) "#maxiter = ", params%maxiter
    write(unit,fmt=*) "#lstd = ", params%lstd
    write(unit,fmt=*) "#ustd = ", params%ustd
    write(unit,fmt=*) "#noise = ", params%noise
    write(unit,fmt=*) "#save_grid = ", params%save_grid

    write(unit,fmt=*) "# "
    
  end subroutine print_parameters_unit


end Module mod_read_parameters
