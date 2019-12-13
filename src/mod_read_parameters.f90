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

    print*, "filename = '",trim(params%filename),"'"

    print*, "fileout = '",trim(params%fileout),"'"
    print*, "timeout = '",trim(params%timeout),"'"
    
    print*, " "
    print*, "______Parameters_____"
    print*, "n_mbb = ", params%n_mbb

    print*, "lambda_tau = ", params%lambda_tau
    print*, "lambda_beta = ", params%lambda_beta
    print*, "lambda_Td = ", params%lambda_Td

    print*, " "
    print*, "lambda_var_tau = ", params%lambda_var_tau
    print*, "lambda_var_beta = ", params%lambda_var_beta
    print*, "lambda_var_Td = ", params%lambda_var_Td
    print*, "lambda_stefan = ", params%lambda_stefan

    print*, " "
    print*, "tau_init = ", params%tau_init
    print*, "beta_init = ", params%beta_init
    print*, "Td_init = ", params%Td_init

    print*, " "
    print*, "tau_init_cib = ", params%tau_init_cib
    print*, "beta_init_cib = ", params%beta_init_cib
    print*, "Td_init_cib = ", params%Td_init_cib

    print*, " "
    print*, "lb_tau = ", params%lb_tau
    print*, "ub_tau = ", params%ub_tau
    print*, "lb_beta = ", params%lb_beta
    print*, "ub_beta = ", params%ub_beta
    print*, "lb_Td = ", params%lb_Td
    print*, "ub_Td = ", params%ub_Td

    print*, " "
    print*, "lb_tau_cib = ", params%lb_tau_cib
    print*, "ub_tau_cib = ", params%ub_tau_cib
    print*, "lb_beta_cib = ", params%lb_beta_cib
    print*, "ub_beta_cib = ", params%ub_beta_cib
    print*, "lb_Td_cib = ", params%lb_Td_cib
    print*, "ub_Td_cib = ", params%ub_Td_cib

    print*, " "
    print*, "l0 = ", params%l0

    print*, " "
    print*, "maxiter_init = ", params%maxiter_init
    print*, "maxiter = ", params%maxiter
    print*, "lstd = ", params%lstd
    print*, "ustd = ", params%ustd
    print*, "noise = ", params%noise
    print*, "save_grid = ", params%save_grid
    print*, "cc = ", params%cc
    print*, "ciba = ", params%ciba
    print*, "radius_tapper = ", params%radius_tapper

    print*, " "
    
  end subroutine print_parameters

  subroutine print_parameters_unit(unit)
    implicit none 

    integer, intent(in) :: unit

    write(unit,fmt=*) "#filename = '",trim(params%filename),"'"

    write(unit,fmt=*) "#fileout = '",trim(params%fileout),"'"
    write(unit,fmt=*) "#timeout = '",trim(params%timeout),"'"
    
    write(unit,fmt=*) "# "
    write(unit,fmt=*) "#______Parameters_____"
    write(unit,fmt=*) "#n_mbb = ", params%n_mbb

    write(unit,fmt=*) "#lambda_tau = ", params%lambda_tau
    write(unit,fmt=*) "#lambda_beta = ", params%lambda_beta
    write(unit,fmt=*) "#lambda_Td = ", params%lambda_Td

    write(unit,fmt=*) "# "
    write(unit,fmt=*) "#lambda_var_tau = ", params%lambda_var_tau
    write(unit,fmt=*) "#lambda_var_beta = ", params%lambda_var_beta
    write(unit,fmt=*) "#lambda_var_Td = ", params%lambda_var_Td
    write(unit,fmt=*) "#lambda_stefan = ", params%lambda_stefan

    write(unit,fmt=*) "# "
    write(unit,fmt=*) "#tau_init = ", params%tau_init
    write(unit,fmt=*) "#beta_init = ", params%beta_init
    write(unit,fmt=*) "#Td_init = ", params%Td_init

    write(unit,fmt=*) " "
    write(unit,fmt=*) "#tau_init_cib = ", params%tau_init_cib
    write(unit,fmt=*) "#beta_init_cib = ", params%beta_init_cib
    write(unit,fmt=*) "#Td_init_cib = ", params%Td_init_cib

    write(unit,fmt=*) "# "
    write(unit,fmt=*) "#lb_tau = ", params%lb_tau
    write(unit,fmt=*) "#ub_tau = ", params%ub_tau
    write(unit,fmt=*) "#lb_beta = ", params%lb_beta
    write(unit,fmt=*) "#ub_beta = ", params%ub_beta
    write(unit,fmt=*) "#lb_Td = ", params%lb_Td
    write(unit,fmt=*) "#ub_Td = ", params%ub_Td

    write(unit,fmt=*) "# "
    write(unit,fmt=*) "#lb_tau_cib = ", params%lb_tau_cib
    write(unit,fmt=*) "#ub_tau_cib = ", params%ub_tau_cib
    write(unit,fmt=*) "#lb_beta_cib = ", params%lb_beta_cib
    write(unit,fmt=*) "#ub_beta_cib = ", params%ub_beta_cib
    write(unit,fmt=*) "#lb_Td_cib = ", params%lb_Td_cib
    write(unit,fmt=*) "#ub_Td_cib = ", params%ub_Td_cib

    write(unit,fmt=*) "# "
    write(unit,fmt=*) "#l0 = ", params%l0

    write(unit,fmt=*) "# "
    write(unit,fmt=*) "#maxiter_init = ", params%maxiter_init
    write(unit,fmt=*) "#maxiter = ", params%maxiter
    write(unit,fmt=*) "#lstd = ", params%lstd
    write(unit,fmt=*) "#ustd = ", params%ustd
    write(unit,fmt=*) "#noise = ", params%noise
    write(unit,fmt=*) "#save_grid = ", params%save_grid
    write(unit,fmt=*) "#cc = ", params%cc
    write(unit,fmt=*) "#ciba = ", params%ciba
    write(unit,fmt=*) "#radius_tapper = ", params%radius_tapper

    write(unit,fmt=*) "# "
    
  end subroutine print_parameters_unit


end Module mod_read_parameters
