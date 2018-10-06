!! This module contains tools to convert fortran object
module mod_convert
  !! This module contains tools to convert fortran object
  implicit none

  private

  public :: str
  
contains

  character(len=20) function str(k)
    !! Convert an integer to string
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
  end function str

end module mod_convert
