module mod_random
  use mod_constants

  implicit none

  private

  public :: c8mat_print_some, c8mat_uniform_01, c8vec_print_part, c8vec_uniform_01, r8mat_print, &
       r8mat_uniform_01, r8vec_print_part, r8vec_uniform_01

contains

  subroutine c8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

    !*****************************************************************************80
    !
    !! C8MAT_PRINT_SOME prints some of a C8MAT.
    !
    !  Discussion:
    !
    !    A C8MAT is a matrix of C8's.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    23 March 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
    !    in the matrix.
    !
    !    Input, complex ( kind = 8 ) A(M,N), the matrix.
    !
    !    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
    !    column, and the last row and column to be printed.
    !
    !    Input, character ( len = * ) TITLE, a title.
    !
    implicit none

    integer ( kind = 4 ), parameter :: incx = 4
    integer ( kind = 4 ) m
    integer ( kind = 4 ) n

    complex ( kind = 8 ) a(m,n)
    character ( len = 20 ) ctemp(incx)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) i2hi
    integer ( kind = 4 ) i2lo
    integer ( kind = 4 ) ihi
    integer ( kind = 4 ) ilo
    integer ( kind = 4 ) inc
    integer ( kind = 4 ) j
    integer ( kind = 4 ) j2
    integer ( kind = 4 ) j2hi
    integer ( kind = 4 ) j2lo
    integer ( kind = 4 ) jhi
    integer ( kind = 4 ) jlo
    character ( len = * )  title
    complex ( kind = 8 ) zero

    zero = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )

    if ( m <= 0 .or. n <= 0 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) '  (None)' 
       return
    end if
    !
    !  Print the columns of the matrix, in strips of INCX.
    !
    do j2lo = jlo, min ( jhi, n ), incx

       j2hi = j2lo + incx - 1
       j2hi = min ( j2hi, n )
       j2hi = min ( j2hi, jhi )

       inc = j2hi + 1 - j2lo

       write ( *, '(a)' ) ' '

       do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i10,10x)' ) j
       end do

       write ( *, '(a,4a20)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
       write ( *, '(a)' ) '  Row'
       write ( *, '(a)' ) '  ---'
       !
       !  Determine the range of the rows in this strip.
       !
       i2lo = max ( ilo, 1 )
       i2hi = min ( ihi, m )

       do i = i2lo, i2hi
          !
          !  Print out (up to) INCX entries in row I, that lie in the current strip.
          !
          do j2 = 1, inc

             j = j2lo - 1 + j2

             if ( a(i,j) == zero ) then
                ctemp(j2) = '       0.0          '
             else if ( imag ( a(i,j) ) == 0.0D+00 ) then
                write ( ctemp(j2), '(g10.3,10x)' ) real ( a(i,j), kind = 8 )
             else
                write ( ctemp(j2), '(2g10.3)' ) a(i,j)
             end if

          end do

          write ( *, '(i5,a1,4a20)' ) i, ':', ( ctemp(j2), j2 = 1, inc )

       end do

    end do

    return
  end subroutine c8mat_print_some
  subroutine c8mat_uniform_01 ( m, n, seed, c )

    !*****************************************************************************80
    !
    !! C8MAT_UNIFORM_01 returns a unit pseudorandom C8MAT.
    !
    !  Discussion:
    !
    !    A C8MAT is a matrix of C8's.
    !
    !    The angles should be uniformly distributed between 0 and 2 * PI,
    !    the square roots of the radius uniformly distributed between 0 and 1.
    !
    !    This results in a uniform distribution of values in the unit circle.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    15 March 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, N, the number of rows and columns
    !    in the matrix.
    !
    !    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
    !    should NOT be 0.  On output, SEED has been updated.
    !
    !    Output, complex ( kind = 8 ) C(M,N), the pseudorandom complex matrix.
    !
    implicit none

    integer ( kind = 4 ) m
    integer ( kind = 4 ) n

    complex ( kind = 8 ) c(m,n)
    integer ( kind = 4 ) i
    integer ( kind = 4 ), parameter :: i4_huge = 2147483647
    integer ( kind = 4 ) j
    integer ( kind = 4 ) k
    real ( kind = 8 ) r
    real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
    integer ( kind = 4 ) seed
    real ( kind = 8 ) theta

    do j = 1, n
       do i = 1, m

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed < 0 ) then
             seed = seed + i4_huge
          end if

          r = sqrt ( real ( seed, kind = 8 ) * 4.656612875D-10 )

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed < 0 ) then
             seed = seed + i4_huge
          end if

          theta = 2.0D+00 * r8_pi * ( real ( seed, kind = 8 ) * 4.656612875D-10 )

          c(i,j) = r * cmplx ( cos ( theta ), sin ( theta ), kind = 8 )

       end do

    end do

    return
  end subroutine c8mat_uniform_01
  subroutine c8vec_print_part ( n, a, max_print, title )

    !*****************************************************************************80
    !
    !! C8VEC_PRINT_PART prints "part" of a C8VEC.
    !
    !  Discussion:
    !
    !    The user specifies MAX_PRINT, the maximum number of lines to print.
    !
    !    If N, the size of the vector, is no more than MAX_PRINT, then
    !    the entire vector is printed, one entry per line.
    !
    !    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
    !    followed by a line of periods suggesting an omission,
    !    and the last entry.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    22 June 2010
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the number of entries of the vector.
    !
    !    Input, complex ( kind = 8 ) A(N), the vector to be printed.
    !
    !    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines
    !    to print.
    !
    !    Input, character ( len = * ) TITLE, a title.
    !
    implicit none

    integer ( kind = 4 ) n

    complex ( kind = 8 ) a(n)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) max_print
    character ( len = * )  title

    if ( max_print <= 0 ) then
       return
    end if

    if ( n <= 0 ) then
       return
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
    write ( *, '(a)' ) ' '

    if ( n <= max_print ) then

       do i = 1, n
          write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
       end do

    else if ( 3 <= max_print ) then

       do i = 1, max_print - 2
          write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
       end do
       write ( *, '(a)' ) '  ........  ..............  ..............'
       i = n
       write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)

    else

       do i = 1, max_print - 1
          write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
       end do
       i = max_print
       write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,a)' ) i, ':', a(i), &
            '...more entries...'

    end if

    return
  end subroutine c8vec_print_part
  subroutine c8vec_uniform_01 ( n, seed, c )

    !*****************************************************************************80
    !
    !! C8VEC_UNIFORM_01 returns a unit pseudorandom C8VEC.
    !
    !  Discussion:
    !
    !    A C8VEC is a vector of C8's.
    !
    !    The angles should be uniformly distributed between 0 and 2 * PI,
    !    the square roots of the radius uniformly distributed between 0 and 1.
    !
    !    This results in a uniform distribution of values in the unit circle.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    15 March 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the number of values to compute.
    !
    !    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should 
    !    NOT be 0.  On output, SEED has been updated.
    !
    !    Output, complex ( kind = 8 ) C(N), the pseudorandom complex vector.
    !
    implicit none

    integer ( kind = 4 ) n

    complex ( kind = 8 ) c(n)
    integer ( kind = 4 ) i
    integer ( kind = 4 ), parameter :: i4_huge = 2147483647
    integer ( kind = 4 ) k
    real    ( kind = 8 ) r
    real    ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
    integer ( kind = 4 ) seed
    real    ( kind = 8 ) theta

    do i = 1, n

       k = seed / 127773

       seed = 16807 * ( seed - k * 127773 ) - k * 2836

       if ( seed < 0 ) then
          seed = seed + i4_huge
       end if

       r = sqrt ( real ( seed, kind = 8 ) * 4.656612875D-10 )

       k = seed / 127773

       seed = 16807 * ( seed - k * 127773 ) - k * 2836

       if ( seed < 0 ) then
          seed = seed + i4_huge
       end if

       theta = 2.0D+00 * r8_pi * ( real ( seed, kind = 8 ) * 4.656612875D-10 )

       c(i) = r * cmplx ( cos ( theta ), sin ( theta ), kind = 8 )

    end do

    return
  end subroutine c8vec_uniform_01
  subroutine r8mat_print ( m, n, a, title )

    !*****************************************************************************80
    !
    !! R8MAT_PRINT prints an R8MAT.
    !
    !  Discussion:
    !
    !    An R8MAT is an array of R8 values.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    12 September 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, the number of rows in A.
    !
    !    Input, integer ( kind = 4 ) N, the number of columns in A.
    !
    !    Input, real ( kind = 8 ) A(M,N), the matrix.
    !
    !    Input, character ( len = * ) TITLE, a title.
    !
    implicit none

    integer ( kind = 4 ) m
    integer ( kind = 4 ) n

    real ( kind = 8 ) a(m,n)
    character ( len = * )  title

    call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

    return
  end subroutine r8mat_print
  subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

    !*****************************************************************************80
    !
    !! R8MAT_PRINT_SOME prints some of an R8MAT.
    !
    !  Discussion:
    !
    !    An R8MAT is an array of R8 values.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    10 September 2009
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
    !
    !    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
    !
    !    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
    !
    !    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
    !
    !    Input, character ( len = * ) TITLE, a title.
    !
    implicit none

    integer ( kind = 4 ), parameter :: incx = 5
    integer ( kind = 4 ) m
    integer ( kind = 4 ) n

    real ( kind = 8 ) a(m,n)
    character ( len = 14 ) ctemp(incx)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) i2hi
    integer ( kind = 4 ) i2lo
    integer ( kind = 4 ) ihi
    integer ( kind = 4 ) ilo
    integer ( kind = 4 ) inc
    integer ( kind = 4 ) j
    integer ( kind = 4 ) j2
    integer ( kind = 4 ) j2hi
    integer ( kind = 4 ) j2lo
    integer ( kind = 4 ) jhi
    integer ( kind = 4 ) jlo
    character ( len = * ) title

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )

    if ( m <= 0 .or. n <= 0 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) '  (None)'
       return
    end if

    do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

       j2hi = j2lo + incx - 1
       j2hi = min ( j2hi, n )
       j2hi = min ( j2hi, jhi )

       inc = j2hi + 1 - j2lo

       write ( *, '(a)' ) ' '

       do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i8,6x)' ) j
       end do

       write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
       write ( *, '(a)' ) '  Row'
       write ( *, '(a)' ) ' '

       i2lo = max ( ilo, 1 )
       i2hi = min ( ihi, m )

       do i = i2lo, i2hi

          do j2 = 1, inc

             j = j2lo - 1 + j2

             if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
                write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
             else
                write ( ctemp(j2), '(g14.6)' ) a(i,j)
             end if

          end do

          write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

       end do

    end do

    return
  end subroutine r8mat_print_some
  subroutine r8mat_uniform_01 ( m, n, seed, r )

    !*****************************************************************************80
    !
    !! R8MAT_UNIFORM_01 fills an R8MAT with unit pseudorandom numbers.
    !
    !  Discussion:
    !
    !    An R8MAT is an array of R8 values.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    11 August 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Paul Bratley, Bennett Fox, Linus Schrage,
    !    A Guide to Simulation,
    !    Springer Verlag, pages 201-202, 1983.
    !
    !    Bennett Fox,
    !    Algorithm 647:
    !    Implementation and Relative Efficiency of Quasirandom
    !    Sequence Generators,
    !    ACM Transactions on Mathematical Software,
    !    Volume 12, Number 4, pages 362-376, 1986.
    !
    !    Peter Lewis, Allen Goodman, James Miller,
    !    A Pseudo-Random Number Generator for the System/360,
    !    IBM Systems Journal,
    !    Volume 8, pages 136-143, 1969.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, N, the number of rows and columns in
    !    the array.
    !
    !    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
    !    should NOT be 0.  On output, SEED has been updated.
    !
    !    Output, real ( kind = 8 ) R(M,N), the array of pseudorandom values.
    !
    implicit none

    integer ( kind = 4 ) m
    integer ( kind = 4 ) n

    integer ( kind = 4 ) i
    integer ( kind = 4 ), parameter :: i4_huge = 2147483647
    integer ( kind = 4 ) j
    integer ( kind = 4 ) k
    integer ( kind = 4 ) seed
    real ( kind = 8 ) r(m,n)

    do j = 1, n

       do i = 1, m

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed < 0 ) then
             seed = seed + i4_huge
          end if

          r(i,j) = real ( seed, kind = 8 ) * 4.656612875D-10

       end do
    end do

    return
  end subroutine r8mat_uniform_01
  subroutine r8vec_print_part ( n, a, max_print, title )

    !*****************************************************************************80
    !
    !! R8VEC_PRINT_PART prints "part" of an R8VEC.
    !
    !  Discussion:
    !
    !    The user specifies MAX_PRINT, the maximum number of lines to print.
    !
    !    If N, the size of the vector, is no more than MAX_PRINT, then
    !    the entire vector is printed, one entry per line.
    !
    !    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
    !    followed by a line of periods suggesting an omission,
    !    and the last entry.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    19 December 2001
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the number of entries of the vector.
    !
    !    Input, real ( kind = 8 ) A(N), the vector to be printed.
    !
    !    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines
    !    to print.
    !
    !    Input, character ( len = * ) TITLE, a title.
    !
    implicit none

    integer ( kind = 4 ) n

    real ( kind = 8 ) a(n)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) max_print
    character ( len = * ) title

    if ( max_print <= 0 ) then
       return
    end if

    if ( n <= 0 ) then
       return
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
    write ( *, '(a)' ) ' '

    if ( n <= max_print ) then

       do i = 1, n
          write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
       end do

    else if ( 3 <= max_print ) then

       do i = 1, max_print - 2
          write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
       end do
       write ( *, '(a)' ) '  ........  ..............'
       i = n
       write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)

    else

       do i = 1, max_print - 1
          write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
       end do
       i = max_print
       write ( *, '(2x,i8,a,1x,g14.6,2x,a)' ) i, ':', a(i), '...more entries...'

    end if

    return
  end subroutine r8vec_print_part
  subroutine r8vec_uniform_01 ( n, seed, r )

    !*****************************************************************************80
    !
    !! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
    !
    !  Discussion:
    !
    !    An R8VEC is a vector of R8's.
    !
    !    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 July 2006
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Paul Bratley, Bennett Fox, Linus Schrage,
    !    A Guide to Simulation,
    !    Springer Verlag, pages 201-202, 1983.
    !
    !    Bennett Fox,
    !    Algorithm 647:
    !    Implementation and Relative Efficiency of Quasirandom
    !    Sequence Generators,
    !    ACM Transactions on Mathematical Software,
    !    Volume 12, Number 4, pages 362-376, 1986.
    !
    !    Peter Lewis, Allen Goodman, James Miller
    !    A Pseudo-Random Number Generator for the System/360,
    !    IBM Systems Journal,
    !    Volume 8, pages 136-143, 1969.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the number of entries in the vector.
    !
    !    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
    !    should NOT be 0.  On output, SEED has been updated.
    !
    !    Output, real ( kind = 8 ) R(N), the vector of pseudorandom values.
    !
    implicit none

    integer ( kind = 4 ) n

    integer ( kind = 4 ) i
    integer ( kind = 4 ) k
    integer ( kind = 4 ) seed
    real ( kind = 8 ) r(n)

    if ( seed == 0 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
       write ( *, '(a)' ) '  Input value of SEED = 0.'
       stop
    end if

    do i = 1, n

       k = seed / 127773

       seed = 16807 * ( seed - k * 127773 ) - k * 2836

       if ( seed < 0 ) then
          seed = seed + 2147483647
       end if

       r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

    end do

    return
  end subroutine r8vec_uniform_01


end module mod_random
