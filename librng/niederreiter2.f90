subroutine calcc2 ( dim_num, cj )

!*****************************************************************************80
!
!! CALCC2 computes the constants C(I,J,R).
!
!  Discussion:
!
!    As far as possible, Niederreiter's notation is used.
!
!    For each value of I, we first calculate all the corresponding
!    values of C.  These are held in the array CI.  All these
!    values are either 0 or 1.  
!
!    Next we pack the values into the
!    array CJ, in such a way that CJ(I,R) holds the values of C
!    for the indicated values of I and R and for every value of
!    J from 1 to NBITS.  The most significant bit of CJ(I,R)
!    (not counting the sign bit) is C(I,1,R) and the least
!    significant bit is C(I,NBITS,R).
!
!  Modified:
!
!    30 March 2003
!
!  Author:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter
!
!  Reference:
!
!    Rudolf Lidl, Harald Niederreiter, 
!    Finite Fields,
!    Second Edition,
!    Cambridge University Press, 1997,
!    ISBN: 0521392314,
!    LC: QA247.3.L53
!
!    Harald Niederreiter,
!    Low-discrepancy and low-dispersion sequences,
!    Journal of Number Theory,
!    Volume 30, 1988, pages 51-70.
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the sequence to be generated.
!
!    Output, integer CJ(DIM_MAX,0:NBITS-1), the packed values of 
!    Niederreiter's C(I,J,R)
!
!  Local Parameters:
!
!    Local, integer MAXDEG, the highest degree of polynomial
!    to be handled. 
!
!    Local, integer DIM_MAX, the maximum dimension that will be used.
!
!    Local, integer MAXE; we need DIM_MAX irreducible polynomials over Z2.
!    MAXE is the highest degree among these.
!
!    Local, integer MAXV, the maximum possible index used in V.
!
!    Local, integer NBITS, the number of bits (not counting the sign) in a
!    fixed-point integer.
!
  implicit none

  integer, parameter :: maxdeg = 50
  integer, parameter :: dim_max = 20
  integer, parameter :: maxe = 6
  integer, parameter :: nbits = 31

  integer, parameter :: maxv = nbits + maxe

  integer add(0:1,0:1)
  integer b(-1:maxdeg)
  integer ci(nbits,0:nbits-1)
  integer cj(dim_max,0:nbits-1)
  integer dim_num
  integer e
  integer i
  integer, dimension ( dim_max, -1:maxe ) :: irred
  integer j
  integer mul(0:1,0:1)
  integer px(-1:maxdeg)
  integer r
  integer sub(0:1,0:1)
  integer term
  integer u
  integer v(0:maxv)
!
!  Here we supply the coefficients and the
!  degrees of the first DIM_MAX irreducible polynomials over Z2.
!
!  They are taken from Lidl and Niederreiter.
!
!  The order of these polynomials is the same as the order in
!  file 'gfplys.dat' used by the general program.
!
!  In this block PX(I, -1) is the degree of the Ith polynomial,
!  and PX(I, N) is the coefficient of x**n in the Ith polynomial.
!
  irred(1:dim_max,-1:maxe) = 0

  irred( 1,-1:1) = (/ 1,0,1 /)
  irred( 2,-1:1) = (/ 1,1,1 /)
  irred( 3,-1:2) = (/ 2,1,1,1 /)
  irred( 4,-1:3) = (/ 3,1,1,0,1 /)
  irred( 5,-1:3) = (/ 3,1,0,1,1 /)
  irred( 6,-1:4) = (/ 4,1,1,0,0,1 /)
  irred( 7,-1:4) = (/ 4,1,0,0,1,1 /)
  irred( 8,-1:4) = (/ 4,1,1,1,1,1 /)
  irred( 9,-1:5) = (/ 5,1,0,1,0,0,1 /)
  irred(10,-1:5) = (/ 5,1,0,0,1,0,1 /)
  irred(11,-1:5) = (/ 5,1,1,1,1,0,1 /)
  irred(12,-1:5) = (/ 5,1,1,1,0,1,1 /)
  irred(13,-1:5) = (/ 5,1,1,0,1,1,1 /)
  irred(14,-1:5) = (/ 5,1,0,1,1,1,1 /)
  irred(15,-1:6) = (/ 6,1,1,0,0,0,0,1 /)
  irred(16,-1:6) = (/ 6,1,0,0,1,0,0,1 /)
  irred(17,-1:6) = (/ 6,1,1,1,0,1,0,1 /)
  irred(18,-1:6) = (/ 6,1,1,0,1,1,0,1 /)
  irred(19,-1:6) = (/ 6,1,0,0,0,0,1,1 /)
  irred(20,-1:6) = (/ 6,1,1,1,0,0,1,1 /)
!
!  Prepare to work in Z2.
!
  call setfld2 ( add, mul, sub )

  do i = 1, dim_num
!
!  For each dimension, we need to calculate powers of an
!  appropriate irreducible polynomial:  see Niederreiter
!  page 65, just below equation (19).
!
!  Copy the appropriate irreducible polynomial into PX,
!  and its degree into E.  Set polynomial B = PX ** 0 = 1.
!  M is the degree of B.  Subsequently B will hold higher
!  powers of PX.
!
    e = irred(i,-1)

    do j = -1, e
      px(j) = irred(i,j)
    end do

    b(-1) = 0
    b(0) = 1
!
!  Niederreiter (page 56, after equation (7), defines two
!  variables Q and U.  We do not need Q explicitly, but we do need U.
!
    u = 0

    do j = 1, nbits
!
!  If U = 0, we need to set B to the next power of PX
!  and recalculate V.  This is done by subroutine CALCV2.
!
      if ( u == 0 ) then
        call calcv2 ( maxv, px, add, mul, sub, b, v )
      end if
!
!  Now C is obtained from V.  Niederreiter obtains A from V (page 65, 
!  near the bottom), and then gets C from A (page 56, equation (7)).  
!  However this can be done in one step.  Here CI(J,R) corresponds to
!  Niederreiter's C(I,J,R).
!
      do r = 0, nbits-1
        ci(j,r) = v(r+u)
      end do
!
!  Increment U.  
!
!  If U = E, then U = 0 and in Niederreiter's
!  paper Q = Q + 1.  Here, however, Q is not used explicitly.
!
      u = u + 1
      if ( u == e ) then
        u = 0
      end if

    end do
!
!  The array CI now holds the values of C(I,J,R) for this value
!  of I.  We pack them into array CJ so that CJ(I,R) holds all
!  the values of C(I,J,R) for J from 1 to NBITS.
!
    do r = 0, nbits-1
      term = 0
      do j = 1, nbits
        term = 2 * term + ci(j,r)
      end do
      cj(i,r) = term
    end do

  end do

  return
end
subroutine calcv2 ( maxv, px, add, mul, sub, b, v )

!*****************************************************************************80
!
!! CALCV2 calculates the constants V(J,R).
!
!  Discussion:
!
!    This program calculates the values of the constants V(J,R) as
!    described in the reference (BFN) section 3.3.  It is called from 
!    either CALCC or CALCC2.  
!
!    Polynomials stored as arrays have the coefficient of degree N 
!    in POLY(N), and the degree of the polynomial in POLY(-1).  
!
!    A polynomial which is identically 0 is given degree -1.
!
!  Modified:
!
!    30 March 2003
!
!  Author:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Algorithm 738: 
!    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 4, 1994, pages 494-495.
!
!  Parameters:
!
!    Input, integer MAXV, the dimension of the array V.
!
!    Input, integer PX(-1:MAXDEG), the appropriate irreducible polynomial 
!    for the dimension currently being considered.  The degree of PX 
!    will be called E.
!
!    Input, integer ADD(0:1,0:1), MUL(0:1,0:1), SUB(0:1,0:1), the 
!    addition, multiplication, and subtraction tables, mod 2.
!
!    Input/output, integer B(-1:MAXDEG).  On input, B is the polynomial 
!    defined in section 2.3 of BFN.  The degree of B implicitly defines 
!    the parameter J of section 3.3, by degree(B) = E*(J-1).  On output,
!    B has been multiplied by PX, so its degree is now E * J.
!
!    Output, integer V(0:MAXV), the computed V array.
!
!  Local Parameters:
!
!    Local, integer ARBIT, indicates where the user can place
!    an arbitrary element of the field of order 2.  This means 
!    0 <= ARBIT < 2.  
!
!    Local, integer BIGM, is the M used in section 3.3.
!    It differs from the [little] m used in section 2.3,
!    denoted here by M.
!
!    Local, integer MAXDEG, the highest degree of polynomial
!    to be handled. 
!
!    Local, integer NONZER, shows where the user must put an arbitrary 
!    non-zero element of the field.  For the code, this means 
!    0 < NONZER < 2.
!
  implicit none

  integer, parameter :: maxdeg = 50

  integer add(0:1,0:1)
  integer, parameter :: arbit = 1
  integer b(-1:maxdeg)
  integer bigm
  integer e
  integer h(-1:maxdeg)
  integer i
  integer j
  integer kj
  integer m
  integer maxv
  integer mul(0:1,0:1)
  integer, parameter :: nonzer = 1
  integer p
  integer px(-1:maxdeg)
  integer q
  integer r
  integer sub(0:1,0:1)
  integer term
  integer v(0:maxv)

  p = 2
  q = 2
  e = px(-1)
!
!  The polynomial H is PX**(J-1), which is the value of B on arrival.
!
!  In section 3.3, the values of Hi are defined with a minus sign:
!  don't forget this if you use them later!
!
  do i = -1, b(-1)
    h(i) = b(i)
  end do

  bigm = h(-1)
!
!  Multiply B by PX so B becomes PX**J.
!  In section 2.3, the values of Bi are defined with a minus sign:
!  don't forget this if you use them later!
!
  call plymul2 ( add, mul, px, b, b )
  m = b(-1)
!
!  We don't use J explicitly anywhere, but here it is just in case.
!
  j = m / e
!
!  Now choose a value of Kj as defined in section 3.3.
!  We must have 0 <= Kj < E*J = M.
!  The limit condition on Kj does not seem very relevant
!  in this program.
!
  kj = bigm
!
!  Choose values of V in accordance with the conditions in section 3.3.
!
  do r = 0, kj-1
    v(r) = 0
  end do
  v(kj) = 1

  if ( kj < bigm ) then

    term = sub ( 0, h(kj) )

    do r = kj+1, bigm-1

      v(r) = arbit
!
!  Check the condition of section 3.3,
!  remembering that the H's have the opposite sign.
!
      term = sub ( term, mul ( h(r), v(r) ) )

    end do
!
!  Now V(BIGM) is anything but TERM.
!
    v(bigm) = add ( nonzer, term )

    do r = bigm+1, m-1
      v(r) = arbit
    end do

  else

    do r = kj+1, m-1
      v(r) = arbit
    end do

  end if
!
!  Calculate the remaining V's using the recursion of section 2.3,
!  remembering that the B's have the opposite sign.
!
  do r = 0, maxv-m
    term = 0
    do i = 0, m-1
      term = sub ( term, mul ( b(i), v(r+i) ) )
    end do
    v(r+m) = term
  end do

  return
end

subroutine niederreiter2 ( dim_num, seed, quasi )

!*****************************************************************************80
!
!! NIEDERREITER2 returns an element of the Niederreiter sequence base 2.
!
!  Modified:
!
!    28 March 2003
!
!  Reference:
!
!    Harald Niederreiter,
!    Low-discrepancy and low-dispersion sequences,
!    Journal of Number Theory,
!    Volume 30, 1988, pages 51-70.
!
!  Parameters:
!
!    Input, integer DIM_NUM, the dimension of the sequence to be generated.
!
!    Input/output, integer SEED, the index of the element entry to
!    compute.  On output, SEED is typically reset by this routine
!    to SEED+1.
!
!    Output, real ( kr8 ) QUASI(DIM_NUM), the next quasirandom vector.
!
!  Local Parameters:
!
!    Local, integer DIM_MAX, the maximum dimension that will be used.
!
!    Local, integer NBITS, the number of bits (not counting the sign) in a
!    fixed-point integer.
!
!    Local, real ( kr8 ) RECIP, the multiplier which changes the
!    integers in NEXTQ into the required real values in QUASI.
!
!    Local, integer CJ(DIM_MAX,0:NBITS-1), the packed values of 
!    Niederreiter's C(I,J,R).
!
!    Local, integer DIM_SAVE, the spatial dimension of the sequence
!    as specified on an initialization call.
!
!    Local, integer NEXTQ(DIM_MAX), the numerators of the next item in the
!    series.  These are like Niederreiter's XI(N) (page 54) except that
!    N is implicit, and the NEXTQ are integers.  To obtain
!    the values of XI(N), multiply by RECIP.
!
  USE kinds
  implicit none

  integer dim_num
  integer, parameter :: dim_max = 20
  integer, parameter :: nbits = 31

  integer, save, dimension (dim_max,0:nbits-1) :: cj
  integer, save :: dim_save = 0
  integer gray
  integer i
  integer, save, dimension ( dim_max ) :: nextq
  real    ( kr8 ) quasi(dim_num)
  integer r
  real    ( kr8 ), parameter :: recip = 2.0D+00**(-nbits)
  integer seed
  integer, save :: seed_save = 0
!
!  Initialization.
!
  if ( dim_save < 1 .or. dim_num /= dim_save .or. seed <= 0 ) then

    if ( dim_num <= 0 .or. dim_max < dim_num ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NIEDERREITER2 - Fatal error!'
      write ( *, '(a)' ) '  Bad spatial dimension.'
      stop
    end if

    dim_save = dim_num

    if ( seed < 0 ) then
      seed = 0
    end if

    seed_save = seed
!
!  Calculate the C array.
!
    call calcc2 ( dim_save, cj )

  end if
!
!  Set up NEXTQ appropriately, depending on the Gray code of SEED.
!
!  You can do this every time, starting NEXTQ back at 0,
!  or you can do it once, and then carry the value of NEXTQ
!  around from the previous computation.
!
  if ( seed /= seed_save + 1 ) then

    gray = ieor ( seed, seed / 2 )

    nextq(1:dim_save) = 0

    r = 0

    do while ( gray /= 0 )

      if ( mod ( gray, 2 ) /= 0 ) then
        do i = 1, dim_save
          nextq(i) = ieor ( nextq(i), cj(i,r) )
        end do
      end if

      gray = gray / 2
      r = r + 1

    end do

  end if
!
!  Multiply the numerators in NEXTQ by RECIP to get the next
!  quasi-random vector.
!
  quasi(1:dim_save) = real ( nextq(1:dim_save), kr8 ) * recip
!
!  Find the position of the right-hand zero in SEED.  This
!  is the bit that changes in the Gray-code representation as
!  we go from SEED to SEED+1.
!
  r = 0
  i = seed

  do while ( mod ( i, 2 ) /= 0 )
    r = r + 1
    i = i / 2
  end do
!
!  Check that we have not passed 2**NBITS calls.
!
  if ( nbits <= r ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NIEDERREITER2 - Fatal error!'
    write ( *, '(a)' ) '  Too many calls!'
    stop
  end if
!
!  Compute the new numerators in vector NEXTQ.
!
  do i = 1, dim_save
    nextq(i) = ieor ( nextq(i), cj(i,r) )
  end do

  seed_save = seed
  seed = seed + 1

  return
end
subroutine niederreiter2_generate ( dim_num, n, seed, r )

!*****************************************************************************80
!
!! NIEDERREITER2_GENERATE generates a set of Niederreiter values.
!
!  Modified:
!
!    19 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer N, the number of points desired.
!
!    Input/output, integer SEED, a seed for the random number generator.
!
!    Output, real ( kr8 ) R(DIM_NUM,N), the points.
!
  USE kinds
  implicit none

  integer dim_num
  integer n

  integer j
  real    ( kr8 ) r(dim_num,n)
  integer seed

  do j = 1, n
    call niederreiter2 ( dim_num, seed, r(1:dim_num,j) )
  end do

  return
end
subroutine niederreiter_write ( dim_num, n, base, skip, r, file_out_name )

!*****************************************************************************80
!
!! NIEDERREITER_WRITE writes a set of Niederreiter values to a file.
!
!  Modified:
!
!    04 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIM_NUM, the spatial dimension.
!
!    Input, integer N, the number of points desired.
!
!    Input, integer BASE, the base.
!
!    Input, integer SKIP, the number of initial points skipped.
!
!    Input, real ( kr8 ) R(DIM_NUM,N), the points.
!
!    Input, character ( len = * ) FILE_OUT_NAME, the name of the
!    file to which the output should be written.
!
  USE kinds
  implicit none

  integer dim_num
  integer n

  integer base
  character ( len = * ) file_out_name
  integer file_out_unit
  integer j
  real    ( kr8 ) r(dim_num,n)
  integer skip
  character ( len = 40 ) string

  call get_unit ( file_out_unit )

  open ( unit = file_out_unit, file = file_out_name, status = 'replace' )

  call timestring ( string )

  write ( file_out_unit, '(a)'       ) '#  ' // trim ( file_out_name )
  write ( file_out_unit, '(a)'       ) &
    '#  created by NIEDERREITER2.F90.'
  write ( file_out_unit, '(a)'       ) '#'
  write ( file_out_unit, '(a)'       ) '#  File generated on ' &
    // trim ( string )
  write ( file_out_unit, '(a)'       ) '#'
  write ( file_out_unit, '(a,i8)'    ) &
    '#  Spatial dimension DIM_NUM = ', dim_num
  write ( file_out_unit, '(a,i8)'    ) '#  Number of points N = ', n
  write ( file_out_unit, '(a,g14.6)' ) '#  EPSILON (unit roundoff) = ', &
    epsilon ( r(1,1) )
  write ( file_out_unit, '(a,i8)'    ) '#  Base: ', base
  write ( file_out_unit, '(a,i8)'    ) '#  Initial values skipped = ', skip
  write ( file_out_unit, '(a)'       ) '#'

  write ( string, '(a,i3,a)' ) '(', dim_num, '(2x,f10.6))'
  do j = 1, n
    write ( file_out_unit, string ) r(1:dim_num,j)
  end do

  close ( unit = file_out_unit )

  return
end
subroutine plymul2 ( add, mul, pa, pb, pc )

!*****************************************************************************80
!
!! PLYMUL2 multiplies two polynomials in the field of order 2.
!
!  Discussion:
!
!    Polynomials stored as arrays have the coefficient of degree N in 
!    POLY(N), and the degree of the polynomial in POLY(-1).  
!
!    A polynomial which is identically 0 is given degree -1.
!
!  Modified:
!
!    30 March 2003
!
!  Author:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Algorithm 738: 
!    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 4, 1994, pages 494-495.
!
!  Parameters:
!
!    Input, integer ADD(0:1,0:1), MUL(0:1,0:1),
!    the addition and multiplication tables, mod 2.
!
!    Input, integer PA(-1:MAXDEG), PB(-1:MAXDEG), two polynomials
!    to be multiplied.
!
!    Output, integer PC(-1:MAXDEG), the product polynomial.
!
!  Local Parameters:
!
!    Local, integer MAXDEG, the highest degree of polynomial
!    to be handled. 
!
  implicit none

  integer, parameter :: maxdeg = 50

  integer add(0:1,0:1)
  integer dega
  integer degb
  integer degc
  integer i
  integer j
  integer mul(0:1,0:1)
  integer pa(-1:maxdeg)
  integer pb(-1:maxdeg)
  integer pc(-1:maxdeg)
  integer pt(-1:maxdeg)
  integer term

  dega = pa(-1)
  degb = pb(-1)

  if ( dega == -1 .or. degb == -1 ) then
    degc = -1
  else
    degc = dega + degb
  end if

  if ( maxdeg < degc ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLYMUL - Fatal error!'
    write ( *, '(a)' ) '  Degree of the product exceeds MAXDEG.'
    stop
  end if

  do i = 0, degc
    term = 0
    do j = max ( 0, i-dega ), min ( degb, i )
      term = add ( term, mul ( pa(i-j), pb(j) ) )
    end do
    pt(i) = term
  end do

  pc(-1) = degc

  do i = 0, degc
    pc(i) = pt(i)
  end do

  do i = degc+1, maxdeg
    pc(i) = 0
  end do

  return
end
subroutine setfld2 ( add, mul, sub )

!*****************************************************************************80
!
!! SETFLD2 sets up arithmetic tables for the finite field of order 2.
!
!  Modified:
!
!    30 March 2003
!
!  Author:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Harald Niederreiter,
!    Algorithm 738: 
!    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
!    ACM Transactions on Mathematical Software,
!    Volume 20, Number 4, 1994, pages 494-495.
!
!  Parameters:
!
!    Output, integer ADD(0:1,0:1), MUL(0:1,0:1), SUB(0:1,0:1), the 
!    addition, multiplication, and subtraction tables, mod 2.
!
  implicit none

  integer add(0:1,0:1)
  integer i
  integer j
  integer mul(0:1,0:1)
  integer p
  integer q
  integer sub(0:1,0:1)

  q = 2

  p = 2

  do i = 0, q-1
    do j = 0, q-1
      add(i,j) = mod ( i + j, p )
    end do
  end do

  do i = 0, q-1
    do j = 0, q-1
      mul(i,j) = mod ( i * j, p )
    end do
  end do

  do i = 0, q-1
    do j = 0, q-1
      sub(add(i,j), i) = j
    end do
  end do

  return
end


