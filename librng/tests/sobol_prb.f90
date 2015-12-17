      program main

      !*****************************************************************************80
      !
      !! SOBOL_PRB calls a set of problems for SOBOL.
      !
      !  Modified:
      !
      !    05 July 2006
      !
      !  Author:
      !
      !    John Burkardt
      !
      
      
      USE qr_generator_mod
      USE cmaes_param_mod
        implicit none
        
        
        INTEGER,PARAMETER                       :: n = 50
        INTEGER :: i
        REAL*8, DIMENSION(n)         :: nums,nnums
        REAL*8				:: r8i8_uniform_01
        INTEGER             :: seed
        integer ( kind = 4 ) tau_sobol
      INTEGER MAXBIT,INIT,k,method,l
      INTEGER N1,N2,DIMEN,OFFSET,TRANSFORM
      PARAMETER (N1=1000,DIMEN=10)
      INTEGER BASE(DIMEN)
      REAL(MK),DIMENSION(DIMEN,N1)                          :: QN1,QN2
      CHARACTER(len=5):: j,mstring,dstring
      seed = 12356
      TRANSFORM = 0
      

      INIT = 1
      OFFSET = 0

      DO method =0, 5
      
     
      write(mstring,'(I3)') method
      DO k = 0, 3
      write(j,'(I3)') k
      CALL test(DIMEN,N1,seed,QN1,method,k)
      !CALL invnor(DIMEN,N1,QN1,QN2,k)
      
      DO l = 1,DIMEN
          write(dstring,'(I3)') l
          OPEN (UNIT=23, FILE='method'//TRIM(adjustl(mstring))//'inversion'//&
          TRIM(adjustl(j))//'dim'//TRIM(adjustl(dstring))&
     &    //'.txt',status='replace',action='write')!
    !      seed = tau_sobol(1)
          DO i = 1, N1
      	    WRITE (23, *)QN1(l,i)
          ENDDO
          CLOSE (23)    
          END DO
      END DO
      END DO





        stop
      end
      subroutine test01

      !*****************************************************************************80
      !
      !! TEST01 tests IEOR on integer ( kind = 4 ) arguments.
      !
      !  Discussion:
      !
      !    The FORTRAN90 function IEOR returns the bitwise exclusive OR 
      !    of two integers.
      !
      !  Modified:
      !
      !    05 July 2006
      !
      !  Author:
      !
      !    John Burkardt
      !
        implicit none

        integer ( kind = 4 ) i
        integer ( kind = 4 ) i4_uniform
        integer ( kind = 4 ) j
        integer ( kind = 4 ) k
        integer ( kind = 4 ) seed
        integer ( kind = 4 ) test

        seed = 123456789

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST01'
        write ( *, '(a)' ) '  IEOR is a FORTRAN90 function which returns the'
        write ( *, '(a)' ) '  bitwise exclusive OR of two integers.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '         I         J    IEOR(I,J)'
        write ( *, '(a)' ) ' '

        do test = 1, 10
          i = i4_uniform ( 0, 100, seed )
          j = i4_uniform ( 0, 100, seed )
          k = ieor ( i, j )
          write ( *, '(2x,i8,2x,i8,2x,i8)' ) i, j, k
        end do

        return
      end
      subroutine test02

      !*****************************************************************************80
      !
      !! TEST02 tests I4_BIT_HI1.
      !
      !  Modified:
      !
      !    05 July 2006
      !
      !  Author:
      !
      !    John Burkardt
      !
        implicit none

        integer ( kind = 4 ) i
        integer ( kind = 4 ) i4_bit_hi1
        integer ( kind = 4 ) i4_uniform
        integer ( kind = 4 ) j
        integer ( kind = 4 ) seed
        integer ( kind = 4 ) test

        seed = 123456789

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST02'
        write ( *, '(a)' ) '  I4_BIT_HI1 returns the location of the high 1 bit.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '         I    I4_BIT_HI1(I)'
        write ( *, '(a)' ) ' '

        do test = 1, 10
          i = i4_uniform ( 0, 100, seed )
          j = i4_bit_hi1 ( i )
          write ( *, '(2x,i8,2x,i8)' ) i, j
        end do

        return
      end
      subroutine test03

      !*****************************************************************************80
      !
      !! TEST03 tests I4_BIT_LO0.
      !
      !  Modified:
      !
      !    05 July 2006
      !
      !  Author:
      !
      !    John Burkardt
      !
        implicit none

        integer ( kind = 4 ) i
        integer ( kind = 4 ) i4_bit_lo0
        integer ( kind = 4 ) i4_uniform
        integer ( kind = 4 ) j
        integer ( kind = 4 ) seed
        integer ( kind = 4 ) test

        seed = 123456789

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST03'
        write ( *, '(a)' ) '  I4_BIT_LO0 returns the location of the lowest 0 bit.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '         I    I4_BIT_LO0(I)'
        write ( *, '(a)' ) ' '

        do test = 1, 10
          i = i4_uniform ( 0, 100, seed )
          j = i4_bit_lo0 ( i )
          write ( *, '(2x,i8,2x,i8)' ) i, j
        end do

        return
      end
      subroutine test04

      !*****************************************************************************80
      !
      !! TEST04 tests I4_SOBOL.
      !
      !  Discussion:
      !
      !    This routine uses the default integer precision, which is
      !    presumed to correspond to a KIND of 4.
      !
      !  Modified:
      !
      !    25 June 2008
      !
      !  Author:
      !
      !    John Burkardt
      !
        implicit none

        integer ( kind = 4 ), parameter :: dim_max = 4

        integer ( kind = 4 ) dim_num
        integer ( kind = 4 ) i
        integer ( kind = 4 ) j
        real    ( kind = 4 ) r(dim_max)
        integer ( kind = 4 ) seed
        integer ( kind = 4 ) seed_in
        integer ( kind = 4 ) seed_out

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST04'
        write ( *, '(a)' ) '  I4_SOBOL computes the next element of a Sobol sequence.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  In this test, we call I4_SOBOL repeatedly.'

        do dim_num = 2, dim_max

          seed = 0

          write ( *, '(a)' ) ' '
          write ( *, '(a,i8)' ) '  Using dimension DIM_NUM =   ', dim_num
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '      Seed      Seed      I4_Sobol'
          write ( *, '(a)' ) '        In       Out'
          write ( *, '(a)' ) ' '
          do i = 0, 110
            seed_in = seed
            call i4_sobol ( dim_num, seed, r )
            seed_out = seed
            if ( i <= 11 .or. 95 <= i ) then
              write ( *, '(2x,i8,2x,i8,2x,4f12.6)' ) seed_in, seed_out, r(1:dim_num)
            else if ( i == 12 ) then
              write ( *, '(a)' ) '......................'
            end if
          end do

        end do

        dim_num = 2

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Repeat the 2D calculation, but start with different seeds.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Using dimension DIM_NUM =   ', dim_num

        do j = 0, 10

          seed = j

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '      Seed      Seed      I4_Sobol'
          write ( *, '(a)' ) '        In       Out'
          write ( *, '(a)' ) ' '
          do i = 0, 110
            seed_in = seed
            call i4_sobol ( dim_num, seed, r )
            seed_out = seed
            if ( i <= 11 .or. 95 <= i ) then
              write ( *, '(2x,i8,2x,i8,2x,4f12.6)' ) seed_in, seed_out, r(1:dim_num)
            else if ( i == 12 ) then
              write ( *, '(a)' ) '......................'
            end if
          end do

        end do

        return
      end
      subroutine test05

      !*****************************************************************************80
      !
      !! TEST05 tests I4_SOBOL.
      !
      !  Discussion:
      !
      !    This routine uses the default integer precision, which is
      !    presumed to correspond to a KIND of 4.
      !
      !  Modified:
      !
      !    28 June 2006
      !
      !  Author:
      !
      !    John Burkardt
      !
        implicit none

        integer ( kind = 4 ), parameter :: dim_num = 3

        integer ( kind = 4 ) i
        real    ( kind = 4 ) r(dim_num)
        integer ( kind = 4 ) seed
        integer ( kind = 4 ) seed_in
        integer ( kind = 4 ) seed_out

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST05'
        write ( *, '(a)' ) '  I4_SOBOL computes the next element of a Sobol sequence.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  In this test, we demonstrate how the SEED can be'
        write ( *, '(a)' ) '  manipulated to skip ahead in the sequence, or'
        write ( *, '(a)' ) '  to come back to any part of the sequence.'

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Using dimension DIM_NUM =   ', dim_num

        seed = 0

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '      Seed      Seed      I4_Sobol'
        write ( *, '(a)' ) '        In       Out'
        write ( *, '(a)' ) ' '
        do i = 0, 10
          seed_in = seed
          call i4_sobol ( dim_num, seed, r )
          seed_out = seed
          write ( *, '(2x,i8,2x,i8,2x,4f12.6)' ) seed_in, seed_out, r(1:dim_num)
        end do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Jump ahead by increasing SEED:'
        write ( *, '(a)' ) ' '

        seed = 100

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '      Seed      Seed      I4_Sobol'
        write ( *, '(a)' ) '        In       Out'
        write ( *, '(a)' ) ' '
        do i = 1, 5
          seed_in = seed
          call i4_sobol ( dim_num, seed, r )
          seed_out = seed
          write ( *, '(2x,i8,2x,i8,2x,4f12.6)' ) seed_in, seed_out, r(1:dim_num)
        end do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Jump back by decreasing SEED:'
        write ( *, '(a)' ) ' '

        seed = 3

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '      Seed      Seed      I4_Sobol'
        write ( *, '(a)' ) '        In       Out'
        write ( *, '(a)' ) ' '
        do i = 0, 10
          seed_in = seed
          call i4_sobol ( dim_num, seed, r )
          seed_out = seed
          write ( *, '(2x,i8,2x,i8,2x,4f12.6)' ) seed_in, seed_out, r(1:dim_num)
        end do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Jump ahead by increasing SEED:'
        write ( *, '(a)' ) ' '

        seed = 98

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '      Seed      Seed      I4_Sobol'
        write ( *, '(a)' ) '        In       Out'
        write ( *, '(a)' ) ' '
        do i = 1, 5
          seed_in = seed
          call i4_sobol ( dim_num, seed, r )
          seed_out = seed
          write ( *, '(2x,i8,2x,i8,2x,4f12.6)' ) seed_in, seed_out, r(1:dim_num)
        end do

        return
      end
      subroutine test06

      !*****************************************************************************80
      !
      !! TEST06 tests I4_SOBOL.
      !
      !  Discussion:
      !
      !    In this test, we get a few samples at high dimensions.
      !
      !  Modified:
      !
      !    28 June 2006
      !
      !  Author:
      !
      !    John Burkardt
      !
        implicit none

        integer ( kind = 4 ), parameter :: test_num = 3

        integer ( kind = 4 ) dim_num
        integer ( kind = 4 ), dimension ( 3 ) :: dim_num_test = (/ 100, 500, 1000 /)
        integer ( kind = 4 ) i
        real    ( kind = 4 ), allocatable, dimension ( : ) :: r
        integer ( kind = 4 ) seed
        integer ( kind = 4 ) seed_in
        integer ( kind = 4 ) seed_out
        integer ( kind = 4 ) test

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST06'
        write ( *, '(a)' ) '  I4_SOBOL computes the next element of a Sobol sequence.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  In this test, we get a few samples at high dimension.'
        write ( *, '(a)' ) '  We only print the first and last 2 entries of each'
        write ( *, '(a)' ) '  sample.'

        do test = 1, test_num

          dim_num = dim_num_test(test)
          allocate ( r(1:dim_num) )

          write ( *, '(a)' ) ' '
          write ( *, '(a,i8)' ) '  Using dimension DIM_NUM =   ', dim_num

          seed = 0

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '      Seed      Seed      I4_Sobol'
          write ( *, '(a)' ) '        In       Out   (First 2, Last 2)'
          write ( *, '(a)' ) ' '

          do i = 0, 10

            seed_in = seed
            call i4_sobol ( dim_num, seed, r )
            seed_out = seed
            write ( *, '(2x,i8,2x,i8,2x,4f12.6)' ) &
              seed_in, seed_out, r(1:2), r(dim_num-1:dim_num)

          end do

          seed = 100000

          do i = 11, 15

            seed_in = seed
            call i4_sobol ( dim_num, seed, r )
            seed_out = seed
            write ( *, '(2x,i8,2x,i8,2x,4f12.6)' ) &
              seed_in, seed_out, r(1:2), r(dim_num-1:dim_num)

          end do

          deallocate ( r )

        end do

        return
      end
      subroutine test07

      !*****************************************************************************80
      !
      !! TEST07 tests I4_XOR.
      !
      !  Discussion:
      !
      !    The function I4_XOR returns the bitwise exclusive OR 
      !    of two integers.
      !
      !  Modified:
      !
      !    05 July 2006
      !
      !  Author:
      !
      !    John Burkardt
      !
        implicit none

        integer ( kind = 4 ) i
        integer ( kind = 4 ) i4_uniform
        integer ( kind = 4 ) i4_xor
        integer ( kind = 4 ) j
        integer ( kind = 4 ) k
        integer ( kind = 4 ) seed
        integer ( kind = 4 ) test

        seed = 123456789

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST07'
        write ( *, '(a)' ) '  I4_XOR is a FORTRAN90 function which returns the'
        write ( *, '(a)' ) '  bitwise exclusive OR of two integers.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '         I         J    I4_XOR(I,J)'
        write ( *, '(a)' ) ' '

        do test = 1, 10
          i = i4_uniform ( 0, 100, seed )
          j = i4_uniform ( 0, 100, seed )
          k = i4_xor ( i, j )
          write ( *, '(2x,i8,2x,i8,2x,i8)' ) i, j, k
        end do

        return
      end
      subroutine test08

      !*****************************************************************************80
      !
      !! TEST08 tests IEOR on integer ( kind = 8 ) arguments.
      !
      !  Discussion:
      !
      !    The FORTRAN90 function IEOR returns the bitwise exclusive OR 
      !    of two integers.
      !
      !  Modified:
      !
      !    05 July 2006
      !
      !  Author:
      !
      !    John Burkardt
      !
        implicit none

        integer ( kind = 8 ) i
        integer ( kind = 8 ) i8_uniform
        integer ( kind = 8 ) j
        integer ( kind = 8 ) k
        integer ( kind = 8 ) seed
        integer ( kind = 8 ) test

        seed = 123456789

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST08'
        write ( *, '(a)' ) '  IEOR is a FORTRAN90 function which returns the'
        write ( *, '(a)' ) '  bitwise exclusive OR of two integers.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '         I         J    IEOR(I,J)'
        write ( *, '(a)' ) ' '

        do test = 1, 10
          i = i8_uniform ( 0, 100, seed )
          j = i8_uniform ( 0, 100, seed )
          k = ieor ( i, j )
          write ( *, '(2x,i8,2x,i8,2x,i8)' ) i, j, k
        end do

        return
      end
      subroutine test09

      !*****************************************************************************80
      !
      !! TEST09 tests I8_BIT_HI1.
      !
      !  Discussion:
      !
      !    This routine uses integer precision KIND of 8.
      !
      !  Modified:
      !
      !    05 July 2006
      !
      !  Author:
      !
      !    John Burkardt
      !
        implicit none

        integer ( kind = 8 ) i
        integer ( kind = 8 ) i8_bit_hi1
        integer ( kind = 8 ) i8_uniform
        integer ( kind = 8 ) j
        integer ( kind = 8 ) seed
        integer ( kind = 8 ) test

        seed = 123456789

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST09'
        write ( *, '(a)' ) '  I8_BIT_HI1 returns the location of the high 1 bit.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '         I    I8_BIT_HI1(I)'
        write ( *, '(a)' ) ' '

        do test = 1, 10
          i = i8_uniform ( 0, 100, seed )
          j = i8_bit_hi1 ( i )
          write ( *, '(2x,i8,2x,i8)' ) i, j
        end do

        return
      end
      subroutine test10

      !*****************************************************************************80
      !
      !! TEST10 tests I8_BIT_LO0.
      !
      !  Discussion:
      !
      !    This routine uses integer precision KIND of 8.
      !
      !  Modified:
      !
      !    05 July 2006
      !
      !  Author:
      !
      !    John Burkardt
      !
        implicit none

        integer ( kind = 8 ) i
        integer ( kind = 8 ) i8_bit_lo0
        integer ( kind = 8 ) i8_uniform
        integer ( kind = 8 ) j
        integer ( kind = 8 ) seed
        integer ( kind = 8 ) test

        seed = 123456789

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST10'
        write ( *, '(a)' ) '  I8_BIT_LO0 returns the location of the lowest 0 bit.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '         I    I8_BIT_LO0(I)'
        write ( *, '(a)' ) ' '

        do test = 1, 10
          i = i8_uniform ( 0, 100, seed )
          j = i8_bit_lo0 ( i )
          write ( *, '(2x,i8,2x,i8)' ) i, j
        end do

        return
      end
      subroutine test11

      !*****************************************************************************80
      !
      !! TEST11 tests I8_SOBOL.
      !
      !  Discussion:
      !
      !    This routine uses integer precision KIND of 8.
      !
      !  Modified:
      !
      !    28 June 2006
      !
      !  Author:
      !
      !    John Burkardt
      !
        implicit none

        integer ( kind = 8 ), parameter :: dim_max = 4

        integer ( kind = 8 ) dim_num
        integer ( kind = 8 ) i
        real    ( kind = 8 ), dimension ( dim_max ) :: r
        integer ( kind = 8 ) seed
        integer ( kind = 8 ) seed_in
        integer ( kind = 8 ) seed_out

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST11'
        write ( *, '(a)' ) '  I8_SOBOL computes the next element of a Sobol sequence.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  In this test, we call I8_SOBOL repeatedly.'

        do dim_num = 2, dim_max

          seed = 0

          write ( *, '(a)' ) ' '
          write ( *, '(a,i8)' ) '  Using dimension DIM_NUM =   ', dim_num
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '      Seed      Seed      I8_Sobol'
          write ( *, '(a)' ) '        In       Out'
          write ( *, '(a)' ) ' '
          do i = 0, 110
            seed_in = seed
            call i8_sobol ( dim_num, seed, r )
            seed_out = seed
            if ( i <= 11 .or. 95 <= i ) then
              write ( *, '(2x,i8,2x,i8,2x,4f12.6)' ) seed_in, seed_out, r(1:dim_num)
            else if ( i == 12 ) then
              write ( *, '(a)' ) '......................'
            end if
          end do

        end do

        return
      end
      subroutine test12

      !*****************************************************************************80
      !
      !! TEST12 tests I8_SOBOL.
      !
      !  Discussion:
      !
      !    This routine uses integer precision KIND of 8.
      !
      !  Modified:
      !
      !    28 June 2006
      !
      !  Author:
      !
      !    John Burkardt
      !
        implicit none

        integer ( kind = 8 ), parameter :: dim_num = 3

        integer ( kind = 8 ) :: i
        real    ( kind = 8 ), dimension ( dim_num ) :: r
        integer ( kind = 8 ) :: seed
        integer ( kind = 8 ) :: seed_in
        integer ( kind = 8 ) :: seed_out

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST12'
        write ( *, '(a)' ) '  I8_SOBOL computes the next element of a Sobol sequence.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  In this test, we demonstrate how the SEED can be'
        write ( *, '(a)' ) '  manipulated to skip ahead in the sequence, or'
        write ( *, '(a)' ) '  to come back to any part of the sequence.'

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Using dimension DIM_NUM =   ', dim_num

        seed = 0

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '      Seed      Seed      I8_Sobol'
        write ( *, '(a)' ) '        In       Out'
        write ( *, '(a)' ) ' '
        do i = 0, 10
          seed_in = seed
          call i8_sobol ( dim_num, seed, r )
          seed_out = seed
          write ( *, '(2x,i8,2x,i8,2x,4f12.6)' ) seed_in, seed_out, r(1:dim_num)
        end do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Jump ahead by increasing SEED:'
        write ( *, '(a)' ) ' '

        seed = 100

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '      Seed      Seed      I8_Sobol'
        write ( *, '(a)' ) '        In       Out'
        write ( *, '(a)' ) ' '
        do i = 1, 5
          seed_in = seed
          call i8_sobol ( dim_num, seed, r )
          seed_out = seed
          write ( *, '(2x,i8,2x,i8,2x,4f12.6)' ) seed_in, seed_out, r(1:dim_num)
        end do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Jump back by decreasing SEED:'
        write ( *, '(a)' ) ' '

        seed = 3

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '      Seed      Seed      I8_Sobol'
        write ( *, '(a)' ) '        In       Out'
        write ( *, '(a)' ) ' '
        do i = 0, 10
          seed_in = seed
          call i8_sobol ( dim_num, seed, r )
          seed_out = seed
          write ( *, '(2x,i8,2x,i8,2x,4f12.6)' ) seed_in, seed_out, r(1:dim_num)
        end do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Jump ahead by increasing SEED:'
        write ( *, '(a)' ) ' '

        seed = 98

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '      Seed      Seed      I8_Sobol'
        write ( *, '(a)' ) '        In       Out'
        write ( *, '(a)' ) ' '
        do i = 1, 5
          seed_in = seed
          call i8_sobol ( dim_num, seed, r )
          seed_out = seed
          write ( *, '(2x,i8,2x,i8,2x,4f12.6)' ) seed_in, seed_out, r(1:dim_num)
        end do

        return
      end
      subroutine test13

      !*****************************************************************************80
      !
      !! TEST13 tests I8_SOBOL.
      !
      !  Discussion:
      !
      !    In this test, we get a few samples at high dimensions.
      !
      !  Modified:
      !
      !    28 June 2006
      !
      !  Author:
      !
      !    John Burkardt
      !
        implicit none

        integer ( kind = 8 ), parameter :: test_num = 3

        integer ( kind = 8 ) dim_num
        integer ( kind = 8 ), dimension ( 3 ) :: dim_num_test = (/ 100, 500, 1000 /)
        integer ( kind = 8 ) i
        real    ( kind = 8 ), allocatable, dimension ( : ) :: r
        integer ( kind = 8 ) seed
        integer ( kind = 8 ) seed_in
        integer ( kind = 8 ) seed_out
        integer ( kind = 8 ) test

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST13'
        write ( *, '(a)' ) '  I8_SOBOL computes the next element of a Sobol sequence.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  In this test, we get a few samples at high dimension.'
        write ( *, '(a)' ) '  We only print the first and last 2 entries of each'
        write ( *, '(a)' ) '  sample.'

        do test = 1, test_num

          dim_num = dim_num_test(test)
          allocate ( r(1:dim_num) )

          write ( *, '(a)' ) ' '
          write ( *, '(a,i8)' ) '  Using dimension DIM_NUM =   ', dim_num

          seed = 0

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '      Seed      Seed      I8_Sobol'
          write ( *, '(a)' ) '        In       Out   (First 2, Last 2)'
          write ( *, '(a)' ) ' '

          do i = 0, 10

            seed_in = seed
            call i8_sobol ( dim_num, seed, r )
            seed_out = seed
            write ( *, '(2x,i8,2x,i8,2x,4f12.6)' ) &
              seed_in, seed_out, r(1:2), r(dim_num-1:dim_num)

          end do

          seed = 100000

          do i = 11, 15

            seed_in = seed
            call i8_sobol ( dim_num, seed, r )
            seed_out = seed
            write ( *, '(2x,i8,2x,i8,2x,4f12.6)' ) &
              seed_in, seed_out, r(1:2), r(dim_num-1:dim_num)

          end do

          deallocate ( r )

        end do

        return
      end
      subroutine test14

      !*****************************************************************************80
      !
      !! TEST14 tests I8_XOR.
      !
      !  Discussion:
      !
      !    The function I8_XOR returns the bitwise exclusive OR 
      !    of two integers.
      !
      !  Modified:
      !
      !    05 July 2006
      !
      !  Author:
      !
      !    John Burkardt
      !
        implicit none

        integer ( kind = 8 ) i
        integer ( kind = 8 ) i8_uniform
        integer ( kind = 8 ) i8_xor
        integer ( kind = 8 ) j
        integer ( kind = 8 ) k
        integer ( kind = 8 ) seed
        integer ( kind = 8 ) test

        seed = 123456789

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST14'
        write ( *, '(a)' ) '  I8_XOR is a FORTRAN90 function which returns the'
        write ( *, '(a)' ) '  bitwise exclusive OR of two integers.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '         I         J    I8_XOR(I,J)'
        write ( *, '(a)' ) ' '

        do test = 1, 100
          i = i8_uniform ( 0, 1000, seed )
          j = i8_uniform ( 0, 1000, seed )
          k = i8_xor ( i, j )
          write ( *, '(2x,i8,2x,i8,2x,i8)' ) i, j, k
        end do

        return
      end
