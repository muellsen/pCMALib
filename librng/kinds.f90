       	  !-------------------------------------------------------------------------
      !  Module       :                   kinds
      !-------------------------------------------------------------------------
      MODULE kinds
      IMPLICIT NONE
          !===========================================================================
          !
          !  Kind type parameters
          !
          !  It is expected that we will have at least the following:
          !
          !    kr4 -> 4 byte real
          !    kr8 -> 8 byte real
          !    ki4 -> 4 byte integer
          !    ki8 -> 8 byte integer
          !
          !===========================================================================
          integer, parameter :: kr4 = selected_real_kind(6, 37)
          integer, parameter :: kr8 = selected_real_kind(9, 99)

          integer, parameter :: ki4 = selected_int_kind(5)
          integer, parameter :: ki8 = selected_int_kind(10) 

      END MODULE kinds
      
