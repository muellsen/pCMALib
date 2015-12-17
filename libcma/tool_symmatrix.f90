	  !-------------------------------------------------------------------------
      !  Subroutine   :          tool_symmatrix
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This Subroutine constructs a symmetric matrix based on
      !					the upper triangular part of the square input Matrix C.
      !					It returns the new symmetric matrix C and also its upper      
      !					triangular part triuC
      !
      !  Input		  : n			(I)	dimension of C, C must be square
	  !
	  !  Input/Output : C(:,:)		(R) Symmetric matrix
	  !					triuC(:,:)	(R) Upper triangular part of sym. matrix
      !                
      !  Remarks      : 
      !
      !  References   : 
	  !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  pCMALib: a parallel fortran 90 library for the evolution strategy with
      !           covariance matrix adaptation
      !  Christian L. Mueller, Benedikt Baumgartner, Georg Ofenbeck
      !  MOSAIC group, ETH Zurich, Switzerland
      !-------------------------------------------------------------------------
#if(__KIND == __SINGLE_PRECISION)
      SUBROUTINE tool_symmatrix_s(C,n,triuC)
#elif(__KIND == __DOUBLE_PRECISION)
      SUBROUTINE tool_symmatrix_d(C,n,triuC)
#endif
      IMPLICIT NONE
      
#if(__KIND == __SINGLE_PRECISION)
	  INTEGER,PARAMETER	:: MK = KIND(1.E0)
#else
	  INTEGER,PARAMETER	:: MK = KIND(1.D0)
#endif       
      !-------------------------------------------------------------------------
      !  Parameters
      !-------------------------------------------------------------------------
      INTEGER, INTENT(in)							:: n
      REAL(MK),DIMENSION(n,n),INTENT(inout)			:: C
      REAL(MK),DIMENSION(n,n),INTENT(inout)			:: triuC

      !-------------------------------------------------------------------------
      !  Local Variables
      !-------------------------------------------------------------------------
      REAL(MK),DIMENSION(n,n)						:: triuC1,triuC2	
      LOGICAL,DIMENSION(n,n)						:: mask
      INTEGER										:: i,k
      
      	!-----------------------------------------------------------------------
        !  Enforce symmetry
      	!-----------------------------------------------------------------------
	  	triuC1 = 0.
	  	triuC2 = 0.
	  	mask = .FALSE.
	  	DO k = 1, n			! Create mask for upper triangular matrix
	  		DO i = 1, n
	  			IF(k .GE. i) THEN
	  				mask(i,k) = .TRUE.
	  			ELSE
	  				mask(i,k) = .FALSE.
	  			END IF
	  		END DO
	  	END DO

	  	WHERE(mask)
	  		triuC1 = C			! Get upper triangular matrix of C
	  	END WHERE				! See Matlab: triu(X)
	  	
	  	DO i = 1, n		! set new mask, see triu(X,k), here: k=1
	  		mask(i,i) = .FALSE.
	  	END DO
	  	
	  	WHERE(mask)
	  		triuC2 = C			! Get second triangular matrix
	  	END WHERE
	  	
	  	C = triuC1 + transpose(triuC2)	! Enforce symmetry
        triuC = triuC1


      RETURN
      END SUBROUTINE