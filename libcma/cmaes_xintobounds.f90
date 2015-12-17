      !-------------------------------------------------------------------------
      !  Subroutine   :          cmaes_xintobounds
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This Subroutine makes sure that a given vector is within
      !					given boundaries (LBounds/UBounds) and sets xout
      !					to LBound/UBound if this is not the case
      !
      !  Input		  : LBounds(:)		(R) Lower bound
      !					UBounds(:)		(R) Upper bound               
      !					N				(I) Dimension of LBounds/NBounds
      !					x(:)			(R) Vector to place within bounds dim(N)
      !
      !  Output		  : xout(:)			(R) Output vector, dim(N)
      !					idx(:)			(L) Index vector, idx(i) == .TRUE. ,
      !										if x(i) was changed, OPTIONAL
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
      SUBROUTINE cmaes_xintobounds(LBounds,UBounds,N,x,xout,idx)
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE cmaes_param_mod
      USE cmaes_mod,only:countOutOfBounds
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Parameters
      !-------------------------------------------------------------------------
	  INTEGER, INTENT(in)			  					:: N
	  REAL(MK),DIMENSION(N),INTENT(in)    				:: LBounds,UBounds
	  REAL(MK),DIMENSION(N),INTENT(in) 					:: x
	  REAL(MK),DIMENSION(N),INTENT(out)					:: xout
	  LOGICAL,DIMENSION(N),INTENT(out),OPTIONAL			:: idx
	  
	  !-------------------------------------------------------------------------
	  !  Local variables
	  !-------------------------------------------------------------------------
	  INTEGER											:: i
	  !-------------------------------------------------------------------------
	  !  Set idx == true if x out of bounds and correct x
	  !-------------------------------------------------------------------------
	  IF(present(idx)) THEN
		  DO i = 1, N
		  	IF(x(i) .LT. LBounds(i)) THEN
		  		xout(i) = LBounds(i)
		  		idx(i) = .TRUE.
		  		countOutOfBounds = countOutOfBounds + 1
		  	ELSE IF(x(i) .GT. UBounds(i)) THEN
		  		xout(i) = UBounds(i)
		  		idx(i) = .TRUE.
		  		countOutOfBounds = countOutOfBounds + 1
		  	ELSE
		  		xout(i) = x(i)
		  		idx(i) = .FALSE.
		  	END IF
		  END DO
	  ELSE
	  	  DO i = 1, N
		  	IF(x(i) .LT. LBounds(i)) THEN
		  		xout(i) = LBounds(i)
		  		countOutOfBounds = countOutOfBounds + 1
		  	ELSE IF(x(i) .GT. UBounds(i)) THEN
		  		xout(i) = UBounds(i)
		  		countOutOfBounds = countOutOfBounds + 1
		  	ELSE 
		  		xout(i) = x(i)
		  	END IF
		  END DO
	  END IF
	  !IF any(idx) WRITE(*,*) 'Returned from xintobounds'
	  RETURN      
      END SUBROUTINE
