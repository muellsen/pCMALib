	  !-------------------------------------------------------------------------
      !  Function     :                  tool_myrange
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Computes range of values of two different vectors. The
      !					function returns maxVal(vectors) - minVal(vectors).
      !					Optionally, a real number can be given which should be
      !					excluded for computations.
      !
      !  Example	  :	x() = [0,1,2,0,4]	y() = [5,8,3]	exc = 0
      !					tool_myrange(x,5,exc,y,3) = 7	(8-1)
      !                
      !	 Input		  : x(:)		(R) Vector to compute range for
      !					n			(I) Dimension of x
      !					y(:)		(R) Second Vector for range computation
      !					m			(I) Dimension of y
      !					exc			(R)	Exceptional value (not considered for
      !									range computation: 0,Inf,NaN, etc.)
      !
      !	 Output		  : tool_myrange		(R) Overall range of x !and! y
      !
      !  Remarks      : 
      !
      !  References   :	this function is a Fortran translation of tool_myrange in
      !					cmaes.m, lines 1613f.
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  pCMALib: a parallel fortran 90 library for the evolution strategy with
      !           covariance matrix adaptation
      !  Christian L. Mueller, Benedikt Baumgartner, Georg Ofenbeck
      !  MOSAIC group, ETH Zurich, Switzerland
      !-------------------------------------------------------------------------
      FUNCTION tool_myrange(x,n,y,m,exc)
      USE cmaes_param_mod
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Parameters
      !-------------------------------------------------------------------------
      real (mk) :: tool_myrange
      INTEGER,INTENT(in)						:: n
      REAL(MK),DIMENSION(n),INTENT(in)			:: x
      INTEGER,INTENT(in)						:: m
      REAL(MK),DIMENSION(m),INTENT(in)			:: y
      REAL(MK),INTENT(in),OPTIONAL				:: exc
      
      !-------------------------------------------------------------------------
      !  Local Variables
      !-------------------------------------------------------------------------
      REAL(MK),DIMENSION(n+m)					:: work
	  LOGICAL,DIMENSION(n+m)					:: mask
	  INTEGER									:: i
	  
	  !-------------------------------------------------------------------------
	  !  Concatenate x and y
	  !-------------------------------------------------------------------------
	  DO i = 1, n
	  	work(i) = x(i)
	  END DO
	  DO i = n+1, n+m
	  	work(i) = y(i-n)
	  END DO
	  
	  !-------------------------------------------------------------------------
	  !  Check for valid values
	  !-------------------------------------------------------------------------
	  mask = .TRUE.
	  IF(present(exc)) THEN
	  	DO i = 1, n+m
	  		IF(work(i) .EQ. exc) mask(i) = .FALSE.
	  	END DO
	  END IF
	  
	  !-------------------------------------------------------------------------
	  !  Compute range
	  !-------------------------------------------------------------------------
	  tool_myrange = maxval(work,mask) - minval(work,mask)
      
      RETURN
      END FUNCTION tool_myrange
