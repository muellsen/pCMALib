	  !-------------------------------------------------------------------------
      !  Subroutine   :          cmaes_checkparams
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This Subroutine checks the consistency of
      !					CMA parameters and variables
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
      SUBROUTINE cmaes_checkparams()
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE cmaes_mod
      USE cmaes_opts_mod
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Local Variables
      !-------------------------------------------------------------------------
      INTEGER	:: i
      
      !-------------------------------------------------------------------------
      !  Check lambda and mu
      !-------------------------------------------------------------------------
	  IF(lambda .LE. 1) STOP 'Lambda must be greater than 1.'
		
	  IF(mu .LT. 1) &
		STOP 'parent number mu must be greater or equal to 1.'
		
	  IF(lambda .LE. mu) &
        STOP 'Lambda must be greater than mu'
      !-------------------------------------------------------------------------	    
      !  Check sigma
      !-------------------------------------------------------------------------
      DO i = 1,input%N
      	IF(input%insigma(i) .LE. 0) &
      		STOP 'Initial search volume (SIGMA) must be greater than zero'
      END DO
      IF(maxval(input%insigma)/minval(input%insigma) .GT. 1.E6) &
      	STOP 'Initial search volume (SIGMA) badly conditioned'
      !-------------------------------------------------------------------------
      !  Check boundaries	
      !-------------------------------------------------------------------------
      IF((size(options%LBounds) .NE. input%N) .OR. &
      	(size(options%UBounds) .NE. input%N)) &
      	STOP 'Bounds have wrong dimensions'
      	
	  DO i = 1, input%N
	  	IF(options%LBounds(i) .GT. options%UBounds(i)) &
	  		STOP 'Lower bound must be smaller than upper bound'
	  END DO
      !-------------------------------------------------------------------------
      !  Check xstart
      !-------------------------------------------------------------------------
      IF(size(input%xstart) .NE. input%N) STOP 'xstart has wrong dimension'
      DO i = 1, input%N
      	IF((input%xstart(i) .LT. options%LBounds(i)) .OR. &
      	(input%xstart(i) .GT. options%UBounds(i))) & 
      		STOP 'xstart not within boundaries'
      END DO

      	
      RETURN
      END SUBROUTINE cmaes_checkparams