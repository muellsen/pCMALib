	  !-------------------------------------------------------------------------
      !  Function     :                  cmaes_myprctile
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This function computes the percentiles in vector perc
      !					from vector inar
      !                
      !	 Input		  : inar(:)		(R) Vector containing data to compute
      !									percentiles from
      !					perc(:)		(R)	Vector containing percentiles to compute
      !					
      !					N,rN		(I) Dimension of vecor inar(N) and perc(rN)
      !
      !	 Output		  : res(:)		(R) Vector containing computed percentiles
      !
      !  Remarks      : 
      !
      !  References   :	this function is a Fortran translation of cmaes_myprctile in
      !					cmaes.m, lines 1618ff.
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  pCMALib: a parallel fortran 90 library for the evolution strategy with
      !           covariance matrix adaptation
      !  Christian L. Mueller, Benedikt Baumgartner, Georg Ofenbeck
      !  MOSAIC group, ETH Zurich, Switzerland
      !-------------------------------------------------------------------------
  	  SUBROUTINE cmaes_myprctile(inar, N, perc, rN, res)
  	  USE cmaes_param_mod
  	  USE tool_mrgrnk_mod
  	  IMPLICIT NONE
  	  !-------------------------------------------------------------------------
  	  !  Arguments
  	  !-------------------------------------------------------------------------
  	  INTEGER,INTENT(in)					:: N
  	  INTEGER,INTENT(in)					:: rN
  	  REAL(MK),INTENT(in)					:: inar(N)
  	  REAL(MK),INTENT(in)					:: perc(rN)
  	  REAL(MK),DIMENSION(rN),INTENT(out) 	:: res
  	  
  	  !-------------------------------------------------------------------------
  	  !  Local Variables
  	  !-------------------------------------------------------------------------
  	  REAL(MK),DIMENSION(N)	:: sar, availablepercentiles
  	  INTEGER,DIMENSION(N)	:: idx	! Used to call SORTC correctly (l.52)
  	  INTEGER 				:: i,k

	  availablepercentiles = 0.0_MK  	  
  	  !-------------------------------------------------------------------------
	  ! Sort inar in ascending order
	  !-------------------------------------------------------------------------
  	  CALL mrgrnk(inar,idx)
  	  DO i = 1, N
  	  	sar(i) = inar(idx(i))
  	  END DO
  	  DO i = 1, rN
  	  	IF(perc(i) .LE. (100.0_MK*0.5_MK/real(N))) THEN
  	  		res(i) = sar(1)
  	  	ELSE IF(perc(i) .GE. (100.0_MK*(real(N)-0.5_MK)/real(N)) ) THEN
  	  		res(i) = sar(N)
  	  	ELSE
  	  		!-------------------------------------------------------------------
  	  		! Find largest index smaller than required percentile
  	  		!-------------------------------------------------------------------
			DO k = 1,N
				availablepercentiles(k) = 100.0 * (real(k)-0.5) / real(N)
			END DO 
  	  		find:DO k = 1, N
  	  			IF(availablepercentiles(k) .GE. perc(i)) EXIT find
  	  		END DO find
  	  		k = k - 1
  	  		!-------------------------------------------------------------------
  	  		! Interpolate linearly
  	  		!-------------------------------------------------------------------
  	  		res(i) = sar(k) + (sar(k+1)-sar(k)) * (perc(i) &
  	  			-availablepercentiles(k)) / (availablepercentiles(k+1) - &
  	  			availablepercentiles(k))
  	  	END IF
  	  END DO

	  RETURN
  	  END SUBROUTINE cmaes_myprctile
  	  