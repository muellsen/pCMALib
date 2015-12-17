	  !-------------------------------------------------------------------------
      !  Subroutine   :          cmaes_handlebounds
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This Subroutine implements boundary handling as in
      !					cmaes.m lines 728-795, it 'penalizes' function fitness
      !					with an additional term if sample was out of bounds.
      !
      !  Input		  : m				(I)	rows of arx/arxvalid
      !					n				(I)	columns of arx/arxvalid
      !					arxvalid(:,:)	(R) matrix of valid sample points
      !					arx(:,:)		(R) matrix of all sample points
      !                
      !  Remarks      : The technical details of boundary handling are not
      !					completely understood yet, so this is just a 1:1
      !					translation of Matlab Code (cmaes.m lines 728ff).
      !					bnd.flgscale as it appears in Matlab code, is not imple-
      !					mented as it is once set to 0 and never changed afterw.
      !
      !  References   : 
      !					@article{Hansen:2007,
	  !					Author = {Hansen, Nikolaus},
	  !					Keywords = {CMA},
	  !					Title = {{The CMA Evolution Strategy: A Tutorial}},
	  !					Year = {2007}}
	  !					
      !  Revisions    :
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  pCMALib: a parallel fortran 90 library for the evolution strategy with
      !           covariance matrix adaptation
      !  Christian L. Mueller, Benedikt Baumgartner, Georg Ofenbeck
      !  MOSAIC group, ETH Zurich, Switzerland
      !-------------------------------------------------------------------------
      SUBROUTINE cmaes_handlebounds(arxvalid,arx,m,n)
      !-------------------------------------------------------------------------
	  !  Modules
	  !-------------------------------------------------------------------------
	  USE cmaes_param_mod
      USE cmaes_mod
      USE cmaes_opts_mod
      USE tool_mrgrnk_mod
      IMPLICIT NONE
      !-------------------------------------------------------------------------
	  !  Parameters
	  !-------------------------------------------------------------------------      
      INTEGER,INTENT(in)								:: m,n
      REAL(MK),DIMENSION(m,n),INTENT(in)				:: arxvalid
      REAL(MK),DIMENSION(m,n),INTENT(in)				:: arx
      !-------------------------------------------------------------------------
	  !  Local Variables
	  !-------------------------------------------------------------------------      
      REAL(MK),DIMENSION(2)								:: val
      REAL(MK)											:: value ! temp
      REAL(MK),DIMENSION(input%N)						:: diag
      REAL(MK)											:: meandiag
      INTEGER											:: i,k,maxI
      LOGICAL,DIMENSION(size(bnd%dfithist))				:: mask
      LOGICAL,DIMENSION(input%N)						:: idx
      REAL(MK),DIMENSION(input%N)						:: tx
      INTEGER,DIMENSION(size(bnd%dfithist))				:: dfitidx
      REAL(MK),DIMENSION(size(bnd%dfithist))			:: dfitsort
      REAL(MK),DIMENSION(2)								:: prct
      
	  dfitsort = 0.0_MK
	  dfitidx = 0
	  prct = (/ 25.0_MK,75.0_MK /)
	  maxI = 1

      !-------------------------------------------------------------------------
	  !  Initialize matrix mask
	  !-------------------------------------------------------------------------
      DO i = 1, size(bnd%dfithist)
      	IF(bnd%dfithist(i) .GE. 0.0_MK) THEN
      		mask(i) = .TRUE.
      	ELSE
      		mask(i) = .FALSE.
      	END IF
      END DO
      
      !-------------------------------------------------------------------------
	  !  Get main diagonal of unscaled covariance matrix C
	  !-------------------------------------------------------------------------
	  DO i = 1, input%N
	  	diag(i) = C(i,i)
	  END DO

	  !-------------------------------------------------------------------------
	  !  ...and calculate mean value
	  !-------------------------------------------------------------------------
      meandiag = sum(diag)/real(input%N)

      !-------------------------------------------------------------------------
  	  !  Get delta fitness values
  	  !-------------------------------------------------------------------------
  	  CALL cmaes_myprctile(fitness%raw, lambda, prct, 2, val)
  	  value = (val(2) - val(1)) / real(input%N) / meandiag / sigma**2

  	  !-------------------------------------------------------------------------
  	  !  Catch non-sensible values
  	  !-------------------------------------------------------------------------
  	  IF(value .GE. posInf) THEN
#ifdef __HAVE_MPI__
	  WRITE(*,*) 'Process ', MY_RANK, ' warning: Non-finite fitness range'
#else		  	  
  	  WRITE(*,*) 'Warning: Non-finite fitness range'
#endif  	
  		value = maxval(bnd%dfithist,mask)
  	  ELSE IF(value .EQ. 0.0_MK) THEN
  	  	! Happens if all points are out of bounds or the same
  		value = minval(bnd%dfithist,mask)
  	  ELSE IF(bnd%validfitval .EQ. 0) THEN	! First sensible val
  		bnd%dfithist = -1.0_MK
  		bnd%validfitval = 1
      END IF
      
      !-------------------------------------------------------------------------
  	  !  Store delta fitness values
  	  !-------------------------------------------------------------------------
      bndd:DO i = 1, size(bnd%dfithist)
      	IF(bnd%dfithist(i) .LT. 0.0_MK) THEN	! Store at first unassigned position
      		bnd%dfithist(i) = value
      		maxI = i ! max index of non-negative elements of bnd%dfithist vector
      		EXIT bndd
      	ELSE IF(i .EQ. size(bnd%dfithist)) THEN		! If all positions are
      		DO k = 1, (size(bnd%dfithist)-1)		! assigned, then shift
      			bnd%dfithist(k) = bnd%dfithist(k+1) ! values to the left
      		END DO					! and store new value at last position
      		bnd%dfithist(i) = value
      		maxI = i ! max index of non-negative elements of bnd%dfithist vector
      	END IF
      END DO bndd

      CALL cmaes_xintobounds(options%LBounds,options%Ubounds,input%N,xmean,tx,idx)

      !-------------------------------------------------------------------------
  	  !  Set initial weights
  	  !-------------------------------------------------------------------------
  	  IF(bnd%iniphase) THEN
  	    !WRITE(*,*) 'Initial dfit val: ', value
  	  	IF(any(idx)) THEN
  	  	  IF(maxI .EQ. 1) THEN
  	  	    value = bnd%dfithist(1)
  	  	  ELSE	! Median is saved in 'value'
			!sort bnd%dfithist in ascending order
			CALL mrgrnk(bnd%dfithist(1:maxI),dfitidx)
			DO k = 1, maxI
				dfitsort(k) = bnd%dfithist(dfitidx(k))
			END DO
			!and get median value (if maxI == 2 it is the mean)
			IF(mod(maxI,2) .EQ. 0) THEN
				value = (dfitsort(maxI/2)+dfitsort(maxI/2+1))/2.0_MK
			ELSE
				value = dfitsort((maxI-1)/2+1)
			END IF
	      END IF
	      !idx
	      WHERE(bnd%isbounded)
	      	diag = diag/meandiag
	      	bnd%weights = 2.0002_MK * value / diag
	      END WHERE
	      IF((bnd%validfitval .EQ. 1) .AND. (countIter-lastRestart) .GT. 2) bnd%iniphase = .FALSE.
	    END IF  
      END IF
      
      !-------------------------------------------------------------------------
  	  !  Increase weights
  	  !-------------------------------------------------------------------------
      IF(any(idx)) THEN

      	!-----------------------------------------------------------------------
  	    !  Judge distance of xmean to boundary
  	    !-----------------------------------------------------------------------
      	tx = xmean - tx
      	DO i = 1, input%N
      		idx(i) = ((idx(i)) .AND. (abs(tx(i)) .GT. 3.0_MK* &
      			max(1.0_MK,sqrt(real(input%N))/mueff) * sigma * sqrt(diag(i))))
      		idx(i) = ( idx(i) .AND. &
      			(sign(1.0_MK,tx(i)) .EQ. sign(1.0_MK,(xmean(i)-xold(i)))) )
      	END DO
      	!-----------------------------------------------------------------------
  	    !  Only increase if xmean is moving away
  	    !-----------------------------------------------------------------------
  	    WHERE(idx)
  	    	bnd%weights = 1.2_MK**(max(1.0_MK,mueff/10.0_MK/real(input%N))) &
  	    		* bnd%weights
  	    END WHERE
      END IF

      !-------------------------------------------------------------------------
  	  !  Calculate scaling biased to unity, product is one
  	  !  (omitted since version 2.54 of cmaes.m)
  	  !-------------------------------------------------------------------------
  	  !mean = sum(log(diag))/real(input%N)
      !bnd%scales = exp(0.9*(log(diag)) - mean)

      !-------------------------------------------------------------------------
  	  !  Assigned penalized fitness
  	  !-------------------------------------------------------------------------

      bnd%arpenalty = matmul( (bnd%weights / bnd%scales), (arxvalid - arx)**2 )
	  fitness%sel = fitness%raw + bnd%arpenalty

	  RETURN
	  END SUBROUTINE cmaes_handlebounds      
	  
	  
	  
	  
	  
	  
