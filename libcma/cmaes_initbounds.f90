	  !-------------------------------------------------------------------------
      !  Subroutine   :          cmaes_initbounds
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This Subroutine checks if bounds were set by the user 
      !					and if sigma was chosen reasonable 
      !					with respect to given boundaries and possibly
      !					adapts sigma
      !                
      !  Remarks      : TODO: maybe we need more bnd parameters as in cmaes.m
      !					lines 521ff.
      !
      !  References   : cmaes_xintobounds (in file cmaes_initbounds.f90 at the bottom)
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  pCMALib: a parallel fortran 90 library for the evolution strategy with
      !           covariance matrix adaptation
      !  Christian L. Mueller, Benedikt Baumgartner, Georg Ofenbeck
      !  MOSAIC group, ETH Zurich, Switzerland
      !-------------------------------------------------------------------------
      SUBROUTINE cmaes_initbounds()
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE cmaes_param_mod
      USE cmaes_mod
      USE cmaes_opts_mod
      IMPLICIT NONE
      
      
      !-------------------------------------------------------------------------
      !  Interfaces
      !-------------------------------------------------------------------------
      INTERFACE
          SUBROUTINE cmaes_xintobounds(LBounds,UBounds,N,x,xout,idx)
          USE cmaes_param_mod
          USE cmaes_mod,only:countOutOfBounds
          IMPLICIT NONE
          INTEGER, INTENT(in)			  					:: N
          REAL(MK),DIMENSION(N),INTENT(in)    				:: LBounds,UBounds
          REAL(MK),DIMENSION(N),INTENT(in) 					:: x
          REAL(MK),DIMENSION(N),INTENT(out)					:: xout
          LOGICAL,DIMENSION(N),INTENT(out),OPTIONAL			:: idx
          END SUBROUTINE
      END INTERFACE
      !-------------------------------------------------------------------------
      !  Local Variables
      !-------------------------------------------------------------------------
      INTEGER							:: i,allocStat
      LOGICAL,DIMENSION(input%N) 		:: idx
      REAL(MK),DIMENSION(input%N)		:: dd
      LOGICAL							:: flgsigma
      CHARACTER(len=10)					:: proc
      
      !-------------------------------------------------------------------------
      !  Initialize Local Variables
      !-------------------------------------------------------------------------
      flgsigma = .FALSE.
#ifdef __HAVE_MPI__
	  WRITE(proc,'(I10)') MY_RANK
#else
	  proc = ''
#endif	       
      
      !-------------------------------------------------------------------------
      !  Check if bounds have right dimension
      !-------------------------------------------------------------------------
      IF((size(options%LBounds) .NE. input%N) .OR. &
      (size(options%UBounds) .NE. input%N)) &
      	STOP 'Bounds have wrong dimensions'
      
      !-------------------------------------------------------------------------
      !  Allocate bnd%isbounded
      !-------------------------------------------------------------------------
      IF (.NOT. ALLOCATED(bnd%isbounded)) THEN
          ALLOCATE(bnd%isbounded(input%N),stat=allocStat)
          IF(allocStat .NE. 0) STOP 'Error allocating bnd.isbounded'
      END IF
      bnd%isbounded = .FALSE.   
      
      !-------------------------------------------------------------------------
      !  Check if bounds were set by the user
      !-------------------------------------------------------------------------
      activeBnd: DO i = 1, input%N
      	IF((options%LBounds(i) .GT. negInf) .OR. &
      	(options%UBounds(i) .LT. posInf)) THEN
      		bnd%isactive = .TRUE.
      		bnd%isbounded(i) = .TRUE.
      	END IF
      END DO activeBnd
	  
	  !-------------------------------------------------------------------------
      ! If bounds were set explicitly...
      !-------------------------------------------------------------------------
      IF(bnd%isactive) THEN
      	!-----------------------------------------------------------------------
      	!  Check if lower bound smaller than upper bound
      	!-----------------------------------------------------------------------
      	DO i = 1, input%N
	  		IF(options%LBounds(i) .GT. options%UBounds(i)) &
	  		STOP 'Lower bound must be smaller than upper bound'
	  	END DO
	  	
	  	!-----------------------------------------------------------------------
	  	! Check if xmean within bounds, correct if not
	  	!-----------------------------------------------------------------------
      	CALL cmaes_xintobounds(options%LBounds,options%UBounds,input%N,&
      														xmean,xmean,idx)
      	findChg: DO i = 1, input%N
      				IF(idx(i)) THEN
#ifdef __HAVE_MPI__
		WRITE(*,*) 'Process ',trim(adjustl(proc)),':',&
#else
		WRITE(*,*) &		
#endif		
      					'Initial point was out of bounds, corrected'
      					EXIT findChg
      				END IF
      			 END DO findChg
      	!-----------------------------------------------------------------------
      	! Allocate bnd%weights
      	!-----------------------------------------------------------------------
      	IF (.NOT. ALLOCATED(bnd%weights)) THEN
      	    ALLOCATE(bnd%weights(input%N),stat=allocStat)
      	    IF(allocStat .NE. 0) STOP 'Error allocating bnd.weights'
      	END IF
      	bnd%weights = 0.0_MK
      	!-----------------------------------------------------------------------
      	! Allocate bnd%arpenalty
      	!-----------------------------------------------------------------------
      	IF (ALLOCATED(bnd%arpenalty)) DEALLOCATE(bnd%arpenalty)
  	    ALLOCATE(bnd%arpenalty(lambda),stat=allocStat)
  	    IF(allocStat .NE. 0) STOP 'Error allocating bnd.arpenalty'
  
      	bnd%arpenalty = 0.0_MK
      	
      	!-----------------------------------------------------------------------
      	! Allocate scaling vector
      	!-----------------------------------------------------------------------
      	IF (.NOT. ALLOCATED(bnd%scales)) THEN
      	    ALLOCATE(bnd%scales(input%N),stat=allocStat)
      	    IF(allocStat .NE. 0) STOP 'Error allocating bnd.scales'
      	END IF
      	bnd%scales = 1.0_MK
      	!-----------------------------------------------------------------------
  		! Assign max. variable change(s) and get main diagonal of C as we are
  		! already looping...
  		!-----------------------------------------------------------------------
      	DO i = 1, input%N
			maxdx(i) = min(maxdx(i),(options%UBounds(i)-options%LBounds(i))/2.)
			dd(i) = C(i,i)
		END DO
		
		!-----------------------------------------------------------------------
		! Check if sigma is too large
		!-----------------------------------------------------------------------
		fact: DO i = 1, input%N
			IF(sigma*sqrt(dd(i)) .GT. maxdx(i)) THEN
				flgsigma = .TRUE.
				EXIT fact
			END IF
		END DO fact
		!-----------------------------------------------------------------------
		! Rescale sigma if necessary
		!-----------------------------------------------------------------------
		IF(flgsigma) THEN
			fac = minval(maxdx/sqrt(dd))/sigma
			sigma = sigma*fac
#ifdef __HAVE_MPI__
		WRITE(*,*) 'Process ',trim(adjustl(proc)),':',&
#else		
		WRITE(*,*) &
#endif				
			'Warning: Initial sigma multiplied by the factor ',fac, &
			' because it was larger than half of one of the boundary intervals'
		END IF
		!-----------------------------------------------------------------------
		! Set idx (indicates if bounds are real numbers)
		!-----------------------------------------------------------------------
		idx = .FALSE.
		DO i = 1, input%N	
			IF((options%LBounds(i) .GT. negInf) .AND. &
				(options%UBounds(i) .LT. posInf)) idx(i) = .TRUE.
		END DO
		!-----------------------------------------------------------------------
		! Check if sigma is too small
		!-----------------------------------------------------------------------
		scal:DO i = 1, input%N
			IF(idx(i)) THEN
				IF(5*sigma*sqrt(dd(i)) .LT. &
					(options%UBounds(i) - options%LBounds(i))) THEN
#ifdef __HAVE_MPI__
		WRITE(*,*) 'Process ',trim(adjustl(proc)),':',&
#else	
		WRITE(*,*)	&
#endif						
				'Initial SIGMA is in at least one coordinate',&
	       		' much smaller than the given boundary intervals.',&
	       		' For reasonable global search performance SIGMA should be',&
	       		' between 0.2 and 0.5 of the bounded interval in each ',&
	       		'coordinate.'
	       		WRITE(*,*)
					EXIT scal
				END IF
			END IF
		END DO scal
		
		!-----------------------------------------------------------------------
		!  Set more boundary variables
		!-----------------------------------------------------------------------
      	IF (ALLOCATED(bnd%dfithist)) DEALLOCATE(bnd%dfithist)
  	    ALLOCATE(bnd%dfithist(20+(3*input%N)/lambda),stat=allocStat)
  	    IF(allocStat .NE. 0) STOP 'Error allocating bnd.dfithist'

      	bnd%dfithist = -1.0_MK
		bnd%dfithist(1) = 1.0_MK
		
		bnd%iniphase = .TRUE.
		bnd%validfitval = 0

      END IF	! bnd%isactive

      
      RETURN
      END SUBROUTINE
      
      
     
