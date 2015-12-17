 	  !-------------------------------------------------------------------------
      !  Subroutine   :          cmaes_init
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This Subroutine initializes the CMA-ES
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
      SUBROUTINE cmaes_init
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE cmaes_param_mod
      USE cmaes_mod
      USE cmaes_opts_mod
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER							:: i
      INTEGER							:: allocStat
      REAL(MK)							:: temp = 0.
      REAL(MK)							:: vecnorm,vecsum
      INTEGER							:: N
      REAL(MK),DIMENSION(input%N,input%N)           :: BinvD
      
      N = input%N
      !-------------------------------------------------------------------------
      !  Allocate and calculate xmean
      !-------------------------------------------------------------------------
      IF (ALLOCATED(xmean)) THEN
         xmean = input%xstart
      ELSE
          ALLOCATE(xmean(N),stat=allocStat)
          IF(allocStat .NE. 0) STOP 'Error allocating xmean'
          xmean = input%xstart
      END IF

      
      
      IF (.NOT. ALLOCATED(xold)) THEN
          ALLOCATE(xold(N),stat=allocStat)
          IF(allocStat .NE. 0) STOP 'Error allocating xold'
      END IF
      xold = 0.0_MK
      
      !-------------------------------------------------------------------------
      !  Strategy internal parameter setting: Selection
      !-------------------------------------------------------------------------
      !  Lambda
      IF(options%PopSize .EQ. 0) THEN
      	lambda = 4+int(3.*log(REAL(N)))
      	options%PopSize = lambda
      ELSE
      	lambda = options%PopSize
      END IF
      IF (initial_Popsize .EQ. 0) THEN
        initial_Popsize = options%Popsize
      END IF
      
      !  mu
      IF(options%ParentNumber .EQ. 0) THEN
      	mu = int(lambda/2)
      	!options%ParentNumber = mu
      ELSE
      	mu = options%ParentNumber
      END IF
      
      !  weights
      IF (ALLOCATED(weights)) DEALLOCATE(weights)
  	  ALLOCATE(weights(mu),stat=allocStat)
  	  IF(allocStat .NE. 0) STOP 'Error allocating weights'
  	  ! Weights equal
  	  IF(options%RecombinationWeights .EQ. 1) THEN
  			weights = 1.
  		! Weights linear
  		ELSE IF(options%RecombinationWeights .EQ. 2) THEN
  			DO i = 1,mu
  				weights(i) = mu + 1 - real(i)
  			END DO
  		! Weights superlinear
  	    ELSE IF(options%RecombinationWeights .EQ. 3) THEN
  			DO i = 1,mu
  				weights(i) = log(REAL(mu)+1.)-log(REAL(i))
  			END DO
  		ELSE 
  			STOP 'RecombinationWeights not set correctly (see cmaes_opts_mod)'
  	  END IF
       
      !  mueff
      vecsum = sum(weights)
      vecnorm = sqrt(sum(weights**2))
      mueff = vecsum**2 / vecnorm**2
      
      
      !-------------------------------------------------------------------------
      !  Evaluate options
      !-------------------------------------------------------------------------
      IF(options%StopMaxIter .EQ. 0.) &
      	options%StopMaxIter = int(1.E3 * (N+5)**2 / sqrt(real(lambda)))

      IF(options%StopTolX .EQ. 0.) &
      	options%StopTolX = 1.E-11 * maxval(input%insigma)

      IF(options%StopTolUpX .EQ. 0.) &
       	options%StopTolUpX = 1.E3 * maxval(input%insigma)


	  IF(.NOT. ASSOCIATED(options%DiffMaxChange)) THEN
      	ALLOCATE(options%DiffMaxChange(N),stat=allocStat)
      	IF(allocStat .NE. 0) STOP 'Error allocating DiffMaxChange'
      	options%DiffMaxChange = posInf
      END IF
	  
	  IF(.NOT. ASSOCIATED(options%DiffMinChange)) THEN
      	ALLOCATE(options%DiffMinChange(N),stat=allocStat)
      	IF(allocStat .NE. 0) STOP 'Error allocating DiffMinChange'
      	options%DiffMinChange = 0
      END IF

	  IF(.NOT. ASSOCIATED(options%LBounds)) THEN 
      	ALLOCATE(options%LBounds(N),stat=allocStat)
      	IF(allocStat .NE. 0) STOP 'Error allocating LBounds'
      	options%LBounds = negInf
      END IF
      
	  IF(.NOT. ASSOCIATED(options%UBounds)) THEN      
      	ALLOCATE(options%UBounds(N),stat=allocStat)
      	IF(allocStat .NE. 0) STOP 'Error allocating UBounds'
      	options%UBounds = posInf
      END IF
    
       
      !-------------------------------------------------------------------------
      !  Strategy internal parameter setting: Adaptation
      !-------------------------------------------------------------------------
      cc = 4.0_MK/(REAL(N)+4.0_MK)
      cs = REAL((mueff+2))/REAL((N+mueff+3))
      mucov = mueff
      ccov = (1._MK/mucov) * 2._MK/(REAL(N)+1.41)**2 + (1._MK - 1._MK/mucov) * &
      		min(1.0_MK,(2.0_MK*mueff-1._MK)/((REAL(N)+2._MK)**2 + mueff))
      damps = (1._MK + 2._MK*max(0.0_MK,sqrt((mueff-1.0_MK)/real(N+1))-1.0_MK)) * &
      		max(0.3_MK,1.0_MK - 1.0_MK*REAL(N)/min(real(options%StopMaxIter, mk), &
      		options%StopMaxFunEvals/lambda)) + cs
      sigma = maxval(input%insigma)
      
      IF (.NOT. ALLOCATED(pc)) THEN
          ALLOCATE(pc(N),stat=allocStat)
          IF(allocStat .NE. 0) STOP 'Error allocating pc'
      END IF
      pc = 0.0_MK
      IF (.NOT. ALLOCATED(ps)) THEN
          ALLOCATE(ps(N),stat=allocStat)
          IF(allocStat .NE. 0) STOP 'Error allocating ps'
      END IF
      ps = 0.0_MK
      
      IF (.NOT. ALLOCATED(B)) THEN
          ALLOCATE(B(N,N),stat=allocStat)
          IF(allocStat .NE. 0) STOP 'Error allocating B'
      END IF
      B = 0.0_MK
      DO i = 1,N
      	B(i,i) = 1.0_MK
      END DO   	
      
      IF (.NOT. ALLOCATED(D)) THEN
  	      ALLOCATE(D(N,N),stat=allocStat)
  	      IF(allocStat .NE. 0) STOP 'Error allocating D'
  	  END IF
  	  D = 0.0_MK
  	  DO i = 1, N
  	  	D(i,i) = input%insigma(i)/maxval(input%insigma)
  	  END DO

      IF (.NOT. ALLOCATED(BD)) THEN 
  	      ALLOCATE(BD(N,N),stat=allocStat)
          IF(allocStat .NE. 0) STOP 'Error allocating internal Matrix BD'
      END IF
      BD = matmul(B,D) 
      
      IF (.NOT. ALLOCATED(C)) THEN
          ALLOCATE(C(N,N),stat=allocStat)
	      IF(allocStat .NE. 0) STOP 'Error allocating internal Matrix C'
	  END IF
      C = matmul(BD,transpose(BD))
      
    
       !------------------------------------------------------------------------
       ! compute inverse of C for bfgs
       !------------------------------------------------------------------------
      IF (options%BFGS_use ) THEN
        IF (.NOT. ALLOCATED(invsC)) THEN
          ALLOCATE(invsC(N,N),stat=allocStat)
	      IF(allocStat .NE. 0) STOP 'Error allocating internal Matrix C'
	    END IF
	  
      
       invsC = 0._MK
        DO i = 1, N
            invsC(i,i) = 1.0_MK/(sigma**2*D(i,i))
        END DO
      ! B * invD*
#ifdef __SP
        CALL SGEMM('N','N',input%N,input%N,input%N,1.0,B,input%N,invsC,input%N,&
                                                        0.0,BinvD,input%N)		
#else
        CALL DGEMM('N','N',input%N,input%N,input%N,1.0d0,B,input%N,invsC,input%N,&
                                                        0.0d0,BinvD,input%N)
#endif
      ! (B * invD*) * B^t
#ifdef __SP
        CALL SGEMM('N','N',input%N,input%N,input%N,1.0,BinvD,input%N,&
        transpose(B),input%N,0.0,invsC,input%N)		
#else
        CALL DGEMM('N','N',input%N,input%N,input%N,1.0d0,BinvD,input%N,&
        transpose(B),input%N,0.0d0,invsC,input%N)
#endif
      END IF

      
!      ALLOCATE(C_PSO(N,N),stat=allocStat)
!      IF(allocStat .NE. 0) STOP 'Error allocating internal Matrix C_PSO'
!      C_PSO = 0._MK
      !-------------------------------------------------------------------------
      !  Allocate Fitness variables
      !-------------------------------------------------------------------------
      
      ! little hack to round up to next integer
      temp = int(real(3*10*N)/real(lambda))
      IF(real(temp) .NE. real(3*10*N)/real(lambda)) &
      	temp = temp + 1.0_MK

      IF (ALLOCATED(fitness%hist)) DEALLOCATE(fitness%hist)
      ALLOCATE(fitness%hist(10+int(temp)),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating fitness.hist'
      fitness%hist = posInf

      IF (ALLOCATED(fitness%histsel)) DEALLOCATE(fitness%histsel)
      ALLOCATE(fitness%histsel(10+int(temp)),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating fitness.histsel'
      fitness%histsel = posInf
      
      IF (ALLOCATED(fitness%raw)) DEALLOCATE(fitness%raw)
      ALLOCATE(fitness%raw(lambda),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating fitness.raw'
      fitness%raw = posInf
      
      IF (ALLOCATED(fitness%sel)) DEALLOCATE(fitness%sel)
	  ALLOCATE(fitness%sel(lambda),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating fitness.sel'
      fitness%sel = 0.0_MK
      
      IF (ALLOCATED(fitness%idx)) DEALLOCATE(fitness%idx)
      ALLOCATE(fitness%idx(lambda),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating fitness.idx'
      DO i = 1, lambda
      	fitness%idx(i) = real(i)
      END DO
      
      IF (ALLOCATED(fitness%idxsel)) DEALLOCATE(fitness%idxsel)
      ALLOCATE(fitness%idxsel(lambda),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating fitness.idxsel'
      DO i = 1, lambda
      	fitness%idxsel(i) = real(i)
      END DO    
     
      !-------------------------------------------------------------------------
      !  Set more variables
      !-------------------------------------------------------------------------
      IF (.NOT. ALLOCATED(maxdx)) THEN
          ALLOCATE(maxdx(N),stat=allocStat)
          IF(allocStat .NE. 0) STOP 'Error allocating maxdx'
      END IF
	  maxdx = options%DiffMaxChange
  	   
	  IF (.NOT. ALLOCATED(mindx)) THEN
	      ALLOCATE(mindx(N),stat=allocStat)
          IF(allocStat .NE. 0) STOP 'Error allocating mindx'
      END IF
	  mindx = options%DiffMinChange
	  
	  !expectation of ||N(0,I)|| == norm(randn(N,1))
	  chiN = sqrt(REAL(N))*((1.0_MK-1.0_MK/(4.0_MK*REAL(N)))+1.0_MK/(21.0_MK*REAL(N*N)))
	  
	  ! normalize recombination weights array
	  weights = weights/sum(weights)
	  
	  stopflag = ''

	  
      RETURN
      END SUBROUTINE cmaes_init
     
