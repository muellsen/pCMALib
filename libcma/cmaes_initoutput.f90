 	  !-------------------------------------------------------------------------
      !  Subroutine   :          cmaes_initoutput
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This Subroutine initializes CMA output variables
      !					in cmaes_out_mod
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
      SUBROUTINE cmaes_initoutput(fitfun)
      USE cmaes_param_mod
      USE cmaes_mod
      USE cmaes_out_mod
      USE cmaes_opts_mod
      IMPLICIT NONE
      EXTERNAL	:: fitfun
      REAL(MK)			:: cmaes_funcwrap
      INTEGER :: allocStat, i
      REAL(MK),DIMENSION(input%N) :: temp
      
      !-------------------------------------------------------------------------
      !  Initialize solution values
      !-------------------------------------------------------------------------
      out%solutions%evals = 0
      ALLOCATE(out%solutions%mean%x(input%N),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating out.solutions.mean.x'
      out%solutions%mean%x = posInf
      out%solutions%mean%f = posInf
      out%solutions%mean%evals = 0
      
      ALLOCATE(out%solutions%recentbest%x(input%N),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating out.solutions.recentbest.x'
      out%solutions%recentbest%x = posInf
      out%solutions%recentbest%f = posInf
      out%solutions%recentbest%evals = 0
      
      ALLOCATE(out%solutions%recentworst%x(input%N),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating out.solutions.recentworst.x'
      out%solutions%recentworst%x = posInf
      out%solutions%recentworst%f = posInf
      out%solutions%recentworst%evals = 0
      
      
      !-------------------------------------------------------------------------
      !  Initialize bestever output values
      !-------------------------------------------------------------------------
      ALLOCATE(out%solutions%bestever%x(input%N),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating out.solutions.bestever.x'
      out%solutions%bestever%x = xmean
      out%solutions%bestever%f = posInf
      out%solutions%bestever%evals = countEval
      bestever = out%solutions%bestever
      
      IF(options%EvalInitialX) THEN
      
      	fitness%hist(1) = cmaes_funcwrap(xmean,input%N,fitfun)
      	fitness%histsel(1) = fitness%hist(1)
      	counteval = counteval + 1
      	IF(fitness%hist(1) .LT. out%solutions%bestever%f) THEN
      		out%solutions%bestever%x = xmean
      		out%solutions%bestever%f = fitness%hist(1)
      		out%solutions%bestever%evals = counteval
      		bestever = out%solutions%bestever
        END IF
      ELSE
      	fitness%histsel(1) = posInf
      	fitness%hist(1) = posInf
      END IF
      
      out%evals = counteval
      out%stopflag = 0.0 !	TODO
      
	  !-------------------------------------------------------------------------
      !  Initialize output history values
      !-------------------------------------------------------------------------
      out%hist%evals = countEval
      out%hist%iterations = countIter
     
      ALLOCATE(out%hist%mean%x(input%N),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating out.hist.mean.x'
      out%hist%mean%x = xmean
      out%hist%mean%f = fitness%hist(1)
      out%hist%mean%evals = countEval
     
      ALLOCATE(out%hist%recentbest%x(input%N),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating out.hist.recentbest.x'
      out%hist%recentbest%x = xmean
      out%hist%recentbest%f = fitness%hist(1)
      out%hist%recentbest%evals = countEval
     
      ALLOCATE(out%hist%recentworst%x(input%N),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating out.hist.recentworst.x'
      out%hist%recentworst%x = xmean
      out%hist%recentworst%f = fitness%hist(1)
      out%hist%recentworst%evals = countEval
      
      out%hist%param%evals = countEval
      out%hist%param%iterations = countIter
      out%hist%param%sigma = sigma
      DO i = 1, input%N
      	temp(i) = C(i,i)
      END DO
      out%hist%param%maxstd = sigma * sqrt(maxval(temp))
      out%hist%param%minstd = sigma * sqrt(minval(temp))
      out%hist%param%maxstdidx = maxval(temp)
      out%hist%param%minstdidx = minval(temp)
      
      ALLOCATE(out%histParamArr%stds(input%N),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating out.histParamArr.stds'
      
      out%histParamArr%stds = sigma * sqrt(temp)
      
            DO i = 1, input%N
      	temp(i) = D(i,i)
      END DO
      out%hist%param%maxD = maxval(temp)
      out%hist%param%minD = minval(temp)
      
      out%histParamArr%evals = countEval
      out%histParamArr%iterations = countIter
      out%histParamArr%sigma = sigma
      ALLOCATE(out%histParamArr%diagD(input%N),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating out.histParamArr.diagD'
      out%histParamArr%diagD = temp
      ALLOCATE(out%histParamArr%Bmax(input%N),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating out.histParamArr.Bmax'  
      ALLOCATE(out%histParamArr%Bmin(input%N),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating out.histParamArr.Bmin'      
      DO i = 1, input%N    
      	out%histParamArr%Bmax(i) = B(i,out%hist%param%maxstdidx)
      	out%histParamArr%Bmin(i) = B(i,out%hist%param%minstdidx)
      END DO
      
      outiter = 0
      
      RETURN
      END SUBROUTINE cmaes_initoutput