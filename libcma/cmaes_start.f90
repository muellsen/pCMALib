 	  !-------------------------------------------------------------------------
      !  Subroutine   :          cmaes_start
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This Subroutine starts the CMA-ES
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
      SUBROUTINE cmaes_start(fitfun,xstart,insigma,varargin)
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE cmaes_param_mod
      USE cmaes_mod, ONLY : input
      USE cmaes_opts_mod
      USE CEC2005
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      EXTERNAL							:: fitfun
      REAL(MK), INTENT(in) 				:: xstart(:)
      REAL(MK), INTENT(in)				:: insigma(:)
      REAL(MK), INTENT(in), OPTIONAL 	:: varargin
      
      !-------------------------------------------------------------------------
      !  Local Variables
      !-------------------------------------------------------------------------
      INTEGER							:: i
      INTEGER							:: allocStat
      REAL(MK)				         	:: cmaes_funcwrap
      REAL                              :: temp
      !-------------------------------------------------------------------------
      !  Get problem dimension from xstart
      !-------------------------------------------------------------------------
      input%N = size(xstart)
      
      !-------------------------------------------------------------------------
      !  Set xstart
      !-------------------------------------------------------------------------
      ALLOCATE(input%xstart(size(xstart)), stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating xstart'
      DO i = 1,input%N
      	input%xstart(i) = xstart(i)
      END DO
      
      !-------------------------------------------------------------------------
      !  Set sigma
      !-------------------------------------------------------------------------
      ALLOCATE(input%insigma(size(xstart)), stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating sigma'
      DO i = 1,input%N
      	input%insigma(i) = insigma(i)
      END DO
      
      
      !-------------------------------------------------------------------------
      !  Initialize CEC2005 parameters if benchmark is used
      !-------------------------------------------------------------------------
      IF (options%use_CEC ) THEN 
          CALL prepare_benchmark(options%Benchfctnr,options%dimensions,1)
      ELSEIF (options%use_BBOB ) THEN 
          i = options%dimensions
#ifdef __BBOB
          CALL init_BBOB(i,options%Benchfctnr,&
      &                   options%output_folder)
#endif
      END IF
      
      !-------------------------------------------------------------------------
      !  Initialize BBOB parameters if benchmark is used
      !-------------------------------------------------------------------------

      
      
      
      
      
      !-------------------------------------------------------------------------
      !  Initialize CMA parameters
      !-------------------------------------------------------------------------
#ifdef __HAVE_MPI__
	  IF(MY_RANK .EQ. 0) THEN
	  WRITE(*,*) '******************************************************'
	  WRITE(*,*) 'Warnings:'
	  END IF
#else
	  WRITE(*,*) '******************************************************'
	  WRITE(*,*) 'Warnings:'	  
#endif	
      CALL cmaes_init()
      
      !-------------------------------------------------------------------------
      !  Initialize Boundary Handling
      !-------------------------------------------------------------------------
      CALL cmaes_initbounds()
      
      !-------------------------------------------------------------------------
      !  Check CMA consistency
      !-------------------------------------------------------------------------
      CALL cmaes_checkparams()
      
      !-------------------------------------------------------------------------
      !  Initialize output and logging variables 
      !-------------------------------------------------------------------------
      CALL cmaes_initoutput(fitfun)
      
      !-------------------------------------------------------------------------
      !  Start CMA run
      !-------------------------------------------------------------------------
      CALL cmaes_run(fitfun)


      !-------------------------------------------------------------------------
      !  Write ouput txt (if flgouttxt == .TRUE.)
      !-------------------------------------------------------------------------
	  IF (options%flgouttxt) THEN
      CALL cmaes_write_output()
	  END IF
      


      
      !-------------------------------------------------------------------------
      !  Create mat file if Matlab is present
      !-------------------------------------------------------------------------
#ifdef __HAVE_MATLAB__
		IF (options%flgGenData) THEN
      		WRITE(*,*) 'Calling Matlab Engine to write .mat File'
      		CALL cma2mat()
      	END IF

      	! Close MATLAB engine if used
        IF (options%use_MATFUNC) THEN 
            temp=cmaes_funcwrap(xstart,-1,fitfun) 
        ENDIF
#endif
	  
	  !-------------------------------------------------------------------------
      !  Free Memory
      !-------------------------------------------------------------------------
	  CALL cmaes_freememory()
	  
	  !-------------------------------------------------------------------------
      !  free CEC2005 resources if benchmark is used
      !-------------------------------------------------------------------------
      IF (options%use_CEC ) THEN 
          CALL close_benchmark()
      ENDIF
      IF (options%use_BBOB ) THEN
#ifdef __BBOB      
         CALL fgeneric_finalize()
#endif      
      ENDIF
	  


      RETURN
      END SUBROUTINE cmaes_start
