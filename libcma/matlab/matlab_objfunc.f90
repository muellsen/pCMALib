#include "fintrf.h"

	   !-------------------------------------------------------------------------
       !  Function	  :          MATLAB_OBJFUNC
       !-------------------------------------------------------------------------
       !
       !  Purpose      : This is a template function how to call
       !				 MATLAB objective functions from Fortran
       !
       !  Input		  :  vars
       !				 m
       !				 n
       !				 lbounds
       !				 ubounds
       !
       !  Output	   : res
       !
       !  Remarks      : Each process starts a MATLAB engine at the beginning of 
       !				 the optimization run and closes it when the process terminates
       !
       !  References   : see www.mathworks.com for further info
       !
       !  Revisions    :
       !-------------------------------------------------------------------------
       !-------------------------------------------------------------------------
      !  pCMALib: a parallel fortran 90 library for the evolution strategy with
      !           covariance matrix adaptation
      !  Christian L. Mueller, Benedikt Baumgartner, Georg Ofenbeck
      !  MOSAIC group, ETH Zurich, Switzerland
      !-------------------------------------------------------------------------
       SUBROUTINE matlab_objfunc(res,vars,m,n,lbounds,ubounds) 
	  
	   !-------------------------------------------------------------------------
	   !  Modules
	   !-------------------------------------------------------------------------
	   USE cmaes_opts_mod
	   
	             
	   ! Parameters
	   REAL(MK), DIMENSION(n),  INTENT(out)           :: res
	   REAL(MK), DIMENSION(m,n),INTENT(in)            :: vars
	   INTEGER, INTENT(in)                            :: m
	   INTEGER, INTENT(in)                            :: n
	   
	   REAL(MK), OPTIONAL                             :: lbounds
	   REAL(MK), OPTIONAL                             :: ubounds
	  
	   
	   ! Variables needed for communication with MATLAB engine
	   CHARACTER(LEN=256)                    		  :: MATLAB_CMD
	   mwPointer 									  :: mxCreateDoubleMatrix
       mwPointer 									  :: mxGetPr
       mwPointer 									  :: engOpen,engClose
	   mwPointer 									  :: engGetVariable
       mwPointer 									  :: engPutVariable, engEvalString
	   mwPointer 									  :: temp, status
	   
	   ! Static variables that control the unique start of engine
	   mwPointer, SAVE 								  :: ep 
	   INTEGER, SAVE 								  :: engineStarted
	   
	   ! Variables for objective function variables and function value
	   mwPointer V, R 
       
       ! Adapt the local path to the MATLAB executable HERE
       ! This is e.g. the path on a Apple MacBook Pro
       MATLAB_CMD = '/Applications/MATLAB_R2009b.app/bin/matlab -nodisplay -nojvm -nosplash'
	   !MATLAB_CMD = '/cluster/apps/matlab/7.6/x86_64/bin/matlab -nodisplay -nojvm -nosplash'

	   ! Terminate MATLAB engine when the FORTRAN process is terminated in cmaesRun.f90 (line 914)		
	   IF (m.EQ.-1) THEN
	   		status=engClose(ep)    	   
	   		IF (status .ne. 0) THEN
          		WRITE(*,*) 'Failed to close MATLAB engine, kill process manually'
          	ELSE
          		WRITE(*,*) 'Closed MATLAB engine'
          	ENDIF
          	res=0.0
          	RETURN
	   ENDIF
	   	
	   ! Start MATLAB engine
       IF (engineStarted .NE. 1) THEN
          ep = engOpen(MATLAB_CMD)
       	  engineStarted = 1
       ENDIF

	   ! Fill variables var into a MATLAB pointer	
       V = mxCreateDoubleMatrix(m, n, 0)
       CALL mxCopyReal8ToPtr(vars, mxGetPr(V), m)

	   ! Place the variable V into the MATLAB workspace.
       status = engPutVariable(ep, 'V', V)

       IF (status .ne. 0) THEN
          WRITE(6,*) 'engPutVariable failed'
          STOP
       ENDIF

       ! Evaluate the objective function in MATLAB
       ! In this example we calculate the sphere function f(x) = sum(x.^2)

       IF (engEvalString(ep, 'funcVal = sum(V.^2);') .ne. 0) THEN
          write(6,*) 'engEvalString failed'
          STOP
       ENDIF

       ! Collect the objective function value from the MATLAB workspace    
       R = engGetVariable(ep, 'funcVal')
       CALL mxCopyPtrToReal8(mxGetPr(R), res, n)

       ! Destroy local variables
       CALL mxDestroyArray(V)
       CALL mxDestroyArray(R)
       
       ! Uncomment the following lines to check whether the sphere function in MATLAB
       ! produces the correct values 
       ! WRITE(*,*) '--------------------------------'
       ! WRITE(*,*) 'Function value in MATLAB: ', res
       ! WRITE(*,*) 'Function value in FORTRAN:', sum(vars**2)
       ! WRITE(*,*) '--------------------------------'

            
       RETURN	   
	   END SUBROUTINE matlab_objfunc
