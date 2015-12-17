#include "fintrf.h"

 	  !-------------------------------------------------------------------------
      !  Subroutine   :          cma2mat
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This Subroutine writes CMA data to a .mat file called
      !					cmaesData.mat .
      !                
      !  Remarks      : Be careful with Logicals! Matlab's fortran api doesnt
      !					allow to write Logicals, so the values written in these
      !					fields need to be transformed afterwards if necessary,
      !					using Matlab's logical() function.
      !
      !  References   : MATLAB C and Fortran API
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  pCMALib: a parallel fortran 90 library for the evolution strategy with
      !           covariance matrix adaptation
      !  Christian L. Mueller, Benedikt Baumgartner, Georg Ofenbeck
      !  MOSAIC group, ETH Zurich, Switzerland
      !-------------------------------------------------------------------------
      SUBROUTINE cma2mat
      !-------------------------------------------------------------------------
	  !  Modules
	  !-------------------------------------------------------------------------
	  USE cmaes_param_mod
      USE cmaes_mod
      USE cmaes_out_mod
      USE cmaes_opts_mod
      USE IFPORT
      IMPLICIT NONE
	  !-------------------------------------------------------------------------
	  !  Pointers (mwpointer defined in fintrf.h)
	  !-------------------------------------------------------------------------
      mwpointer mxCreateDoubleMatrix, mxCreateDoubleScalar, mxCreateString, &
     			mxCreateStructArray
      mwpointer matGetVariable, mxGetPr, matOpen
      mwpointer mp
      mwpointer pInputN, pInputXstart, pInputInsigma
      mwpointer pSelLambda, pSelMu, pSelWeights, pSelMueff
      mwpointer pCountEval, pCountEvalNaN, pCountIter, pCountOutOfBounds
      mwpointer pOutBest, pOutBestX, pOutBestF, pOutBestEvals
      mwpointer pOpts,pOptsLBounds,pOptsUBounds, pOptsStopFitness, &
      			pOptsStopMaxFunEvals, pOptsStopMaxIter, pOptsStopFunEvals, &
      			pOptsStopIter, pOptsStopTolX, pOptsStopTolUpX, pOptsStopTolFun,&
      			pOptsStopTolHistFun, pOptsStopOnWarnings, pOptsDiffMaxChange, &
      			pOptsDiffMinChange, pOptsWarnOnEqualFunctionValues, &
      			pOptsEvalInitialX, pOptsPopSize, pOptsParentNumber, &
      			pOptsRecombinationWeights, pbe3,pbe4,pbe5,pacc_evals, &
      			pbesthistory
      mwpointer pStopFlag, pFuncName
      
	  !-------------------------------------------------------------------------
	  !  Functions
	  !-------------------------------------------------------------------------
      INTEGER				:: matClose, mxIsFromGlobalWS
      INTEGER matPutVariable, matPutVariableAsGlobal, matDeleteVariable
	
	  !-------------------------------------------------------------------------
	  !  Local Variables
	  !-------------------------------------------------------------------------
	  ! write status (0 on success)
	  INTEGER							:: status,i
	  ! filename ending (depends on current process)
	  CHARACTER(len=10)					:: proc
	  ! variables used for structs
	  mwsize,DIMENSION(2)				:: dims
      mwsize 							:: mInd
      mwsize 							:: nInd 
	  INTEGER*4,DIMENSION(19) 			:: nFieldInds
      INTEGER*4 						:: nFieldInd
      mwIndex 							:: fieldIndex = 1
          
      CHARACTER(len=30),DIMENSION(19)	:: fieldnames
	  CHARACTER(len=30),DIMENSION(3)	:: fieldnames2
	  CHARACTER(len=300)                    :: errmsg
	  CHARACTER(len=200)                :: tmp_dir,cur_dir
      
      ! input
	  REAL(MK)							:: matInputN
	  REAL(MK),DIMENSION(input%N)		:: matInputXstart
	  REAL(MK),DIMENSION(input%N)		:: matInputInsigma
	  ! selection
	  REAL(MK)							:: matSelLambda
	  REAL(MK)							:: matSelMu
	  REAL(MK),DIMENSION(mu)			:: matSelWeights
	  REAL(MK)							:: matSelMueff
	  ! counts
	  REAL(MK)							:: matCountEval
	  REAL(MK)							:: matCountEvalNaN
	  REAL(MK)							:: matCountIter
	  REAL(MK)							:: matCountOutOfBounds
	  ! output
	  REAL(MK),DIMENSION(input%N)		:: matOutBestX
	  REAL(MK)							:: matOutBestF
	  REAL(MK)							:: matOutBestEvals
	  REAL(MK),DIMENSION(options%StopMaxFunEvals/options%record_modulo) :: &
	                                       mathist
	  ! options
	  REAL(MK),DIMENSION(input%N)		:: matOptsLBounds,matOptsUBounds
	  REAL(MK)							:: matOptsStopFitness
	  REAL(MK)							:: matOptsStopMaxFunEvals
	  REAL(MK)							:: matOptsStopMaxIter
	  REAL(MK)							:: matOptsStopFunEvals
	  REAL(MK)							:: matOptsStopIter
	  REAL(MK)							:: matOptsStopTolX
	  REAL(MK)							:: matOptsStopTolUpX
	  REAL(MK)							:: matOptsStopTolFun
	  REAL(MK)							:: matOptsStopTolHistFun
	  REAL(MK)							:: matOptsStopOnWarnings
	  REAL(MK),DIMENSION(input%N)		:: matOptsDiffMaxChange
	  REAL(MK),DIMENSION(input%N)		:: matOptsDiffMinChange
	  REAL(MK)							:: matOptsWarnOnEqualFunctionValues
	  REAL(MK)							:: matOptsEvalInitialX
	  REAL(MK)							:: matOptsPopSize
	  REAL(MK)							:: matOptsParentNumber
	  REAL(MK)							:: matOptsRecombinationWeights
	  REAL(MK)                          :: matbe3,matbe4,matbe5
      REAL(MK)                          :: matacc_evals
	  
	  !-------------------------------------------------------------------------
	  !  Initialize local variables
	  !-------------------------------------------------------------------------
	  dims(1) = 1
	  dims(2) = 1
	  
	  matInputN = real(input%N)
	  matInputXstart = input%xstart
	  matInputInsigma = input%insigma
	  
	  matSelLambda = real(lambda)
	  matSelMu = real(mu)
	  matSelWeights = weights
	  matSelMueff = mueff
	  
	  matCountEval = real(countEval)
	  matCountEvalNaN = real(countEvalNaN)
	  matCountIter = real(countIter)
	  matCountOutOfBounds = real(countOutOfBounds)
	  
	  matOutBestX = bestever%x
	  matOutBestF = bestever%f
	  matOutBestEvals = real(bestever%evals)
	  
	  matOptsLBounds = options%LBounds
	  matOptsUBounds = options%UBounds
	  matOptsStopFitness = options%StopFitness
	  matOptsStopMaxFunEvals = options%StopMaxFunEvals
	  matOptsStopMaxIter = real(options%StopMaxIter)
	  matOptsStopFunEvals = options%StopFunEvals
	  matOptsStopIter = options%StopIter
	  matOptsStopTolX = options%StopTolX
	  matOptsStopTolUpX = options%StopTolUpX
	  matOptsStopTolFun = options%StopTolFun
	  matOptsStopTolHistFun = options%StopTolHistFun
	  
	  matbe3 = be3
	  matbe4 = be4
	  matbe5 = be5
	  matacc_evals = acc_evals
	  
	  IF (options%record_besthist ) THEN
	    mathist = besthist
	  END IF
	  
	  
	  !Little work-around as logical values are not supported by
	  !Matlabs Fortran API
	  IF(options%StopOnWarnings) THEN
	  	matOptsStopOnWarnings = 1
	  ELSE
	  	matOptsStopOnWarnings = 0
	  END IF
	  IF(options%WarnOnEqualFunctionValues) THEN
	  	matOptsWarnOnEqualFunctionValues = 1
	  ELSE
	  	matOptsWarnOnEqualFunctionValues = 0
	  END IF
	  IF(options%EvalInitialX) THEN
	  	matOptsEvalInitialX = 1
	  ELSE
	  	matOptsEvalInitialX = 0
	  END IF
	  
	  matOptsDiffMaxChange = options%DiffMaxChange
	  matOptsDiffMinChange = options%DiffMinChange
	  matOptsPopSize = options%PopSize
	  matOptsParentNumber = options%ParentNumber
	  matOptsRecombinationWeights = options%RecombinationWeights

	  !-------------------------------------------------------------------------
	  !  Open MAT-file for writing
	  !-------------------------------------------------------------------------	
#ifdef __HAVE_MPI__
	  WRITE(proc,'(I10)') MY_RANK
	  proc = '_' // trim(adjustl(proc)) 
#else
	  proc = ''
#endif	
      status = GETCWD( cur_dir )
      tmp_dir=trim(adjustl(options%output_folder))
      
      
      !workaround cause matOpen cant handle long strings
      
      status = CHDIR(tmp_dir)
      select case (status) 
      case (2)  ! ENOENT        
       errmsg = 'The directory '//TRIM(tmp_dir)//' does not exist'
      case (20) ! ENOTDIR   
       errmsg = TRIM(tmp_dir)//' is not a directory'
      case default 
        write (errmsg,*) 'Error with code ', status     
      end select
      tmp_dir = 'cmaesData' &
     & // trim(adjustl(proc))  // '.mat'
      mp = matOpen(tmp_dir, 'w')
      IF(mp .eq. 0) THEN
         WRITE(6,*) 'Can''t open ''cmaesData.mat'' for writing.'
         WRITE(6,*) '(Do you have write permission in this directory?)'
         STOP
      END IF
      tmp_dir = cur_dir
      status = CHDIR(tmp_dir)
      select case (status) 
      case (2)  ! ENOENT        
       errmsg = 'The directory '//TRIM(tmp_dir)//' does not exist'
      case (20) ! ENOTDIR   
       errmsg = TRIM(tmp_dir)//' is not a directory'
      case default 
        write (errmsg,*) 'Error with code ', status     
      end select
      

		!	  status = matPutVariable(mp, 'phonebook', pbook)
		!	  IF(status .NE. 0) THEN
		!         WRITE(6,*) 'matPutVariable ''phonebook'' failed'
		!         STOP
		!      END IF

	  !-------------------------------------------------------------------------
	  !  WRITE OPTIONS STRUCT
	  !-------------------------------------------------------------------------
	  ! Create struct
	  fieldnames(1) = 'LBounds'
	  fieldnames(2) = 'UBounds'
	  fieldnames(3) = 'StopFitness'
	  fieldnames(4) = 'StopMaxFunEvals'
	  fieldnames(5) = 'StopMaxIter'
	  fieldnames(6) = 'StopFunEvals'
	  fieldnames(7) = 'StopIter'
	  fieldnames(8) = 'StopTolX'
	  fieldnames(9) = 'StopTolUpX'
	  fieldnames(10) = 'StopTolFun'
	  fieldnames(11) = 'StopTolHistFun'
	  fieldnames(12) = 'StopOnWarnings'
	  fieldnames(13) = 'WarnOnEqualFunctionValues'
	  fieldnames(14) = 'DiffMaxChange'
	  fieldnames(15) = 'DiffMinChange'
	  fieldnames(16) = 'EvalInitialX'
	  fieldnames(17) = 'PopSize'
	  fieldnames(18) = 'ParentNumber'
	  fieldnames(19) = 'RecombinationWeights'


      mInd = 2
      nFieldInd = 19
      pOpts = mxCreateStructArray(mInd,dims,nFieldInd,fieldnames)
	  ! Create pointers to LBounds/UBounds arrays
	  ! Important to use mInd,nInd of type mwsize as arguments	
      mInd = input%N
      nInd = 1
	  pOptsLBounds = mxCreateDoubleMatrix(mInd,nInd,0)
	  CALL mxCopyReal8ToPtr(matOptsLBounds,mxGetPr(pOptsLBounds),mInd)
	  pOptsUBounds = mxCreateDoubleMatrix(mInd,nInd,0)
	  CALL mxCopyReal8ToPtr(matOptsUBounds,mxGetPr(pOptsUBounds),mInd)
	  ! Create pointers to StopFlags
	  pOptsStopFitness = mxCreateDoubleScalar(matOptsStopFitness)
      pOptsStopMaxFunEvals = mxCreateDoubleScalar(matOptsStopMaxFunEvals)
      pOptsStopMaxIter = mxCreateDoubleScalar(matOptsStopMaxIter)
      pOptsStopFunEvals = mxCreateDoubleScalar(matOptsStopFunEvals)
      pOptsStopIter = mxCreateDoubleScalar(matOptsStopIter)
      pOptsStopTolX = mxCreateDoubleScalar(matOptsStopTolX)
      pOptsStopTolUpX = mxCreateDoubleScalar(matOptsStopTolUpX)
      pOptsStopTolFun = mxCreateDoubleScalar(matOptsStopTolFun)
      pOptsStopTolHistFun = mxCreateDoubleScalar(matOptsStopTolHistFun)
      pOptsStopOnWarnings = mxCreateDoubleScalar(matOptsStopOnWarnings)
      ! Create pointers to other options
      pOptsWarnOnEqualFunctionValues = mxCreateDoubleScalar( &
            									matOptsWarnOnEqualFunctionValues)
      ! Important to use mInd,nInd of type mwsize as arguments
      mInd = input%N
      nInd = 1
      pOptsDiffMaxChange = mxCreateDoubleMatrix(mInd,nInd,0)
      CALL mxCopyReal8ToPtr(matOptsDiffMaxChange,mxGetPr(pOptsDiffMaxChange),&
      mInd)
      pOptsDiffMinChange = mxCreateDoubleMatrix(mInd,nInd,0)
      CALL mxCopyReal8ToPtr(matOptsDiffMinChange,mxGetPr(pOptsDiffMinChange),&
      mInd)


      pOptsEvalInitialX = mxCreateDoubleScalar(matOptsEvalInitialX)
      pOptsPopSize = mxCreateDoubleScalar(matOptsPopSize)
      pOptsParentNumber = mxCreateDoubleScalar(matOptsParentNumber)
      pOptsRecombinationWeights = &
      					mxCreateDoubleScalar(matOptsRecombinationWeights)
      mInd=1
	  ! Fill struct field indices
      DO i=1,19 
      	nFieldInds(i) = i
      END DO
	  ! Important to use nFieldInds of type INTEGER*4 as arguments
	  ! Important to use fieldIndex of type mwIndex as arguments	
	  CALL mxSetFieldByNumber(pOpts,fieldIndex,nFieldInds(1),pOptsLBounds)
      CALL mxSetFieldByNumber(pOpts,fieldIndex,nFieldInds(2),pOptsUBounds)
	  CALL mxSetFieldByNumber(pOpts,fieldIndex,nFieldInds(3),pOptsStopFitness)
	  CALL mxSetFieldByNumber(pOpts,fieldIndex,nFieldInds(4),pOptsStopMaxFunEvals)
	  CALL mxSetFieldByNumber(pOpts,fieldIndex,nFieldInds(5),pOptsStopMaxIter)
	  CALL mxSetFieldByNumber(pOpts,fieldIndex,nFieldInds(6),pOptsStopFunEvals)
	  CALL mxSetFieldByNumber(pOpts,fieldIndex,nFieldInds(7),pOptsStopIter)
	  CALL mxSetFieldByNumber(pOpts,fieldIndex,nFieldInds(8),pOptsStopTolX)
	  CALL mxSetFieldByNumber(pOpts,fieldIndex,nFieldInds(9),pOptsStopTolUpX)
	  CALL mxSetFieldByNumber(pOpts,fieldIndex,nFieldInds(10),pOptsStopTolFun)
	  CALL mxSetFieldByNumber(pOpts,fieldIndex,nFieldInds(11),pOptsStopTolHistFun)
	  CALL mxSetFieldByNumber(pOpts,fieldIndex,nFieldInds(12),pOptsStopOnWarnings)
	  CALL mxSetFieldByNumber(pOpts,fieldIndex,nFieldInds(13),pOptsWarnOnEqualFunctionValues)
	  CALL mxSetFieldByNumber(pOpts,fieldIndex,nFieldInds(14),pOptsDiffMaxChange)
	  CALL mxSetFieldByNumber(pOpts,fieldIndex,nFieldInds(15),pOptsDiffMinChange)
	  CALL mxSetFieldByNumber(pOpts,fieldIndex,nFieldInds(16),pOptsEvalInitialX)
	  CALL mxSetFieldByNumber(pOpts,fieldIndex,nFieldInds(17),pOptsPopSize)
	  CALL mxSetFieldByNumber(pOpts,fieldIndex,nFieldInds(18),pOptsParentNumber)
	  CALL mxSetFieldByNumber(pOpts,fieldIndex,nFieldInds(19),pOptsRecombinationWeights)

	  
	  status = matPutVariable(mp, 'opts'//trim(adjustl(proc)), pOpts)
	  IF(status .NE. 0) THEN
         WRITE(6,*) 'matPutVariable ''opts'' failed'
         STOP
      END IF
       
      
	  !-------------------------------------------------------------------------
	  !  WRITE INPUT VARIABLES
	  !-------------------------------------------------------------------------
	  
	  !-------------------------------------------------------------------------
	  !  write input%N
	  !-------------------------------------------------------------------------
	  pInputN = mxCreateDoubleScalar(matInputN)
	  status = matPutVariable(mp, 'N'//trim(adjustl(proc)), pInputN)
	  IF(status .NE. 0) THEN
         WRITE(6,*) 'matPutVariable ''N'' failed'
         STOP
      END IF
     ! WRITE(*,*)'pInputN'
      
	  !-------------------------------------------------------------------------
	  !  write input%xstart
	  !-------------------------------------------------------------------------
	  mInd = input%N
	  nInd = 1
	  pInputXstart = mxCreateDoubleMatrix(mInd,nInd,0)
	  CALL mxCopyReal8ToPtr(matInputXstart, mxGetPr(pInputXstart), mInd)
	  status = matPutVariable(mp, 'xstart'//trim(adjustl(proc)), pInputXstart)
	  IF(status .NE. 0) THEN
         WRITE(6,*) 'matPutVariable ''xstart'' failed'
         STOP
      END IF  

      !-------------------------------------------------------------------------
	  !  write input%insigma
	  !-------------------------------------------------------------------------
	  pInputInsigma = mxCreateDoubleMatrix(input%N,mInd,0)
	  CALL mxCopyReal8ToPtr(matInputInsigma, mxGetPr(pInputInsigma), input%N)
	  status = matPutVariable(mp, 'insigma'//trim(adjustl(proc)), pInputInsigma)
	  IF(status .NE. 0) THEN
	  	WRITE(6,*) 'matPutVariable ''insigma'' failed'
	  	STOP
	  END IF
	  
	  !-------------------------------------------------------------------------
	  !  WRITE STRATEGY SELECTION PARAMETERS
	  !-------------------------------------------------------------------------
	
      !-------------------------------------------------------------------------
	  !  write lambda
	  !-------------------------------------------------------------------------	  
	  pSelLambda = mxCreateDoubleScalar(matSelLambda)
	  status = matPutVariable(mp, 'lambda'//trim(adjustl(proc)), pSelLambda)
	  IF(status .NE. 0) THEN
         WRITE(6,*) 'matPutVariable ''lambda'' failed'
         STOP
      END IF
      
      !-------------------------------------------------------------------------
	  !  write mu
	  !-------------------------------------------------------------------------	      
	  pSelMu = mxCreateDoubleScalar(matSelMu)
	  status = matPutVariable(mp, 'mu'//trim(adjustl(proc)), pSelMu)
	  IF(status .NE. 0) THEN
         WRITE(6,*) 'matPutVariable ''mu'' failed'
         STOP
      END IF      
      
      
      !-------------------------------------------------------------------------
	  !  write mueff
	  !-------------------------------------------------------------------------	      
	  pSelMueff = mxCreateDoubleScalar(matSelMueff)
	  status = matPutVariable(mp, 'mueff'//trim(adjustl(proc)), pSelMueff)
	  IF(status .NE. 0) THEN
         WRITE(6,*) 'matPutVariable ''mueff'' failed'
         STOP
      END IF      
    

      !-------------------------------------------------------------------------
	  !  write weights
	  !-------------------------------------------------------------------------
	  mInd = mu
      nInd = 1
         
      pSelWeights = mxCreateDoubleMatrix(mInd,nInd,0)
         
	  CALL mxCopyReal8ToPtr(matSelWeights, mxGetPr(pSelWeights), mInd)
	  status = matPutVariable(mp, 'weights'//trim(adjustl(proc)), pSelWeights)
	  IF(status .NE. 0) THEN
	  	WRITE(6,*) 'matPutVariable ''weights'' failed'
	  	STOP
	  END IF
	  
	  !-------------------------------------------------------------------------
	  !  write bestever history
	  !-------------------------------------------------------------------------
	  IF (options%record_besthist ) THEN
	  mInd = options%StopMaxFunEvals / options%record_modulo
	  nInd = 1
	  pbesthistory = mxCreateDoubleMatrix(mInd,nInd,0)
	  CALL mxCopyReal8ToPtr(mathist, mxGetPr(pbesthistory), mInd)
	  status = matPutVariable(mp, 'besthist'//trim(adjustl(proc)), pbesthistory)
	  IF(status .NE. 0) THEN
         WRITE(6,*) 'matPutVariable ''besthist'' failed'
         STOP
      END IF  
      END IF
   

      !-------------------------------------------------------------------------
	  !  WRITE COUNTS
	  !-------------------------------------------------------------------------
	  
      !-------------------------------------------------------------------------
	  !  write countEval
	  !-------------------------------------------------------------------------	  
	  pCountEval = mxCreateDoubleScalar(matCountEval)
	  status = matPutVariable(mp, 'countEval'//trim(adjustl(proc)), pCountEval)
	  IF(status .NE. 0) THEN
         WRITE(6,*) 'matPutVariable ''countEval'' failed'
         STOP
      END IF      
   
      !-------------------------------------------------------------------------
	  !  write countEvalNaN
	  !-------------------------------------------------------------------------	  
	  pCountEvalNaN = mxCreateDoubleScalar(matCountEvalNaN)
	  status = matPutVariable(mp, 'countEvalNaN'//trim(adjustl(proc)),&
	  															 pCountEvalNaN)
	  IF(status .NE. 0) THEN
         WRITE(6,*) 'matPutVariable ''countEvalNaN'' failed'
         STOP
      END IF      
      
      !-------------------------------------------------------------------------
	  !  write countIter
	  !-------------------------------------------------------------------------	  
	  pCountIter = mxCreateDoubleScalar(matCountIter)
	  status = matPutVariable(mp, 'countIter'//trim(adjustl(proc)), pCountIter)
	  IF(status .NE. 0) THEN
         WRITE(6,*) 'matPutVariable ''countIter'' failed'
         STOP
      END IF     

      !-------------------------------------------------------------------------
	  !  write countOutOfBounds
	  !-------------------------------------------------------------------------	  
	  pCountOutOfBounds = mxCreateDoubleScalar(matCountOutOfBounds)
	  status = matPutVariable(mp, 'countOutOfBounds'//trim(adjustl(proc)), &
	  											pCountOutOfBounds)
	  IF(status .NE. 0) THEN
         WRITE(6,*) 'matPutVariable ''countOutOfBounds'' failed'
         STOP
      END IF   
   
   
      !-------------------------------------------------------------------------
	  !  write Be4,Be5,Be3,acc_evals
	  !------------------------------------------------------------------------
	  IF (options%benchmark ) THEN
	  pbe3 = mxCreateDoubleScalar(matbe3)
	  status = matPutVariable(mp, 'be3'//trim(adjustl(proc)), &
	  											pbe3)
	  IF(status .NE. 0) THEN
         WRITE(6,*) 'matPutVariable ''be3'' failed'
         STOP
      END IF   
      pbe4 = mxCreateDoubleScalar(matbe4)
	  status = matPutVariable(mp, 'be4'//trim(adjustl(proc)), &
	  											pbe4)
	  IF(status .NE. 0) THEN
         WRITE(6,*) 'matPutVariable ''be4'' failed'
         STOP
      END IF   
      pbe5 = mxCreateDoubleScalar(matbe5)
	  status = matPutVariable(mp, 'be5'//trim(adjustl(proc)), &
	  											pbe5)
	  IF(status .NE. 0) THEN
         WRITE(6,*) 'matPutVariable ''be5'' failed'
         STOP
      END IF
      pacc_evals = mxCreateDoubleScalar(matacc_evals)
	  status = matPutVariable(mp, 'acc_evals'//trim(adjustl(proc)), &
	  											pacc_evals)
	  IF(status .NE. 0) THEN
         WRITE(6,*) 'matPutVariable ''acc_evals'' failed'
         STOP
      END IF
      
      END IF
      
      !-------------------------------------------------------------------------
	  !  WRITE OUTPUT (bestever struct)
	  !-------------------------------------------------------------------------
	  ! Important to use fieldnames2 of type CHARACTER(len=30),DIMENSION(3)
	  ! Important to use nFieldInd of type INTEGER*4 as argument
	  	   
	  fieldnames2(1) = 'x'
	  fieldnames2(2) = 'f'
	  fieldnames2(3) = 'evals'
      mInd = 2
      nFieldInd = 3
	  pOutBest = mxCreateStructArray(mInd,dims,nFieldInd,fieldnames2)
    
      mInd = input%N
      nInd = 1
      pOutBestX = mxCreateDoubleMatrix(mInd,nInd,0)	            
      CALL mxCopyReal8ToPtr(matOutBestX,mxGetPr(pOutBestX), mInd)
	  pOutBestF = mxCreateDoubleScalar(matOutBestF)
	  pOutBestEvals = mxCreateDoubleScalar(matOutBestEvals)
      
      nFieldInd = 1
	  CALL mxSetFieldByNumber(pOutBest,fieldIndex,nFieldInd,pOutBestX)
	  nFieldInd = 2
      CALL mxSetFieldByNumber(pOutBest,fieldIndex,nFieldInd,pOutBestF)
	  nFieldInd = 3
      CALL mxSetFieldByNumber(pOutBest,fieldIndex,nFieldInd,pOutBestEvals)
	  
	  status = matPutVariable(mp, 'bestever'//trim(adjustl(proc)), pOutBest)
	  IF(status .NE. 0) THEN
	  	WRITE(6,*) 'matPutVariable ''bestever'' failed'
	  	STOP
	  END IF
	  
	  !-------------------------------------------------------------------------
	  !  WRITE STOPFLAG
	  !------------------------------------------------------------------------- 
	  pStopFlag = mxCreateString(stopflag)
	  status = matPutVariable(mp, 'stopflag'//trim(adjustl(proc)), pStopFlag)
	  IF(status .NE. 0) THEN
	  	WRITE(6,*) 'matPutVariable ''stopflag'' failed'
	  	STOP
	  END IF
	  
	  !-------------------------------------------------------------------------
	  !  WRITE FUNCTION NAME
	  !------------------------------------------------------------------------- 
	  pFuncName = mxCreateString(options%funcName)
	  status = matPutVariable(mp, 'funcName'//trim(adjustl(proc)), pFuncName)
	  IF(status .NE. 0) THEN
	  	WRITE(6,*) 'matPutVariable ''funcName'' failed'
	  	STOP
	  END IF	   

      status = matClose(mp)
      IF(status .NE. 0) THEN
         WRITE(6,*) 'Error closing MAT-file'
         STOP
      END IF
      
      !-------------------------------------------------------------------------
	  !  Clean up memory
	  !-------------------------------------------------------------------------
	  ! Description mxDestroyArray (C and Fortran)
      ! mxDestroyArray deallocates the memory occupied by the specified mxArray.  
      ! mxDestroyArray not only deallocates the memory occupied by the mxArray's 
      ! characteristics fields (such as m and n), but also deallocates all the 
      ! mxArray's associated data arrays, such as pr and pi for complex arrays, 
      ! ir and jc for sparse arrays, fields of structure arrays, and cells of cell arrays. 
      ! Do not call mxDestroyArray on an mxArray you are returning on the left-hand side.
	  
      CALL mxDestroyArray(pInputN)
      CALL mxDestroyArray(pInputXstart)
      CALL mxDestroyArray(pInputInsigma)
      CALL mxDestroyArray(pSelLambda)
      CALL mxDestroyArray(pSelMu)
      CALL mxDestroyArray(pSelWeights)
      CALL mxDestroyArray(pSelMueff)
      CALL mxDestroyArray(pCountEval)
      CALL mxDestroyArray(pCountEvalNaN)
      CALL mxDestroyArray(pCountIter)
      CALL mxDestroyArray(pCountOutOfBounds)
      CALL mxDestroyArray(pOutBest)  
      CALL mxDestroyArray(pOpts)
      
      CALL mxDestroyArray(pStopFlag)
      CALL mxDestroyArray(pFuncName)
      	  
	  ! Destruction of these arrays leads to segmentation faults under linux 
      !CALL mxDestroyArray(pOutBestX)
      !CALL mxDestroyArray(pOutBestF)
      !CALL mxDestroyArray(pOutBestEvals)
      !CALL mxDestroyArray(pOptsLBounds)
      !CALL mxDestroyArray(pOptsUBounds)
      !CALL mxDestroyArray(pOptsStopFitness)
      !CALL mxDestroyArray(pOptsStopMaxFunEvals)
      !CALL mxDestroyArray(pOptsStopMaxIter)
      !CALL mxDestroyArray(pOptsStopFunEvals)
      !CALL mxDestroyArray(pOptsStopIter)
      !CALL mxDestroyArray(pOptsStopTolX)
      !CALL mxDestroyArray(pOptsStopTolUpX) 
      !CALL mxDestroyArray(pOptsStopTolFun)
      !CALL mxDestroyArray(pOptsStopTolHistFun)
      !CALL mxDestroyArray(pOptsStopOnWarnings)
      !CALL mxDestroyArray(pOptsDiffMaxChange)
      !CALL mxDestroyArray(pOptsDiffMinChange)
      !CALL mxDestroyArray(pOptsWarnOnEqualFunctionValues)
      !CALL mxDestroyArray(pOptsEvalInitialX)
      !CALL mxDestroyArray(pOptsPopSize)
      !CALL mxDestroyArray(pOptsParentNumber)
      !CALL mxDestroyArray(pOptsRecombinationWeights)


      
	  WRITE(6,*)
#ifdef __HAVE_MPI__	  
      tmp_dir='Process ' // trim(adjustl(proc))  // ': '// &
     & 'data saved to ' // trim(adjustl(options%output_folder)) &
     & // '/cmaesData' // trim(adjustl(proc)) // '.mat'
      WRITE(6,*) 
#else  
	  WRITE(6,*) 'Data saved to ' // trim(adjustl(options%output_folder)) &
	 & // '/cmaesData.mat' 
#endif	     	
      
      RETURN
      END SUBROUTINE
