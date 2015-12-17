      !-------------------------------------------------------------------------
      !  Subroutine   :          cmaes_write_output
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This Subroutine writes CMA data to multiple txt files
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
      SUBROUTINE cmaes_write_output
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE cmaes_param_mod
      USE cmaes_mod
      USE cmaes_out_mod
      USE cmaes_opts_mod
#ifdef __HAVE_MPI__    
      USE mpi_mod
      USE pso_mod
#endif
      
      
      !-------------------------------------------------------------------------
      !  Local Variables
      !-------------------------------------------------------------------------
      INTEGER                                            :: uNum
      CHARACTER(len=200)                                 :: tmp_dir
      INTEGER                                            :: ioError
      CHARACTER(len=20),SAVE                             :: formN
#ifndef __HAVE_MPI__
      CHARACTER(len=10)                                  :: proc
#endif      
      !-------------------------------------------------------------------------
      !  Initialize some variables
      !-------------------------------------------------------------------------           
      uNum = 26+22

#ifndef __HAVE_MPI__    
      ! variable from mpi_mod
      proc = ''    
#endif  

      CALL tool_formatarrays(formN,input%N)
      !-------------------------------------------------------------------------
      !  bestever%x
      !-------------------------------------------------------------------------            
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/out_bestever_x'& 
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)
      WRITE(uNum,formN) bestever%x
      CLOSE(unit=uNum)
      !-------------------------------------------------------------------------
      !  bestever%f
      !-------------------------------------------------------------------------            
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/out_bestever_f'& 
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)
      WRITE(uNum,*) bestever%f
      CLOSE(unit=uNum)      
      !-------------------------------------------------------------------------
      !  bestever%evals
      !-------------------------------------------------------------------------            
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/out_bestever_evals'& 
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)
      WRITE(uNum,*) bestever%evals
      CLOSE(unit=uNum)      
      
    
      !-------------------------------------------------------------------------
      !  N
      !-------------------------------------------------------------------------            
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/out_N'& 
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)
      WRITE(uNum,*) real(input%N)
      CLOSE(unit=uNum)
      !-------------------------------------------------------------------------
      !  countEval
      !-------------------------------------------------------------------------            
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/out_countEval'& 
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)
      WRITE(uNum,*) real(countEval)
      CLOSE(unit=uNum)
      !-------------------------------------------------------------------------
      !  countEvalNaN
      !-------------------------------------------------------------------------            
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/out_countEvalNaN'& 
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)
      WRITE(uNum,*) real(countEvalNaN)
      CLOSE(unit=uNum)    
      !-------------------------------------------------------------------------
      !  countIter
      !-------------------------------------------------------------------------            
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/out_countIter'& 
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)
      WRITE(uNum,*)  real(countIter)
      CLOSE(unit=uNum)          
      !-------------------------------------------------------------------------
      !  countOutOfBounds
      !-------------------------------------------------------------------------            
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/out_countOutOfBounds'& 
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)
      WRITE(uNum,*)  real(countOutOfBounds)
      CLOSE(unit=uNum)          
      !-------------------------------------------------------------------------
      !  funcName
      !-------------------------------------------------------------------------            
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/out_funcName'& 
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)
      WRITE(uNum,*)  options%funcName
      CLOSE(unit=uNum)  
      !-------------------------------------------------------------------------
      !  insigma
      !-------------------------------------------------------------------------            
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/out_insigma'& 
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)
      WRITE(uNum,formN)  input%insigma
      CLOSE(unit=uNum)              
      !-------------------------------------------------------------------------
      !  lambda
      !-------------------------------------------------------------------------            
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/out_lambda'& 
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)
      WRITE(uNum,*)  real(lambda)
      CLOSE(unit=uNum)
      !-------------------------------------------------------------------------
      !  mu
      !-------------------------------------------------------------------------            
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/out_mu'& 
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)
      WRITE(uNum,*)  real(mu)
      CLOSE(unit=uNum)     
      !-------------------------------------------------------------------------
      !  mueff
      !-------------------------------------------------------------------------            
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/out_mueff'& 
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)
      WRITE(uNum,*) mueff
      CLOSE(unit=uNum)    
      !-------------------------------------------------------------------------
      ! stopflag
      !-------------------------------------------------------------------------            
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/out_stopflag'& 
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)
      WRITE(uNum,*) stopflag
      CLOSE(unit=uNum)          
      
      
       
      !-------------------------------------------------------------------------
      !  weights
      !-------------------------------------------------------------------------            
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/out_weights'& 
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)
      WRITE(uNum,formN) weights
      CLOSE(unit=uNum)     
      !-------------------------------------------------------------------------
      !  xstart
      !-------------------------------------------------------------------------            
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/out_xstart'& 
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)
      WRITE(uNum,formN) input%xstart
      CLOSE(unit=uNum)    
      
      !-------------------------------------------------------------------------
      !  BFGS Counts if bfgs is used
      !-------------------------------------------------------------------------            
      IF (options%BFGS_use ) THEN
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/out_BFGSexitcounter'& 
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)
      WRITE(uNum,*) BFGSexitcounter(1)
      WRITE(uNum,*) BFGSexitcounter(2)
      WRITE(uNum,*) BFGSexitcounter(3)
      WRITE(uNum,*) BFGSexitcounter(4)
      WRITE(uNum,*) BFGSexitcounter(5)
      WRITE(uNum,*) BFGSexitcounter(6)
      CLOSE(unit=uNum)
      END IF
      
      !-----------------------------------------------------------------------
      !  Parse and read input data
      !-----------------------------------------------------------------------
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/out_settings'& 
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)
      
      
      
      
      WRITE(uNum,*) 'REL_SIGMA'
        WRITE(uNum,*)  options%rel_sigma
      WRITE(uNum,*) 'ABS_SIGMA'
        WRITE(uNum,*)  options%abs_sigma
      WRITE(uNum,*) 'ALLDIM_LBOUNDS'
        WRITE(uNum,*)  options%alldim_lBounds
      WRITE(uNum,*) 'ALLDIM_UBOUNDS'
        WRITE(uNum,*)  options%alldim_uBounds
      WRITE(uNum,*) 'DIMENSIONS'
        WRITE(uNum,*)  options%dimensions
      WRITE(uNum,*) 'POPSIZE'
        WRITE(uNum,*)  options%PopSize
      WRITE(uNum,*) 'USE_CEC'
        WRITE(uNum,*)  options%use_CEC
      WRITE(uNum,*) 'USE_LJ'
        WRITE(uNum,*)  options%use_LJ
      WRITE(uNum,*) 'USE_BBOB'
        WRITE(uNum,*)  options%use_BBOB
      WRITE(uNum,*) 'USE_TIP'
        WRITE(uNum,*)  options%use_TIP
      WRITE(uNum,*) 'BENCHFCTNR'
        WRITE(uNum,*)  options%Benchfctnr
      WRITE(uNum,*) 'OUTPUT_FOLDER'
          WRITE(uNum,*)  options%output_folder
      WRITE(uNum,*) 'STOPFITNESS'
        WRITE(uNum,*)  options%StopFitness
      WRITE(uNum,*) 'STOPMAXFUNEVALS'
        WRITE(uNum,*)  options%StopMaxFunEvals
      WRITE(uNum,*) 'STOPMAXITER'
        WRITE(uNum,*)  options%StopMaxIter
      WRITE(uNum,*) 'STOPTOLX'
        WRITE(uNum,*)  options%StopTolX
      WRITE(uNum,*) 'STOPTOLUPX'
        WRITE(uNum,*)  options%StopTolUpX
      WRITE(uNum,*) 'STOPTOLFUN'
        WRITE(uNum,*)  options%StopTolFun
      WRITE(uNum,*) 'STOPTOLHISTFUN'
        WRITE(uNum,*)  options%StopTolHistFun
      WRITE(uNum,*) 'STOPONWARNINGS'
        WRITE(uNum,*)  options%StopOnWarnings
       WRITE(uNum,*) 'WARNONEQUALFUNCTIONVALUES'
         WRITE(uNum,*)  options%WarnOnEqualFunctionValues
      WRITE(uNum,*) 'EVALINITIALX'
        WRITE(uNum,*)  options%EvalInitialX
      WRITE(uNum,*) 'PARENTNUMBER'
        WRITE(uNum,*)  options%ParentNumber
      WRITE(uNum,*) 'RECOMBINATIONWEIGHTS'
        WRITE(uNum,*)  options%RecombinationWeights
      WRITE(uNum,*) 'VERBOSEMODULO'
        WRITE(uNum,*)  options%VerboseModulo
      WRITE(uNum,*) 'FLGGENDATA'
        WRITE(uNum,*)  options%flgGenData
      WRITE(uNum,*) 'INTGENDATA'
        WRITE(uNum,*)  options%intGenData
      WRITE(uNum,*) 'FUNCNAME'
        WRITE(uNum,*)  options%funcName
      WRITE(uNum,*) 'PSCMA'
        WRITE(uNum,*)  options%pscma
      WRITE(uNum,*) 'PSOWEIGHT'
        WRITE(uNum,*)  options%psoWeight
      WRITE(uNum,*) 'PSOFREQ'
        WRITE(uNum,*)  options%psoFreq
      WRITE(uNum,*) 'ACCURACY'
        WRITE(uNum,*)  options%accuracy
      WRITE(uNum,*) 'GLOBAL_MIN'
        WRITE(uNum,*)  options%global_min
      WRITE(uNum,*) 'USE_SEED'
        WRITE(uNum,*)  options%use_seed
      WRITE(uNum,*) 'SILENTMPI'
        WRITE(uNum,*)  options%silentMPI
      WRITE(uNum,*) 'QR_SAMPLING'
        WRITE(uNum,*)  options%qr_sampling
      WRITE(uNum,*) 'QR_SAMPLER'
        WRITE(uNum,*)  options%qr_sampler
      WRITE(uNum,*) 'QR_SCRAMBLING'
        WRITE(uNum,*)  options%qr_scrambling
      WRITE(uNum,*) 'QR_INVERTER'
        WRITE(uNum,*)  options%qr_inverter
      WRITE(uNum,*) 'USE_INIT_BOUNDS'
        WRITE(uNum,*)  options%use_init_bounds
      WRITE(uNum,*) 'INIT_LBOUNDS'
        WRITE(uNum,*)  options%init_lBounds
      WRITE(uNum,*) 'INIT_UBOUNDS'
        WRITE(uNum,*)  options%init_uBounds
      WRITE(uNum,*) 'SEED_FOLDER'
          WRITE(uNum,*) options%seed_folder
      WRITE(uNum,*) 'BENCHMARK'
        WRITE(uNum,*)  options%benchmark
      WRITE(uNum,*) 'RECORD_ACCURACY'
        WRITE(uNum,*)  options%record_accuracy
      WRITE(uNum,*) 'RESTART_CMA'
        WRITE(uNum,*)  options%restart_cma
      WRITE(uNum,*) 'RESTART_TYPE'
        WRITE(uNum,*)  options%restart_type
      WRITE(uNum,*) 'RESTARTS'
        WRITE(uNum,*)  options%Restarts
      WRITE(uNum,*) 'INCPOPSIZE'
        WRITE(uNum,*)  options%IncPopSize
      WRITE(uNum,*) 'RECORD_BESTHIST'
        WRITE(uNum,*)  options%record_besthist
      WRITE(uNum,*) 'RECORD_MODULO'
        WRITE(uNum,*)  options%record_modulo
      WRITE(uNum,*) 'MAXINCFAC'
        WRITE(uNum,*)  options%MaxIncFac
      WRITE(uNum,*) 'BFGS_USE'
        WRITE(uNum,*)  options%BFGS_use
      WRITE(uNum,*) 'BFGS_FACTR'
        WRITE(uNum,*)  options%BFGS_factr
      WRITE(uNum,*) 'BFGS_PGTOL'
        WRITE(uNum,*)  options%BFGS_pgtol
      WRITE(uNum,*) 'BFGS_GRAD_STEPSIZE'
        WRITE(uNum,*)  options%BFGS_grad_stepsize
      WRITE(uNum,*) 'BFGS_POSITION'
        WRITE(uNum,*)  options%BFGS_position
      WRITE(uNum,*) 'BFGS_CENTRAL_DIFFERENCE'
        WRITE(uNum,*)  options%BFGS_central_difference
      WRITE(uNum,*) 'BFGS_SIGMASTEP_MIN'
        WRITE(uNum,*)  options%BFGS_sigmastep_min
      WRITE(uNum,*) 'BFGS_CONFPROP'
        WRITE(uNum,*)  options%BFGS_confprop      
      WRITE(uNum,*) 'STOPTIME'
        WRITE(uNum,*)  options%StopTime
      WRITE(uNum,*) 'STOPTIMEHH'
        WRITE(uNum,*)  options%StopTimeHH
      WRITE(uNum,*) 'STOPTIMEMM'
        WRITE(uNum,*)  options%StopTimeMM
      WRITE(uNum,*) 'STOPTIMESS'
        WRITE(uNum,*)  options%StopTimeSS
      WRITE(uNum,*) 'CECFOLDERS'
          WRITE(uNum,*) options%CECFolders
      WRITE(uNum,*) 'WRITE_PDB'
        WRITE(uNum,*)  options%write_pdb            
      WRITE(uNum,*) 'LJ_COMP'
        WRITE(uNum,*)  options%LJ_comp
      WRITE(uNum,*) 'USE_LJ_COMP'
        WRITE(uNum,*)  options%use_LJ_comp
      WRITE(uNum,*) 'USE_DF'
        WRITE(uNum,*)  options%use_DF
      WRITE(uNum,*) 'USE_MATFUNC'
        WRITE(uNum,*)  options%use_MATFUNC
      WRITE(uNum,*) 'DF_S'
        WRITE(uNum,*)  options%DF_s
      WRITE(uNum,*) 'DF_RAST'
        WRITE(uNum,*)  options%DF_rast
      WRITE(uNum,*) 'USE_RANDOM_LANDSCAPE'
        WRITE(uNum,*)  options%use_RANDOM_LANDSCAPE
      WRITE(uNum,*) 'FLGOUTTXT'
        WRITE(uNum,*)  options%flgouttxt
      WRITE(uNum,*) 'FLGGENTRACE'
        WRITE(uNum,*)  options%flgGenTrace
      CLOSE(unit=uNum)        
      
      
      END SUBROUTINE cmaes_write_output
