      !-------------------------------------------------------------------------
      !  Module       :                   cmaes_opts_mod
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains all CMA options the user can set
      !                
      !  Remarks      : User options can be set by including this module and
      !					assigning values to the options-type, i.e.
      !					options%StopFitness = ...
      !					PLEASE MIND: some options are arrays and need to be
      !					ALLOCATED EXPLICITLY before usage!!!
      !					Naming and default values taken from cmaes.m, variables
      !					are initialized in cmaesInit.f90, taking into account
      !					the user input and options settings
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
      MODULE cmaes_opts_mod
      USE cmaes_param_mod
      IMPLICIT NONE
      SAVE
      REAL(MK), PARAMETER :: posInf = huge(posInf)
      REAL(MK), PARAMETER :: negInf = -huge(negInf)
      
      
      TYPE :: opts
        !-----------------------------------------------------------------------
        !  Cructial Parameters you might wanna set
        !-----------------------------------------------------------------------
        !defines the start sigma as absolut value (doesnt need boundaries - used if non zero
        !rel_sigma will be ignored if abs_sigma is set)
        REAL(MK)    :: abs_sigma = 0.0_MK
        !defines the start sigma relative to the boundaries (to init bounds if set)
        REAL(MK)    :: rel_sigma = 0.2_MK

        !Lower and Upper Boundary in all Dimensions
        REAL(MK)    :: alldim_lBounds = negInf
        REAL(MK)    :: alldim_uBounds = posInf

        !if set true CMA will start within init Bounds
        LOGICAL     :: use_init_bounds = .FALSE.
        REAL(MK)    :: init_lBounds = negInf
        REAL(MK)    :: init_uBounds = posInf

        !Number of Dimension the Optimization Problem has
        INTEGER     :: dimensions = 2
        
        !-----------------------------------------------------------------------
        !  Output Parameters
        !-----------------------------------------------------------------------
        !folder in which the output should be written (relative to program)
        CHARACTER(len=200)       :: output_folder = 'out'

      
        !-----------------------------------------------------------------------
        !  Stopping criteria (value of stop flag)
        !-----------------------------------------------------------------------
        != '-Inf stop if f(xmin) < stopfitness, minimization'
        REAL(MK) 	:: StopFitness = negInf
        != 'Inf  maximal number of fevals'
        REAL(MK) 	:: StopMaxFunEvals = posInf 
        != '1e3*(N+5)^2/sqrt(popsize) maximal number of iterations'
        INTEGER :: StopMaxIter = 0    
        != 'Inf  stop after resp. evaluation'
        REAL(MK) 	:: StopFunEvals = posInf
        != 'Inf  stop after resp. iteration'
        REAL(MK) 	:: StopIter = posInf 
        != '1e-11*max(insigma) stop if x-change smaller TolX' 
        REAL(MK) 	:: StopTolX = 0. 
        != '1e3*max(insigma) stop if x-changes larger TolUpX'     
        REAL(MK) 	:: StopTolUpX = 0. 
        != '1e-12 stop if fun-changes smaller TolFun'     
        REAL(MK) 	:: StopTolFun = 1.E-12  
        != '1e-13 stop if back fun-changes smaller TolHistFun'    
        REAL(MK) 	:: StopTolHistFun = 1.E-13 
        ! = 'yes' 
        LOGICAL 	:: StopOnWarnings = .TRUE.
        LOGICAL     :: StopTime       = .FALSE.
        REAL(MK)    :: StopTimeHH     = 00_MK
        REAL(MK)    :: StopTimeMM     = 00_MK
        REAL(MK)    :: StopTimeSS     = 00_MK
        !-----------------------------------------------------------------------
        !  Others
        !-----------------------------------------------------------------------
        != 'Inf  maximal variable change(s)
        REAL(MK),POINTER,DIMENSION(:) 		:: DiffMaxChange 
        != '0    minimal variable change(s)
        REAL(MK),POINTER,DIMENSION(:) 		:: DiffMinChange 
        != 'on'
        LOGICAL 					  	:: WarnOnEqualFunctionValues = .TRUE. 
        != '-Inf lower bounds
        REAL(MK),POINTER,DIMENSION(:) 		:: LBounds 
        != 'Inf  upper bounds
        REAL(MK),POINTER,DIMENSION(:) 		:: UBounds 
        != 'yes  evaluation of initial solution'
        LOGICAL 					  		:: EvalInitialX = .TRUE.
!		!= '(4 + floor(3*log(N)))  population size, lambda'
        INTEGER 					 		:: PopSize = 0    
        != 'floor(popsize/2)     popsize equals lambda'  
        INTEGER 							:: ParentNumber = 0
        != 'superlinear decrease(3) or linear(2), or equal(1)'
        INTEGER 							:: RecombinationWeights = 3 
        ! >=0, messaging after every i-th iteration
        INTEGER 							:: VerboseModulo = 100 
        
        ! write output data into text files        
#ifdef __HAVE_MATLAB__
        LOGICAL                             :: flgouttxt = .FALSE.
#else
        LOGICAL                             :: flgouttxt = .TRUE.
#endif        
        ! log data from each generation (C,D,sigma,etc)
        LOGICAL								:: flgGenData = .FALSE.
        ! log data from the trace (bestever.f and bestever.x) in 
        !intervals given by intGenData, this is included in flgGenData
        LOGICAL                             :: flgGenTrace = .FALSE.
        ! generation interval to log data, default is interval = 1
        INTEGER								:: intGenData = 1
        !Objective Function Name
        CHARACTER(len=200)					:: funcName = ''
        ! only one process generates output
        LOGICAL                             :: silentMPI = .TRUE.

        
        !-----------------------------------------------------------------------
        !  Particle Swarm Parameters
        !-----------------------------------------------------------------------
        !Switch PS-CMA on or off
        LOGICAL 							:: pscma = .FALSE.
        !Set weight for coupling C(cma) with C(pso), any value between 0 and 1.
        !Example: psoWeight = 0.3 --> C = 0.3*C(cma) + 0.7*C(pso)
        REAL(MK)							:: psoWeight = 0.0_MK
        !Specify how often a PSO update is performed. For example psoFreq = 5 
        !means that every 5th iteration an update will be done
        INTEGER								:: psoFreq = 200
        
        
        !-----------------------------------------------------------------------
        !  Random Number Generator Parameters
        !-----------------------------------------------------------------------
        
        !Use Quasi Random Sampling on/off
        LOGICAL 							:: qr_sampling = .TRUE.
        !which sampler to use
      !                             0.............Sobol
      !                    (default)1.............Sobol with scrambling 
      !                             2.............Halton
      !                             3.............Halton R implementation
      !                             4.............Faure (buggy!)
      !                             5.............Niederreiter
        INTEGER                            :: qr_sampler = 1

      !inverter (I) which Inverter to use (Uniform -> Normal(0,1))
      !                             -1............dont use any inverter (uniform)
      !                             0.............Moros Inverse
      !                    (default)1             Peter J. Acklam’s Inverter
      !                             2.............Inverter from the R Implementation
        INTEGER                            :: qr_inverter = 1
 			!which type of scrambling to use when using Sobol with scrambling

		!                   0 - NO SCRAMBLING
		!                   1 - OWEN TYPE SCRAMBLING
		!                   2 - FAURE-TEZUKA TYPE SCRAMBLING
		!                   3 - OWEN + FAURE-TEZUKA TYPE SCRAMBLING
        INTEGER                            :: qr_scrambling = 0
        
        !----------------------------------------------------------------------
        ! Restart / Inc Popsize CMA Parameters
        !----------------------------------------------------------------------
        ! If true CMA will restart after converging 
        LOGICAL :: restart_cma = .FALSE.
        ! Type of Restart CMA should make after converging
        INTEGER :: restart_type = 0
        !                         0 restart randomly within bounds (default)
        !                         1 restart from point of convergence
        !                         2 restart from same startpoint all the time
         !% number of restarts '; 
        INTEGER :: Restarts  = 0   ! if set to zero CMA will run till max. Iteration or Funevals are reached
                                   ! 
        !% multiplier for population size before each restart';
        REAL(MK) :: IncPopSize   = 1.25_MK
        !% maximum size multiplier in comparison to the orignal population size (ignoring negative values)
        ! this setting is important because for big problem sizes one easily gets into memory troubles
        INTEGER :: MaxIncFac = 100
        
        
        !-----------------------------------------------------------------------
        !  Benchmark variables
        !  Based on the CEC 2005 Benchmark, the following variables can be
        !  assigned to measure performance
        !-----------------------------------------------------------------------
        
        INTEGER             :: Benchfctnr = 1     ! also used for BBOB
        LOGICAL				:: benchmark = .FALSE. ! logs when certain f-val are reaced etc.
        REAL(MK)			:: global_min = 0.0_MK !0 mean no global min!!!!
        REAL(MK)			:: accuracy = 0.0_MK
        REAL(MK)            :: record_accuracy = 0.0_MK !used to check when optimization reaches gmin with this accuracy
        LOGICAL             :: use_seed = .FALSE.   !read in seed from file
        CHARACTER(len=200)  :: seed_folder = 'false'          !folder containing seed.txt
        LOGICAL             :: record_besthist = .FALSE.
        INTEGER             :: record_modulo = 100  !how many recordings (equal spaced) of the bestever should be made
        CHARACTER(len=200)  :: CECFolders = '' !path of the support Data Folder relative to the working dir (or absolut)
        LOGICAL             :: write_pdb = .FALSE. !if pdb should be written if using LJ/Water
        LOGICAL             :: use_BBOB = .FALSE.
        LOGICAL             :: use_CEC   = .FALSE.
        LOGICAL             :: use_LJ = .FALSE.
        LOGICAL             :: use_TIP = .FALSE.
        LOGICAL             :: use_LJ_comp = .FALSE.
        LOGICAL             :: use_DF = .FALSE.
        LOGICAL             :: use_RANDOM_LANDSCAPE = .FALSE.
        LOGICAL             :: use_MATFUNC = .FALSE.
        
        ! varialbes for the LJ Compression case
        REAL(MK)            :: LJ_comp = 1_MK
        
		! variables for the DoubleFunnel Test case
        REAL(MK)            :: DF_d = 0_MK
        REAL(MK)            :: DF_s = 0_MK
        LOGICAL             :: DF_rast = .TRUE.
        
        
        
        !-----------------------------------------------------------------------
        !  BFGS variables
        !  BFGS is a gradient search to assist CMA
        !-----------------------------------------------------------------------
        
        LOGICAL             :: BFGS_use = .FALSE.
        
        REAL(MK)            :: BFGS_grad_stepsize = 10d-5
        LOGICAL             :: BFGS_central_difference = .FALSE.
        INTEGER(MK)         :: BFGS_position = 2
                               ! 1 = replace X Values by Local Minimum X
                               ! 2 = replace F Values with F Values at Local Minimum
        
        
        ! BFGS_confprop is the upper limit, on how far BFGS is allowed to travel within the searchspace
        ! relative to the current density of the multivariant gauss of CMA
        ! adjusted in terms of probability 
        ! (e.g. 0.95 means that BFGS is only allowed to travel within the 0.95% probability volumn of the multivariant gauss)
        REAL(MK)            :: BFGS_confprop = 0.95d0
        ! BFGS_sigmastep_min is the lower limit, on how far BFGS has to travel within each iteration relative
        ! to the current sigma
        REAL(MK)            :: BFGS_sigmastep_min = 0.1d0
                               
        REAL(MK)            :: BFGS_factr  = 0_MK !if zero StopTolFun/machine eps is used
        REAL(MK)            :: BFGS_pgtol  = 1.0d-5
        
        INTEGER            :: BFGS_print = -1
        !It controls the frequency and type of output generated:
        !c        iprint<0    no output is generated;
        !c        iprint=0    print only one line at the last iteration;
        !c        0<iprint<99 print also f and |proj g| every iprint iterations;
        !c        iprint=99   print details of every iteration except n-vectors;
        !c        iprint=100  print also the changes of active set and final x;
        !c        iprint>100  print details of every iteration including x and g;
        !c       When iprint > 0, the file iterate.dat will be created to
        !c                        summarize the iteration..
        INTEGER             :: BFGS_maxiter = 300
	   
	   
        
                
        !-----------------------------------------------------------------------
        !  internal used parameters (calculated from the parameters set above)
        !-----------------------------------------------------------------------
        REAL(MK),POINTER,DIMENSION(:) :: xstart
        !the final used sigma
        REAL(MK),POINTER,DIMENSION(:)   :: insigma
        INTEGER                         :: seed
        
        
        
        
        
        
        
        

        ! The following lines state variables from cmaes.m that haven't been
        ! implemented (yet)
        
        !% objective function FUN accepts NxM matrix, with M>1?';
        !INTEGER :: EvalParallel != 'no  
       
        !% display messages like initial and final message';
        !INTEGER :: Display  != 'on   
        !
        !INTEGER :: Plotting != 'on   % plot while running';
        !% resume former run from SaveFile';
        !INTEGER :: Resume  != 'no     
        !% off==do some additional (minor) problem capturing';
        !INTEGER :: Science != 'off  
        !% [on|final|off][-v6] save data to file';
        !INTEGER :: Saving  !=    'on   
        ! % if >1 record data less frequently after gen=100';
        !INTEGER :: SaveModulo != '1   
        !% max. percentage of time for recording data';
        !INTEGER :: SaveTime   != '25   
        
        
      END TYPE opts
      
      TYPE(opts) :: options
      


      
      END MODULE cmaes_opts_mod
