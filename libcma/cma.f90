      !-------------------------------------------------------------------------
      !  Program	  :          cma
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This is the main entry point for pcmalib
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
      !  Christian L. Mueller, Benedikt Baumgartner, Georg Ofenbeck, Christian L. Mueller
      !  MOSAIC group, ETH Zurich, Switzerland
      !-------------------------------------------------------------------------
      
      PROGRAM cma
      USE cmaes_mod
      USE cmaes_param_mod
      USE cmaes_opts_mod
      USE CEC2005
      IMPLICIT NONE
#ifdef __HAVE_MPI__      
      INCLUDE "mpif.h"
#endif      
      INTERFACE cmaes_start
      	SUBROUTINE cmaes_start(fitfun,xstart,insigma,varargin)
	      USE cmaes_param_mod
    	  USE cmaes_mod, ONLY : input
    	  IMPLICIT NONE
	      EXTERNAL							:: fitfun	      
    	  REAL(MK), INTENT(in) 				:: xstart(:)
    	  REAL(MK), INTENT(in)				:: insigma(:)
      	  REAL(MK), INTENT(in), OPTIONAL 	:: varargin
      	END SUBROUTINE cmaes_start
      END INTERFACE cmaes_start
   
      INTEGER 				:: ierr,i
!      REAL(MK),DIMENSION(d) :: xstart
      
      CHARACTER(len=10)		:: proc
	  ! Random number generator
	  INTEGER				:: seed
      REAL(MK)				:: ZBQLUAB
	  CHARACTER(len=256)    :: ctrlfile
      
      
#ifdef __BBOB
      EXTERNAL                          :: benchmark_bbob
#endif      
      EXTERNAL                          :: LJ
      EXTERNAL                          :: TIP4P
      EXTERNAL                          :: LJ_POT_COMP
      EXTERNAL                          :: DoubleFunnel
      EXTERNAL                          :: random_landscape
      EXTERNAL                          :: user_function
#ifdef  __HAVE_MATLAB__
      EXTERNAL                          :: matlab_objfunc
#endif
      CALL GETARG(1 , ctrlfile)
 !     WRITE(*,*) ctrlfile

   



#ifdef __HAVE_MPI__
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MY_RANK, ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NUM_CMA_RUNS, ierr)
      WRITE(proc,'(I10)') MY_RANK
#endif
      
      CALL CPU_TIME(time_start)
      CALL cmaes_readparams(ctrlfile,ierr)



      ! Start CMA-ES/PS-CMA-ES
      IF (options%use_TIP ) THEN
        CALL cmaes_start(TIP4P,options%xstart,options%insigma)
      ELSEIF (options%use_LJ ) THEN      
        CALL cmaes_start(LJ,options%xstart,options%insigma)
      ELSEIF (options%use_LJ_comp ) THEN
        CALL cmaes_start(LJ_POT_COMP,options%xstart,options%insigma)
      ELSEIF (options%use_DF ) THEN
        CALL cmaes_start(DoubleFunnel,options%xstart,options%insigma)        
      ELSEIF (options%use_CEC ) THEN 
	    CALL cmaes_start(benchmark,options%xstart,options%insigma)
	  ELSEIF (options%use_BBOB ) THEN
#ifdef __BBOB	  
	    CALL cmaes_start(benchmark_bbob,options%xstart,options%insigma)
#else
      WRITE(*,*) 'You compiled CMA without BBOB and therefore cant call it'	    
#endif
      ELSEIF (options%use_MATFUNC ) THEN
#ifdef __HAVE_MATLAB__ 
      CALL cmaes_start(matlab_objfunc,options%xstart,options%insigma)
#else
      WRITE(*,*) 'You compiled CMA without MATLAB support and therefore cant call it'
#endif
	    
	  ELSEIF (options%use_RANDOM_LANDSCAPE ) THEN
	    CALL ZBQLINI(1)
	    CALL cmaes_start(random_landscape,options%xstart,options%insigma) !this is mainly to test CMA-ES behaviour
	  ELSE	  	    
	    CALL cmaes_start(user_function,options%xstart,options%insigma)
	    !   will be called with 
	    !   CALL fitfun(f,x(N),N,1,options%LBounds,options%UBounds)
	  END IF

	  
#ifdef __HAVE_MPI__
	  CALL MPI_FINALIZE(ierr)
#endif

	  
      END PROGRAM cma
