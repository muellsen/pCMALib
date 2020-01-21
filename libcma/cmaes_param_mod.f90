      !-------------------------------------------------------------------------
      !  Module       :                   cmaes_param_mod
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains parameters used in CMA
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
      MODULE cmaes_param_mod
      IMPLICIT NONE
      SAVE
#ifdef __SP
	  INTEGER,PARAMETER		:: MK = KIND(1.E0)	  
#else
	  INTEGER,PARAMETER		:: MK = KIND(1.D0)	  
#endif	  	        
        
        
      !------------------------------------------------------------------
      ! Notice - the Random Number Part uses its own kind types since most
      ! functions are provided for both precision types
      ! Can be found in kinds.f90 within the librng folder
      !------------------------------------------------------------------

	  
#ifdef __HAVE_MPI__      
      INTEGER				:: MY_RANK
      INTEGER				:: NUM_CMA_RUNS
#endif


      END MODULE cmaes_param_mod
