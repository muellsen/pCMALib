 	  !-------------------------------------------------------------------------
      !  Module       :                   cmaes_out_mod
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains routines and variables
      !					for CMA output
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
	  MODULE cmaes_out_mod
	  USE cmaes_param_mod
	  IMPLICIT NONE
	  SAVE
	  
	    !-----------------------------------------------------------------------
      	!  bestever solution struct (accessed via solution struct)
      	!-----------------------------------------------------------------------	
		TYPE best
			REAL(MK),ALLOCATABLE,DIMENSION(:)	:: x
			REAL(MK)							:: f
			INTEGER								:: evals
		END TYPE best
		
		!-----------------------------------------------------------------------
      	!  Mean struct (accessed via history struct)
      	!-----------------------------------------------------------------------	
		TYPE mean
			REAL(MK),ALLOCATABLE,DIMENSION(:)	:: x
			REAL(MK)							:: f
			INTEGER								:: evals
		END TYPE mean
		
		!-----------------------------------------------------------------------
      	!  Recentbest struct (accessed via history struct)
      	!-----------------------------------------------------------------------	
		TYPE recBest
			REAL(MK),ALLOCATABLE,DIMENSION(:)	:: x
			REAL(MK)							:: f
			INTEGER								:: evals
		END TYPE recBest
		
		!-----------------------------------------------------------------------
      	!  Recentworst struct (accessed via history struct)
      	!-----------------------------------------------------------------------	
		TYPE recWorst
			REAL(MK),ALLOCATABLE,DIMENSION(:)	:: x
			REAL(MK)							:: f
			INTEGER								:: evals
		END TYPE recWorst
		
		!-----------------------------------------------------------------------
      	!  Single parameters struct (accessed via history struct)
      	!-----------------------------------------------------------------------	
		TYPE param
			INTEGER								:: evals
			INTEGER								:: iterations
			REAL(MK)							:: sigma
			REAL(MK)							:: maxstd
			REAL(MK)							:: minstd
			REAL(MK)							:: maxD
			REAL(MK)							:: minD
			INTEGER								:: maxstdidx
			INTEGER								:: minstdidx
			!Comment see line 602
		END TYPE param
		
		!-----------------------------------------------------------------------
      	!  History struct (accessed via output struct)
      	!-----------------------------------------------------------------------	
		TYPE history
			TYPE(recBest) 	:: recentbest
			TYPE(recWorst) 	:: recentworst
			TYPE(mean) 		:: mean
			TYPE(param) 	:: param
			INTEGER 		:: evals
			INTEGER 		:: iterations
		END TYPE history

		!-----------------------------------------------------------------------
      	!  Parameter array struct (accessed via output struct)
      	!-----------------------------------------------------------------------				
		TYPE paramArr
			INTEGER 							:: evals
			INTEGER 							:: iterations
			REAl(MK)							:: sigma
			REAL(MK),ALLOCATABLE,DIMENSION(:)	:: diagD
			REAL(MK),ALLOCATABLE,DIMENSION(:)	:: stds
			REAL(MK),ALLOCATABLE,DIMENSION(:)	:: Bmax
			REAL(MK),ALLOCATABLE,DIMENSION(:)	:: Bmin
			! comment line 614
		END TYPE paramArr
		
		!-----------------------------------------------------------------------
      	!  Soultion struct (accessed via output struct)
      	!-----------------------------------------------------------------------	
		TYPE sol
			INTEGER	   		:: evals
			TYPE(best) 		:: bestever
			TYPE(mean) 		:: mean
			TYPE(recBest)	:: recentbest
			TYPE(recWorst)	:: recentworst
		END TYPE sol
		
		!-----------------------------------------------------------------------
      	!  Output struct
      	!-----------------------------------------------------------------------	
		TYPE :: output
			TYPE(sol) 		:: solutions
			TYPE(history) 	:: hist
			TYPE(paramArr) 	:: histParamArr
			INTEGER	  		:: evals
			REAL(MK)	  	:: stopflag
		END TYPE output
		
		!-----------------------------------------------------------------------
		!  The following variables are directly accessible to program routines
		!-----------------------------------------------------------------------
		TYPE(output) :: out ! gives access to all the other variables defined
							! in this module
		TYPE(best)	 :: bestever
		INTEGER		 :: outiter

        !-----------------------------------------------------------------------
      	!  Variables used to save information to be written into .mat files
      	!----------------------------------------------------------------------
        REAL(MK)     :: be3 = -1d0 
        REAL(MK)     :: be4 = -1d0
        REAL(MK)     :: be5 = -1d0
        REAL(MK),ALLOCATABLE,DIMENSION(:)     :: besthist
        
        INTEGER     :: acc_evals = -1
		  
	  END MODULE cmaes_out_mod