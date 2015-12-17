 	  !-------------------------------------------------------------------------
      !  Module       :                   cmaes_mod
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains structures and variables 
      !					used in CMA
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
      MODULE cmaes_mod
      	USE cmaes_param_mod
      	USE tool_eigendecomp_mod 
      	IMPLICIT NONE
      	SAVE
      	
      	!-----------------------------------------------------------------------
      	!  CMA input parameters
      	!-----------------------------------------------------------------------
      	TYPE :: inp
      		REAL(MK),ALLOCATABLE,DIMENSION(:) :: xstart
      		REAL(MK),ALLOCATABLE,DIMENSION(:) :: insigma
      		INTEGER						  :: N
      	END TYPE inp
      	!  Create instance of type inp to use in routines
      	TYPE(inp) :: input
      	
      	!-----------------------------------------------------------------------
      	!  CMA counting variables
      	!-----------------------------------------------------------------------
      	INTEGER :: countEval = 0
      	INTEGER :: countEvalNaN = 0
      	INTEGER :: countIter = 0
      	INTEGER	:: countOutOfBounds = 0
      	INTEGER :: lastRestart = 0
      	INTEGER :: countRestarts = 0
      	INTEGER :: initial_Popsize = 0
      	INTEGER :: countBFGSEval = 0
      	!this is used to output the reason for BFGS iteration exits at the end an CMA run
	    INTEGER,DIMENSION(6) :: BFGSexitcounter = 0
	    !(1) regular exits
	    !(2) regular exits cause of pgtol
	    !(3) error exits
	    !(4) BFGS_sigmastep_min exits
	    !(5) Confprop exit
	    !(6) MaxIter exit
      	!-----------------------------------------------------------------------
      	!  CMA timing
      	!-----------------------------------------------------------------------
      	
      	REAL(MK) :: time_start
      	REAL(MK) :: time_end
      	
      	!-----------------------------------------------------------------------
      	!  Miscellaneous variables
      	!-----------------------------------------------------------------------
      	REAL(MK),ALLOCATABLE,DIMENSION(:) :: xmean
      	REAL(MK),ALLOCATABLE,DIMENSION(:) :: xold
      	REAL(MK),ALLOCATABLE,DIMENSION(:) :: maxdx
      	REAL(MK),ALLOCATABLE,DIMENSION(:) :: mindx
      	REAL(MK)						  :: fac
      	REAL(MK)						  :: fmean
      	CHARACTER(len=20)			  :: stopflag
      	
		!expectation of ||N(0,I)|| == norm(randn(N,1))
      	REAL(MK)						  :: chiN
      	
      	!-----------------------------------------------------------------------
      	!  Strategy parameter setting: Selection
      	!-----------------------------------------------------------------------
      	INTEGER							:: lambda = 0
      	INTEGER							:: mu = 0
      	REAL(MK),ALLOCATABLE,DIMENSION(:) 	:: weights
      	REAL(MK)							:: mueff
      	
      	!-----------------------------------------------------------------------
      	!  Strategy parameter setting: Adaption
      	!-----------------------------------------------------------------------
      	! Time constant for cumulation of covariance matrix
      	REAL(MK)							:: cc 
      	! Time constant for cumulation of step size control
      	REAL(MK)							:: cs
      	! Size of mu used for calculating learning rate ccov
      	REAL(MK)							:: mucov
      	! Learning rate for covariance matrix
      	REAL(MK)							:: ccov
      	! Damping for step size control, usually close to one
      	REAL(MK)							:: damps
      	! Overall standard deviation
      	REAL(MK) 							:: sigma
      	! Evolution paths for C and sigma 
      	REAL(MK),ALLOCATABLE,DIMENSION(:)	:: pc,ps
      	! B defines the coordinate system
      	REAL(MK),ALLOCATABLE,DIMENSION(:,:) :: B
      	! Diagonal matrix D defines the scaling
      	REAL(MK),ALLOCATABLE,DIMENSION(:,:) :: D
      	! For speed up only
      	REAL(MK),ALLOCATABLE,DIMENSION(:,:) :: BD
      	! Covariance matrix
      	REAL(MK),ALLOCATABLE,DIMENSION(:,:) :: C
      	REAL(MK),ALLOCATABLE,DIMENSION(:,:) :: invsC
      	
      	
      	REAL(MK)                            :: chi
      	!REAL(MK),ALLOCATABLE,DIMENSION(:,:) :: C_PSO
!REAL(MK),DIMENSION(2)						:: gb_vec
		
		!-----------------------------------------------------------------------
      	!  Fitness values
      	!-----------------------------------------------------------------------
      	TYPE :: fit
			REAL(MK),ALLOCATABLE,DIMENSION(:)		:: hist
			REAL(MK),ALLOCATABLE,DIMENSION(:)		:: histsel
			REAL(MK),ALLOCATABLE,DIMENSION(:)		:: raw
			REAL(MK),ALLOCATABLE,DIMENSION(:)		:: sel
			INTEGER,ALLOCATABLE,DIMENSION(:)		:: idx
			INTEGER,ALLOCATABLE,DIMENSION(:)		:: idxsel
		END TYPE fit
		TYPE(fit) :: fitness
		
		!-----------------------------------------------------------------------
      	!  Boundary handling
      	!-----------------------------------------------------------------------
		TYPE :: bound
			LOGICAL 							:: isactive = .FALSE.
			REAL(MK),ALLOCATABLE,DIMENSION(:)		:: weights
			REAL(MK),ALLOCATABLE,DIMENSION(:)		:: scales
			LOGICAL,ALLOCATABLE,DIMENSION(:)	:: isbounded
			REAL(MK),ALLOCATABLE,DIMENSION(:)		:: dfithist
			REAL(MK),ALLOCATABLE,DIMENSION(:)		:: arpenalty
			LOGICAL								:: iniphase		
			INTEGER								:: validfitval	
		END TYPE bound
		TYPE(bound) :: bnd
		
		!-----------------------------------------------------------------------
      	!  Module Procedures
      	!-----------------------------------------------------------------------
!		CONTAINS
!#include "cmaes_xintobounds.f90"
		
      END MODULE cmaes_mod
