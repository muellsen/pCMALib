#include "define.h"
 	  !-------------------------------------------------------------------------
      !  Module       :                   pso_mod
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains variables used for Particle Swarm
      !					Update
      !                
      !  Remarks      : Variables are initialized in Subroutine pso_init()
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

      MODULE pso_mod
      USE cmaes_param_mod
	  IMPLICIT NONE
	  SAVE
#ifdef __HAVE_MPI__
	  !-------------------------------------------------------------------------
	  !  Variables
	  !-------------------------------------------------------------------------	  
	  REAL(MK),ALLOCATABLE,DIMENSION(:,:)			:: C_PSO
	  REAL(MK),ALLOCATABLE,DIMENSION(:,:)			:: B_rot
	  REAL(MK),ALLOCATABLE,DIMENSION(:)				:: b_main
	  REAL(MK),ALLOCATABLE,DIMENSION(:)				:: bias
	  REAL(MK),ALLOCATABLE,DIMENSION(:)				:: gb_vec
	  REAL(MK),ALLOCATABLE,DIMENSION(:)				:: gb_vec_xold
	  REAL(MK),ALLOCATABLE,DIMENSION(:)				:: local_best_vec


	  !-------------------------------------------------------------------------
	  !  Module Procedures
	  !-------------------------------------------------------------------------
	  CONTAINS
	  !-------------------------------------------------------------------------
      !  Subroutine   :          pso_init
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This Subroutine initializes the variables of pso_mod
      !
      !  Input		  : 
	  !
      !  Remarks      : 
      !
      !  References   : 
	  !					
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: pso_init.f90,v $
      !  Revision 1.4  2008/04/18 13:02:56  paulb
      !  repository update
      !
      !  Revision 1.1  2008/03/14 18:49:03  paulb
      !  Corrected wrong datatype in MPI_BCAST, split up variables of cmaesRun 
      !	 into modules to keep code readable
      !
      !
      !-------------------------------------------------------------------------
      !  Covariance Matrix Adaption Library (LIBCMA)
      !  Benedikt Baumgartner
      !  Computational Biophysics Lab, ETH Zurich, Switzerland
      !-------------------------------------------------------------------------
      SUBROUTINE pso_init(N)
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Parameters
      !-------------------------------------------------------------------------
      INTEGER, INTENT(in)							:: N
      
      !-------------------------------------------------------------------------
      !  Local Variables
      !-------------------------------------------------------------------------
      INTEGER										:: allocStat
      
      
      ALLOCATE(C_PSO(N,N),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating C_PSO'
      C_PSO = 0.0_MK
      
      ALLOCATE(B_rot(N,N),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating B_rot'
      B_rot = 0.0_MK
      
      ALLOCATE(b_main(N),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating b_main'
      b_main = 0.0_MK
      
      ALLOCATE(bias(N),stat=allocStat)      
      IF(allocStat .NE. 0) STOP 'Error allocating bias'
      bias = 0.0_MK
      
      ALLOCATE(gb_vec(N),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating gb_vec'
      gb_vec = 0.0_MK
      
      ALLOCATE(gb_vec_xold(N),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating gb_vec_xold'
      gb_vec_xold = 0.0_MK
      
      ALLOCATE(local_best_vec(N),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating local_best_vec'
      local_best_vec = 0.0_MK
      
      RETURN
      END SUBROUTINE pso_init
!#include "psoUpdate.f90"
 	  !-------------------------------------------------------------------------
      !  Subroutine   :                   cmaes_reinitialize
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Reinitializes CMA run at new position xmean_new
      !
      !  Input		  :	xmean_new(:)		(R)	array of length N. Specifies the
      !											reinitialization position
      !					N					(I)	Length of xmean_new
      !
      !	 Input/Output :
      !                
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: cmaes_reinitialize.f90,v $
      !  Revision 1.3  2008/04/18 13:02:56  paulb
      !  repository update
      !
      !
      !-------------------------------------------------------------------------
      !  pCMALib: a parallel fortran 90 library for the evolution strategy with
      !           covariance matrix adaptation
      !  Christian L. Mueller, Benedikt Baumgartner, Georg Ofenbeck, Christian L. Mueller
      !  MOSAIC group, ETH Zurich, Switzerland
      !-------------------------------------------------------------------------

	  SUBROUTINE cmaes_reinitialize(xmean_new,N)
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE cmaes_param_mod
      USE cmaes_mod
      USE cmaes_opts_mod
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Parameters
      !-------------------------------------------------------------------------
      INTEGER,INTENT(in)							:: N
      REAL(MK),DIMENSION(N),INTENT(in)				:: xmean_new
      
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER										:: i
      
      ! Reinitialize CMA variables
      pc = 0.0_MK
      ps = 0.0_MK
      C = 0.0_MK
      B = 0.0_MK
      D = 0.0_MK
      BD = 0.0_MK
      DO i = 1,N
      	C(i,i) = 1._MK
      	B(i,i) = 1._MK
      	D(i,i) = 1._MK
      	BD(i,i) = 1._MK
      END DO 
      
      xmean = xmean_new
      xold = xmean
      cc = 4./(REAL(N)+4.)
      cs = (mueff+2)/(N+mueff+3)
      mucov = mueff
      ccov = (1./mucov) * 2/(N+1.41)**2 + (1 - 1/mucov) * &
      		min(1.,(2*mueff-1)/((N+2)**2 + mueff))
      damps = (1 + 2*max(0.,sqrt((mueff-1.)/real(N+1))-1.)) * &
      		max(0.3,1 - N/min(real(options%StopMaxIter), &
      		options%StopMaxFunEvals/lambda)) + cs
      sigma = maxval(input%insigma)
      
      RETURN
      END SUBROUTINE  	
#endif      
      END MODULE pso_mod
