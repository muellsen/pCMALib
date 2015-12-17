 	  !-------------------------------------------------------------------------
      !  Module       :                   cmaes_run_mod
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains variables used in cmaes_run()
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
      MODULE cmaes_run_mod
      USE cmaes_param_mod
	  USE cmaes_opts_mod
      USE cmaes_mod
      USE cmaes_out_mod
      USE tool_mrgrnk_mod
      USE tool_symmatrix_mod
#ifdef __HAVE_MPI__    
      USE mpi_mod
      USE pso_mod
#endif
      USE qr_gen_mod
      	
      	
	  !-------------------------------------------------------------------------
	  !  Module Procedures
	  !-------------------------------------------------------------------------
      CONTAINS
#include "cmaes_handlebounds.f90"
#include "cmaes_writegen.f90"

      END MODULE cmaes_run_mod