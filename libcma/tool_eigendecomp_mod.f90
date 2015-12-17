#include "define.h"
 	  !-------------------------------------------------------------------------
      !  Module       :                   tool_eigendecomp_mod
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module handles the call of the tool_eigendecomp
      !					Subroutine for single or double precision
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
	  MODULE tool_eigendecomp_mod
	  
	  INTERFACE tool_eigendecomp
	  	MODULE PROCEDURE tool_eigendecomp_s
	  	MODULE PROCEDURE tool_eigendecomp_d
	  END INTERFACE
	  
	  
	  CONTAINS
#define __KIND __SINGLE_PRECISION
#include "tool_eigendecomp.f90"
#undef __KIND

#define __KIND __DOUBLE_PRECISION	  
#include "tool_eigendecomp.f90"	  
#undef __KIND	  
	  
	  
	  END MODULE tool_eigendecomp_mod
	  
	  