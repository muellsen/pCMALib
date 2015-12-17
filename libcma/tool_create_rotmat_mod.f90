#include "define.h"
 	  !-------------------------------------------------------------------------
      !  Module       :                   tool_create_rotmat_mod
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module handles the call of the tool_create_rotmat
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
	  MODULE tool_create_rotmat_mod
	  
	  INTERFACE tool_create_rotmat
	  	MODULE PROCEDURE tool_create_rotmat_s
	  	MODULE PROCEDURE tool_create_rotmat_d
	  END INTERFACE
	  
	  
	  CONTAINS
#define __KIND __SINGLE_PRECISION
#include "tool_create_rotmat.f90"
#undef __KIND

#define __KIND __DOUBLE_PRECISION	  
#include "tool_create_rotmat.f90"	  
#undef __KIND	  
	  
	  
	  END MODULE tool_create_rotmat_mod