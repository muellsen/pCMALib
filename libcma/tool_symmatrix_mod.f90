#include "define.h"
 	  !-------------------------------------------------------------------------
      !  Module       :                   tool_symmatrix_mod
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module handles the call of the tool_symmatrix
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
	  MODULE tool_symmatrix_mod
	  
	  INTERFACE tool_symmatrix
	  	MODULE PROCEDURE tool_symmatrix_s
	  	MODULE PROCEDURE tool_symmatrix_d
	  END INTERFACE
	  
	  
	  CONTAINS
#define __KIND __SINGLE_PRECISION
#include "tool_symmatrix.f90"
#undef __KIND

#define __KIND __DOUBLE_PRECISION	  
#include "tool_symmatrix.f90"	  
#undef __KIND	  
	  
	  
	  END MODULE tool_symmatrix_mod