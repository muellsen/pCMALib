	  !-------------------------------------------------------------------------
      !  Function	  :          cmaes_funcwrap
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This is a wrapper function
      !
      !  Input		  : x(:)	(R)	Vector to evaluate fitfun for
      !                 N		(I) Dimension of x(:)
      !					fitfun		External Subroutine defined in cma.f90
      !
      !	 Output		  : Function value at x(:)
      !
      !  Remarks      : 
      !
      !  References   : cma.f90
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  pCMALib: a parallel fortran 90 library for the evolution strategy with
      !           covariance matrix adaptation
      !  Christian L. Mueller, Benedikt Baumgartner, Georg Ofenbeck
      !  MOSAIC group, ETH Zurich, Switzerland
      !-------------------------------------------------------------------------
#ifdef __SP
	  REAL FUNCTION cmaes_funcwrap(x,N,fitfun)
#else      
	  DOUBLE PRECISION FUNCTION cmaes_funcwrap(x,N,fitfun)
#endif	  
	  USE cmaes_param_mod
	  USE cmaes_opts_mod
	  USE cmaes_mod
	  IMPLICIT NONE
	  !-------------------------------------------------------------------------
	  !  Parameters
	  !-------------------------------------------------------------------------
	  INTEGER,INTENT(in)					:: N
	  REAL(MK),DIMENSION(N),INTENT(in)		:: x
	  !-------------------------------------------------------------------------
	  !  Externals
	  !-------------------------------------------------------------------------
	  EXTERNAL							:: fitfun
	  !-------------------------------------------------------------------------
	  !  Local Variables
	  !-------------------------------------------------------------------------
	  REAL(MK),DIMENSION(1)				    :: temp1
	  REAL(MK),DIMENSION(N,1)				:: temp2
	  
	  temp2(:,1) = x
	  temp1(1) = 0.0_MK
	  
	  countEval = countEval + 1
	  CALL fitfun(temp1,temp2,N,1,options%LBounds,options%UBounds)
	  !CALL fitfun(temp1,temp2,N,1)
	  cmaes_funcwrap = temp1(1)
	  RETURN
	  END FUNCTION