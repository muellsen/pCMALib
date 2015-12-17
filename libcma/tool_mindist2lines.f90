	  !-------------------------------------------------------------------------
      !  Subroutine   :                   tool_mindist2lines
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This Routine calculates the minimum distance between
      !					two lines in n-dimensional space. Return values are the
      !					minimum distance and the corresponding points on the two
      !					lines.
      !
      !  Input		  :	N		(I)	Number of dimensions
      !					X0(:)	(R) Starting point of first vector
      !					X(:)	(R) First vector
      !					Y0(:)	(R) Starting point of second vector
      !					Y(:)	(R) Second vector
      !
      !	 Output 	  : DIST	(R) minimum distance between X and Y
      !					XMIN(:)	(R) corresponding coordinates of DIST on X
      !					YMIN(:)	(R) corresponding coordinates of DIST on Y
      !                
      !  Remarks      : DIST = ||XMIN - YMIN|| with XMIN = X0 + Xt and 
      !					YMIN = Y0 + Ys,	s and t being scalars
      !
      !	 References	  :	@article{Bard:2001,
	  !					Author = {Bard, Michael; Himel, Denny},
	  !					Keywords = {Minimum Distance},
	  !					Title = {{The Minimum Distance between Two Lines in
	  !								n-Space}},
	  !					Year = {2001}}
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  pCMALib: a parallel fortran 90 library for the evolution strategy with
      !           covariance matrix adaptation
      !  Christian L. Mueller, Benedikt Baumgartner, Georg Ofenbeck
      !  MOSAIC group, ETH Zurich, Switzerland
      !-------------------------------------------------------------------------
      SUBROUTINE tool_mindist2lines(N,X0,X,Y0,Y,DIST,XMIN,YMIN)
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE cmaes_param_mod
      IMPLICIT NONE
      
	  !-------------------------------------------------------------------------
      !  Parameters
      !-------------------------------------------------------------------------
      INTEGER,INTENT(in)						:: N
      REAL(MK),DIMENSION(N),INTENT(in)			:: X0,X,Y0,Y
      REAL(MK),INTENT(out)						:: DIST
      REAL(MK),DIMENSION(N),INTENT(out)			:: Ymin,Xmin						
      
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK)									:: A,B,C,D,E,F
      REAL(MK)									:: s,t
      
      
      A = sum(X**2)
	  B = 2.0_MK*(dot_product(X,X0) - dot_product(X,Y0))
	  C = 2.0_MK*(dot_product(X,Y))
	  D = 2.0_MK*(dot_product(Y,Y0) - dot_product(Y,X0))
	  E = sum(Y**2)
	  F = sum(X0**2) + sum(Y0**2)
	  
	  ! A unique solution only exists if
	  ! 	| 2A  -C |
	  ! det |		 | != 0
	  !		| -C  2E |
	  ! If no unique solution exists, xmin=x0, ymin=y0, dist = ||x0-y0||
	  IF(((2*A + 2*E) - ((-C)*(-C))) .EQ. 0.0_MK) THEN
	  	!WRITE(*,*) MY_RANK,': Warning: no unique distance in tool_mindist2lines()'
	  	XMIN = X0
	  	YMIN = Y0
	  	DIST = sqrt(sum((X0 - Y0)**2))
	  	RETURN
	  END IF
	  
	  s = (2.0_MK*A*D + B*C) / (C**2 - 4.0_MK*A*E)
	  t = (C*s - B) / 2.0_MK*A
	  
	  DIST = sqrt(( (B*C*D + B**2*E + C**2*F + A*(D**2-4.0_MK*E*F)) / &
	  													(C**2-4.0_MK*A*E))**2)

	  XMIN = X0 + X*t
	  YMIN = Y0 + Y*s	  													
      !DIST = sqrt(sum((YMIN-XMIN)**2))
      
      RETURN
      END SUBROUTINE