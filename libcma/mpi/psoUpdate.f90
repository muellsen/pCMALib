 	  !-------------------------------------------------------------------------
      !  Subroutine   :                   psoUpdate
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Perform Particle Swarm Update in current CMA generation.
      !					Biases xmean and rotates Covariance Matrix.
      !
      !  Input		  :	N				(I) Number of dimensions
      !					local_best(:)	(R) bestever position of current process
      !
      !	 Input/Output : triuC(:,:)		(R) upper triangular part of C
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

	  SUBROUTINE psoUpdate(N,triuC,local_best)
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE tool_create_rotmat_mod
      USE cmaes_mod
      USE cmaes_opts_mod
      USE tool_symmatrix_mod
      USE tool_eigendecomp_mod
      USE mpi_mod
      USE pso_mod
      IMPLICIT NONE
	  INCLUDE "mpif.h"
	  
	  !-------------------------------------------------------------------------
      !  Interfaces
      !-------------------------------------------------------------------------
      INTERFACE
	      SUBROUTINE xIntoBounds(LBounds,UBounds,N,x,xout,idx)
	      USE cmaes_param_mod
	      USE cmaes_mod,only:countOutOfBounds
	      IMPLICIT NONE
		  REAL(MK),DIMENSION(N),INTENT(in)    				:: LBounds,UBounds
		  INTEGER, INTENT(in)			  					:: N
		  REAL(MK),DIMENSION(N),INTENT(in) 					:: x
		  REAL(MK),DIMENSION(N),INTENT(out)					:: xout
		  LOGICAL,DIMENSION(N),INTENT(out),OPTIONAL			:: idx
		  END SUBROUTINE
	  END INTERFACE	  

	  !-------------------------------------------------------------------------
      !  Parameters
      !-------------------------------------------------------------------------
	  INTEGER, INTENT(in)					:: N
	  REAL(MK),DIMENSION(N),INTENT(in)		:: local_best
	  REAL(MK),DIMENSION(N,N),INTENT(inout)	:: triuC

	  !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
	  REAL(MK),DIMENSION(N,N)				:: R
	  INTEGER								:: i
	  REAL(MK),DIMENSION(N) 				:: multmp_vN
	  REAL(MK),DIMENSION(N,N)				:: multmp_mN
	  REAL(MK)								:: gb_vec_length,bias_weight
	  REAL(MK)								:: DUMMY,r1,r2
	  REAL(MK)								:: ZBQLU01
	  REAL(MK),DIMENSION(N)					:: xtmp
	  REAL(MK),DIMENSION(N,N)				:: D_square
	  

      	!-----------------------------------------------------------------------
		! Bias xmean
		!-----------------------------------------------------------------------
		bias_weight = options%psoWeight
		gb_vec = GLOBAL_X_BEST-xmean
		gb_vec_length = sqrt(sum(gb_vec**2))
		local_best_vec = local_best-xmean
		bias = 0.0_MK
		xtmp = 0.0_MK
		DUMMY = 0.0_MK
		multmp_vN = 0.0_MK
		multmp_mN = 0.0_MK
		r1 = ZBQLU01(DUMMY)
		r2 = ZBQLU01(DUMMY)	
		

		!IF(sqrt(sum(gb_vec**2)) .GE. 2._MK*sigma) THEN
!		IF(gb_vec_length .GE. 0.0_MK) THEN
!			IF(sigma .LT. gb_vec_length) THEN
!				!gb_vec = gb_vec/gb_vec_length ! norm of gb_vec
!				bias = sigma*gb_vec/gb_vec_length
!				IF(sigma/gb_vec_length .LE. 0.3_MK*gb_vec_length) THEN
!					bias = 0.5_MK*gb_vec	
!				END IF
!			!
!			ELSE
!				!IF(sigma .GE. gb_vec_length) THEN
!					bias = gb_vec
!				!ELSE
!				!	bias = sigma*gb_vec/gb_vec_length
!			END IF		
!			!END IF					
!			xmean = xmean + bias
!		END IF
	 	
	 	
	 	IF(gb_vec_length .GT. 0.0_MK) THEN
	 		IF(sigma .LT. gb_vec_length) THEN
	 			IF(sigma/gb_vec_length .LE. 0.1_MK*gb_vec_length) THEN
	 			!IF(D(N,N) .LE. 1.E-4) THEN
	 			    !WRITE(*,*) 'Process ', MY_RANK,'with dist',gb_vec_length,' and ', sigma,' applies bias'
	 				bias = 0.5_MK*gb_vec
	 				!bias = 3.0_MK*sigma*gb_vec/gb_vec_length
	 			ELSE
	 			!WRITE(*,*) MY_RANK,'sigmabias'	
	 				bias = sigma*gb_vec/gb_vec_length
	 			END IF	
	 		ELSE
	 			bias = 0.0_MK!r1*gb_vec
	 		END IF
	 	END IF		
		
		!bias = 2.0_MK*r1*(local_best_vec) + 2.0_MK*r2*(gb_vec)
		xmean = xmean + bias
	
		!-----------------------------------------------------------------------
		! Rotate Eigenvectors and construct C_PSO if psoWeight < 1
		!-----------------------------------------------------------------------

		IF (options%psoWeight .LT. 1.0_MK) THEN

		b_main = B(:,N)
		gb_vec_xold = GLOBAL_X_BEST-xold
	
		CALL tool_create_rotmat(b_main,gb_vec_xold,N,R)
		!R = 0.0_MK
        !DO i = 1,N
      	!  R(i,i) = 1.0_MK
        !END DO
			
		DO i = 1, N
		! B_rot(:,i) = R * B(:,i)
		
#ifdef __SP
			CALL SGEMV('N',N,N,1.0,R,N,B(:,i),1,0.0,multmp_vN,1)
#else      		
      		CALL dgemv('N',N,N,1.0d0,R,N,B(:,i),1,0.0d0,multmp_vN,1)
#endif			
			!B_rot(:,i) = matmul(R,B(:,i))
			B_rot(:,i) = multmp_vN
		END DO
		
		! C_PSO = B_rot * D^2 * B_rot'
		
		D_square = D**2
#ifdef __SP
		CALL SGEMM('N','T',N,N,N,1.0,D_square,N,B_rot,N,0.0,multmp_mN,N)
		CALL SGEMM('N','N',N,N,N,1.0,B_rot,N,multmp_mN,N,0.0,C_PSO,N)		
#else		
		CALL dgemm('N','T',N,N,N,1.0d0,D_square,N,B_rot,N,0.0d0,multmp_mN,N)
		CALL dgemm('N','N',N,N,N,1.0d0,B_rot,N,multmp_mN,N,0.0d0,C_PSO,N)
#endif

		!C_PSO = matmul(B_rot,matmul(D**2,transpose(B_rot)))
	

	  	!-----------------------------------------------------------------------
        !  Enforce symmetry
      	!-----------------------------------------------------------------------
		CALL tool_symmatrix(C_PSO,N,triuC)
	  	
	  	
!! TODO !!	  	
	  	!-----------------------------------------------------------------------
        !  Eigen decomposition, D=diagonal matrix of eigenvalues, 
        !  B=normalized eigenvectors
      	!-----------------------------------------------------------------------
	  	!CALL eigenDecomp(triuC,N,D,B)
  		END IF
  		
      RETURN
      END SUBROUTINE
