	  !-------------------------------------------------------------------------
      !  Subroutine   :          tool_create_rotmat
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Create a n-dimensional Rotation Matrix 'R' that rotates 
      !					a vector 'S' onto another vector 'T', such that
      !					R*S = a*T with 'a' being a scalar.
      !
      !  Input		  : S(:)		(R) Start vector
      !					T(:)		(R) Target vector
      !					N			(I)	Length of S and T
      !	
	  !  Output		  : R(:,:)		(R) Rotation Matrix that rotates S onto T
      !                
      !  Remarks      : Instead of the intrinsic matmul function the BLAS
      !					routines XGEMM,XGEMV are used. This gives a speedup of
      !					roughly 90%. For Givens Rotation we use BLAS routine
      !					XROTG
      !  References   : 
      !  Revisions    :
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  pCMALib: a parallel fortran 90 library for the evolution strategy with
      !           covariance matrix adaptation
      !  Christian L. Mueller, Benedikt Baumgartner, Georg Ofenbeck
      !  MOSAIC group, ETH Zurich, Switzerland
      !-------------------------------------------------------------------------
#if(__KIND == __SINGLE_PRECISION)
      SUBROUTINE tool_create_rotmat_s(S,T,N,R)
#elif(__KIND == __DOUBLE_PRECISION)
      SUBROUTINE tool_create_rotmat_d(S,T,N,R)
#endif
      IMPLICIT NONE
      
#if(__KIND == __SINGLE_PRECISION)
	  INTEGER,PARAMETER	:: MK = KIND(1.E0)
#else
	  INTEGER,PARAMETER	:: MK = KIND(1.D0)
#endif       
      !-------------------------------------------------------------------------
      !  Parameters
      !-------------------------------------------------------------------------
      INTEGER, INTENT(in)							:: N
      REAL(MK),DIMENSION(N),INTENT(in)				:: S
      REAL(MK),DIMENSION(N),INTENT(in)				:: T
      REAL(MK),DIMENSION(N,N),INTENT(out)			:: R

   
      !-------------------------------------------------------------------------
      !  Local Variables
      !------------------------------------------------------------------------
	  REAL(MK),DIMENSION(N)							:: S_work,T_work,S_sup,T_sup
	  REAL(MK),DIMENSION(N,N)						:: R_tar,R_tmp, R_sup
	  REAL(MK)										:: G_S,G_C
	  REAL(MK),DIMENSION(2)							:: S_tmp,T_tmp
	  INTEGER										:: p,q,i
     
      !-------------------------------------------------------------------------
      !  Check input
      !-------------------------------------------------------------------------
      IF(size(S) .NE. size(T)) WRITE(*,*) &
      		'tool_create_rotmat: input vectors have different size'
      IF(size(S) .NE. N) WRITE(*,*) &
      		'tool_create_rotmat: input vectors are not of size N'
      		
      !-------------------------------------------------------------------------
      !  Initialize variables
      !-------------------------------------------------------------------------
      S_work = S
      T_work = T

      R = 0.0_MK
      DO i = 1,N
      	R(i,i) = 1.0_MK
      END DO	  
	  R_tar = R
	  R_tmp = R
  
	  !-------------------------------------------------------------------------
      !  CREATE ROTATION MATRIX R
      !-------------------------------------------------------------------------
	  DO p = N-1,1,-1
	  	DO q = N,p+1,-1
	  		! Get vectors in current plane
	  		T_tmp(1) = T_work(p)
      		T_tmp(2) = T_work(q)
      		S_tmp(1) = S_work(p)
      		S_tmp(2) = S_work(q)
      		
      		!-------------------------------------------------------------------
      		! Perform Givens Rotation on start vector
      		!-------------------------------------------------------------------
      		! SROTG/drotg compute G_C and G_S such that
      		! | G_C -G_S |   | A |   | r |
      		! |          | * |   | = |   |
      		! | G_S  G_C |   | B |   | 0 |
      		!
#if(__KIND == __SINGLE_PRECISION)
	  		CALL SROTG(S_tmp(1),S_tmp(2),G_C,G_S)
#else      		
      		CALL drotg(S_tmp(1),S_tmp(2),G_C,G_S)
#endif      
			! Check direction of rotation
			IF(S_tmp(1) .LT. 0) THEN
				G_C = -G_C
				G_S = -G_S
			END IF
				
			! Build a Rotation Matrix out of G_C and G_S		
      		R_tmp = 0.0_MK
      		DO i = 1,N
      			R_tmp(i,i) = 1.0_MK
      		END DO
      		R_tmp(p,p) = G_C
      		R_tmp(q,q) = G_C
      		R_tmp(p,q) = G_S
      		R_tmp(q,p) = -G_S
	
      		! Rotate start vector and update R      		
#if(__KIND == __SINGLE_PRECISION)      		
      		! S_work = R_tmp*S_work
      		CALL SGEMV('N',N,N,1.0,R_tmp,N,S_work,1,0.0,S_sup,1)
      		S_work = S_sup
      		! R = R_tmp*R
      		CALL SGEMM('N','N',N,N,N,1.0,R_tmp,N,R,N,0.0,R_sup,N)
      		R = R_sup
#else
      		! S_work = R_tmp*S_work
      		CALL dgemv('N',N,N,1.0d0,R_tmp,N,S_work,1,0.0d0,S_sup,1)
      		S_work = S_sup
      		! R = R_tmp*R
      		CALL dgemm('N','N',N,N,N,1.0d0,R_tmp,N,R,N,0.0d0,R_sup,N)
      		R = R_sup
#endif      		
      		
      		      		
      		!-------------------------------------------------------------------
      		! Perform Givens Rotation on target vector
      		!-------------------------------------------------------------------
#if(__KIND == __SINGLE_PRECISION)
	  		CALL SROTG(T_tmp(1),T_tmp(2),G_C,G_S)
#else      		
      		CALL drotg(T_tmp(1),T_tmp(2),G_C,G_S)
#endif         
			IF(T_tmp(1) .LT. 0) THEN
				G_C = -G_C
				G_S = -G_S
			END IF
				
	     	R_tmp = 0.0_MK	
			DO i = 1,N
      			R_tmp(i,i) = 1.0_MK
      		END DO	
      		R_tmp(p,p) = G_C
      		R_tmp(q,q) = G_C
      		R_tmp(p,q) = G_S
      		R_tmp(q,p) = -G_S
      		
      		! Rotate target vector and update R_tar
#if(__KIND == __SINGLE_PRECISION)       		
			! T_work = R_tmp*T_work  		
      		CALL SGEMV('N',N,N,1.0,R_tmp,N,T_work,1,0.0,T_sup,1)
      		T_work = T_sup
      		!R_tar = R_tmp*R_tar
      		CALL SGEMM('N','N',N,N,N,1.0,R_tmp,N,R_tar,N,0.0,R_sup,N)
      		R_tar = R_sup
#else     
			! T_work = R_tmp*T_work  		
      		CALL dgemv('N',N,N,1.0d0,R_tmp,N,T_work,1,0.0d0,T_sup,1)
      		T_work = T_sup
      		!R_tar = R_tmp*R_tar
      		CALL dgemm('N','N',N,N,N,1.0d0,R_tmp,N,R_tar,N,0.0d0,R_sup,N)
      		R_tar = R_sup
#endif      		  		
      		
	  	END DO
	  END DO
	  !print*,'blas',R_tar
  
	  !-------------------------------------------------------------------------
      !  Construct final Rotation Matrix: R = inv(R_tar)*R
      !-------------------------------------------------------------------------
#if(__KIND == __SINGLE_PRECISION)   
	  CALL SGEMM('T','N',N,N,N,1.0,R_tar,N,R,N,0.0,R_sup,N)
#else
	  CALL dgemm('T','N',N,N,N,1.0d0,R_tar,N,R,N,0.0d0,R_sup,N)
#endif
      R = R_sup
	  
      RETURN
      END SUBROUTINE		