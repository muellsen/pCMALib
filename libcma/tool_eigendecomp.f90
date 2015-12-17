	  !-------------------------------------------------------------------------
      !  Subroutine   :          tool_eigendecomp
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This Subroutine performs an eigenvalue decomposition of
      !					a symmetric matrix. It computes the matrix's eigenvalues      
      !					and optionally its eigenvectors
      !					it is mainly a wrapper function for the corresponding
      !					LAPACK routine
      !
      !  Input		  : A(:,:)		(R) Upper triangular part of sym. matrix
      !					n			(I)	Leading dimension of A
	  !
	  !  Output		  : eVals(:,:)	(R) Diagonal matrix with eigenvalues of A
	  !									on the main diagonal in ascending order
	  !					eVecs(:,:)	(R)	optional, contains orthonormal
	  !									eigenvectors of A
	  !                 info        (I) The Info Result Var return from Lapack
      !                
      !  Remarks      : TODO, calling the lapack routines causes the
      !					compiler to produce multiple errors like:
      !					multiple definitions of symbol _cos, _exp, _log, etc.
      !
      !  References   : 
      !					LAPACK User Guide ( http://www.netlib.org/lapack/lug/ )
	  !					Subroutine SSYEV, dsyev
      !  Revisions    :
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  pCMALib: a parallel fortran 90 library for the evolution strategy with
      !           covariance matrix adaptation
      !  Christian L. Mueller, Benedikt Baumgartner, Georg Ofenbeck
      !  MOSAIC group, ETH Zurich, Switzerland
      !-------------------------------------------------------------------------
#if(__KIND == __SINGLE_PRECISION)
      SUBROUTINE tool_eigendecomp_s(A,n,eVals,eVecs,info)
#elif(__KIND == __DOUBLE_PRECISION)
      SUBROUTINE tool_eigendecomp_d(A,n,eVals,eVecs,info)
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
      INTEGER, INTENT(in)							:: n
      REAL(MK),DIMENSION(n,n),INTENT(in)			:: A
      REAL(MK),DIMENSION(n,n),INTENT(out)			:: eVals
      REAL(MK),DIMENSION(n,n),OPTIONAL,INTENT(out)	:: eVecs
      INTEGER                                       :: info
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------
#if(__KIND == __SINGLE_PRECISION)      
      EXTERNAL 									:: SSYEV
#else
	  EXTERNAL 									:: dsyev 
#endif   
      !-------------------------------------------------------------------------
      !  Local Variables
      !-------------------------------------------------------------------------
      CHARACTER(len=1)							:: job
      INTEGER									:: lWork
      REAL(MK),ALLOCATABLE,DIMENSION(:)			:: work
      INTEGER									:: i
      REAL(MK),DIMENSION(n)						:: eValWork
      INTEGER                                           :: allocStat
     
      
      
      !-------------------------------------------------------------------------
      !  Specify job (only compute Eigenvalues('N') or also Eigenvectors('V') )
      !-------------------------------------------------------------------------
      IF(present(eVecs)) THEN
      	job = 'V'
      ELSE
      	job = 'N'
      END IF
      
      !-------------------------------------------------------------------------
      !  Make a copy of A (so it isnt changed unintentionally)
      !-------------------------------------------------------------------------
      eVecs = A
      
      !-------------------------------------------------------------------------
      !  Set length of work array
      !-------------------------------------------------------------------------
      !lWork = 2*n*n+6*n+1
      lWork = -1
      ALLOCATE(work(1),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating work(1)'
#if(__KIND == __SINGLE_PRECISION)
      CALL SSYEV(job,'U',n,eVecs,n,eValWork,work,-1,info)
#else
      CALL DSYEV(job,'U',n,eVecs,n,eValWork,work,-1,info)
#endif
      i = work(1)
      DEALLOCATE(work)
      
      ALLOCATE(work(i),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating work'
      !-------------------------------------------------------------------------
      !  Call LAPACK Routine
      !-------------------------------------------------------------------------
#if(__KIND == __SINGLE_PRECISION)
      CALL SSYEV(job,'U',n,eVecs,n,eValWork,work,lWork,info)
#else
      CALL DSYEV(job,'U',n,eVecs,n,eValWork,work,i,info)
#endif
      !-------------------------------------------------------------------------
      !  Error management
      !-------------------------------------------------------------------------
      IF(info .GT. 0) THEN
        WRITE(*,*) 'EVD failed to converge'
      END IF
      IF(info .LT. 0) WRITE(*,*) 'EVD: Illegal value'
      
      !-------------------------------------------------------------------------
      !  Set output variables
      !-------------------------------------------------------------------------
      DO i = 1,n
      	eVals(i,i) = eValWork(i) ! Convert Vector to Matrix
      END DO


      RETURN
      END SUBROUTINE