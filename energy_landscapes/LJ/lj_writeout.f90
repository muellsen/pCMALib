 	  !-------------------------------------------------------------------------
      !  SUBROUTINE   :                   LJ_write_out
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Writes a pdb file representing a LJ cluster
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

      SUBROUTINE LJ_write_out(filename,X,m,n,Fmin) !,V,ELJ,GTEST,SECT)
      USE cmaes_param_mod
      USE cmaes_mod
      IMPLICIT NONE
      
      CHARACTER(LEN = 200),INTENT(IN)       :: filename
      INTEGER, INTENT(in)                            :: m ! (3*LJ_N-6) dim vector
      INTEGER, INTENT(in)                            :: n
      REAL(MK), DIMENSION(m,n),INTENT(in)            :: X ! Position
      REAL(MK),INTENT(in)			 	 			 :: Fmin !Energy difference to min
      
      !local vars
      REAL(MK), DIMENSION(m+6)                         :: X_full
      REAL(MK), DIMENSION(m+6)                         :: V_full
      INTEGER                                          :: i
      LOGICAL                                          :: GTEST
      
        ! local variables
      INTEGER,DIMENSION(:),ALLOCATABLE      :: Npart_vec
      REAL(MK),DIMENSION(:,:),ALLOCATABLE   :: xp_iproc
      INTEGER                               :: ip,count,iproc
      CHARACTER                             :: altloc
      CHARACTER (LEN =  4 )                 :: atom_name
      CHARACTER                             :: chains
      CHARACTER (LEN =  2 )                 :: charge
      CHARACTER (LEN =  2 )                 :: element
      CHARACTER                             :: icode
      INTEGER                               :: resno
      REAL(MK)                              :: occ,f_write
      CHARACTER (LEN =  3 )                 :: resname
      CHARACTER (LEN =  4 )                 :: segid
      INTEGER                               :: info
      
      !====================================================================!
      ! initialize
      
      atom_name = 'AR  '
      altloc = ' '
      resname = '   '
      chains = ' '
      resno = 0
      icode = ' '
      occ = 1.0D+00
      segid = '    '
      element = '  '
      charge = '  '
      
      X_full = 0
      X_full(4) = X(1,1)
      X_full(7) = X(2,1)
      X_full(8) = X(3,1)
      
      
      DO i = 4,m
        X_full(i+6) = X(i,1)
      END DO
      
      f_write = min(Fmin,9.9999_MK)

      !ATOM   1607  N   LYS A 233     -26.180   4.759 -18.385  1.00 41.53           N  

      !====================================================================!
      ! open file
      IF (countIter .EQ. 1) THEN
          OPEN(99,FILE=filename,FORM='FORMATTED',STATUS='REPLACE',IOSTAT=info)
          IF (info .NE. 0) THEN
             WRITE(*,*)'ERROR in lj_writeout: opening file failed.'
             GOTO 9999
          ENDIF          
      END IF
      
      WRITE (99,'(a6,i5)')'MODEL ',countRestarts
      DO i = 1,m+6,3
      WRITE (99, &
      '(a6,i5,1x,a4,a1,a3,a1,1x,i4,a1,3x,3f8.3,f6.2,f6.3,6x,a4,a2,a2)' ) &
      'ATOM  ', i, atom_name, altloc, resname, chains, &
     & resno,icode,X_full(i),X_full(i+1),X_full(i+2),occ,f_write,segid,  &
     & element, charge
      END DO
      
      !====================================================================!
      ! terminate pdb-file and close file
      WRITE (99,'(a6)')'ENDMDL'
      
         !WRITE (99, &
         !     '(a6,i5,1x,4x,a1,a3,1x,a1,i4,a1)' ) &
         !     'TER   ', ip+1, altloc, resname, chains, resno, icode
         !CLOSE(99,IOSTAT=info)
         !IF (info .NE. 0) THEN
         !   WRITE(*,*)'ERROR in lj_writeout: closing file failed.'
         !   GOTO 9999
         !ENDIF
      
            
9999  CONTINUE ! jump here upon error

      END SUBROUTINE LJ_write_out
