 	  !-------------------------------------------------------------------------
      !  SUBROUTINE   :                   water_write_out
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Writes a pdb file representing a water cluster
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
      

      SUBROUTINE water_write_out(filename,X,m,n) !,V,ELJ,GTEST,SECT)
      USE cmaes_param_mod
      USE cmaes_mod
      IMPLICIT NONE
      
      CHARACTER(LEN = 200),INTENT(IN)                :: filename
      INTEGER, INTENT(in)                            :: m ! (3*LJ_N-6) dim vector
      INTEGER, INTENT(in)                            :: n
      REAL(MK), DIMENSION(m,n),INTENT(in)            :: X
      
      !local vars
      REAL(MK), DIMENSION(m/2*3)                       :: X_full
      
      !REAL(MK), DIMENSION(m+6)                         :: V_full
      INTEGER                                          :: i,j
      LOGICAL                                          :: GTEST
      
        ! local variables
      INTEGER,DIMENSION(:),ALLOCATABLE      :: Npart_vec
      REAL(MK),DIMENSION(:,:),ALLOCATABLE :: xp_iproc
      INTEGER                               :: ip,count,iproc
      CHARACTER                             :: altloc
      CHARACTER (LEN =  4 )                 :: atom_name
      CHARACTER                             :: chains
      CHARACTER (LEN =  2 )                 :: charge
      CHARACTER (LEN =  2 )                 :: element
      CHARACTER                             :: icode
      INTEGER                               :: resno
      REAL(mk)                               :: occ
      CHARACTER (LEN =  3 )                 :: resname
      CHARACTER (LEN =  4 )                 :: segid
      INTEGER                               :: info
      REAL(MK),DIMENSION(9)                 :: Cords
      !====================================================================!
      ! initialize
      
      
      altloc = ' '
      resname = '   '
      chains = ' '
      resno = 0
      icode = ' '
      occ = 1.0D+00
      segid = '    '
      element = '  '
      charge = '  '
      
      
      
      
      
      DO i = 1,m/2,3
        CALL TIPIO( X(i,1), X(i+1,1), X(i+2,1),&
        X(m/2+i,1),X(m/2+1+i,1),X(m/2+2+i,1), Cords)
        DO j = 1,9
            X_full((i-1)/3*9+j) = Cords(j)
        END DO
      END DO
      
      !ATOM   1607  N   LYS A 233     -26.180   4.759 -18.385  1.00 41.53           N  

      !====================================================================!
      ! open file
      IF (countIter .EQ. 1) THEN
          OPEN(99,FILE=filename,FORM='FORMATTED',STATUS='REPLACE',IOSTAT=info)
          IF (info .NE. 0) THEN
             WRITE(*,*)'ERROR in water_writeout: opening file failed.'
             GOTO 9999
          ENDIF          
      END IF
      
      WRITE (99,'(a6,i5)')'MODEL ',countIter + 1
      DO i = 1,m/2*3,9
      atom_name = 'O   '
      WRITE (99, &
      '(a6,i5,1x,a4,a1,a3,a1,1x,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)' ) &
      'ATOM  ', i, atom_name, altloc, resname, chains, &
     & resno,icode,X_full(i),X_full(i+1),X_full(i+2),occ,0d0,segid,  &
     & element, charge
      atom_name = 'H   '
      WRITE (99, &
      '(a6,i5,1x,a4,a1,a3,a1,1x,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)' ) &
      'ATOM  ', i+1, atom_name, altloc, resname, chains, &
     & resno,icode,X_full(i+3),X_full(i+4),X_full(i+5),occ,0d0,segid,  &
     & element, charge
      WRITE (99, &
      '(a6,i5,1x,a4,a1,a3,a1,1x,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)' ) &
      'ATOM  ', i+2, atom_name, altloc, resname, chains, &
     & resno,icode,X_full(i+6),X_full(i+7),X_full(i+8),occ,0d0,segid,  &
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

      END SUBROUTINE water_write_out

      
      
      
      
      SUBROUTINE TIPIO(X1, Y1, Z1, L1, M1, N1, COORDS)
      IMPLICIT NONE
      DOUBLE PRECISION X1, Y1, Z1, COORDS(*), C2A1, M1, L1, N1, ALPHA1, CA1, S1, C3A1, L12, M12, N12

      L12=L1**2
      M12=M1**2
      N12=N1**2
      ALPHA1=SQRT(L12+M12+N12)
      CA1=COS(ALPHA1)
      C2A1=CA1
      IF (ALPHA1.LT.0.0001D0) THEN
!        C3A1=(-ALPHA1/2+ALPHA1**3/24)
         C3A1=(-0.5D0+ALPHA1**2/24.0D0)
         S1=(1.0D0-ALPHA1**2/6)
      ELSE
         C3A1=(CA1-1.0D0)/ALPHA1**2
         S1=SIN(ALPHA1)/ALPHA1
      ENDIF

      COORDS(1) = X1
      COORDS(2) = Y1
      COORDS(3) = Z1    
      COORDS(4) = 0.756950327*c2a1 - c3a1*l1*(0.756950327*l1 - 0.585882276*n1) + 0.585882276*m1*s1 + X1
      COORDS(5) = -(c3a1*m1*(0.756950327*l1 - 0.585882276*n1)) + (-0.585882276*l1 - 0.756950327*n1)*s1 + Y1
      COORDS(6) = -0.585882276*c2a1 - c3a1*(0.756950327*l1 - 0.585882276*n1)*n1 + 0.756950327*m1*s1 + Z1
      COORDS(7) = -0.756950327*c2a1 + c3a1*l1*(0.756950327*l1 + 0.585882276*n1) + 0.585882276*m1*s1 + X1
      COORDS(8) = c3a1*m1*(0.756950327*l1 + 0.585882276*n1) + (-0.585882276*l1 + 0.756950327*n1)*s1 + Y1
      COORDS(9) = -0.585882276*c2a1 + c3a1*(0.756950327*l1 + 0.585882276*n1)*n1 - 0.756950327*m1*s1 + Z1

      RETURN
      END
