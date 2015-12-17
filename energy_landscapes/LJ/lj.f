
C   GMIN: A program for finding global minima
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of GMIN.
C
C   GMIN is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   GMIN is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C
C  Energy and gradient for LJ.

      SUBROUTINE LJ(ELJ,X,m,n,lbounds,ubounds) !,V,ELJ,GTEST,SECT)
      USE cmaes_param_mod
      IMPLICIT NONE
      
      INTEGER, INTENT(in)                            :: m ! (3*LJ_N-6) dim vector
      INTEGER, INTENT(in)                            :: n
      REAL(MK), DIMENSION(n),  INTENT(out)           :: ELJ
      REAL(MK), DIMENSION(m,n),INTENT(in)            :: X

      REAL(MK), OPTIONAL                             :: lbounds
      REAL(MK), OPTIONAL                             :: ubounds
      
      !local vars
      REAL(MK), DIMENSION(m+6)                         :: X_full
      REAL(MK), DIMENSION(m+6)                         :: V_full
      INTEGER                                          :: i,j
      LOGICAL                                          :: GTEST
      REAL(MK)                                         :: ELJn
      
      DO j = 1, n
      GTEST = .TRUE.
      X_full = 0
      X_full(4) = X(1,1)
      X_full(7) = X(2,1)
      X_full(8) = X(3,1)
      
      DO i = 4,m
        X_full(i+6) = X(i,n)
      END DO
      
      
      GTEST = .FALSE.
      !CALL LJDIFF((m+6)/3,X_full)
      CALL LJ_1(ELJn,X_full,m+6,n,V_full,GTEST)
      ELJ(n) = ELJn
      END DO
      
      END SUBROUTINE LJ

      SUBROUTINE LJ_GRAD(ELJ,X,V,m,n) !,V,ELJ,GTEST,SECT)
      USE cmaes_param_mod
      IMPLICIT NONE
      
      INTEGER, INTENT(in)                            :: m ! (3*LJ_N-6) dim vector
      INTEGER, INTENT(in)                            :: n
      REAL(MK), DIMENSION(n),  INTENT(out)           :: ELJ
      REAL(MK), DIMENSION(m,n),INTENT(in)            :: X
      REAL(MK), DIMENSION(m),INTENT(OUT)             :: V

      !REAL(MK), OPTIONAL                             :: lbounds
      !REAL(MK), OPTIONAL                             :: ubounds
      
      !local vars
      REAL(MK), DIMENSION(m+6)                       :: X_full
      REAL(MK), DIMENSION(m+6)                         :: V_full
      INTEGER                                          :: i,j
      LOGICAL                                          :: GTEST
      REAL(MK)                                         :: ELJn
      
      
      DO j = 1, n
      GTEST = .TRUE.
      X_full = 0
      X_full(4) = X(1,1)
      X_full(7) = X(2,1)
      X_full(8) = X(3,1)
      

      DO i = 4,m
        X_full(i+6) = X(i,1)
      END DO
      GTEST = .TRUE.
      !CALL LJDIFF((m+6)/3,X_full)
      CALL LJ_1(ELJn,X_full,m+6,n,V_full,GTEST)
      ELJ(n) = ELJn
      
      
      DO i = 4,m
        V(i) = V_full(i+6)
      END DO
      V(1) = V_full(4)
      V(2) = V_full(7)
      V(3) = V_full(8)
      END DO
      
      END SUBROUTINE LJ_GRAD



C
      SUBROUTINE LJ_1(ELJ,X,m,n,V,GTEST)
      USE cmaes_param_mod
      IMPLICIT NONE
      
      INTEGER, INTENT(in)                            :: m ! (3*LJ_N-6) dim vector
      INTEGER, INTENT(in)                            :: n
      REAL(MK),  INTENT(out)                         :: ELJ
      REAL(MK), DIMENSION(m),INTENT(in)              :: X

      !REAL(MK), OPTIONAL                             :: lbounds
      !REAL(MK), OPTIONAL                             :: ubounds
      
      !Locals
      LOGICAL GTEST,SECT
      INTEGER J1, J2, J3, J4
      INTEGER NATOMS
      DOUBLE PRECISION Xt(3*m), DIST, V(3*m), G(m,m), VT(m),
     1                 R6, ELJt, DUMMYX, DUMMYY, DUMMYZ, XMUL2, DUMMY
      LOGICAL EVAP, evapreject
      COMMON /EV/ EVAP, evapreject
      
      
      
      SECT = .FALSE.
      NATOMS = m/3
!  
      IF (SECT) THEN
         CALL LJDIFF(NATOMS, X)
         RETURN
      ENDIF
C      EVAP=.FALSE.
      ELJ=0.0D0
      DO J1=1,NATOMS 
         VT(J1)=0.0D0
      ENDDO
      IF (GTEST) THEN
         DO J1=1,NATOMS
            J3=3*J1
            DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
C            IF ((DIST.GT.RADIUS).AND.(.NOT.LJMFT)) THEN
CC              IF (DEBUG) WRITE(*,'(A,I3,A,3G15.5,A,G15.5,A,G15.5)') 
CC    1                    'Atom ',J1,' at ',X(J3-2),X(J3-1),X(J3),' DIST=',DIST,' RADIUS=',RADIUS
C               EVAP=.TRUE.
C               ELJ=ELJ+10.0D2*(DIST-RADIUS)**2
C            ENDIF
            G(J1,J1)=0.0D0
            DO J2=J1+1,NATOMS
               J4=3*J2
        DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               DIST=1.0D0/DIST
               R6=DIST**3
               DUMMY=R6*(R6-1.0D0)
               VT(J1)=VT(J1)+DUMMY
               VT(J2)=VT(J2)+DUMMY
               ELJ=ELJ+DUMMY
               DIST=DIST*R6
               G(J2,J1)=-24.0D0*(2.0D0*R6-1.0D0)*DIST
               G(J1,J2)=G(J2,J1)
            ENDDO
         ENDDO
      ELSE
         DO J1=1,NATOMS
            J3=3*J1
            DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
C            IF ((DIST.GT.RADIUS).AND.(.NOT.LJMFT)) THEN
C               EVAP=.TRUE.
CC              IF (DEBUG) WRITE(*,'(A,I3,A,3G15.5,A,G15.5,A,G15.5)') 
CC    1                    'Atom ',J1,' at ',X(J3-2),X(J3-1),X(J3),' DIST=',DIST,' RADIUS=',RADIUS
C               ELJ=ELJ+10.0D2*(DIST-RADIUS)**2
C            ENDIF
            DO J2=J1+1,NATOMS
               J4=3*J2
        DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               DIST=1.0D0/DIST
               R6=DIST**3
               DUMMY=R6*(R6-1.0D0)
               VT(J1)=VT(J1)+DUMMY
               VT(J2)=VT(J2)+DUMMY
               ELJ=ELJ+DUMMY
            ENDDO
         ENDDO
      ENDIF
C      IF (DEBUG.AND.EVAP) THEN
CC        WRITE(*,'(A)') 'An atom has evaporated - dumping coordinates'
CC        WRITE(40,'(I4)') NATOMS
CC        WRITE(40,'(A,I4,A,F15.5)') 'energy after evap=',ELJ
CC        WRITE(40,'(A2,3F20.10)') ('LJ ',X(3*(J1-1)+1),X(3*(J1-1)+2),X(3*(J1-1)+3),J1=1,NATOMS)
C      ENDIF
      ELJ=ELJ*4.0D0

      IF (.NOT.GTEST) RETURN

      DO J1=1,NATOMS
         J3=3*J1
         IF (.FALSE. ) THEN!SEEDT.AND.(J1.GT.NATOMS-NSEED).AND.FREEZECORE) THEN
            V(J3-2)=0.0D0
            V(J3-1)=0.0D0
            V(J3)=0.0D0
         ELSE
            DUMMYX=0.0D0
            DUMMYY=0.0D0
            DUMMYZ=0.0D0
            DO J4=1,NATOMS
               J2=3*J4
               XMUL2=G(J4,J1)
               DUMMYX=DUMMYX+XMUL2*(X(J3-2)-X(J2-2))
               DUMMYY=DUMMYY+XMUL2*(X(J3-1)-X(J2-1))
               DUMMYZ=DUMMYZ+XMUL2*(X(J3)  -X(J2))
            ENDDO
            DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
C           IF ((DIST.GT.RADIUS).AND.(.NOT.LJMFT)) THEN
C              DUMMYX=DUMMYX+10.0D2*16.0D0*(DIST-RADIUS)*X(J3-2)
C              DUMMYY=DUMMYY+10.0D2*16.0D0*(DIST-RADIUS)*X(J3-1)
C              DUMMYZ=DUMMYZ+10.0D2*16.0D0*(DIST-RADIUS)*X(J3)
C           ENDIF
            V(J3-2)=DUMMYX
            V(J3-1)=DUMMYY
            V(J3)=DUMMYZ
         ENDIF
      ENDDO
!     WRITE(MYUNIT,'(A)') 'lj> coords and VT for atom 192:'
!     WRITE(MYUNIT,'(I6,4F20.10)') (J1,VT(J1),X(3*(J1-1)+1),X(3*(J1-1)+2),X(3*(J1-1)+3),J1=192,192)

      RETURN
      END
C
C*************************************************************************
C
C  Subroutine LJDIFF calculates the cartesian second
C  derivative matrix analytically. Reduced units.
C
C*************************************************************************
C
      SUBROUTINE LJDIFF(N, X)
      IMPLICIT NONE
      INTEGER N, J1, J2,natoms
      DOUBLE PRECISION X(3*N), 
     1                 R2(N,N), 
     2                 R8(N,N), G(N,N),
     3                 R14(N,N), F(N,N)
C 
C  Store distance matrices.
C
      
      DO J1=1,N
         R2(J1,J1)=0.0D0
         R8(J1,J1)=0.0D0
         R14(J1,J1)=0.0D0
         DO J2=J1+1,N
            R2(J2,J1)=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1               +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2               +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
            R2(J2,J1)=1.0D0/R2(J2,J1)
            R8(J2,J1)=R2(J2,J1)**4
            R14(J2,J1)=R8(J2,J1)*R8(J2,J1)/R2(J2,J1)
            R2(J1,J2)=R2(J2,J1)
         ENDDO
      ENDDO

      CALL LJS(G,F,R2,R14,R8,X,N,N)

      RETURN
      END

C*****************************************************************************

      SUBROUTINE LJS(G,F,R2,R14,R8,X,N,NATOMS)
      !USE commons
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6,NATOMS
      DOUBLE PRECISION G(NATOMS,NATOMS), R14(NATOMS,NATOMS),
     1                 R8(NATOMS,NATOMS),
     1                 R2(NATOMS,NATOMS), F(NATOMS,NATOMS), 
     2                 X(3*NATOMS),DUMMY,HESS(3*NATOMS,3*NATOMS)

      DO J1=1,N
         G(J1,J1)=0.0D0
         F(J1,J1)=0.0D0
         DO J2=J1+1,N 
            F(J2,J1)=672.0D0*R14(J2,J1)-192.0D0*R8(J2,J1)
            F(J1,J2)=F(J2,J1)
            G(J2,J1)=-24.0D0*(2.0D0*R14(J2,J1)-R8(J2,J1))
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO
C
C  Now do the hessian. First are the entirely diagonal terms.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0
            DO J4=1,N
               DUMMY=DUMMY+F(J4,J1)*R2(J4,J1)*
     1                 (X(J3)-X(3*(J4-1)+J2))**2 + G(J4,J1)   
            ENDDO
            HESS(J3,J3)=DUMMY
         ENDDO
      ENDDO
C
C  Next are the terms where x_i and x_j are on the same atom
C  but are different, e.g. y and z.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J2+1,3
               DUMMY=0.0D0
               DO J5=1,N
                  DUMMY=DUMMY + F(J5,J1)*R2(J5,J1)* 
     1           (X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4)) 
               ENDDO
               HESS(3*(J1-1)+J4,J3)=DUMMY
            ENDDO
         ENDDO
      ENDDO
C
C  Case III, different atoms, same cartesian coordinate.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*
     1                           (X(J3)-X(3*(J4-1)+J2))**2-G(J4,J1) 
            ENDDO
         ENDDO
      ENDDO
C
C  Case IV: different atoms and different cartesian coordinates.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               DO J5=1,J2-1
                  J6=3*(J4-1)+J5
                  HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)
     1                    *(X(J3)-X(3*(J4-1)+J2))
     2                    *(X(3*(J1-1)+J5)-X(J6))
                  HESS(3*(J4-1)+J2,3*(J1-1)+J5)=HESS(J6,J3)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C  Symmetrise Hessian
C
      DO J1=1,3*N
         DO J2=J1+1,3*N
            HESS(J1,J2)=HESS(J2,J1)
         ENDDO
      ENDDO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Energy and gradient for LJ with a cutoff
C
!      SUBROUTINE LJCUT(X,V,ELJ,GTEST,SECT)
!      USE commons
!      IMPLICIT NONE
!      LOGICAL GTEST,SECT
!      INTEGER J1, J2, J3, J4
!      DOUBLE PRECISION X(3*NATOMS), DIST, V(3*NATOMS), G(NATOMS,NATOMS), 
!     1                 R6, ELJ, DUMMYX, DUMMYY, DUMMYZ, XMUL2, DUMMY, CUTOFF2
!      LOGICAL EVAP, evapreject
!      COMMON /EV/ EVAP, evapreject
!
!      CUTOFF2=CUTOFF**2
!      IF (SECT) THEN
!         CALL LJDIFF(NATOMS, X)
!         RETURN
!      ENDIF
!C     EVAP=.FALSE.
!      ELJ=0.0D0
!      DO J1=1,NATOMS 
!         VT(J1)=0.0D0
!      ENDDO
!      IF (GTEST) THEN
!         DO J1=1,NATOMS
!            J3=3*J1
!            DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
!C            IF ((DIST.GT.RADIUS).AND.(.NOT.LJMFT)) THEN
!CC              IF (DEBUG) WRITE(*,'(A,I3,A,3G15.5,A,G15.5,A,G15.5)') 
!CC    1                    'Atom ',J1,' at ',X(J3-2),X(J3-1),X(J3),' DIST=',DIST,' RADIUS=',RADIUS
!C               EVAP=.TRUE.
!C               ELJ=ELJ+10.0D2*(DIST-RADIUS)**2
!C            ENDIF
!            G(J1,J1)=0.0D0
!            atom2: DO J2=J1+1,NATOMS
!               J4=3*J2
!               DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
!               IF (DIST.GT.CUTOFF2) CYCLE atom2
!               DIST=1.0D0/DIST
!               R6=DIST**3
!               DUMMY=R6*(R6-1.0D0)
!               VT(J1)=VT(J1)+DUMMY
!               VT(J2)=VT(J2)+DUMMY
!               ELJ=ELJ+DUMMY
!               DIST=DIST*R6
!               G(J2,J1)=-24.0D0*(2.0D0*R6-1.0D0)*DIST
!               G(J1,J2)=G(J2,J1)
!            ENDDO atom2
!         ENDDO
!      ELSE
!         DO J1=1,NATOMS
!            J3=3*J1
!            DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
!C            IF ((DIST.GT.RADIUS).AND.(.NOT.LJMFT)) THEN
!C               EVAP=.TRUE.
!CC              IF (DEBUG) WRITE(*,'(A,I3,A,3G15.5,A,G15.5,A,G15.5)') 
!CC    1                    'Atom ',J1,' at ',X(J3-2),X(J3-1),X(J3),' DIST=',DIST,' RADIUS=',RADIUS
!C               ELJ=ELJ+10.0D2*(DIST-RADIUS)**2
!C            ENDIF
!            atom22: DO J2=J1+1,NATOMS
!               J4=3*J2
!               DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
!               IF (DIST.GT.CUTOFF2) CYCLE atom22
!               DIST=1.0D0/DIST
!               R6=DIST**3
!               DUMMY=R6*(R6-1.0D0)
!               VT(J1)=VT(J1)+DUMMY
!               VT(J2)=VT(J2)+DUMMY
!               ELJ=ELJ+DUMMY
!            ENDDO atom22
!         ENDDO
!      ENDIF
!C      IF (DEBUG.AND.EVAP) THEN
!CC        WRITE(*,'(A)') 'An atom has evaporated - dumping coordinates'
!CC        WRITE(40,'(I4)') NATOMS
!CC        WRITE(40,'(A,I4,A,F15.5)') 'energy after evap=',ELJ
!CC        WRITE(40,'(A2,3F20.10)') ('LJ ',X(3*(J1-1)+1),X(3*(J1-1)+2),X(3*(J1-1)+3),J1=1,NATOMS)
!C      ENDIF
!      ELJ=ELJ*4.0D0
!
!      IF (.NOT.GTEST) RETURN
!
!      DO J1=1,NATOMS
!         J3=3*J1
!         IF (SEEDT.AND.(J1.GT.NATOMS-NSEED).AND.FREEZECORE) THEN
!            V(J3-2)=0.0D0
!            V(J3-1)=0.0D0
!            V(J3)=0.0D0
!         ELSE
!            DUMMYX=0.0D0
!            DUMMYY=0.0D0
!            DUMMYZ=0.0D0
!            DO J4=1,NATOMS
!               J2=3*J4
!               XMUL2=G(J4,J1)
!               DUMMYX=DUMMYX+XMUL2*(X(J3-2)-X(J2-2))
!               DUMMYY=DUMMYY+XMUL2*(X(J3-1)-X(J2-1))
!               DUMMYZ=DUMMYZ+XMUL2*(X(J3)  -X(J2))
!            ENDDO
!C           DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
!C           IF ((DIST.GT.RADIUS).AND.(.NOT.LJMFT)) THEN
!C              DUMMYX=DUMMYX+10.0D2*16.0D0*(DIST-RADIUS)*X(J3-2)
!C              DUMMYY=DUMMYY+10.0D2*16.0D0*(DIST-RADIUS)*X(J3-1)
!C              DUMMYZ=DUMMYZ+10.0D2*16.0D0*(DIST-RADIUS)*X(J3)
!C           ENDIF
!            V(J3-2)=DUMMYX
!            V(J3-1)=DUMMYY
!            V(J3)=DUMMYZ
!         ENDIF
!      ENDDO
!
!      RETURN
!      END
