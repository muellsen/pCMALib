 !-------------------------------------------------------------------------
      !  Module       :            Double Funnel
      !-------------------------------------------------------------------------
      !
      !  Purpose      : a fortran implementation of the Double Funnel Test
      !
      !  Remarks      : 
      !
      !  References   :  The Impact of Global Structure on Search
      !                  Monte Lunacek, Darrell Whitley, and Andrew Sutton
      !                  
      !  Revisions    :
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  pCMALib: a parallel fortran 90 library for the evolution strategy with
      !           covariance matrix adaptation
      !  Christian L. Mueller, Benedikt Baumgartner, Georg Ofenbeck
      !  MOSAIC group, ETH Zurich, Switzerland
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! doubleFunnel
      !-------------------------------------------------------------------------
      !  Purpose      : runs benchmark on a given functionnr
      ! 
      !  Input       :  funcnr  (I) the number of the function to be benchmarked 
      !                 m       (I) m Dimension of the matrix (Rows)
      !                               = Dimension of Vectors
      !                 n       (I) n Dimension of the matrix (Columns)
      !                               = # of Vectors to process
      !  Input/Output:  vars    (R) the matrix with the input values of size m*n
      !
      !  Output:        res(R) the vector with the results
      !
      !-------------------------------------------------------------------------
      
      
      !biases,shiftmatrix,rotatematrix
      SUBROUTINE DoubleFunnel(res,vars,m,n,lbounds,ubounds) 
#ifdef __SP  
#define __CMA           
#endif
#ifdef __DP      
#define __CMA    
#endif


#ifdef __CMA
      USE cmaes_param_mod
      USE cmaes_opts_mod
#endif


#ifndef __CMA    
       INTEGER,PARAMETER                             :: MK = 8
#endif

      
      
      REAL(MK)                                   :: Pi = 3.141592653589793&
     &23846264338327950288419716939937510d0                  
      REAL(MK),DIMENSION(n),INTENT(out)            :: res
      REAL(MK),DIMENSION(m,n),INTENT(in)           :: vars
      INTEGER,INTENT(in)                           :: m
      INTEGER,INTENT(in)                           :: n
      REAL(MK),DIMENSION(m),OPTIONAL              :: lbounds
      REAL(MK),DIMENSION(m),OPTIONAL              :: ubounds
      
       
      ! Locally used variables
      REAL(MK),DIMENSION(m)                        :: ub,lb
      REAL(MK)                                     :: s1,s2,s,d
      REAL(MK)                                     :: u1,u2
      LOGICAL,DIMENSION(n)                         :: boundviolation
      LOGICAL                                      :: rast
      
      IF (.NOT. PRESENT(lbounds)) THEN
        lb = -5_MK
      ELSE
        lb = lbounds
      END IF
      
      IF (.NOT. PRESENT(ubounds)) THEN
        ub = 5_MK
      ELSE
        ub = ubounds
      END IF
      
      
#ifdef __CMA
      s = options%DF_s
      d = options%DF_d
      rast = options%DF_rast
#else
      s = 0.7d0
      d = 0
      rast = .TRUE.
#endif
      
      
      u1 = 2.5d0
      u2 = -2.5d0
      
      IF (s .NE. 0d0) THEN
        d = u1**2-s*u2**2
      ELSE
        s = -(d - u1**2)/u2**2
      END IF
      
      DO i = 1,n
            s1 = 0
            s2 = 0
            DO j = 1,m
                s1 = s1 + (vars(j,i) - u1)**2
                s2 = s2 + (vars(j,i) - u2)**2
            END DO
            s2 = d*m+s*s2
            res(i) = MIN(s1,s2)
            IF (rast ) THEN
                s1 = 0
                DO j = 1,m
                s1 = s1 + (1-COS(2*PI*(vars(j,i) - u1)))
                END DO
                s1 = 10*s1
                res(i) = res(i) + s1
            END IF
        
      END DO    
      
      
      
      
      
      END SUBROUTINE DoubleFunnel
