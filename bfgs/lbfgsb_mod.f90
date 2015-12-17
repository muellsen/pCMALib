      !-------------------------------------------------------------------------
      !  Module       :            lbfgsb_mod
      !-------------------------------------------------------------------------
      
      !References:

      !   [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
      !   memory algorithm for bound constrained optimization'',
      !   SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.

      !   [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
      !   Subroutines for Large Scale Bound Constrained Optimization''
      !   Tech. Report, NAM-11, EECS Department, Northwestern University,
      !   1994.
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  pCMALib: a parallel fortran 90 library for the evolution strategy with
      !           covariance matrix adaptation
      !  Christian L. Mueller, Benedikt Baumgartner, Georg Ofenbeck
      !  MOSAIC group, ETH Zurich, Switzerland
      !-------------------------------------------------------------------------
      
      
      MODULE lbfgsb_mod
      
      USE cmaes_param_mod
      USE cmaes_opts_mod
      USE cmaes_mod
      USE CEC2005 
      
      IMPLICIT NONE
      
      
      !This common Block is used in the lbfgs_nocedal line search
      !see description there
      !----------------
      COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
      DOUBLE PRECISION GTOL,STPMIN,STPMAX
      INTEGER MP,LP
      !----------------
      

      CONTAINS
      
      !-------------------------------------------------------------------------
      ! SUBROUTINE:     lbfgsb
      !-------------------------------------------------------------------------
      !  Purpose      : runs lbfgsb algorithm on test problem
      ! 
      !  Input       :  fcn     (E) External Function to call for Fitness Value
      !                 n       (I) Dimensions of problem
      !                 m       (I) Number of of limited memory corrections stored
      !                 x(n)    (R) Starting point of dimension n
      !
      ! Input/Output    f       (R) Function Value (at start/end of bfgs)
      !
      !
      !  Input(optional): 
      !                  
      !-------------------------------------------------------------------------
      
      
      SUBROUTINE lbfgsb(fcn,n,m,x,f)
      !Parameters
      INTEGER                                   :: n,m,i,j
      
      REAL(MK),DIMENSION(n)                     :: x,x1,x2,l,u
      REAL(MK),DIMENSION(n)                     :: g
      REAL(MK)                                  :: f
      REAL(MK),DIMENSION(n)                     :: DIAG
      EXTERNAL                                  :: fcn
      INTEGER                                   :: iflag,nbd(n)
	  INTEGER									:: iwa(3*n),isave(44)
      REAL(MK),DIMENSION(n*(2*m+1)+2*m)         :: W
      !INTEGER,DIMENSION(2)                      :: iprint
      INTEGER                                   :: iprint
      REAL(MK)                                  :: cmaes_funcwrap
	  REAL(MK),DIMENSION(n,n)                   :: A 
      REAL(MK)                                  :: DNRM2,wa(2*m*n+4*n+12*m*m+12*m)
      REAL(MK)                                  :: srad, distx
      REAL(MK), DIMENSION(n)                    :: x_org, x_diff, x_g, x_y, arrr
      REAL(MK)                                  :: tmpreal,factr,pgtol,dsave(29)     
      CHARACTER*60								:: task,csave
	  LOGICAL									:: lsave(4)
	  REAL(MK)                                  :: dpmeps
      INTEGER                                   :: iteration_counter
      iprint = options%BFGS_print

      ! nbd is a status variable for boundaries (see lbfgs-b.f)
      nbd = 0
      l = options%LBounds
      u = options%UBounds
      iteration_counter = 0
	  
	  tmpreal = dpmeps()
	  IF(options%BFGS_factr > 0_MK) THEN
	    factr = options%BFGS_factr
	  ELSE
	    factr = options%StopTolFun/tmpreal
	  END IF      
      pgtol= options%BFGS_pgtol
	  
	  task = 'START'
	  iprint = -1 
      x_org = 0.0d0
      
      
      !radius for the minimum step size
      srad = options%BFGS_sigmastep_min * sigma
      distx = 1_MK
      
      DO
         !saving the origin point
          x_org = x
          
           CALL setulb(n,m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint,&
                csave,lsave,isave,dsave)
          iteration_counter = iteration_counter + 1
      
      IF (iteration_counter .GT. options%BFGS_maxiter) THEN
        BFGSexitcounter(6) = BFGSexitcounter(6) + 1
        !EXIT
      END IF
      

      IF (.NOT. (task(1:2) .eq. 'FG' .OR. task(1:5) .eq. 'NEW_X')) THEN
       IF (task(1:11) .eq. 'CONVERGENCE') THEN
	    IF (task(1:48) .eq. 'CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL')&
	    THEN
	       BFGSexitcounter(2) = BFGSexitcounter(2) + 1
	    END IF
	   BFGSexitcounter(1) = BFGSexitcounter(1) + 1
       ENDIF
       IF (task(1:4) .eq. 'ABNO' .OR. task(1:5) .eq. 'ERROR') THEN
        BFGSexitcounter(3) = BFGSexitcounter(3) + 1
       ENDIF
       EXIT

      END IF


          
          IF (options%use_LJ) THEN
          !only get new x/g if we moved
           !IF (task(1:2) .eq. 'FG') THEN
            CALL LJ_GRAD(f,x,g,n,1) !For LJ optimization use the analytical gradient
            countBFGSEval = countBFGSEval + 1  
           !END IF
          ELSEIF (options%use_TIP) THEN
            CALL TIP4P_GRAD(f,x,g,n,1) !For TIP4P optimization use the analytical gradient
            countBFGSEval = countBFGSEval + 1
          ELSE
          

            IF (options%Benchfctnr .EQ. 9) THEN
                !only get new x/g if we moved
                !IF (task(1:2) .eq. 'FG') THEN
                    CALL F09_grad(g,x,n,n,options%LBounds,options%UBounds)!gradient
                !CALL my_grad_approx(n,x,fcn,f,g1) !otherwise approximate the Gradient
                    f = cmaes_funcwrap(x,n,fcn)
                    countBFGSEval = countBFGSEval + 1 !+ n
                !END IF
            ELSE
            
            !only get new x/g if we moved
                !IF (task(1:2) .eq. 'FG') THEN
                  f = cmaes_funcwrap(x,n,fcn)
                  countBFGSEval = countBFGSEval + 1 
                  CALL my_grad_approx(n,x,fcn,f,g) !otherwise approximate the Gradient
                !END IF
                
                
            END IF            
          END IF
          
!------------------------------------------------------------------------
     
      IF (task(1:5) .eq. 'NEW_X') THEN
          !checking if we moved from the origin point
          CALL VDSUB( n, x, x_org, x_diff )
          distx = DNRM2(n, x_diff, 1)    
         !----------------------------------------------------------------------
         !check if BFGS did a step that travel at least BFGS_sigmastep_min
          IF (tmpreal .GT. chi ) THEN    
           !IF (distx .GT. 0 .AND. distx .LT. srad) THEN
             BFGSexitcounter(4) = BFGSexitcounter(4) + 1 
             !EXIT
           END IF      
         !----------------------------------------------------------------------
         !check if BFGS is still within the multivariant gaussian with
         !prob. set in BFGS_confprop
          x1 = x-xmean
#ifdef __SP
            CALL SGEMV('N',input%N,input%N,1.0,invsC,input%N,x1,&
                                                            1,0.0,x2,1)
#else      		
            CALL DGEMV('N',input%N,input%N,1.0d0,invsC,input%N,x1,&
                                                            1,0.0d0,x2,1)
#endif
            tmpreal = DOT_PRODUCT(x1,x2)

 		  IF (tmpreal .GT. chi) THEN    
 		      BFGSexitcounter(5) = BFGSexitcounter(5) + 1
 		      !EXIT 
 		  END IF     
      END IF
      END DO
      
      END SUBROUTINE lbfgsb

      !-------------------------------------------------------------------------
      !  FUNCTION   :          my_grad_approx
      !-------------------------------------------------------------------------
      !
      !  Purpose      : this function should estimate gradients by finite differences
      !
      !  Input		  : 
      !
      !  Output       : 
      !  Remarks      : 
      !
      !-------------------------------------------------------------------------
      SUBROUTINE my_grad_approx(n,x,fcn,f,g)
      
      INTEGER                                   :: n
      REAL(MK),DIMENSION(n)                     :: x,g,xtmp
      REAL(MK)                                  :: f
      EXTERNAL                                  :: fcn
      
      
      !local vars
      REAL(MK)                                  :: cmaes_funcwrap
      INTEGER                                   :: i
      REAL(MK)                                  :: f_plus,f_minus
      
      
      IF (options%BFGS_central_difference ) THEN
          DO i = 1, n
            xtmp = x
            xtmp(i) = xtmp(i) + options%BFGS_grad_stepsize/2
            f_plus =  cmaes_funcwrap(xtmp,n,fcn)
            countBFGSEval = countBFGSEval + 1 
            xtmp(i) = xtmp(i) - options%BFGS_grad_stepsize
            f_minus =  cmaes_funcwrap(xtmp,n,fcn)
            countBFGSEval = countBFGSEval + 1 
            g(i) = (f_plus - f_minus)/options%BFGS_grad_stepsize
          END DO
      ELSE
           DO i = 1, n
            xtmp = x
            xtmp(i) = xtmp(i) + options%BFGS_grad_stepsize
            countBFGSEval = countBFGSEval + 1 
            f_plus =  cmaes_funcwrap(xtmp,n,fcn)
            !countBFGSEval = countBFGSEval + 1 
            g(i) = (f_plus - f)/options%BFGS_grad_stepsize
          END DO
      END IF
      
      
      
      
      
      
      
      END SUBROUTINE my_grad_approx
      
      
      
      !the gradient solver needs specific function call
      SUBROUTINE gradient_wrap (n,x,f, fcn)
      INTEGER                   :: n
      REAL(MK),DIMENSION(n)     :: x
      REAL(MK)                  :: f
      REAL(MK)                  :: cmaes_funcwrap
      EXTERNAL                  :: fcn
      
      
      
      
      f = cmaes_funcwrap(x,n,fcn)
      
      
      
      
      END SUBROUTINE gradient_wrap
      
      
      

      
      
      
      
      END MODULE lbfgsb_mod
