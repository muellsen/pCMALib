      !-------------------------------------------------------------------------
      !  SUBROUTINE   :          invnor
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This subroutine generates a normally distributed sequence
      !                 from a given sequence with values between [0, 1).
      !
      !  Input		  : m        		    (I)  spatial dimension
      !                 n                   (I)  number of points to generate
      !                 seq(m, n)           (R) sequence to be transformed
      !                 inverter            (I) which inverter to use (1-3)
      !
      !  Input/Output : iseq(m,n)           (R) normally distributed sequence
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
#ifdef __SP  
#define __CMA           
#endif
#ifdef __DP      
#define __CMA    
#endif
      
	  SUBROUTINE invnor(m, n, seq, iseq, inverter)

      use, intrinsic :: ieee_arithmetic
      
#ifndef __CMA    
      INTEGER,PARAMETER                             :: MK = kind(1.d0)
#else
      USE cmaes_param_mod
#endif
      IMPLICIT NONE

      INTEGER, INTENT(in) 					           :: m, n
	  REAL(MK), DIMENSION(m, n), INTENT(IN)    :: seq
	  REAL(MK), DIMENSION(m, n), INTENT(INOUT) :: iseq
	  INTEGER                                          :: i, j,inverter
	  REAL(MK)                                         :: morosinv
	  REAL                                             :: dinvnorm
	  REAL                                             :: HQNORM
	  
	
	  IF (inverter .EQ. 0) THEN
	      DO i = 1, m
      	    DO j = 1, n
      	        IF (seq(i,j) .EQ. 0.0_MK) THEN
      	            iseq(i,j) = 0.0_MK
      	        ELSE
      		        iseq(i, j) = morosinv(seq(i, j))
      		    IF (ieee_is_nan(iseq(i,j)) ) THEN
      		        iseq(i,j) = 0.0_MK
      		    END IF
      		    END IF
      	    ENDDO
          ENDDO
      ELSE
        IF (inverter .EQ. 1) THEN
	      DO i = 1, m
      	    DO j = 1, n
      	     IF (seq(i,j) .EQ. 0.0_MK) THEN
      	            iseq(i,j) = 0.0_MK
      	        ELSE
      		    iseq(i, j) = dinvnorm(seq(i, j))
      		    IF (ieee_is_nan(iseq(i,j)) ) THEN
      		        iseq(i,j) = 0.0_MK
      		    END IF
      		    END IF
      	    ENDDO
          ENDDO
        ELSE
         DO i = 1, m
      	    DO j = 1, n
      	     IF (seq(i,j) .EQ. 0.0_MK) THEN
      	            iseq(i,j) = 0.0_MK
      	        ELSE
      		    iseq(i, j) = HQNORM(seq(i, j)) !from prng_halton_sobol_R
      		    IF (ieee_is_nan(iseq(i,j)) ) THEN
      		        iseq(i,j) = 0.0_MK
      		    END IF
      		    END IF
      	    ENDDO
          ENDDO
        END IF
      END IF
    
      END SUBROUTINE invnor
    
      !-------------------------------------------------------------------------
      !  FUNCTION   :          morosinv
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This function converts a number in [0, 1) range into
      !                 a normal number.
      !
      !  Input		  : prob                (R) number in [0, 1) range
      !
      !  Output       : morosinv            (R) normal number
      !
      !  Remarks      : 
      !
      !  References   : http://www.wwz.unibas.ch/witheo/yvan/software/
      !                 modInverseNormalCDF.bas
      !-------------------------------------------------------------------------

      FUNCTION morosinv(u)

      
#ifndef __CMA    
      INTEGER,PARAMETER                             :: MK = kind(1.d0)
#else
      USE cmaes_param_mod
#endif    
      IMPLICIT NONE  
      REAL(MK) :: u
	  REAL(MK) :: x,r
	  REAL(MK) :: morosinv
      REAL(MK) :: a0 = 2.50662823884D0
      REAL(MK) :: a1 = -18.61500062529D0
	  REAL(MK) :: a2 = 41.39119773534D0
	  REAL(MK) :: a3 = -25.44106049637D0
      REAL(MK) :: b0 = -8.4735109309D0
      REAL(MK) :: b1 = 23.08336743743D0
      REAL(MK) :: b2 = -21.06224101826D0
      REAL(MK) :: b3 = 3.13082909833D0
      REAL(MK) :: c0 = 0.337475482272615D0
      REAL(MK) :: c1 = 0.976169019091719D0
      REAL(MK) :: c2 = 0.1607979714918209D0
      REAL(MK) :: c3 = 2.76438810333863D-2
      REAL(MK) :: c4 = 3.8405729373609D-3
      REAL(MK) :: c5 = 3.951896511919D-4
      REAL(MK) :: c6 = 3.21767881768D-5
      REAL(MK) :: c7 = 2.888167364D-7
	  REAL(MK) :: c8 = 3.960315187D-7


      ! check if argument is in legal range; throw an error if not.
      IF ((u .LE. -1.0D0) .OR. (u .GE. 1.0D0)) THEN
	  	PRINT *, 'Invalid probablity value'
      ENDIF
    
      x = u - 0.5D0

      IF (abs(x) < 0.42D0) THEN
      	r = x*x
        r = x * (((a3 * r + a2) * r + a1) * r + a0) &
     & / ((((b3 * r + b2) * r + b1) * r + b0) * r + 1.0D0)
	  ELSE
	    r = u;
    	IF (x > 0.0D0) THEN
    	    r = 1.0D0 - u
    	END IF
    	
        
        r = log(-log(r))    
        r=c0+r*(c1+r*(c2+r*(c3+r*(c4+r*(c5+r*(c6+r*(c7+r*c8)))))))
		IF (x < 0.0D0) r = -r
	  END IF
      
    
      morosinv = r

      END FUNCTION morosinv

      !-------------------------------------------------------------------------
      !  FUNCTION   :          dinvnorm
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This function converts a number in [0, 1) range into
      !                 a normal number.
      !
      !  Input		  : p                   (R) number in [0, 1) range
      !
      !  Output       : dinvnorm            (R) normal number
      !
      !  Remarks      : 
      !
      !  References   :         ren-raw chen, rutgers business school
      !                         normal inverse
      !                         translate from
      !                         http://home.online.no/~pjacklam/notes/invnorm
      !                         a routine written by john herrero
      !-------------------------------------------------------------------------

      REAL function dinvnorm(p)
#ifndef __CMA    
      INTEGER,PARAMETER                             :: MK = kind(1.d0)
#else
      USE cmaes_param_mod
#endif
      IMPLICIT NONE
      REAL(MK) p,p_low,p_high
      REAL(MK) a1,a2,a3,a4,a5,a6
      REAL(MK) b1,b2,b3,b4,b5
      REAL(MK) c1,c2,c3,c4,c5,c6
      REAL(MK) d1,d2,d3,d4
      REAL(MK) z,q,r
      a1=-39.6968302866538D0
      a2=220.946098424521D0
      a3=-275.928510446969D0
      a4=138.357751867269D0
      a5=-30.6647980661472D0
      a6=2.50662827745924D0
      b1=-54.4760987982241D0
      b2=161.585836858041D0
      b3=-155.698979859887D0
      b4=66.8013118877197D0
      b5=-13.2806815528857D0
      c1=-0.00778489400243029D0
      c2=-0.322396458041136D0
      c3=-2.40075827716184D0
      c4=-2.54973253934373D0
      c5=4.37466414146497D0
      c6=2.93816398269878D0
      d1=0.00778469570904146D0
      d2=0.32246712907004D0
      d3=2.445134137143D0
      d4=3.75440866190742D0
      p_low=0.02425D0
      p_high=1-p_low
      if(p.lt.p_low) goto 201
      if(p.ge.p_low) goto 301
201   q=SQRT(-2*log(p))
      z=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/&
     &((((d1*q+d2)*q+d3)*q+d4)*q+1)
      goto 204
301   if((p.ge.p_low).and.(p.le.p_high)) goto 202
      if(p.gt.p_high) goto 302
202   q=p-0.5
      r=q*q
      z=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/&
     &(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
      goto 204
302   if((p.gt.p_high).and.(p.lt.1)) goto 203
203   q=SQRT(-2*LOG(1-p))
      z=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/&
     &((((d1*q+d2)*q+d3)*q+d4)*q+1)
204   dinvnorm=z
      return
      end



!  
