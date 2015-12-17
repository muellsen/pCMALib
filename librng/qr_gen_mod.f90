      !-------------------------------------------------------------------------
      !  Module       :            qr_generator_mod
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! SUBROUTINE:     qr_generator
      !-------------------------------------------------------------------------
      !  Purpose      : returns a Qausi Random Number Series
      ! 
      !  Input       :
      !                 n       (I) Dimensions to generate
      !                 m       (I) Number of points per Dimension
      !
      !  Input(optional): 
      !                   sampler (I) which QR Sampler to use. Choises are:
      !                             0.............Sobol
      !                    (default)1.............Sobol with scrambeling 
      !                             2.............Halton
      !                             3.............Halton R implementation
      !                             4.............Faure (buggy!)
      !                             5.............Niederreiter
      !
      !                   inverter (I) which Inverter to use (Uniform -> Normal(0,1))
      !                             -1............dont use any inverter (uniform)
      !                             0.............Moros Inverse
      !                    (default)1             Peter J. Acklam's Inverter
      !                             2.............Inverter from the R Implementation
      !                 
      !
      !  Input/Output:  
      !                 seed1   (I) Seed that is used (and changed! - save new seed!)
      !
      !
      !  Output:        vars    (R) n,m dimensional array containing the Quasi Random Sample
      !
      !
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
      MODULE qr_gen_mod

#ifndef __CMA    
      INTEGER,PARAMETER                             :: MK = kind(1.d0)
#else
      USE cmaes_param_mod
      USE cmaes_opts_mod
#endif
      IMPLICIT NONE
      
      INTEGER, PARAMETER                              :: maxbit=30
      REAL(MK),DIMENSION(:),POINTER                   :: sobol_quasi
      INTEGER                                         :: ll
      INTEGER,DIMENSION(:,:),POINTER                  :: sv
      INTEGER                                         :: init_sobol = 1
      INTEGER                                         :: init_faure = 1
      INTEGER                                         :: count
      CONTAINS
      
  
      SUBROUTINE qr_generator(n,m,seed1,vars,sampler,inverter)
      
      !Parameters
      INTEGER                                         :: n,m
      INTEGER,INTENT(INOUT)                           :: seed1
      REAL(MK),DIMENSION(n,m),INTENT(OUT)             :: vars
      INTEGER, OPTIONAL                               :: sampler, inverter

      !local vars
      
      INTEGER                                         :: m_sampler,m_inverter
      INTEGER                                         :: i, prime
      INTEGER                                         :: iflag,step,qs,prime_ge
      
      INTEGER,DIMENSION(n)                            :: base,off,init,trans
      INTEGER                                         :: trans_s
      REAL(MK),DIMENSION(n,m)                         :: quasi
      INTEGER*8                                       :: seed8
      REAL(MK),DIMENSION(m,n)                         :: vars_transpose
      REAL(MK),DIMENSION(n,m)                         :: vars8
      INTEGER                                         :: alloc_status
      REAL                                            :: tmp_real
      
      IF (.NOT. PRESENT(sampler)) THEN
        m_sampler = 1
      ELSE
        m_sampler = sampler
      END IF
      
      IF (.NOT. PRESENT(inverter)) THEN
        m_inverter = 1
      ELSE
        m_inverter = inverter
      END IF

      
      IF (m_sampler .EQ. 0) THEN
         IF (n .GT. 1111) THEN
           WRITE(*,*) 'Error - SOBOL called for more than 1111 dimensions!'
           STOP
         END IF
          !seed8 = SNGL(seed1)   !seed8 cause of INTERGER*8 comp. problems
	  seed8 = seed1
          DO i = 1, m
            Call i8_sobol ( n, seed8, vars8(:,i))
            vars(:,i) = vars8(:,i)
          END DO
          seed1 = seed8
            tmp_real = MAXVAL(vars)
            IF (tmp_real .GT. 1 ) THEN
                WRITE (*,*) 'INVALID SOBOL VAL'
                STOP
            END IF
      ELSEIF (m_sampler .EQ. 1) THEN
        trans_s = 0
        iflag = 0        
        IF (init_sobol .EQ. 1) THEN
            ALLOCATE(sobol_quasi(n),stat=alloc_status)
            ALLOCATE(sv(n,maxbit),stat=alloc_status)            
        END IF
#ifdef __CMA
        CALL SOBOL(vars_transpose, m, n, sobol_quasi ,ll, count, sv,&
     &     options%qr_scrambling, seed1, init_sobol, trans_s)
#else
        CALL SOBOL(vars_transpose, m, n, sobol_quasi ,ll, count, sv,&
     &     2, seed1, init_sobol, trans_s)
#endif
        init_sobol = 0
        vars = TRANSPOSE (vars_transpose)
        tmp_real = MAXVAL(vars)
         IF (tmp_real .GT. 1) THEN
           tmp_real = MAXVAL(vars)
         END IF
        
        
      ELSEIF (m_sampler .EQ. 2) THEN        
        CALL halton_sequence ( n, m, vars8 ) 
        vars = vars8
      ELSEIF (m_sampler .EQ. 3) THEN
        init = 1
        off = 0
        trans = 0
        CALL HALTON1(vars,m,n,base,off,init,trans)
      ELSEIF (m_sampler .EQ. 4) THEN
         IF (n .GT. 10) THEN
         WRITE(*,*) 'Error - Faure called for more than 10 dimensions!'
         STOP
         END IF
       qs = prime_ge(n)
       IF (init_faure .EQ. 1) THEN
        seed1 = -1
        init_sobol = 0
       END IF
       DO i = 1, m
        !seed8 = SNGL(seed1)
        seed8 = seed1
        CALL faure ( n, seed8, vars(:,i) ) 
        seed1 = seed8
       END DO
      ELSEIF (m_sampler .EQ. 5) THEN
         IF (n .GT. 20) THEN
         WRITE(*,*) 'Error - Niederreiter called for more than 20 dimensions!'
         STOP
         END IF
        DO i = 1,m
            CALL niederreiter2 ( n, seed1, vars(:,i) )
        END DO
      END IF
      
      IF(m_inverter .NE. -1) THEN
        CALL invnor(n,m,vars,vars,m_inverter)
      END IF

      END SUBROUTINE qr_generator
      END MODULE qr_gen_mod
