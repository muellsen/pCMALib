      !-------------------------------------------------------------------------
      !  Routine       :            random_landscape
      !-------------------------------------------------------------------------
      !
      !  Purpose      :  returns random fitness values
      !
      !  Remarks      : 
      !
      !                  
      !  Revisions    :
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  pCMALib: a parallel fortran 90 library for the evolution strategy with
      !           covariance matrix adaptation
      !  Christian L. Mueller, Benedikt Baumgartner, Georg Ofenbeck
      !  MOSAIC group, ETH Zurich, Switzerland
      !-------------------------------------------------------------------------




      SUBROUTINE random_landscape(res,vars,m,n,lbounds,ubounds) 
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

      REAL(MK),DIMENSION(n),INTENT(out)            :: res
      REAL(MK),DIMENSION(m,n),INTENT(in)           :: vars
      INTEGER,INTENT(in)                           :: m
      INTEGER,INTENT(in)                           :: n
      REAL(MK),DIMENSION(m),OPTIONAL              :: lbounds
      REAL(MK),DIMENSION(m),OPTIONAL              :: ubounds
      REAL(MK)                                    :: dummy
      REAL(MK)								:: ZBQLU01
      ! Locally used variables
      REAL(MK),DIMENSION(m)                        :: ub,lb
      LOGICAL,DIMENSION(n)                         :: boundviolation

      
      res(1) = ZBQLU01(dummy) * 300
      
      END SUBROUTINE random_landscape
