      
#ifdef __SP  
#define __CMA           
#endif
#ifdef __DP      
#define __CMA    
#endif
      !-------------------------------------------------------------------------
      ! Subroutine      :            benchmark_bbob
      !-------------------------------------------------------------------------
      !
      !  Purpose      : a fortran interface to the C implementation of the BBOB
      !                 
      !
      !  Remarks      : 
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


      

      
      SUBROUTINE benchmark_bbob(res,vars_in,m,n,boundarray_down,boundarray_up)
      !parameters

#ifdef __CMA
      USE cmaes_param_mod
      USE cmaes_opts_mod
#endif
      IMPLICIT NONE      


      REAL(MK),DIMENSION(n),INTENT(out)             :: res
      REAL(MK),DIMENSION(m,n),INTENT(IN)            :: vars_in
      
      INTEGER,INTENT(in)                             :: m
      INTEGER,INTENT(in)                             :: n
      REAL(MK),OPTIONAL,DIMENSION(m)      :: boundarray_up
      REAL(MK),OPTIONAL,DIMENSION(m)      :: boundarray_down
      REAL                                :: fgeneric_evaluate
      EXTERNAL                            :: fgeneric_evaluate
      
#ifdef __BBOB      
      res = fgeneric_evaluate(vars_in);
      !CALL benchmark_c(res,vars_in,m,n)
#endif
      
      
      END SUBROUTINE benchmark_bbob
      