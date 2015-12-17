	  !-------------------------------------------------------------------------
      !  Subroutine   :          cmaes_freememory
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This Subroutine deallocates Memory
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
      SUBROUTINE cmaes_freememory()
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE cmaes_opts_mod
      USE cmaes_mod
      USE cmaes_out_mod
      USE cmaes_run_mod
      USE CEC2005
      IMPLICIT NONE
      INTEGER		:: deallstat
      
      !-------------------------------------------------------------------------
      !  Deallocate CEC2005 memory
      !-------------------------------------------------------------------------
      IF (options%use_CEC ) THEN 
          CALL close_benchmark()
      ENDIF      
      !-------------------------------------------------------------------------
      !  Deallocate cmaes_opts_mod memory
      !-------------------------------------------------------------------------
      DEALLOCATE(options%DiffMaxChange,stat=deallstat)
      DEALLOCATE(options%DiffMinChange,stat=deallstat)
      DEALLOCATE(options%LBounds,stat=deallstat)
      DEALLOCATE(options%UBounds,stat=deallstat)
      DEALLOCATE(options%insigma,stat=deallstat)
      DEALLOCATE(options%xstart,stat=deallstat)
      IF(deallstat .NE. 0) WRITE(*,*) 'Error deallocating cmaes_opts_mod memory'

      !-------------------------------------------------------------------------
      !  Deallocate cmaes_mod memory
      !-------------------------------------------------------------------------      
      DEALLOCATE(input%xstart,stat=deallstat)
      DEALLOCATE(input%insigma,stat=deallstat)
      DEALLOCATE(xmean,stat=deallstat)
      DEALLOCATE(xold,stat=deallstat)
      DEALLOCATE(maxdx,stat=deallstat)
	  DEALLOCATE(mindx,stat=deallstat)      
	  DEALLOCATE(weights,stat=deallstat)
	  DEALLOCATE(pc,stat=deallstat)
	  DEALLOCATE(ps,stat=deallstat)
	  DEALLOCATE(B,stat=deallstat)
	  DEALLOCATE(D,stat=deallstat)
	  DEALLOCATE(BD,stat=deallstat)
	  DEALLOCATE(C,stat=deallstat)
	  DEALLOCATE(fitness%hist,stat=deallstat)
	  DEALLOCATE(fitness%histsel,stat=deallstat)
	  DEALLOCATE(fitness%raw,stat=deallstat)
	  DEALLOCATE(fitness%sel,stat=deallstat)
	  DEALLOCATE(fitness%idx,stat=deallstat)
	  DEALLOCATE(fitness%idxsel,stat=deallstat)
	  DEALLOCATE(bnd%weights,stat=deallstat)
	  DEALLOCATE(bnd%scales,stat=deallstat)
	  DEALLOCATE(bnd%isbounded,stat=deallstat)
	  DEALLOCATE(bnd%dfithist,stat=deallstat)
	  DEALLOCATE(bnd%arpenalty,stat=deallstat)
	  
	  !-------------------------------------------------------------------------
      !  Deallocate cmaes_out_mod memory
      !-------------------------------------------------------------------------
	  DEALLOCATE(bestever%x,stat=deallstat)
	  DEALLOCATE(out%solutions%bestever%x,stat=deallstat)
	  DEALLOCATE(out%solutions%mean%x,stat=deallstat)
	  DEALLOCATE(out%solutions%recentworst%x,stat=deallstat)
	  DEALLOCATE(out%solutions%recentbest%x,stat=deallstat)
	  DEALLOCATE(out%hist%recentbest%x,stat=deallstat)
	  DEALLOCATE(out%hist%recentworst%x,stat=deallstat)
	  DEALLOCATE(out%hist%mean%x,stat=deallstat)
	  DEALLOCATE(out%histParamArr%diagD,stat=deallstat)
	  DEALLOCATE(out%histParamArr%stds,stat=deallstat)
	  DEALLOCATE(out%histParamArr%Bmax,stat=deallstat)
	  DEALLOCATE(out%histParamArr%Bmin,stat=deallstat)

	  !-------------------------------------------------------------------------
      !  Deallocate pso_mod memory
      !-------------------------------------------------------------------------	
#ifdef __HAVE_MPI__        
	  DEALLOCATE(C_PSO,stat=deallstat)
	  DEALLOCATE(B_rot,stat=deallstat)
	  DEALLOCATE(b_main,stat=deallstat)
      DEALLOCATE(bias,stat=deallstat)
      DEALLOCATE(gb_vec,stat=deallstat)
      DEALLOCATE(gb_vec_xold,stat=deallstat)

	  !-------------------------------------------------------------------------
      !  Deallocate mpi_mod memory
      !-------------------------------------------------------------------------      
      DEALLOCATE(GLOBAL_X_BEST,stat=deallstat)
      DEALLOCATE(stop_proc_mask,stat=deallstat)
	  !DEALLOCATE(stop_proc)
#endif
      RETURN
      END SUBROUTINE cmaes_freememory
