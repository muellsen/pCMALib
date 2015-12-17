	  !-------------------------------------------------------------------------
	  !  Subroutine	  :					cmaes_writegen
	  !-------------------------------------------------------------------------
	  ! 
	  !  Purpose	  : writes current generation variables to .txt files in
	  !					directory /out named like the data written.
	  !					Data from subsequent generation is appended at eof.
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
#ifdef __HAVE_MPI__	  
	  SUBROUTINE cmaes_writegen(arx,arxvalid,GLOBAL_X_BEST,m,n,unit)
#else
	  SUBROUTINE cmaes_writegen(arx,arxvalid,m,n,unit)
#endif	  	  
	  !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE cmaes_param_mod
	  USE cmaes_mod
	  USE cmaes_opts_mod
	  USE cmaes_out_mod
	  IMPLICIT NONE
	  
	  !-------------------------------------------------------------------------
      !  Parameters
      !-------------------------------------------------------------------------
	  INTEGER,INTENT(in)							:: m
	  INTEGER,INTENT(in)							:: n
	  REAL(MK),DIMENSION(m,n),INTENT(in)			:: arx
	  REAL(MK),DIMENSION(m,n),INTENT(in)			:: arxvalid
#ifdef __HAVE_MPI__	  
	  REAL(MK),DIMENSION(input%N),INTENT(in)		:: GLOBAL_X_BEST
#endif	  
	  INTEGER,INTENT(in)							:: unit
	  
	  !-------------------------------------------------------------------------
      !  Local Variables
      !-------------------------------------------------------------------------
	  CHARACTER(len=20),SAVE 				:: formN,formHist,formInt,formFit
	  INTEGER								:: i
	  INTEGER,SAVE			 				:: calls = 0
	  INTEGER                               :: histsize
	  
	  !-------------------------------------------------------------------------
      !  Get Array Formatter (only once)
      !-------------------------------------------------------------------------
	  IF(calls .EQ. 0) THEN
	  	CALL tool_formatarrays(formN,input%N)
	  	CALL tool_formatarrays(formFit,lambda)	
	  	CALL tool_formatarrays(formHist,size(fitness%hist)) 
	  	WRITE(formInt,'(I10)') lambda ! Convert lambda to String
	  	formInt = '(3X,' // trim(adjustl(formInt)) // '(I5))'
	  	IF(options%EvalInitialX .AND. options%flgGenData) THEN
	  		DO i = 1, input%N
			  	WRITE(unit,formN) C(i,1:input%N)
			  	WRITE(unit+1,formN) B(i,1:input%N)
			  	WRITE(unit+2,formN) D(i,1:input%N)
			END DO 
			WRITE(unit+3,*) sigma
			WRITE(unit+6,formN) xmean
			calls = calls + 1
			RETURN
	  	END IF
	  END IF
	  
	  IF (options%flgGenData) THEN
	  !-------------------------------------------------------------------------
      !  Write C, B and D
      !-------------------------------------------------------------------------
	  DO i = 1, input%N
	  	WRITE(unit,formN) C(i,1:input%N)
	  	WRITE(unit+1,formN) B(i,1:input%N)
	  	WRITE(unit+2,formN) D(i,1:input%N)
	  END DO 
	  
	  !-------------------------------------------------------------------------
      !  Write some other CMA Parameters
      !-------------------------------------------------------------------------
	  WRITE(unit+3,*) sigma
	  WRITE(unit+4,formN) pc
	  WRITE(unit+5,formN) ps
	  WRITE(unit+6,formN) xmean
	  
	  !-------------------------------------------------------------------------
      !  Write Sampling Coordinates
      !-------------------------------------------------------------------------
	  !DO i = 1, n
	  	!WRITE(unit+7,formN) arx(1:m,i)
	  	!WRITE(unit+8,formN) arxvalid(1:m,i)
	  !END DO
	  
	  DO i = 1, n	  
		WRITE(unit+7,formN) arx(1:m,i)
		WRITE(unit+8,formN) arxvalid(1:m,i)
	  END DO
	  
	  !-------------------------------------------------------------------------
      !  Write Fitness Values
      !-------------------------------------------------------------------------
	  DO i = 1, n
	  WRITE(unit+11,*)   fitness%sel(i)
	  WRITE(unit+9, *)   fitness%raw(i)
	  WRITE(unit+10,*)   fitness%idx(i)  
	  WRITE(unit+12, *)   fitness%idxsel(i)
	  END DO
	  histsize = size(fitness%hist)
	  DO i = 1, histsize
	  WRITE(unit+13,*) fitness%hist(i)
	  WRITE(unit+14,*) fitness%histsel(i)
	  END DO
	  
	  WRITE(unit+15,*) countEval
	  ! Benis code
	  !WRITE(unit+16,*) abs(options%global_min - bestever%f)
	  WRITE(unit+16,*) bestever%f
	  
	  WRITE(unit+17,formN) bestever%x 
#ifdef __HAVE_MPI__	  
	  WRITE(unit+18,formN) GLOBAL_X_BEST
#endif
      WRITE(unit+20,*) n
      WRITE(unit+21,*) histsize
	  ELSEIF (options%flgGenTrace) THEN

      WRITE(unit+16,*) bestever%f
	  WRITE(unit+17,formN) bestever%x
	  
	  

      END IF

	  calls = calls + 1
	  
	  RETURN
	  END SUBROUTINE