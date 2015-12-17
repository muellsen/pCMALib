      !-------------------------------------------------------------------------
      !  Function       :            cmaes_startpoint
      !-------------------------------------------------------------------------
      !
      !  Purpose      : sets the inital Startpoint for CMA
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

      SUBROUTINE cmaes_startpoint(new_start)
      USE cmaes_opts_mod
      USE cmaes_mod
      USE cmaes_param_mod
      USE cmaes_out_mod
#ifdef __HAVE_MPI__      
      INCLUDE "mpif.h"
#endif
      REAL(MK),OPTIONAL,DIMENSION(options%dimensions),INTENT(IN) :: new_start

      !local vars
      CHARACTER(len=10)		:: proc
	  ! Random number generator
      REAL(MK)				:: ZBQLUAB
      !status if its inital initalization
      LOGICAL               :: first_init
      INTEGER               :: pop_mult
      INTEGER							:: allocStat
      INTEGER               :: tmpint


      IF (.NOT. ASSOCIATED(options%xstart)) THEN
        ALLOCATE(options%xstart(options%dimensions),stat=allocStat)
      	IF(allocStat .NE. 0) STOP 'Error allocating options%xstart'
        first_init = .TRUE.
      ELSE
        first_init = .FALSE.
      END IF
      
      
      IF (first_init  .OR. options%restart_type .EQ. 0) THEN
      ! Initialize MPI
#ifdef __HAVE_MPI__      
	  ! Init random number generator
	  IF (first_init) THEN
	  CALL ZBQLINI(options%seed)
	  END IF
	  !CALL SYSTEM_CLOCK(options%seed)
	  
      DO i = 1,options%dimensions
            IF (options%use_init_bounds ) THEN
                options%xstart(i) = ZBQLUAB(options%init_lBounds,&
     &          options%init_uBounds)
            ELSE
                options%xstart(i) = ZBQLUAB(options%LBounds(i),&
     &          options%UBounds(i))
            END IF
      END DO

#else
	  ! Init random number generator	
	  IF (first_init) THEN
	  CALL ZBQLINI(options%seed)
	  END IF
      DO i = 1,options%dimensions
            IF (options%use_init_bounds ) THEN
                options%xstart(i) = ZBQLUAB(options%init_lBounds,&
     &          options%init_uBounds)
            ELSE
                options%xstart(i) = ZBQLUAB(options%LBounds(i),&
     &          options%UBounds(i))
            END IF
      END DO 
#endif
   
      
      
      ELSEIF (options%restart_type .EQ. 1) THEN
                ! We want to restart same location again
                !just changing Seed
                tmpint = options%seed + 1
	            options%seed = tmpint
	            DO i = 1,input%N
      	            input%xstart(i) = options%xstart(i)
                END DO
      ELSEIF (options%restart_type .EQ. 2) THEN
                ! We want to restart at the location where we stopped
                tmpint = options%seed + 1
	            options%seed = tmpint	            
                input%xstart = bestever%x
                
      END IF
      
      
       
      
      IF (.NOT. first_init) THEN        
        !CALL cmaes_reinitialize(options%xstart,options%dimensions)
        
        !if options%MaxIncFac is negative then it can always increase
        IF (options%MaxIncFac .LT. 0 .OR. &
     &           initial_Popsize * options%MaxIncFac & 
     &          .GT. options%Popsize * options%IncPopSize) THEN
            options%Popsize = options%Popsize * options%IncPopSize
        END IF
      DO i = 1,input%N
      	input%insigma(i) = options%insigma(i)
      END DO
        CALL cmaes_init()
        CALL cmaes_initbounds()        
      ELSE
         initial_Popsize = 0
      END IF
      
      END SUBROUTINE cmaes_startpoint
