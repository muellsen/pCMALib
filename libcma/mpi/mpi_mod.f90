#include "define.h"
 	  !-------------------------------------------------------------------------
      !  Module       :                   mpi_mod
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains variables used in MPI communication
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
      !  Christian L. Mueller, Benedikt Baumgartner, Georg Ofenbeck, Christian L. Mueller
      !  MOSAIC group, ETH Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE mpi_mod
      USE cmaes_param_mod
      IMPLICIT NONE
      SAVE
	  !-------------------------------------------------------------------------
	  !  Interfaces
	  !-------------------------------------------------------------------------	

	  !-------------------------------------------------------------------------
	  !  Variables
	  !-------------------------------------------------------------------------	  
	  !  Variables to be exchanged and updated between runs
	  REAL(MK),DIMENSION(2)						:: fbest
	  REAL(MK),ALLOCATABLE,DIMENSION(:)			:: GLOBAL_X_BEST
	  REAL(MK),DIMENSION(2)						:: GLOBAL_F_BEST
	  REAL(MK),DIMENSION(2)						:: GLOBAL_F_BEST_LAST
	  INTEGER									:: GLOBAL_FUN_EVALS
	  INTEGER,ALLOCATABLE,DIMENSION(:)			:: F_BEST_RANK
	  !REAL(MK),ALLOCATABLE,DIMENSION(:,:)		:: GLOBAL_BEST_EVs
	  
	  !  Variables used to stop runs and update communicator
	  LOGICAL									:: STOP_FLAG, stop_me
	  LOGICAL,ALLOCATABLE,DIMENSION(:)			:: stop_proc_mask
	  INTEGER,ALLOCATABLE,DIMENSION(:) 			:: stop_proc
	  INTEGER									:: stop_count
	  INTEGER									:: j, new_rank, group_world  
	  INTEGER									:: group_worker, comm_worker
	  INTEGER									:: comm_size
	  
	  !  Used for output after 10^3,10^4 and 10^5 FES
	  INTEGER									:: writeCount
      CHARACTER(len=10)                                     :: proc
	  

	  
	  
	  !-------------------------------------------------------------------------
	  !  Module Procedures
	  !-------------------------------------------------------------------------
	  CONTAINS
      SUBROUTINE mpi_comm_init(N,bestF,bestX)
      INCLUDE "mpif.h" 
	  INTEGER,PARAMETER	:: MK = KIND(1.D0)
	  !-------------------------------------------------------------------------
      !  Parameters
      !-------------------------------------------------------------------------	
	  INTEGER,INTENT(in)					:: N
	  REAL(MK),INTENT(in)					:: bestF
	  REAL(MK),DIMENSION(N),INTENT(in)		:: bestX
	  	
	  !-------------------------------------------------------------------------
      !  Local Variables
      !-------------------------------------------------------------------------
	  INTEGER								:: allocStat, mpiError
	  	
	  	
	  ! Initialize GLOBAL_X_BEST and local fbest
	  ALLOCATE(GLOBAL_X_BEST(N),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating GLOBAL_X_BEST'
	  GLOBAL_X_BEST = bestX
	  fbest(1) = bestF
	  fbest(2) = MY_RANK
	  
	  ! Initialize F_BEST_RANK
	  ALLOCATE(F_BEST_RANK(NUM_CMA_RUNS),stat=allocStat)
	  IF(allocStat .NE. 0) STOP 'Error allocating F_BEST_RANK'
	  F_BEST_RANK = huge(allocStat)
	  
	  ! Initialize GLOBAL_BEST_EVs
!	  ALLOCATE(GLOBAL_BEST_EVs(N,2),stat=allocStat)
!	  IF(allocStat .NE. 0) STOP 'Error allocating GLOBAL_BEST_EVs'
!	  GLOBAL_F_BEST = 0.0_MK
   	  
   	  ! Initialize variables for stopping processes
	  ALLOCATE(stop_proc_mask(NUM_CMA_RUNS),stat=allocStat)
	  IF(allocStat .NE. 0) STOP 'Error allocating stop_proc_mask'
	  stop_proc_mask = .FALSE.
	  STOP_FLAG = .FALSE.
	  	
	  ! writeCount is used for output after 10^3,10^4 and 10^5 FES
	  writeCount = 0
	  	
	  ! convert process number to string
	  WRITE(proc,'(I10)') MY_RANK
	  proc = '_' // trim(adjustl(proc)) 
	  	
	  
	  ! initialize communication group
	  CALL MPI_COMM_RANK(MPI_COMM_WORLD, new_rank, mpiError)
	  CALL MPI_COMM_GROUP(MPI_COMM_WORLD, group_world, mpiError)
	  CALL MPI_COMM_DUP(MPI_COMM_WORLD, comm_worker, mpiError)
	  group_worker = group_world
	  comm_size = NUM_CMA_RUNS
	  
	  RETURN	
	  END SUBROUTINE
 	  !-------------------------------------------------------------------------
      !  Subroutine   :                   mpi_comm_adapt
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Subroutine adapts the MPI communicator if any CMA run
      !					wants to stop
      !
      !	 Input/Output : stop_me	(L) true if current process was excluded
      !									successfully from communications
      !                
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: mpi_comm_adapt.f90,v $
      !  Revision 1.3  2008/04/18 13:02:56  paulb
      !  repository update
      !
	  !
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  pCMALib: a parallel fortran 90 library for the evolution strategy with
      !           covariance matrix adaptation
      !  Christian L. Mueller, Benedikt Baumgartner, Georg Ofenbeck, Christian L. Mueller
      !  MOSAIC group, ETH Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE mpi_comm_adapt()
      INCLUDE "mpif.h"
      INTEGER					:: ioError,i,j
      
      CALL MPI_ALLGATHER(STOP_FLAG,1,MPI_LOGICAL,stop_proc_mask,1,MPI_LOGICAL,&
	  														comm_worker,ioError)
	  														
	  IF(any(stop_proc_mask)) THEN	! If process wants to stop...

  	  	stop_count = 0
		DO j = 1,NUM_CMA_RUNS	! ...get number of stopping processes...
			IF(stop_proc_mask(j)) THEN
		  		stop_count = stop_count + 1
		  	END IF
		END DO
		IF(stop_count .NE. 0) THEN	! ...and get process number
			ALLOCATE(stop_proc(stop_count))
	  		i = 1
	  		DO j = 1,NUM_CMA_RUNS
	  			IF(stop_proc_mask(j)) THEN
	  				stop_proc(i) = j-1
	  				i = i + 1
	  			END IF
	  		END DO
	  		! stop_proc is now an array containing the ranks of processes to be
	  		! stopped
			! Due to the MPI specifications, excluding a process from
			! communications is done via creating a sub-group first,
			! and create a new communicator, based upon that group,
			! afterwards.
	  		CALL MPI_COMM_GROUP(comm_worker, group_worker, ioError)
			CALL MPI_GROUP_EXCL(group_worker, stop_count, stop_proc,&
														group_worker,ioError)
			CALL MPI_GROUP_RANK(group_worker, new_rank, ioError)
			! Processes not within the group will be assigned the value
			! of MPI_UNDEFINED in new_rank
			CALL MPI_COMM_CREATE(comm_worker, group_worker,&
														comm_worker, ioError)
			
			IF(new_rank .NE. MPI_UNDEFINED) THEN
				CALL MPI_COMM_RANK(comm_worker,new_rank,ioError)
				CALL MPI_COMM_SIZE(comm_worker,comm_size,ioError)		
				DEALLOCATE(F_BEST_RANK)
				ALLOCATE(F_BEST_RANK(comm_size))												
			ELSE
				!EXIT genLoop
				stop_me = .TRUE.
				DEALLOCATE(F_BEST_RANK)
			END IF
		 	DEALLOCATE(stop_proc)
  	  	  END IF
		  	stop_proc_mask = .FALSE.
		END IF
		
		RETURN  	
		END SUBROUTINE  	
 	  !-------------------------------------------------------------------------
      !  Subroutine   :                   mpi_rank
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This Routine ranks the value of vec(1) in the current
      !					communicator and returns the corresponding indices 
      !					(vec(2)) in ascending order in idx
      !
      !  Input		  :	vec(:)		(R) Vector of size 2 with the value to rank 
      !									(vec(1)) and the corresponding mpi 
      !									process rank (vec(2))
      !					comm		(I) MPI communicator
      !					comm_size	(I) Size of communicator (#processes)
      !
      !	 Input/Output : idx(:)		(I) Integer array of size comm_size
      !									containing process ranks in ascending
      !									order
      !                
      !  Remarks      : Example: #proc=2,
      !						proc1: vec(1) = 8.0, vec(2) = my_rank = 0,
      !						proc2: vec(1) = 5.0, vec(2) = my_rank = 1
      !					=>	idx = [1,0]
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: mpi_rank.f90,v $
      !  Revision 1.3  2008/04/18 13:02:56  paulb
      !  repository update
      !
      !
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  pCMALib: a parallel fortran 90 library for the evolution strategy with
      !           covariance matrix adaptation
      !  Christian L. Mueller, Benedikt Baumgartner, Georg Ofenbeck, Christian L. Mueller
      !  MOSAIC group, ETH Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE mpi_rank(vec,idx,comm,comm_size)
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE cmaes_param_mod
      USE tool_mrgrnk_mod
      INCLUDE "mpif.h"
      !-------------------------------------------------------------------------
      !  Parameters
      !-------------------------------------------------------------------------
      INTEGER, PARAMETER							:: msg_size = 2
      REAL(MK),DIMENSION(msg_size),INTENT(in)		:: vec
      INTEGER,INTENT(in)							:: comm,comm_size
      INTEGER,DIMENSION(comm_size),INTENT(inout)	:: idx
      
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER										:: ierr,allocStat,i
      REAL(MK),ALLOCATABLE,DIMENSION(:,:)			:: recv
      
      
	  ! Initialize receive buffer
      ALLOCATE(recv(msg_size,comm_size),stat=allocStat)
      IF(allocStat .NE. 0) &
      			WRITE(*,*) 'Rank ', MY_RANK, ': Error allocating receive buffer'
      
      !IF(comm_size .GT. 1) THEN
      	! gather values in recv
      	CALL MPI_ALLGATHER(vec,msg_size,MPI_DOUBLE_PRECISION,recv,msg_size,&
      											MPI_DOUBLE_PRECISION,comm,ierr)
      											
	  	! rank values of recv in ascending order and get corresponding index      											
      	CALL mrgrnk(recv(1,:),idx)
      
      	!IF(recv(1,idx(1)) .LT. GLOBAL_F_BEST(1)) THEN
      		GLOBAL_F_BEST(1) = recv(1,idx(1))
      		DO i = 1, comm_size
      			idx(i) = idx(i) - 1
      		END DO
      
      		GLOBAL_F_BEST(2) = idx(1)
      	!END IF	
      
!      ELSE	! do this, if only single process is running
!      		IF(vec(1) .LE. GLOBAL_F_BEST(1)) THEN
!      			GLOBAL_F_BEST(1) = vec(1)
!      			GLOBAL_F_BEST(2) = MY_RANK
!      			idx(1) = MY_RANK
!      		END IF
!      END IF		
      
      ! Deallocate memory
      DEALLOCATE(recv,stat=allocStat)
      IF(allocStat .NE. 0) &
      		WRITE(*,*) 'Rank ', MY_RANK, ': Error deallocating receive buffer'      
      
      RETURN
      END SUBROUTINE mpi_rank
!#include "mpi_shift_worst.f90"
      END MODULE mpi_mod
