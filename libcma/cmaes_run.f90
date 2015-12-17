      !-------------------------------------------------------------------------
      !  Subroutine   :          cmaes_run
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This Subroutine implements the CMA Generation loop
      !                
      !  Remarks      :
      !
      !  References   : 
      !					@article{Hansen:2007,
      !					Author = {Hansen, Nikolaus},
      !					Keywords = {CMA},
      !					Title = {{The CMA Evolution Strategy: A Tutorial}},
      !					Year = {2007}}
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  pCMALib: a parallel fortran 90 library for the evolution strategy with
      !           covariance matrix adaptation
      !  Christian L. Mueller, Benedikt Baumgartner, Georg Ofenbeck
      !  MOSAIC group, ETH Zurich, Switzerland
      !-------------------------------------------------------------------------
      SUBROUTINE cmaes_run(fitfun)
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE cmaes_run_mod
      USE lbfgsb_mod
      IMPLICIT NONE
      
#ifdef __HAVE_MPI__
      INCLUDE "mpif.h"
#endif	  

      !-------------------------------------------------------------------------
      !  Interfaces
      !-------------------------------------------------------------------------
      INTERFACE
          SUBROUTINE cmaes_xintobounds(LBounds,UBounds,N,x,xout,idx)
          USE cmaes_param_mod
          USE cmaes_mod,only:countOutOfBounds
          IMPLICIT NONE
          INTEGER, INTENT(in)			  					:: N
          REAL(MK),DIMENSION(N),INTENT(in)    				:: LBounds,UBounds
          REAL(MK),DIMENSION(N),INTENT(in) 					:: x
          REAL(MK),DIMENSION(N),INTENT(out)					:: xout
          LOGICAL,DIMENSION(N),INTENT(out),OPTIONAL			:: idx
          END SUBROUTINE
      END INTERFACE
      
            
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      EXTERNAL 											:: fitfun
      !-------------------------------------------------------------------------
      !  Functions
      !-------------------------------------------------------------------------
      REAL(MK)											:: ZBQLNOR 	!RNG
      !REAL(MK)											:: ZBQLU01
      REAL(MK)											:: cmaes_funcwrap
      REAL(MK)											:: tool_myrange
      !-------------------------------------------------------------------------
      !  Local Variables
      !-------------------------------------------------------------------------
      INTEGER                                           :: info
      REAL(MK),ALLOCATABLE,DIMENSION(:,:)    :: arz, arx, arxvalid
      REAl(MK),DIMENSION(:,:),ALLOCATABLE               :: darz
      REAL(MK),DIMENSION(:,:),ALLOCATABLE               :: tempMat
      REAL(MK),DIMENSION(input%N)						:: diag
      INTEGER											:: i, k, tries
      REAL(MK),DIMENSION(input%N)						:: zmean
      REAL(MK)											:: hsig
      REAL(MK),DIMENSION(input%N,1)						:: colVec
      REAL(MK),DIMENSION(1,input%N)						:: rowVec
      REAL(MK),ALLOCATABLE,DIMENSION(:,:)               :: weightsMat
      REAL(MK),DIMENSION(input%N,input%N)				:: triuC
      REAL(MK),DIMENSION(input%N,input%N)				:: tempMat2
      LOGICAL,DIMENSION(input%N,input%N)				:: mask
      REAL(MK)											:: temp
      REAL(MK),DIMENSION(input%N)						:: xmin
      REAL(MK)											:: fmin
      INTEGER											:: ioError
      INTEGER											:: uNum
      CHARACTER(len=20)									:: outFormat
      CHARACTER(len=30)									:: tmp_char,tmp_char2
      CHARACTER(len=200)                                :: tmp_dir 
      REAL(MK),DIMENSION(:),ALLOCATABLE					:: newSort
      LOGICAL											:: file_write=.TRUE.
      LOGICAL                                           :: dir_exist
      INTEGER                                           :: allocStat
      INTEGER                                           :: tmpint
      REAL(MK)                                          :: tmpreal  
      REAL(KIND(1.E0))                                  :: probability = 0e0
      REAL(KIND(1.E0))		   							:: chi_out = 0e0
      
#ifndef __HAVE_MPI__
      CHARACTER(len=10)                                     :: proc
#endif
      

      ! Temporary variables used for matrix multiplication
      REAL(MK),DIMENSION(input%N) 						:: multmp_vN	  
      REAL(MK),DIMENSION(input%N,input%N)				:: multmp_mN, multmp_mN2	
      REAL(MK),DIMENSION(:,:),ALLOCATABLE :: multmp_mMU	
      REAL(MK),DIMENSION(:,:),ALLOCATABLE :: mutmp
      REAL(MK),DIMENSION(1)                             :: posInf_array 

      !------------------------------------------------------------------------
      ! Allocate local variables
      !------------------------------------------------------------------------
      ALLOCATE(arz(input%N,lambda),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating arz'
      ALLOCATE(arx(input%N,lambda),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating arx'
      ALLOCATE(arxvalid(input%N,lambda),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating arxvalid'
      ALLOCATE(darz(input%N,lambda),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating darz'
      ALLOCATE(tempMat(input%N,mu),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating tempMat'
      ALLOCATE(weightsMat(mu,mu),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating weightsMat'
      ALLOCATE(newSort(lambda),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating newSort'
      ALLOCATE(multmp_mMU(mu,input%N),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating multmp_mMU'
      ALLOCATE(mutmp(input%N,mu),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating mutmp'
      !-------------------------------------------------------------------------
      !  Initialize some variables
      !-------------------------------------------------------------------------	  
      xold = xmean
      arx = 0.0_MK
      arxvalid = 0.0_MK
      xmin = posInf
      fmin = posInf
      uNum = 26
      posInf_array(1) = posInf   
      probability = SNGL(options%BFGS_confprop)
      CALL CHSPPF(probability,input%N,chi_out)   
      chi = REAL (chi_out,MK)
	  !-------------------------------------------------------------------------
      !  Initialize array for the Convergence History
      !-------------------------------------------------------------------------	  
      tmpint = options%StopMaxFunEvals / options%record_modulo
      IF (options%record_besthist ) THEN
         ALLOCATE(besthist(tmpint),stat=allocStat)
         besthist(:) = 0.0_MK
      END IF
      
      

#ifdef __HAVE_MPI__
      !-------------------------------------------------------------------------
      !  Initialize MPI Communication variables (see mpi_mod)
      !-------------------------------------------------------------------------	  
      CALL mpi_comm_init(input%N,bestever%f,bestever%X)

      !-------------------------------------------------------------------------
      !  Initialize Particle Swarm Update variables
      !-------------------------------------------------------------------------
      CALL pso_init(input%N)	
#else
      ! variable from mpi_mod
      proc = ''

#endif	  	  	  

      !-------------------------------------------------------------------------
      !  Create Output Folder if not existing
      !-------------------------------------------------------------------------
      
#ifdef __HAVE_MPI__
      IF (MY_RANK .EQ. 0) THEN
#endif
      inquire(FILE=trim(adjustl(options%output_folder)),exist=dir_exist)
      tmp_dir = 'mkdir '  //  trim(adjustl(options%output_folder))
      if (dir_exist.eqv.(.FALSE.)) then
             CALL system(tmp_dir)
      endif
#ifdef __HAVE_MPI__
      END IF
      CALL MPI_BARRIER(MPI_COMM_WORLD, ioError)
#endif
  
      
       
       
      !-------------------------------------------------------------------------
      !  Save the used seed if other output is generated
      !-------------------------------------------------------------------------
#ifdef __HAVE_MATLAB__
      IF (.TRUE.) THEN !always if .mat is generated (seed not in .mat)
#else
      IF (options%flgouttxt .OR. options%flgGenData .OR. options%flgGenTrace) &
     &THEN !if any textfiles are written
#endif
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/seed'& 
#ifdef __HAVE_MPI__
     & // '_' &
#endif
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
#ifdef __HAVE_MPI__
      OPEN(unit=uNum+20+MY_RANK,file=tmp_dir,status='replace',&
#else
      OPEN(unit=uNum+20,file=tmp_dir,status='replace',&
#endif
     & action='write', iostat=ioError)

      
      IF(ioError .NE. 0) THEN
            WRITE(*,*) 'I/O-Error writing Seed Data.'
            STOP
      END IF		
      
#ifdef __HAVE_MPI__
            WRITE(uNum+20+MY_RANK,*) options%seed
#else
      WRITE(uNum+20,*) options%seed
#endif
	  IF(ioError .NE. 0) THEN
            WRITE(*,*) 'I/O-Error writing Seed Data.'
            STOP
      END IF

#ifdef __HAVE_MPI__
      CLOSE(unit=uNum+20+MY_RANK)
#else
      CLOSE(unit=uNum+20)
#endif
      END IF
      !-----------------------------------------------------------------------
      
      
     
      IF(options%flgGenTrace .AND. .NOT. options%flgGenData) THEN
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/besteverF'&
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum+16,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/besteverX'&
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum+17,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)
        IF(ioError .NE. 0) THEN
            WRITE(*,*) 'I/O-Error writing Generation Trace.'
            STOP
        END IF
      END IF
       
      !-------------------------------------------------------------------------
      !  Open output files (only if flgGenData == .TRUE.)
      !-------------------------------------------------------------------------	              
      
      

      IF(options%flgGenData) THEN
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/C'  // trim(adjustl(proc))&
     & // '.txt'
      OPEN(unit=uNum,file=tmp_dir,status='replace',&
     &action='write', iostat=ioError)
      tmp_dir=trim(adjustl(options%output_folder)) &
     & // '/B'&
     & // trim(adjustl(proc)) //  '.txt'
      OPEN(unit=uNum+1,file=tmp_dir,status='replace',&
     &action='write', iostat=ioError)
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/D' &
     &      // trim(adjustl(proc))   //  '.txt'
      OPEN(unit=uNum+2,file=tmp_dir,&
     &     status='replace', action='write', iostat=ioError)	
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/sigma' &
     &   // trim(adjustl(proc))   //  '.txt'
      OPEN(unit=uNum+3,file=tmp_dir,&
     & status='replace', action='write', iostat=ioError)	
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/pc' &
     &      // trim(adjustl(proc))   //  '.txt'
      OPEN(unit=uNum+4,file=tmp_dir,&
     &     status='replace', action='write', iostat=ioError)
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/ps'&
     &      // trim(adjustl(proc))   //  '.txt'
      OPEN(unit=uNum+5,file=tmp_dir,&
     &     status='replace', action='write', iostat=ioError)	
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/xmean' &
     &      // trim(adjustl(proc))   //  '.txt'
      OPEN(unit=uNum+6,file=tmp_dir,&
     &     status='replace', action='write', iostat=ioError)	
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/arx'&
     &      // trim(adjustl(proc))   //  '.txt'
      OPEN(unit=uNum+7,file=tmp_dir,&
     &     status='replace', action='write', iostat=ioError)	
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/arxvalid'&
     &      // trim(adjustl(proc))   //  '.txt'
      OPEN(unit=uNum+8,file=tmp_dir,&
     &     status='replace', action='write', iostat=ioError)	
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/fitraw'&
     &      // trim(adjustl(proc))   //  '.txt'
      OPEN(unit=uNum+9,file=tmp_dir,&
     &     status='replace', action='write', iostat=ioError)	
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/fitidx'&
     &      // trim(adjustl(proc))   //  '.txt'
      OPEN(unit=uNum+10,file=tmp_dir,&
     &     status='replace', action='write', iostat=ioError)
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/fitsel'&
     &      // trim(adjustl(proc))   //  '.txt'
      OPEN(unit=uNum+11,file=tmp_dir,&
     &     status='replace', action='write', iostat=ioError)	      											
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/fitidxsel' &
     &       // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum+12,file=tmp_dir,status='replace',&
     &action='write', iostat=ioError)	      											
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/fithist'&
     &      // trim(adjustl(proc))   //  '.txt'
      OPEN(unit=uNum+13,file=tmp_dir,&
     &     status='replace', action='write', iostat=ioError)	      											
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/fithistsel'& 
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum+14,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)	
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/Fevals'&
     &      // trim(adjustl(proc))   //  '.txt'
      OPEN(unit=uNum+15,file=tmp_dir,&
     &     status='replace', action='write', iostat=ioError)
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/besteverF'& 
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum+16,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/besteverX'& 
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum+17,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/GLOBAL_X_BEST'&
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum+18,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)	
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/beste345'& 
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum+19,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/popsize'& 
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum+20,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)
      tmp_dir=trim(adjustl(options%output_folder))&
     & // '/histsize'& 
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(unit=uNum+21,file=tmp_dir,status='replace',&
     & action='write', iostat=ioError)
        IF(ioError .NE. 0) THEN
            WRITE(*,*) 'I/O-Error writing Generation Data.'
            STOP
        END IF

            
      
      END IF
      
      
      !-------------------------------------------------------------------------
      !  Write output files for generation 0
      !-------------------------------------------------------------------------
      IF(options%flgGenData) THEN
#ifdef __HAVE_MPI__      
        CALL cmaes_writegen(arx,arxvalid,GLOBAL_X_BEST,input%N,lambda,uNum)
#else
        CALL cmaes_writegen(arx,arxvalid,input%N,lambda,uNum)
#endif      	
      END IF

     
      
      

		
      countRestarts = 0
      
      !-------------------------------------------------------------------------
      !  START CMAES-RUNS (only done once if run without restart flag set)
      !-------------------------------------------------------------------------
      restartLoop: DO WHILE(stopflag .EQ. '')


      !-------------------------------------------------------------------------
      !  START GENERATION LOOP
      !-------------------------------------------------------------------------
      genLoop: DO WHILE(stopflag .EQ. '')      
      
      countIter = countIter + 1
      
      fitness%raw = posInf

      !-------------------------------------------------------------------------
      !  Loop as long as there are NaN values in fitness%raw (cmaes.m 689ff)
      !-------------------------------------------------------------------------
      NaNexists: DO WHILE(1 .EQ. 1)	! is intercepted if k == 0 (line 85)
        !-----------------------------------------------------------------------
        !  Find NaN values in fitness%raw
        !-----------------------------------------------------------------------
        findNaN:DO i = 1, lambda
            k = 0
            IF(fitness%raw(i) .EQ. huge(fitness%raw(1))) THEN
                k = i
                EXIT findNaN
            END IF
        END DO findNaN
        
        !-----------------------------------------------------------------------
        !  Exit outer loop if no more NaN values were found in fitness%raw
        !-----------------------------------------------------------------------
        IF(k .EQ. 0) EXIT NaNexists
          
        !-----------------------------------------------------------------------
        !  Resample until fitness is not NaN
        !-----------------------------------------------------------------------
        tries = 0
        DO WHILE(fitness%raw(k) .EQ. posInf)
            !-------------------------------------------------------------------
            ! Create random number
            !-------------------------------------------------------------------
                        
            IF (options%qr_sampling ) THEN
                CALL qr_generator(input%N,1,options%seed,darz(:,k),&
     &                                  options%qr_sampler,options%qr_inverter)
            ELSE
                DO i = 1, input%N
                darz(i,k) = ZBQLNOR(0.0d0,1.0d0)
!arz(i,k) = matlabran(i,noel)
                END DO
!noel = noel + 1      	  	
            END IF
            
            DO i = 1, input%N
#ifdef __SP      	  		! Little work-around as
                arz(i,k) = sngl(darz(i,k))		! ZBQLNOR is double precision
#else
                arz(i,k) = darz(i,k)
#endif
            END DO
            
            !-------------------------------------------------------------------
            ! Create sample point
            !-------------------------------------------------------------------
            !arx(:,k) = xmean + sigma * matmul(BD,arz(:,k))
#ifdef __SP
            CALL SGEMV('N',input%N,input%N,1.0,BD,input%N,arz(:,k),&
                                                            1,0.0,multmp_vN,1)
#else      		
            CALL DGEMV('N',input%N,input%N,1.0d0,BD,input%N,arz(:,k),&
                                                            1,0.0d0,multmp_vN,1)
#endif
            arx(:,k) = xmean + sigma * multmp_vN
            
            !-------------------------------------------------------------------
            ! Check if sample point is valid (within given bounds)
            !-------------------------------------------------------------------
            IF(bnd%isactive) THEN
              CALL cmaes_xintobounds(options%LBounds,options%UBounds,&
                    input%N,arx(:,k),arxvalid(:,k))
            ELSE
                arxvalid(:,k) = arx(:,k)
            END IF	  
            
            !-------------------------------------------------------------------
            ! Replace Sample point / Sample Fitness via BFGS
            !-------------------------------------------------------------------
           IF (options%BFGS_use) THEN
              !.AND. &
              !((mod(countIter-lastRestart,options%VerboseModulo).EQ.0) .OR. &
              !((countIter-lastRestart).EQ.1))) THEN

            IF (options%BFGS_position .EQ. 1) THEN
              CALL lbfgsb(fitfun,input%N,17,arxvalid(:,k),fitness%raw(k))
              arx(:,k) = arxvalid(:,k)
            ELSE
              arx(:,k) = arxvalid(:,k)
              CALL lbfgsb(fitfun,input%N,17,arx(:,k),fitness%raw(k))
              arx(:,k) = arxvalid(:,k)
            END IF
            !countEval = countEval + 1
            countEval = countEval + 1!countBFGSEval
           !END IF
           ELSE
            !-------------------------------------------------------------------
            ! Otherwise Evaluate function at valid sample point
            !-------------------------------------------------------------------
            fitness%raw(k) = cmaes_funcwrap(arxvalid(:,k),input%N,fitfun)
           END IF
            tries = tries + 1
            !-------------------------------------------------------------------
            ! Check if function value NaN
            !-------------------------------------------------------------------
            IF(fitness%raw(k) .EQ. huge(fitness%raw(1))) &
                countEvalNaN = countEvalNaN + 1
            IF(mod(tries,1000) .EQ. 0) &
                WRITE(*,*) 'NaN objective function values at evaluation ',&
                            countEval
            !WRITE(*,*) 'tries: ', tries
        END DO	! End resampling until != NaN
                        
        
        
        !-----------------------------------------------------------------------
        !  Save bestever%f after 1.e3,1.e4 and 1.e5 iterations.
        !  Save #FES to reach given accuracy
        !  See Suganthan Benchmark Paper
        !-----------------------------------------------------------------------
        IF(options%benchmark) THEN
#ifdef __HAVE_MPI__
        CALL MPI_ALLREDUCE(countEval,GLOBAL_FUN_EVALS,1,MPI_INTEGER,&
                                             MPI_SUM,comm_worker,ioError)
        IF((writeCount .EQ. 0) .AND. (GLOBAL_FUN_EVALS .GE. 1.E3)) THEN
            WRITE(uNum+18,*) abs(options%global_min - bestever%f)
            writeCount = 1
        ELSE IF((writeCount .EQ. 1) .AND. (GLOBAL_FUN_EVALS .GE. 1.E4)) THEN
            WRITE(uNum+18,*) abs(options%global_min - bestever%f)
            writeCount = 2
        ELSE IF((writeCount .EQ. 2) .AND. (GLOBAL_FUN_EVALS .GE. 1.E5)) THEN
            WRITE(uNum+18,*) abs(options%global_min - bestever%f)
            writeCount = 3
        END IF	
        
        IF((abs(options%global_min-bestever%f) .LE. options%record_accuracy) &
     &      .AND. (file_write)) THEN
            tmp_dir=trim(adjustl(options%output_folder)) & 
     &      // '/acc'  //  trim(adjustl(proc)) //  '.txt'
            OPEN(unit=73,file=tmp_dir,status='replace', action='write',&
     &      iostat=ioError)
            WRITE(73,*) GLOBAL_FUN_EVALS
            CLOSE(unit=73)
            file_write = .FALSE.
        END IF								 
#else      	
        IF((countEval .EQ. 1.E3) .OR. &
           (countEval .EQ. 1.E4) .OR. &
           (countEval .EQ. 1.E5)) THEN
#ifdef __HAVE_MATLAB__
            IF (countEval .EQ. 1.E3) THEN
                be3 = abs(options%global_min - bestever%f)
            ELSEIF (countEval .EQ. 1.E4) THEN
                be4 = abs(options%global_min - bestever%f)
            ELSE
                be5 = abs(options%global_min - bestever%f)
            END IF 
#endif
            IF (options%flgGenData ) THEN
                WRITE(uNum+18,*) abs(options%global_min - bestever%f)
            END IF
        END IF
        
        !-----------------------------------------------------------------------
        !  Save bestever history for the convergence graph        
        !  See Suganthan Benchmark Paper
        !-----------------------------------------------------------------------
        IF (options%record_besthist ) THEN
            tmpint = countEval / options%record_modulo 
            IF ( tmpint .GT. 0) THEN
                IF ( besthist(tmpint) .EQ. 0 ) THEN
                   besthist(tmpint) = bestever%f
                END IF
            END IF
        END IF

        !-----------------------------------------------------------------------
        !  
        !  Save #FES to reach given accuracy
        !  See Suganthan Benchmark Paper
        !-----------------------------------------------------------------------

       IF( abs(options%global_min-bestever%f) .LE. options%record_accuracy &
     &     .AND. (file_write)) THEN
         file_write = .FALSE.
#ifdef __HAVE_MATLAB__
            acc_evals = countEval
#endif       
            IF (options%flgGenData ) THEN
                tmp_dir=trim(adjustl(options%output_folder)) & 
     &          // '/acc.txt'
                OPEN(unit=73,file=tmp_dir,status='replace', action='write',&
     &          iostat=ioError)
                WRITE(73,*) countEval
                CLOSE(unit=73)
            END IF
        END IF
#endif
        

    
        
        END IF	! End benchmark
        
      END DO NaNexists	! End find NaN


      fitness%sel = fitness%raw

      !------------------------------------------------------------------------
      !  Handle boundaries (add penalties to fitness%raw)
      !------------------------------------------------------------------------
      IF(bnd%isactive) CALL cmaes_handlebounds(arxvalid,arx,input%N,lambda)
 

      !-------------------------------------------------------------------------
      !  Sort by fitness (idx/idxsel being indices indicating orig. positions)
      !-------------------------------------------------------------------------
      newSort = fitness%sel
      CALL mrgrnk(newSort,fitness%idxsel)
      DO i = 1, lambda
        fitness%sel(i) = newSort(fitness%idxsel(i))
      END DO
      newSort = fitness%raw
      CALL mrgrnk(newSort,fitness%idx)
      DO i = 1, lambda
        fitness%raw(i) = newSort(fitness%idx(i))
      END DO
      
      !-------------------------------------------------------------------------
      ! Record short history of best fitness values
      !-------------------------------------------------------------------------
      DO i = size(fitness%hist), 2, -1	! Move array values to the right
        fitness%hist(i) = fitness%hist(i-1)
        fitness%histsel(i) = fitness%histsel(i-1)
      END DO
      fitness%hist(1) = fitness%raw(1) 	! And save latest best fitness value
      fitness%histsel(1) = fitness%sel(1) ! at position 1
    
      !-------------------------------------------------------------------------
      !  SELECTION AND RECOMBINATION (calculate new xmean)
      !-------------------------------------------------------------------------
      xold = xmean
      DO i = 1, mu		! Extract mu best columns of arx and save it to tempMat
        tempMat(:,i) = arx(:,int(fitness%idxsel(i)))
      END DO
      !xmean = matmul(tempMat,weights)
      
      ! xmean = tempMat * weights	
#ifdef __SP
      CALL SGEMV('N',input%N,mu,1.0,tempMat,input%N,weights,1,0.0,xmean,1)
#else      		
      CALL dgemv('N',input%N,mu,1.0d0,tempMat,input%N,weights,1,0.0d0,xmean,1)
#endif

      
      DO i = 1, mu		! Do the same for arz/zmean
        tempMat(:,i) = arz(:,int(fitness%idxsel(i)))
      END DO
      !zmean = matmul(tempMat,weights)
      
      ! xmean = tempMat * weights
#ifdef __SP
      CALL SGEMV('N',input%N,mu,1.0,tempMat,input%N,weights,1,0.0,zmean,1)
#else      		
      CALL dgemv('N',input%N,mu,1.0d0,tempMat,input%N,weights,1,0.0d0,zmean,1)
#endif
      
      IF(mu .EQ. 1) THEN
        fmean = fitness%sel(1)
      ELSE
        fmean = posInf
      END IF
      
      !-------------------------------------------------------------------------
      !  UPDATE EVOLUTION PATHS (cumulation)
      !-------------------------------------------------------------------------
      !multmp_vN = B * zmean
#ifdef __SP
      CALL SGEMV('N',input%N,input%N,1.0,B,input%N,zmean,1,0.0,multmp_vN,1)
#else      		
      CALL dgemv('N',input%N,input%N,1.0d0,B,input%N,zmean,1,0.0d0,multmp_vN,1)
#endif  
      ps = (1._MK-cs)*ps + (sqrt(cs*(2._MK-cs)*mueff)) * multmp_vN! Eq.40 
                                                                ! in Hansen:2007
      IF(sqrt(sum(ps**2))/sqrt(1._MK-(1.-cs)**(2._MK*real((countIter-lastRestart))))/chiN < &
        1.4_MK + 2._MK/real(input%N+1)) THEN
        hsig = 1._MK
      ELSE
        hsig = 0._MK
      END IF
                                                                ! Eq.42
      pc = (1-cc)*pc + hsig*(sqrt(cc*(2._MK-cc)*mueff)/sigma) * (xmean-xold)
      

#ifdef __HAVE_MPI__
      !-------------------------------------------------------------------------
      !  COMMUNICATE GLOBAL BEST
      !-------------------------------------------------------------------------
      IF(fitness%hist(1) .LT. out%solutions%bestever%f) THEN
        bestever%x = arxvalid(:,int(fitness%idx(1)))
        bestever%f = fitness%hist(1)
      END IF
      fbest(1) = bestever%f
      fbest(2) = new_rank
      ! Find GLOBAL_F_BEST of running processes
      CALL mpi_rank(fbest,F_BEST_RANK,comm_worker,comm_size)
      !CALL MPI_ALLREDUCE(fbest,GLOBAL_F_BEST,1,MPI_2DOUBLE_PRECISION,&
      !										 MPI_MINLOC,comm_worker,ioError)
      IF(countIter .EQ. 1) THEN	  	
        GLOBAL_F_BEST_LAST = GLOBAL_F_BEST
        CALL MPI_BCAST(GLOBAL_X_BEST,input%N,MPI_DOUBLE_PRECISION,&
                                INT(GLOBAL_F_BEST(2)),comm_worker,ioError)
      ELSE	
        ! Broadcast GLOBAL_X_BEST, if GLOBAL_F_BEST has changed
        
        !! Possible TODO, not sure !!
        IF((GLOBAL_F_BEST(1) .NE. GLOBAL_F_BEST_LAST(1))) THEN
            GLOBAL_X_BEST = bestever%x
            CALL MPI_BCAST(GLOBAL_X_BEST,input%N,MPI_DOUBLE_PRECISION,&
                                INT(GLOBAL_F_BEST(2)),comm_worker,ioError)
            GLOBAL_F_BEST_LAST = GLOBAL_F_BEST
        END IF								
      END IF	 
      
      !-------------------------------------------------------------------------
      !  PARTICLE SWARM UPDATE
      !-------------------------------------------------------------------------    	
      IF(options%pscma) THEN  
        IF((countIter .GT. 1).AND.(mod(countIter,options%psoFreq) .EQ. 0)) THEN
            CALL psoUpdate(input%N,triuC,bestever%x)
        END IF
      END IF	  
#endif	  
      !-------------------------------------------------------------------------
      !  ADAPT COVARIANCE MATRIX
      !-------------------------------------------------------------------------
      DO i = 1, input%N
        colVec(i,1) = pc(i)
        rowVec(1,i) = pc(i)
      END DO
      weightsMat = 0.
      DO i = 1, mu
          tempMat(:,i) = xold(:)
          weightsMat(i,i) = weights(i)
      END DO
        
      ! The following lines implement Eq. 43. BLAS routines are used for matrix
      ! multiplication. Unfortunately this decreases readability. An 
      ! implementation with 'matmul' would look like:
      !	C = (1._MK-ccov+(1._MK-hsig)*ccov*cc*(2._MK-cc)/mucov)*C & ! old matrix
      !		+ ccov*(1._MK/mucov)*matmul(colVec,rowVec) &  ! plus rank one update
      !		+ ccov*(1._MK-1._MK/mucov) &				  ! plus rank mu update
      !		*sigma**(-2) * &
      !		matmul((arx(:,int(fitness%idxsel(1:mu))) - tempMat ), &
      !		matmul(weightsMat, &
      !			transpose(arx(:,int(fitness%idxsel(1:mu))) - tempMat)))
       
      IF(ccov .GT. 0.) THEN										! Eq.43
        mutmp = (arx(:,int(fitness%idxsel(1:mu))) - tempMat)		
#ifdef __SP
        CALL SGEMM('N','T',input%N,input%N,1,1.0,colVec,input%N,&
                                        rowVec,input%N,0.0,multmp_mN,input%N)	 
        
        CALL SGEMM('N','T',mu,input%N,mu,1.0,weightsMat,mu,&
                                            mutmp,input%N,0.0,multmp_mMU,mu)	
            
        CALL SGEMM('N','N',input%N,input%N,mu,1.0,mutmp,input%N,&
                                        multmp_mMU,mu,0.0,multmp_mN2,input%N)	
#else	  
        CALL dgemm('N','T',input%N,input%N,1,1.0d0,colVec,input%N,&
                                        rowVec,input%N,0.0d0,multmp_mN,input%N)	 
                                                
        CALL dgemm('N','T',mu,input%N,mu,1.0d0,weightsMat,mu,&
                                            mutmp,input%N,0.0d0,multmp_mMU,mu)
        
        CALL dgemm('N','N',input%N,input%N,mu,1.0d0,mutmp,input%N,&
                                        multmp_mMU,mu,0.0d0,multmp_mN2,input%N)									
#endif
    
        C = (1._MK-ccov+(1._MK-hsig)*ccov*cc*(2._MK-cc)/mucov)*C & ! old matrix
            + ccov*(1._MK/mucov)*multmp_mN &  ! plus rank one update
            + ccov*(1._MK-1._MK/mucov) &				  ! plus rank mu update
            *sigma**(-2) * multmp_mN2  		
!
      END IF	  
      
      

#ifdef __HAVE_MPI__	  
      !-------------------------------------------------------------------------
      ! Couple CMA and PSO
      !-------------------------------------------------------------------------
      IF(options%pscma) THEN
          IF((countIter .GT. 1).AND.(mod(countIter,options%psoFreq) .EQ. 0) &
            .AND. (options%psoWeight .LT. 1.0_MK)) THEN
            C = options%psoWeight*C + (1._MK - options%psoWeight)*C_PSO
          END IF
      END IF
#endif	  
      !-------------------------------------------------------------------------
      !  Remove momentum
      !-------------------------------------------------------------------------
      IF (sum(ps**2)/real(input%N)  .GT. &
          1.5_MK + 10.0_MK*(2.0_MK/real(input%N)) .AND. &
          (fitness%histsel(1) .GT. &
          max(fitness%histsel(2),fitness%histsel(3)))) THEN

      ps = ps * &
           sqrt(real(input%N)*(1+max(0.0_MK,log(sum(ps**2)/real(input%N)))) &
           /sum(ps**2))

	  WRITE(*,*) 'Momentum in ps removed at [',countIter,',',countEval,']'

	  END IF

      !-------------------------------------------------------------------------
      !  ADAPT SIGMA
      !-------------------------------------------------------------------------	  
      sigma =  sigma*exp((sqrt(sum(ps**2))/chiN-1._MK)*cs/damps)		! Eq.41

      !-------------------------------------------------------------------------
      !  Update B and D from C
      !-------------------------------------------------------------------------
      IF((ccov .GT. 0.) .AND. &
            (mod(real((countIter-lastRestart)), real(1._MK/ccov/input%N/10._MK)) .LT. 1._MK)) THEN
        !-----------------------------------------------------------------------
        !  Enforce symmetry
        !-----------------------------------------------------------------------
        CALL tool_symmatrix(C,input%N,triuC)	  	
        
        
        
        !-----------------------------------------------------------------------
        !  Eigen decomposition, D=diagonal matrix of eigenvalues, 
        !  B=normalized eigenvectors
        !-----------------------------------------------------------------------
        CALL tool_eigendecomp(triuC,input%N,D,B,info)
      
        !if tool_eigendecomp returned Error go to the Restart Loop
        IF (info .NE. 0) THEN
            stopflag = 'DeComp Error'
            !STOP
            GOTO 100
        END IF
        !-----------------------------------------------------------------------
        !  Limit condition of C to 1e14 + 1
        !-----------------------------------------------------------------------	  	
        mask = .FALSE.
        DO i = 1, input%N
            mask(i,i) = .TRUE.		! only work on main diagonal
        END DO	  	
        
        IF(minval(D,mask) .LE. 0.) THEN		! If any eigenvalue < 0
            IF(options%StopOnWarnings) THEN
                stopflag = 'warnconditioncov'
            ELSE
                WRITE(*,*) 'Warning: Iteration ', countIter, &
                            ': Eigenvalue (smaller) zero'
                WHERE(D .LT. 0.)
                    D = 0.
                END WHERE
                
                temp = maxval(D,mask)/1.E14
                tempMat2 = 0.
                DO i = 1, input%N
                    tempMat2(i,i) = 1.
                END DO
                
                C = C + temp*tempMat2
                D = D + temp*tempMat2
            END IF	
        END IF	!minval
        !temp = maxval(D,mask)/1.E14_MK - minval(D,mask)
        !  DO i = 1, input%N
  	    !    WRITE (23, *) temp
  	    !END DO
      
        IF(maxval(D,mask) .GT. 1.E14*minval(D,mask)) THEN	! If at upper limit
            IF(options%StopOnWarnings) THEN
                stopflag = 'warnconditioncov'
            ELSE
                WRITE(*,*) 'Warning: Iteration ', countIter, &
                            ': Condition of C at upper limit'
                temp = maxval(D,mask)/1.E14_MK - minval(D,mask)
                C = C + temp*tempMat2
                D = D + temp*tempMat2
            END IF
        END IF
       !------------------------------------------------------------------------
       ! compute inverse of C for bfgs
       !------------------------------------------------------------------------
      IF (options%BFGS_use ) THEN
       invsC = 0._MK
        WHERE(mask)
            invsC = 1_MK/(sigma**2*D)
        END WHERE
      !abusing BD her as temp. matrix since it will be overriden anyway in next
      !block
      ! B * invD*
#ifdef __SP
        CALL SGEMM('N','N',input%N,input%N,input%N,1.0,B,input%N,invsC,input%N,&
                                                        0.0,BD,input%N)		
#else
        CALL DGEMM('N','N',input%N,input%N,input%N,1.0d0,B,input%N,invsC,input%N,&
                                                        0.0d0,BD,input%N)
#endif
      ! (B * invD*) * B^t
#ifdef __SP
        CALL SGEMM('N','N',input%N,input%N,input%N,1.0,BD,input%N,&
        transpose(B),input%N,0.0,invsC,input%N)		
#else
        CALL DGEMM('N','N',input%N,input%N,input%N,1.0d0,BD,input%N,&
        transpose(B),input%N,0.0d0,invsC,input%N)
#endif
      END IF
        
        !-----------------------------------------------------------------------
        !  D contains standard deviations now
        !-----------------------------------------------------------------------
        WHERE(mask)
            D = sqrt(D)
        END WHERE
        

        
#ifdef __SP
        CALL SGEMM('N','N',input%N,input%N,input%N,1.0,B,input%N,D,input%N,&
                                                        0.0,BD,input%N)		
#else
        CALL DGEMM('N','N',input%N,input%N,input%N,1.0d0,B,input%N,D,input%N,&
                                                        0.0d0,BD,input%N)
#endif
        !BD = matmul(B,D)

      END IF	!Update B,D

!IF(comm_size .GE. 3) THEN	 
!IF(mod(countIter,50) .EQ. 0) THEN
!	CALL mpi_shift_worst(input%N)
!END IF	
!END IF	 
!STOP 
!IF(MY_RANK .EQ. 0) THEN
!IF(mod(countIter,2) .EQ. 0) THEN
!	xmean = 0.0_MK
!	CALL cmaes_reinitialize(xmean,input%N)	  
!END IF
!END IF  	  
      
      !-------------------------------------------------------------------------
      !  NUMERICAL ERROR MANAGEMENT
      !-------------------------------------------------------------------------	  
      DO i = 1, input%N
        diag(i) = C(i,i)		! Get diagonal vector of C
      END DO
      
      !-------------------------------------------------------------------------
      !  Adjust maximal coordinate axis deviations
      !-------------------------------------------------------------------------
      IF(any(sigma*sqrt(diag) .GT. maxdx)) sigma = minval(maxdx/sqrt(diag))
      
      !-------------------------------------------------------------------------
      !  Adjust minimal coordinate axis deviations
      !-------------------------------------------------------------------------
      IF(any(sigma*sqrt(diag) .LT. mindx)) &
        sigma = maxval(mindx/sqrt(diag)) * exp(0.05+cs/damps)
      
      !-------------------------------------------------------------------------
      !  Adjust too low coordinate axis deviations
      !-------------------------------------------------------------------------
      mask = .FALSE.
      tempMat2 = 0._MK
      DO i = 1, input%N
        IF(xmean(i) .EQ. xmean(i) + 0.2_MK*sigma*sqrt(diag(i))) &
                                                        mask(i,i) = .TRUE.
        tempMat2(i,i) = diag(i)
      END DO
      
      IF(any(mask)) THEN
        IF(options%stopOnWarnings) THEN
            stopflag = 'warnnoeffectcoord'
        ELSE
            WRITE(*,*) 'Warning: Iteration ', countIter, &
                            ': Coordinate axis std deviation too low'
            WHERE(.NOT. mask)
                tempMat2 = 0._MK
            END WHERE
            C = C + ccov*tempMat2
            sigma = sigma*exp(0.05_MK+cs/damps)
        END IF
      END IF
      
      !-------------------------------------------------------------------------
      !  Adjust step size in case of (numerical) precision problem
      !-------------------------------------------------------------------------
      IF(all(xmean .EQ. xmean + &
        0.1_MK*sigma*BD(:,1+int(mod((countIter-lastRestart),input%N))))) THEN
        i = 1+int(mod((countIter-lastRestart),input%N))
        IF(options%stopOnWarnings) THEN
            stopflag = 'warnnoeffectaxis'
        ELSE
            WRITE(*,*) 'Warning: Iteration ', countIter, &
                            ': Main axis std deviation ', sigma*D(i,i), &
                            ' has no effect'
            sigma = sigma*exp(0.2_MK+cs/damps)
        END IF
            
      END IF
      
      !-------------------------------------------------------------------------
      !  Adjust step size in case of equal function values (flat fitness)
      !-------------------------------------------------------------------------
      temp = int(0.1_MK+lambda/4._MK)				! Little hack to round up
      IF(real(temp) .NE. (0.1_MK+lambda/4._MK)) &
        temp = temp + 1._MK
      IF(fitness%sel(1) .EQ. fitness%sel(1+int(temp))) THEN
        IF(options%WarnOnEqualFunctionValues .AND. options%stopOnWarnings) THEN
            stopflag = 'warnequalfunvals'
        ELSE
            IF(options%WarnOnEqualFunctionValues) THEN
            DO i = 1, input%N
                diag(i) = D(i,i)
            END DO
            WRITE(*,*) &
                'Warning: Iteration ', countIter, ': Equal function values f=',&
                fitness%sel(1), ' at maximal main axis sigma ', &
                sigma*maxval(diag)
            END IF
            sigma = sigma*exp(0.2_MK+cs/damps)
        END IF
      END IF
      
      !-------------------------------------------------------------------------
      !  Adjust step size in case of equal function values
      !-------------------------------------------------------------------------
      temp = tool_myrange(fitness%hist,size(fitness%hist),fitness%sel(1),1,posInf)
      IF((countIter - lastRestart .GT. 2) .AND. (temp .EQ. 0.)) THEN
        IF(options%stopOnWarnings) THEN
            stopflag = 'warnequalfunvalhist'
        ELSE
            DO i = 1, input%N
                diag(i) = D(i,i)
            END DO
            WRITE(*,*) &
                'Warning: Iteration ', countIter, ': Equal function values', &
                'at maximal main axis sigma ', sigma*maxval(diag)	  		
            sigma = sigma*exp(0.2_MK+cs/damps)
        END IF
      END IF
      
      
      !-------------------------------------------------------------------------
      !  KEEP OVERALL BEST SOLUTION
      !-------------------------------------------------------------------------
      out%evals = countEval
      out%solutions%evals = countEval
      out%solutions%mean%x = xmean
      out%solutions%mean%f = fmean
      out%solutions%mean%evals = countEval
      out%solutions%recentbest%x = arxvalid(:,int(fitness%idx(1)))
      out%solutions%recentbest%f = fitness%raw(1)
      out%solutions%recentbest%evals = countEval + fitness%idx(1) - lambda
      out%solutions%recentworst%x = &
        arxvalid(:,int(fitness%idx(size(fitness%idx))))
      out%solutions%recentworst%f = fitness%raw(size(fitness%raw))
      out%solutions%recentworst%evals = countEval &
        + fitness%idx(size(fitness%idx)) - lambda
      IF(fitness%hist(1) .LT. out%solutions%bestever%f) THEN
        out%solutions%bestever%x = arxvalid(:,int(fitness%idx(1)))
        out%solutions%bestever%f = fitness%hist(1)
        out%solutions%bestever%evals = countEval + fitness%idx(1) - lambda
        bestever = out%solutions%bestever
      END IF
      
      !-------------------------------------------------------------------------
      !  Set stop flag
      !-------------------------------------------------------------------------
      IF(sigma*maxval(diag) .EQ. 0.) stopflag = 'bug'	!should never happen ;)
      
      
      IF(countIter .GE. options%StopMaxIter) stopflag = 'maxiter'
      
      DO i = 1, input%N
        diag(i) = C(i,i)
      END DO

      IF(all(sigma*(max(abs(pc),sqrt(diag))) .LT. options%StopTolX)) &
        stopflag = 'tolx'
      
      IF(any(sigma*sqrt(diag) .GT. options%StopTolUpX)) stopflag = 'tolupx'
      
      IF((countIter-lastRestart .GT. 2) .AND. &
        (tool_myrange(fitness%sel,size(fitness%sel),fitness%hist,&
            size(fitness%hist),posInf) .LE. options%StopTolFun)) &
                stopflag = 'tolfun'

      IF((countIter-lastRestart .GE. size(fitness%hist)) .AND. &
        (tool_myrange(fitness%hist,size(fitness%hist),posInf_array,1,posInf) &
            .LE. options%StopTolHistFun)) THEN 
            stopflag = 'tolhistfun'
      END IF
      
      IF((countEval .GE. options%StopFunEvals) .OR. &
            (countIter .GE.options%StopIter)) stopflag = 'stopFunEvals'			
      
      IF (options%StopTime ) THEN
        CALL CPU_TIME(time_end)
        tmpreal = time_end - time_start
        IF (tmpreal .GT. (options%StopTimeHH*60+options%StopTimeMM)*60+ &
     &                    options%StopTimeSS) THEN
            stopflag = 'stopTime'
        END IF
      END IF

#ifdef __HAVE_MPI__
      CALL MPI_ALLREDUCE(countEval,GLOBAL_FUN_EVALS,1,MPI_INTEGER,&
                                             MPI_SUM,comm_worker,ioError)
      IF(GLOBAL_FUN_EVALS .GE. options%StopMaxFunEvals) stopflag = 'maxfunevals'
      IF(GLOBAL_F_BEST(1) .LE. options%StopFitness) stopflag = 'fitness'
#else
      IF(countEval .GE. options%StopMaxFunEvals) stopflag = 'maxfunevals'
      IF(fitness%raw(1) .LE. options%StopFitness) stopflag = 'fitness'
#endif
      !-------------------------------------------------------------------------
      !  Output Generation
      !-------------------------------------------------------------------------
      IF((options%flgGenData .OR.  options%flgGenTrace) .AND. &
     &(mod(countIter,options%intGenData)) .EQ. 0) THEN
#ifdef __HAVE_MPI__      
        CALL cmaes_writegen(arx,arxvalid,GLOBAL_X_BEST,input%N,lambda,uNum)
#else
        CALL cmaes_writegen(arx,arxvalid,input%N,lambda,uNum)
        
#endif 
      END IF
      
      ! Minimal (valid) candiate solution and function value in current population
      xmin = arxvalid(:,int(fitness%idx(1)))
      fmin = fitness%raw(1)
      
      IF(options%VerboseModulo > 0) THEN
        IF(countIter .EQ. 1) THEN
#ifdef __HAVE_MPI__
        IF(MY_RANK .EQ. 0 .OR. (.NOT. options%silentMPI )) THEN
            WRITE(*,*) '******************************************************'
            WRITE(*,*) 'Started MPI-CMA' 
            IF (options%silentMPI ) THEN
                WRITE(*,*) 'To guarantee a decent console output',&
     &              ' only Process 0 is shown'
            END IF
            WRITE(*,*) ' All output data is saved to ',&
                'folder ' // options%output_folder
            WRITE(*,*) '******************************************************'	
      IF(options%pscma) THEN
          IF(MY_RANK .EQ. 0 .OR. (.NOT. options%silentMPI)) THEN
            WRITE(*,*) '---------------------------------'
            WRITE(*,*) 'PSO configuration:'
            WRITE(*,*) 'Weight: ', options%psoWeight
            WRITE(*,*) 'Frequency:', options%psoFreq
            WRITE(*,*) '---------------------------------'
            WRITE(*,*) 'Initial sigma', options%insigma(1)
            WRITE(*,*) '---------------------------------'
          END IF 	
      END IF	
#endif	  	
            ! Do some formatting for a nice console output	
            WRITE(tmp_char,'(I10)') input%N
            tmp_char2 = 'n = '  //  trim(adjustl(tmp_char))  //  ': ('
            WRITE(tmp_char,'(I10)') mu
            tmp_char2 = trim(adjustl(tmp_char2)) // ' ' // trim(adjustl(tmp_char))&
                 // ' ,'
            WRITE(tmp_char,'(I10)') lambda
            tmp_char2 = trim(adjustl(tmp_char2)) //  ' ' // trim(adjustl(tmp_char))
            WRITE(*,*) trim(adjustl(tmp_char2)),' )-CMA-ES on function '&
     &       // trim(adjustl(options%funcName))
            WRITE(*,*) '	Iterat,     #Fevals:	Function Value'
#ifdef __HAVE_MPI__
        END IF
#endif			  		
        END IF
        
#ifdef __HAVE_MPI__
        IF(MY_RANK .EQ. 0 .OR. .NOT. options%silentMPI) THEN
#endif			  	
        IF((mod((countIter-lastRestart), options%VerboseModulo) .EQ. 0) &
            .OR. (stopflag .NE. '') .OR. (countIter .EQ. 1)) THEN
            WRITE(*,*) countIter, ',', countEval, ':', fmin
        IF ((options%use_LJ .OR. options%use_LJ_COMP)  .AND. options%write_pdb ) THEN
                   WRITE(tmp_char,'(I10)') countIter
                   tmp_dir=trim(adjustl(options%output_folder))&
     & // '/LJ' &
     & // trim(adjustl(proc))&
     & // '.pdb'
      tmp_dir = trim(adjustl(tmp_dir))
                  CALL LJ_write_out(tmp_dir,xmin,input%N,1,(fmin-options%StopFitness))
                  END IF
            
      IF (options%use_TIP  .AND. options%write_pdb ) THEN
                   WRITE(tmp_char,'(I10)') countIter
                   tmp_dir=trim(adjustl(options%output_folder))&
     & // '/WATER' &
     & // trim(adjustl(proc))&
     & // '.pdb'
      tmp_dir = trim(adjustl(tmp_dir))
      CALL water_write_out(tmp_dir,xmin,input%N,1)
      END IF
                  
                  
        END IF
#ifdef __HAVE_MPI__
        END IF
#endif			  	
      END IF
      
      
      !-------------------------------------------------------------------------
      !  Reset fitness indices
      !-------------------------------------------------------------------------
      DO i = 1, lambda
        fitness%idx(i) = real(i)
        fitness%idxsel(i) = real(i)
      END DO
      
      
#ifdef __HAVE_MPI__	 
      !-------------------------------------------------------------------------
      !  ADAPT COMMUNICATOR 
      !	 (exclude processes from communications,that are about to stop)
      !------------------------------------------------------------------------- 
      ! Check if any process wants to stop (collect data in stop_proc_mask)
      IF(stopflag .NE. '' .AND. .NOT.(options%restart_cma)) STOP_FLAG = .TRUE.
	  CALL mpi_comm_adapt()
      IF(stop_me) THEN        
        EXIT genLoop
      ENDIF
#endif

      END DO genLoop
      
100   IF (options%restart_cma  .AND. stopflag .NE. 'bug' .AND. &
     &stopflag .NE. 'fitness' &
     & .AND. stopflag .NE. 'maxfunevals' .AND. stopflag .NE. 'maxiter' .AND. &
     & stopflag .NE. 'stopTime') THEN
#ifdef __HAVE_MPI__
        IF(MY_RANK .EQ. 0 .OR. .NOT. options%silentMPI) THEN
#endif
      WRITE(*,*) countIter, ',', countEval, ':', fitness%hist(1)
#ifdef __HAVE_MPI__
        END IF
#endif
      
      IF (options%use_LJ  .AND. options%write_pdb ) THEN
                   WRITE(tmp_char,'(I10)') countIter
                   tmp_dir=trim(adjustl(options%output_folder))&
     & // '/LJ' &
     & // trim(adjustl(proc))&
     & // '.pdb'
      
      
      CALL LJ_write_out(tmp_dir,xmin,input%N,1,(fmin-options%StopFitness))
      END IF
      
      IF (options%use_TIP  .AND. options%write_pdb ) THEN
                   WRITE(tmp_char,'(I10)') countIter
                   tmp_dir=trim(adjustl(options%output_folder))&
     & // '/WATER' &
     & // trim(adjustl(proc))&
     & // '.pdb'
      CALL water_write_out(tmp_dir,xmin,input%N,1)
      END IF
      
      
     
      countRestarts = countRestarts + 1
#ifdef __HAVE_MPI__
        IF(MY_RANK .EQ. 0 .OR. .NOT. options%silentMPI) THEN
#endif
      WRITE(*,*) '------- Restart #', countRestarts, ' Reason: ', stopflag
#ifdef __HAVE_MPI__
        END IF
#endif
      stopflag = ''
      CALL cmaes_startpoint()

#ifdef __HAVE_MPI__
        IF(MY_RANK .EQ. 0 .OR. .NOT. options%silentMPI) THEN
#endif
        ! Do some formatting for a nice console output	
        WRITE(tmp_char,'(I10)') input%N
        tmp_char2 = 'n = '  //  trim(adjustl(tmp_char))  //  ': ('
        WRITE(tmp_char,'(I10)') mu
        tmp_char2 = trim(adjustl(tmp_char2)) // ' ' // trim(adjustl(tmp_char))&
             // ' ,'
        WRITE(tmp_char,'(I10)') lambda
        tmp_char2 = trim(adjustl(tmp_char2)) //  ' ' // trim(adjustl(tmp_char))
        WRITE(*,*) trim(adjustl(tmp_char2)),' )-CMA-ES on function '&
        &       // trim(adjustl(options%funcName))
        WRITE(*,*) '---------------------------------'
        WRITE(*,*) 'Initial sigma', options%insigma(1)
        WRITE(*,*) '---------------------------------'
        WRITE(*,*) '	Iterat,     #Fevals:	Function Value'
#ifdef __HAVE_MPI__
        END IF
#endif
      lastRestart = countIter
      
      
      !------------------------------------------------------------------------
      ! Re-allocate local variables
      !------------------------------------------------------------------------
      
      DEALLOCATE(arz,stat=allocStat)
      DEALLOCATE(arx,stat=allocStat)
      DEALLOCATE(arxvalid,stat=allocStat)
      DEALLOCATE(darz,stat=allocStat)
      DEALLOCATE(tempMat,stat=allocStat)
      DEALLOCATE(weightsMat,stat=allocStat)
      DEALLOCATE(newSort,stat=allocStat)
      DEALLOCATE(multmp_mMU,stat=allocStat)
      DEALLOCATE(mutmp,stat=allocStat)
      ALLOCATE(arz(input%N,lambda),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating arz'
      ALLOCATE(arx(input%N,lambda),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating arx'
      ALLOCATE(arxvalid(input%N,lambda),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating arxvalid'
      ALLOCATE(darz(input%N,lambda),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating darz'
      ALLOCATE(tempMat(input%N,mu),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating tempMat'
      ALLOCATE(weightsMat(mu,mu),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating weightsMat'
      ALLOCATE(newSort(lambda),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating newSort'
      ALLOCATE(multmp_mMU(mu,input%N),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating multmp_mMU'
      ALLOCATE(mutmp(input%N,mu),stat=allocStat)
      IF(allocStat .NE. 0) STOP 'Error allocating mutmp'
      
      
      xold = xmean
      arx = 0.0_MK
      arxvalid = 0.0_MK
      xmin = posInf
      fmin = posInf
      uNum = 26
      posInf_array(1) = posInf
      
      END IF
      END DO restartLoop
      

#ifdef __HAVE_MPI__
      !-------------------------------------------------------------------------
      !  ADAPT COMMUNICATOR IN CASE RESTART IS ACTIVATED
      !	 (exclude processes from communications,that are about to stop)
      !-------------------------------------------------------------------------
      ! Check if any process wants to stop (collect data in stop_proc_mask)
      IF(stopflag .NE. '' .AND. options%restart_cma) THEN
      	STOP_FLAG = .TRUE.
      	CALL mpi_comm_adapt()
      END IF
#endif



      !-------------------------------------------------------------------------
      !  Close Generation Data Files
      !-------------------------------------------------------------------------
      IF(options%flgGenData) THEN
        DO i = 0, 21
            CLOSE(unit=uNum+i)
        END DO
      END IF
      
      
      CLOSE(99,IOSTAT=info)
       CLOSE (23)   
      
      
      !-------------------------------------------------------------------------
      !  Display Result
      !-------------------------------------------------------------------------
#ifdef __HAVE_MPI__
      
      IF(MY_RANK .EQ. 0) THEN
#endif	        
      CALL tool_formatarrays(outFormat,1)
      WRITE(*,*)
      WRITE(*,*) 'Process Bestever.f:',bestever%f !GLOBAL_F_BEST(1)
      WRITE(*,*) 'Process Bestever.x:'
      WRITE(*,outFormat) bestever%x !GLOBAL_X_BEST
      WRITE(*,*) 'Stopflag:',stopflag
      IF (options%BFGS_use ) THEN
	  WRITE(*,*) '------------------BFGS-------------------' 
	  WRITE(*,*) 'BFGS Evals: ', countBFGSEval
	  WRITE(*,*) 'Exit', BFGSexitcounter(1),'Grad exit',BFGSexitcounter(2)  
	  WRITE(*,*) 'Min. Step exit',BFGSexitcounter(4), 'Confprop exit',&
	  BFGSexitcounter(5)
	  WRITE(*,*) 'Error Exit', BFGSexitcounter(3)
	  WRITE(*,*) 'Steps Exit', BFGSexitcounter(6)
	  END IF
      

#ifdef __HAVE_MPI__
      END IF
#endif	        

      RETURN
      END SUBROUTINE cmaes_run
