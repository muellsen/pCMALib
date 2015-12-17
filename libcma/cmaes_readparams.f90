      !-------------------------------------------------------------------------
      ! File       :            cmaes_readparams
      !-------------------------------------------------------------------------
      !
      !  Purpose      : reads all parameters nessecary for CMA optimizer from
      !                 a text file and by that avoids recompling every time
      !
      !  Remarks      : modified from PPM library
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
      SUBROUTINE cmaes_readparams(ctrlfile, info)
      USE cmaes_opts_mod
      USE CEC2005
      IMPLICIT NONE
#ifdef __HAVE_MPI__      
      INCLUDE "mpif.h"
#endif

      INTEGER, INTENT(inout)                :: info
      CHARACTER(len=256), INTENT(in)        :: ctrlfile
      CHARACTER(len=10)		:: proc, tmp
	  ! Random number generator
      REAL(MK)				:: ZBQLUAB

      !-----------------------------------------------------------------------------
      ! Declaration of local variables
      !-----------------------------------------------------------------------------
      
      INTEGER                               :: i,j
      INTEGER                               :: idx,i1,i2
      CHARACTER(LEN=256)                    :: cbuf
      INTEGER                               :: iUnit
      INTEGER                               :: ILEN,ios
      INTEGER                               :: ilenctrl
      INTEGER                               :: ibc,iline
      CHARACTER(LEN=256)                    :: cvalue,carg
      CHARACTER(LEN=256)                    :: bcloc,tvalue
      LOGICAL                               :: lExist
      INTEGER                               :: ierr
      INTEGER                               :: MPTYPE

      
      
#ifdef __HAVE_MPI__      
	  WRITE(proc,'(I10)') MY_RANK
	  tmp = '_' //  trim(adjustl(proc))
	  proc = tmp
#else
      proc = ''
#endif

#ifdef __HAVE_MPI__
      IF (MY_RANK .EQ. 0) THEN  !only rank 0 reads the file
#endif      
      !-----------------------------------------------------------------------------
      ! Definition of file unit
      !-----------------------------------------------------------------------------
      
      iUnit = 23
      info  = 0
      !-----------------------------------------------------------------------------
      ! Check that the parameter file exists
      !-----------------------------------------------------------------------------
      
      ilenctrl = LEN_TRIM(ctrlfile)         ! length of filename string
      INQUIRE(FILE=ctrlfile(1:ilenctrl), EXIST=lExist)
      IF(.NOT. lExist) THEN
         WRITE(*,'(2A)')'No such file: ',ctrlfile(1:ilenctrl)
         Info = 1
         STOP
      END IF
      
      !-----------------------------------------------------------------------------
      ! open the file
      !-----------------------------------------------------------------------------
      
      OPEN(iUnit, FILE=ctrlfile(1:ilenctrl), IOSTAT=ios, ACTION='READ')
      IF(ios .NE. 0) THEN
         WRITE(*,'(2A)')'Failed to open file: ',ctrlfile(1:ilenctrl)
         Info = 1
         STOP
      END IF
      
      !-----------------------------------------------------------------------------
      ! scan file
      !-----------------------------------------------------------------------------
  
      iline = 0
      DO 
         iline = iline + 1        ! increment line 
         READ(iUnit,'(A)',END=100,ERR=200) cbuf 
         ILEN = LEN_TRIM(cbuf)
         
         !--------------------------------------------------------------------------
         !  Skip comment or empty lines 
         !--------------------------------------------------------------------------
         IF(ILEN .GT. 0 .AND. cbuf(1:1) .NE. '#') THEN
            !-----------------------------------------------------------------------
            !  Remove all spaces
            !-----------------------------------------------------------------------
            j = 0
            DO i=1,ILEN
               IF(cbuf(i:i) .NE. ' ') THEN
                  j = j + 1
                  cbuf(j:j) = cbuf(i:i)
               END IF
            END DO
            ILEN = j     ! update length of string
            
            !-----------------------------------------------------------------------
            !  Find position of =
            !-----------------------------------------------------------------------
            idx = INDEX(cbuf,'=')
            
            !-----------------------------------------------------------------------
            !  Exit if = is missing
            !-----------------------------------------------------------------------
            IF(idx .LT. 0) THEN
               WRITE(*,'(A,I5)')'Incorrect line: ',iline
               Info = 1
               STOP
            END IF
            
            !-----------------------------------------------------------------------
            !  Get argument and value
            !-----------------------------------------------------------------------
            carg   = ADJUSTL(cbuf(1:idx-1))
            cvalue = ADJUSTL(cbuf(idx+1:ILEN))
            
            !-----------------------------------------------------------------------
            !  Convert to upper case
            !-----------------------------------------------------------------------
            !CALL UpperCase(carg,idx-1)

            !-----------------------------------------------------------------------
            !  Parse and read input data
            !-----------------------------------------------------------------------
            IF(carg(1:(idx-1)) .EQ. 'REL_SIGMA') THEN
               READ(cvalue,*,iostat=ios,err=200) options%rel_sigma
            ELSEIF(carg(1:(idx-1)) .EQ. 'ABS_SIGMA') THEN
               READ(cvalue,*,iostat=ios,err=200) options%abs_sigma
!            ELSEIF              
!               iruntag = LEN_TRIM(runtag)
            ELSEIF (carg(1:(idx-1)).EQ.'ALLDIM_LBOUNDS') THEN
               READ(cvalue,*,iostat=ios,err=200) options%alldim_lBounds
            ELSEIF (carg(1:(idx-1)).EQ.'ALLDIM_UBOUNDS') THEN
               READ(cvalue,*,iostat=ios,err=200) options%alldim_uBounds
            ELSEIF (carg(1:(idx-1)).EQ.'DIMENSIONS') THEN
               READ(cvalue,*,iostat=ios,err=200) options%dimensions
            ELSEIF (carg(1:(idx-1)).EQ.'POPSIZE') THEN
               READ(cvalue,*,iostat=ios,err=200) options%PopSize
            ELSEIF (carg(1:(idx-1)).EQ.'USE_CEC') THEN
               READ(cvalue,*,iostat=ios,err=200) options%use_CEC
            ELSEIF (carg(1:(idx-1)).EQ.'USE_LJ') THEN
               READ(cvalue,*,iostat=ios,err=200) options%use_LJ
            ELSEIF (carg(1:(idx-1)).EQ.'USE_BBOB') THEN
               READ(cvalue,*,iostat=ios,err=200) options%use_BBOB
            ELSEIF (carg(1:(idx-1)).EQ.'USE_TIP') THEN
               READ(cvalue,*,iostat=ios,err=200) options%use_TIP
            ELSEIF (carg(1:(idx-1)).EQ.'BENCHFCTNR') THEN
               READ(cvalue,*,iostat=ios,err=200) options%Benchfctnr
            ELSEIF (carg(1:(idx-1)).EQ.'OUTPUT_FOLDER') THEN
               READ(cvalue,'(A)',iostat=ios,err=200) options%output_folder
            ELSEIF (carg(1:(idx-1)).EQ.'STOPFITNESS') THEN
               READ(cvalue,*,iostat=ios,err=200) options%StopFitness
            ELSEIF (carg(1:(idx-1)).EQ.'STOPMAXFUNEVALS') THEN
               READ(cvalue,*,iostat=ios,err=200) options%StopMaxFunEvals
            ELSEIF (carg(1:(idx-1)).EQ.'STOPMAXITER') THEN
               READ(cvalue,*,iostat=ios,err=200) options%StopMaxIter
            ELSEIF (carg(1:(idx-1)).EQ.'STOPTOLX') THEN
               READ(cvalue,*,iostat=ios,err=200) options%StopTolX
            ELSEIF (carg(1:(idx-1)).EQ.'STOPTOLUPX') THEN
               READ(cvalue,*,iostat=ios,err=200) options%StopTolUpX
            ELSEIF (carg(1:(idx-1)).EQ.'STOPTOLFUN') THEN
               READ(cvalue,*,iostat=ios,err=200) options%StopTolFun
            ELSEIF (carg(1:(idx-1)).EQ.'STOPTOLHISTFUN') THEN
               READ(cvalue,*,iostat=ios,err=200) options%StopTolHistFun
            ELSEIF (carg(1:(idx-1)).EQ.'STOPONWARNINGS') THEN
               READ(cvalue,*,iostat=ios,err=200) options%StopOnWarnings
        ELSEIF (carg(1:(idx-1)).EQ.'WARNONEQUALFUNCTIONVALUES') THEN
           READ(cvalue,*,iostat=ios,err=200) options%WarnOnEqualFunctionValues
            ELSEIF (carg(1:(idx-1)).EQ.'EVALINITIALX') THEN
               READ(cvalue,*,iostat=ios,err=200) options%EvalInitialX
            ELSEIF (carg(1:(idx-1)).EQ.'PARENTNUMBER') THEN
               READ(cvalue,*,iostat=ios,err=200) options%ParentNumber
            ELSEIF (carg(1:(idx-1)).EQ.'RECOMBINATIONWEIGHTS') THEN
               READ(cvalue,*,iostat=ios,err=200) options%RecombinationWeights
            ELSEIF (carg(1:(idx-1)).EQ.'VERBOSEMODULO') THEN
               READ(cvalue,*,iostat=ios,err=200) options%VerboseModulo
            ELSEIF (carg(1:(idx-1)).EQ.'FLGGENDATA') THEN
               READ(cvalue,*,iostat=ios,err=200) options%flgGenData
            ELSEIF (carg(1:(idx-1)).EQ.'INTGENDATA') THEN
               READ(cvalue,*,iostat=ios,err=200) options%intGenData
            ELSEIF (carg(1:(idx-1)).EQ.'FLGOUTTXT') THEN
               READ(cvalue,*,iostat=ios,err=200) options%flgouttxt
            ELSEIF (carg(1:(idx-1)).EQ.'FLGGENTRACE') THEN
               READ(cvalue,*,iostat=ios,err=200) options%flgGenTrace
            ELSEIF (carg(1:(idx-1)).EQ.'FUNCNAME') THEN
               READ(cvalue,*,iostat=ios,err=200) options%funcName
            ELSEIF (carg(1:(idx-1)).EQ.'PSCMA') THEN
               READ(cvalue,*,iostat=ios,err=200) options%pscma
            ELSEIF (carg(1:(idx-1)).EQ.'PSOWEIGHT') THEN
               READ(cvalue,*,iostat=ios,err=200) options%psoWeight
            ELSEIF (carg(1:(idx-1)).EQ.'PSOFREQ') THEN
               READ(cvalue,*,iostat=ios,err=200) options%psoFreq
            ELSEIF (carg(1:(idx-1)).EQ.'ACCURACY') THEN
               READ(cvalue,*,iostat=ios,err=200) options%accuracy
            ELSEIF (carg(1:(idx-1)).EQ.'GLOBAL_MIN') THEN
               READ(cvalue,*,iostat=ios,err=200) options%global_min
            ELSEIF (carg(1:(idx-1)).EQ.'USE_SEED') THEN
               READ(cvalue,*,iostat=ios,err=200) options%use_seed
            ELSEIF (carg(1:(idx-1)).EQ.'SILENTMPI') THEN
               READ(cvalue,*,iostat=ios,err=200) options%silentMPI
            ELSEIF (carg(1:(idx-1)).EQ.'QR_SAMPLING') THEN
               READ(cvalue,*,iostat=ios,err=200) options%qr_sampling
            ELSEIF (carg(1:(idx-1)).EQ.'QR_SAMPLER') THEN
               READ(cvalue,*,iostat=ios,err=200) options%qr_sampler
            ELSEIF (carg(1:(idx-1)).EQ.'QR_SCRAMBLING') THEN
               READ(cvalue,*,iostat=ios,err=200) options%qr_scrambling
            ELSEIF (carg(1:(idx-1)).EQ.'QR_INVERTER') THEN
               READ(cvalue,*,iostat=ios,err=200) options%qr_inverter
            ELSEIF (carg(1:(idx-1)).EQ.'USE_INIT_BOUNDS') THEN
               READ(cvalue,*,iostat=ios,err=200) options%use_init_bounds
            ELSEIF (carg(1:(idx-1)).EQ.'INIT_LBOUNDS') THEN
               READ(cvalue,*,iostat=ios,err=200) options%init_lBounds
            ELSEIF (carg(1:(idx-1)).EQ.'INIT_UBOUNDS') THEN
               READ(cvalue,*,iostat=ios,err=200) options%init_uBounds
            ELSEIF (carg(1:(idx-1)).EQ.'SEED_FOLDER') THEN
               READ(cvalue,'(A)',iostat=ios,err=200) options%seed_folder
            ELSEIF (carg(1:(idx-1)).EQ.'BENCHMARK') THEN
               READ(cvalue,*,iostat=ios,err=200) options%benchmark
            ELSEIF (carg(1:(idx-1)).EQ.'RECORD_ACCURACY') THEN
               READ(cvalue,*,iostat=ios,err=200) options%record_accuracy
            ELSEIF (carg(1:(idx-1)).EQ.'RESTART_CMA') THEN
               READ(cvalue,*,iostat=ios,err=200) options%restart_cma
            ELSEIF (carg(1:(idx-1)).EQ.'RESTART_TYPE') THEN
               READ(cvalue,*,iostat=ios,err=200) options%restart_type
            ELSEIF (carg(1:(idx-1)).EQ.'RESTARTS') THEN
               READ(cvalue,*,iostat=ios,err=200) options%Restarts
            ELSEIF (carg(1:(idx-1)).EQ.'INCPOPSIZE') THEN
               READ(cvalue,*,iostat=ios,err=200) options%IncPopSize
            ELSEIF (carg(1:(idx-1)).EQ.'RECORD_BESTHIST') THEN
               READ(cvalue,*,iostat=ios,err=200) options%record_besthist
            ELSEIF (carg(1:(idx-1)).EQ.'RECORD_MODULO') THEN
               READ(cvalue,*,iostat=ios,err=200) options%record_modulo
            ELSEIF (carg(1:(idx-1)).EQ.'MAXINCFAC') THEN
               READ(cvalue,*,iostat=ios,err=200) options%MaxIncFac
            ELSEIF (carg(1:(idx-1)).EQ.'BFGS_USE') THEN
               READ(cvalue,*,iostat=ios,err=200) options%BFGS_use
            ELSEIF (carg(1:(idx-1)).EQ.'BFGS_FACTR') THEN
               READ(cvalue,*,iostat=ios,err=200) options%BFGS_factr
            ELSEIF (carg(1:(idx-1)).EQ.'BFGS_PGTOL') THEN
               READ(cvalue,*,iostat=ios,err=200) options%BFGS_pgtol            
            ELSEIF (carg(1:(idx-1)).EQ.'BFGS_POSITION') THEN
               READ(cvalue,*,iostat=ios,err=200) options%BFGS_position
            ELSEIF (carg(1:(idx-1)).EQ.'BFGS_CENTRAL_DIFFERENCE') THEN
               READ(cvalue,*,iostat=ios,err=200) options%BFGS_central_difference            
            ELSEIF (carg(1:(idx-1)).EQ.'BFGS_CONFPROP') THEN
               READ(cvalue,*,iostat=ios,err=200) options%BFGS_confprop
            ELSEIF (carg(1:(idx-1)).EQ.'BFGS_SIGMASTEP_MIN') THEN
               READ(cvalue,*,iostat=ios,err=200) options%BFGS_sigmastep_min
            ELSEIF (carg(1:(idx-1)).EQ.'STOPTIME') THEN
               READ(cvalue,*,iostat=ios,err=200) options%StopTime
            ELSEIF (carg(1:(idx-1)).EQ.'STOPTIMEHH') THEN
               READ(cvalue,*,iostat=ios,err=200) options%StopTimeHH
            ELSEIF (carg(1:(idx-1)).EQ.'STOPTIMEMM') THEN
               READ(cvalue,*,iostat=ios,err=200) options%StopTimeMM
            ELSEIF (carg(1:(idx-1)).EQ.'STOPTIMESS') THEN
               READ(cvalue,*,iostat=ios,err=200) options%StopTimeSS
            ELSEIF (carg(1:(idx-1)).EQ.'CECFOLDERS') THEN
               READ(cvalue,'(A)',iostat=ios,err=200) options%CECFolders
            ELSEIF (carg(1:(idx-1)).EQ.'WRITE_PDB') THEN
               READ(cvalue,*,iostat=ios,err=200) options%write_pdb            
            ELSEIF (carg(1:(idx-1)).EQ.'LJ_COMP') THEN
               READ(cvalue,*,iostat=ios,err=200) options%LJ_comp
            ELSEIF (carg(1:(idx-1)).EQ.'USE_LJ_COMP') THEN
               READ(cvalue,*,iostat=ios,err=200) options%use_LJ_comp
            ELSEIF (carg(1:(idx-1)).EQ.'USE_DF') THEN
               READ(cvalue,*,iostat=ios,err=200) options%use_DF
            ELSEIF (carg(1:(idx-1)).EQ.'USE_MATFUNC') THEN
               READ(cvalue,*,iostat=ios,err=200) options%use_MATFUNC
            ELSEIF (carg(1:(idx-1)).EQ.'DF_S') THEN
               READ(cvalue,*,iostat=ios,err=200) options%DF_s
            ELSEIF (carg(1:(idx-1)).EQ.'DF_RAST') THEN
               READ(cvalue,*,iostat=ios,err=200) options%DF_rast
            ELSEIF (carg(1:(idx-1)).EQ.'USE_RANDOM_LANDSCAPE') THEN
               READ(cvalue,*,iostat=ios,err=200) options%use_RANDOM_LANDSCAPE
            ELSE
               WRITE(*,*) 'ERROR reading input File - unkown Parameter ', carg
               STOP
            END IF
            
                
         END IF
      END DO


      !-----------------------------------------------------------------------------
      !  Something went wrong if we got here
      !-----------------------------------------------------------------------------
  200 CONTINUE
      WRITE(*,'(A,I5,2A)') 'Error reading line: ',iline,&
     &                 ' of file: ',ctrlfile(1:ilenctrl)
      ILEN = LEN_TRIM(cbuf)
      WRITE(*,'(A)') cbuf(1:ILEN)
      Info = 1
      STOP
      
      !-----------------------------------------------------------------------------
      !  End of file
      !-----------------------------------------------------------------------------
  100 Info = 0
      
      !-----------------------------------------------------------------------------
      !  Close file
      !-----------------------------------------------------------------------------
      CLOSE(iUnit)
      
      

      
      
      
#ifdef __HAVE_MPI__
      END IF
      
      !---------------------------------------------------------------------
      !  Determine MPI data type
      !---------------------------------------------------------------------
#ifdef   __SP
          MPTYPE = MPI_REAL
#else
          MPTYPE = MPI_DOUBLE_PRECISION
#endif
      
      
      !-----------------------------------------------------------------------------
      !  Communicate readin file with all processes
      !-----------------------------------------------------------------------------
      
      CALL MPI_Bcast(options%dimensions,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%rel_sigma,1,MPTYPE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%abs_sigma,1,MPTYPE,0,MPI_COMM_WORLD,ierr)      
      CALL MPI_Bcast(options%alldim_lBounds,1,MPTYPE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%alldim_uBounds,1,MPTYPE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%PopSize,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%use_CEC,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%use_LJ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%use_TIP,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%use_BBOB,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%Benchfctnr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%output_folder,200,MPI_CHARACTER,0,&
     &MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%StopFitness,1,MPTYPE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%StopMaxFunEvals,1,MPTYPE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%StopMaxIter,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%StopTolX,1,MPTYPE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%StopTolUpX,1,MPTYPE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%StopTolFun,1,MPTYPE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%StopTolHistFun,1,MPTYPE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%StopOnWarnings,1,MPI_LOGICAL,0,MPI_COMM_WORLD,&
     &ierr)
      CALL MPI_Bcast(options%WarnOnEqualFunctionValues,1,MPI_LOGICAL,0,&
     &MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%EvalInitialX,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%ParentNumber,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%RecombinationWeights,1,MPI_INTEGER,0,&
     &MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%VerboseModulo,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%flgGenData,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%intGenData,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%funcName,200,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%CECFolders,200,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%silentMPI,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%pscma,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%psoWeight,1,MPTYPE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%psoFreq,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%accuracy,1,MPTYPE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%global_min,1,MPTYPE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%use_seed,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%seed_folder,200,MPI_CHARACTER,0,MPI_COMM_WORLD,&
     &ierr)
      CALL MPI_Bcast(options%qr_sampling,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%qr_sampler,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%qr_inverter,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%use_init_bounds,1,MPI_LOGICAL,0,MPI_COMM_WORLD,&
     &ierr)
      CALL MPI_Bcast(options%init_uBounds,1,MPTYPE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%init_lBounds,1,MPTYPE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%benchmark,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%record_accuracy,1,MPTYPE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%restart_cma,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%restart_type,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%Restarts,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%IncPopSize,1,MPTYPE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%record_besthist,1,MPI_LOGICAL,0,MPI_COMM_WORLD,&
     &ierr)
      CALL MPI_Bcast(options%record_modulo,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%MaxIncFac,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%BFGS_use,1,MPI_LOGICAL,0,MPI_COMM_WORLD,&
     &ierr)
      CALL MPI_Bcast(options%BFGS_position,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%BFGS_factr,1,MPTYPE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%BFGS_pgtol,1,MPTYPE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%BFGS_central_difference,1,MPI_LOGICAL,0,&
     &MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%BFGS_confprop,1,MPTYPE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%BFGS_sigmastep_min,1,MPTYPE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%StopTime,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%StopTimeHH,1,MPTYPE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%StopTimeMM,1,MPTYPE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%StopTimeSS,1,MPTYPE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%write_pdb,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%LJ_comp,1,MPTYPE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%use_LJ_comp,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%use_MATFUNC,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%flgouttxt,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%flgGenTrace,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%use_DF,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%DF_rast,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%DF_s,1,MPTYPE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast(options%use_RANDOM_LANDSCAPE,1,MPI_LOGICAL,0&
     &,MPI_COMM_WORLD,ierr)

      
      
      CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
#endif      


      !-------------------------------------------------------------------------
      ! Initalize the Seed for the Programm
      !------------------------------------------------------------------------
      
      IF (.NOT. options%use_seed) THEN
        CALL SYSTEM_CLOCK(options%seed)
#ifdef __HAVE_MPI__    
	  WRITE(proc,'(I10)') MY_RANK
      options%seed  = options%seed + MY_RANK**2 !making sure that every proc has different seed
#endif
      ELSE 
      !-----------------------------------------------------------------------------
      ! open the file
      !-----------------------------------------------------------------------------
      
      options%seed_folder = TRIM(adjustl(options%seed_folder)) &      
     & // '/seed'& 
     &      // trim(adjustl(proc))  &
     &      //  '.txt'
      OPEN(iUnit, FILE=options%seed_folder, IOSTAT=ios, ACTION='READ')
      IF(ios .NE. 0) THEN
         WRITE(*,'(2A)')'Failed to open file: ',options%seed_folder
         Info = 1
         STOP
      END IF
      
      READ(iUnit,*) options%seed
      !WRITE(*,*)'Seed : ',proc, ' ' ,options%seed
      
      CLOSE(iUnit)
      
      !    carg = TRIM(options%seed_folder)
      !    OPEN(40, FILE=carg,IOSTAT=ios, ACTION='READ')
      !    IF(ios .NE. 0) THEN
      !       WRITE(*,'(2A)')'Failed to open file: ',carg
      !    END IF
      !    READ (40,*,iostat=ios) options%seed
      !     IF(ios .NE. 0) THEN
      !       WRITE(*,'(2A)')'Failed to open file: ',carg
      !    END IF
      END IF


      !-----------------------------------------------------------------------------
      ! Process readin values
      !-----------------------------------------------------------------------------
      
      
      options%output_folder = trim(adjustl(options%output_folder))
      IF(options%dimensions .NE. 0) THEN    !Dimension must be set!!!!
          ALLOCATE(options%LBounds(options%dimensions))
          ALLOCATE(options%UBounds(options%dimensions))
          

          options%LBounds = options%alldim_lBounds
          options%UBounds = options%alldim_uBounds
          
          ALLOCATE(options%insigma(options%dimensions))
        IF (options%alldim_lBounds .EQ. negInf .OR. &
      &     options%alldim_uBounds .EQ. posInf) THEN
            IF (.NOT. options%use_init_bounds) THEN
                		WRITE(*,*)	&
				'Bounds not set and no init Bounds given',&
	       		' Therefore Init Bounds are set to [0,1] to get',&
	       		' a reasonable startpoint'
	       		WRITE(*,*)
	        options%use_init_bounds = .TRUE.
	        options%init_lBounds = 0.0_MK
	        options%init_uBounds = 1.0_MK
            END IF
        END IF
        
        IF (options%abs_sigma .NE. 0.0_MK) THEN
            options%insigma(:) = options%abs_sigma
        ELSEIF(options%rel_sigma .NE. 0) THEN
         IF (options%use_init_bounds ) THEN
            options%insigma(:) = options%rel_sigma*&
      &                       (options%init_uBounds-options%init_lBounds)
         ELSE
            options%insigma(:) = options%rel_sigma*&
      &                       MAXVAL(options%UBounds-options%LBounds)
         END IF
        END IF
        
      IF (options%global_min .NE. 0d0 .AND. options%accuracy .NE. 0d0) THEN
        options%StopFitness = options%global_min + options%accuracy
      END IF
        


        
      !-----------------------------------------------------------------------------
      ! Set the inital Starting Point
      !-----------------------------------------------------------------------------
       CALL cmaes_startpoint()
        
        
        
      
      
      END IF
      
      
      
      
      
      !-----------------------------------------------------------------------------
      !  Return 
      !-----------------------------------------------------------------------------
      IF (.FALSE.) THEN
 9999 STOP
      END IF
      END SUBROUTINE cmaes_readparams




