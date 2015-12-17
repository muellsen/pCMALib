      !-------------------------------------------------------------------------
      !  Module       :            CEC2005
      !-------------------------------------------------------------------------
      !
      !  Purpose      : a fortran implementation of the CEC05 Test Functions
      !
      !  Remarks      : you always need to call init_benchmark at start and 
      !                 close_benchmark at the end
      !                 use "benchmark" to call a specific benchmark function
      !
      !  References   : Problem Definitions and Evaluation Criteria for the 
      !                 CEC 2005 Special Session on Real-Parameter Optimization
      !                 P. N. Suganthan,N. Hansen,J. J. Liang,K. Deb,
      !                 Y. -P. Chen,A. Auger,S. Tiwari
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  pCMALib: a parallel fortran 90 library for the evolution strategy with
      !           covariance matrix adaptation
      !  Christian L. Mueller, Benedikt Baumgartner, Georg Ofenbeck
      !  MOSAIC group, ETH Zurich, Switzerland
      !-------------------------------------------------------------------------
      
      MODULE CEC2005      
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

      IMPLICIT none      

      

      LOGICAL,PARAMETER                             :: noise = .TRUE.


       
#ifndef __CMA    
       INTEGER,PARAMETER                             :: MK = 8
#endif

     
      REAL(MK)                          :: testf_NAN  = huge(testf_NAN)
      REAL(MK),DIMENSION(:),POINTER     :: biases
      REAL(MK),DIMENSION(:),POINTER     :: shiftmatrix

      REAL(MK),DIMENSION(:,:),POINTER   :: rotatematrix
      REAL(MK),DIMENSION(:,:),POINTER   :: schwefel_m_a,schwefel_m_b  
      REAL(MK),DIMENSION(:),POINTER     :: schwefel_mA
      CHARACTER(len=200)                :: function_name
      CHARACTER(len=200)                :: path
      INTEGER                           :: prepared_function_nr
      REAL(MK)                          :: Pi = 3.141592653589793&
     &23846264338327950288419716939937510d0                  
      
      !-------------------------------------------------------------------------
      !  Fields for Function Nr 15
      ! 
      !-------------------------------------------------------------------------
      REAL(MK),DIMENSION(:,:),POINTER                :: hybrid_shiftmatrix
      REAL(MK),DIMENSION(:,:,:),POINTER              :: hybrid_identity


      INTEGER  :: job_num_func = 10
      REAL(MK),DIMENSION(10)    :: job_sigma = &
     &               (/1.0_MK,1.0_MK,1.0_MK,1.0_MK,1.0_MK,&
     &            1.0_MK,1.0_MK,1.0_MK,1.0_MK,1.0_MK /)
      REAL(MK),DIMENSION(10)    :: job_lambda = &
     &          (/1.0_MK,1.0_MK,10.0_MK,10.0_MK,&
     &            5.0_MK/60.0_MK,5.0_MK/60.0_MK,5.0_MK/32.0_MK,5.0_MK/32.0_MK,&
     &            5.0_MK/100.0_MK,5.0_MK/100.0_MK /)
      REAL(MK),DIMENSION(10)    :: job_biases = &
     &        (/0.0_MK,100.0_MK,200.0_MK,300.0_MK,400.0_MK,&
     &            500.0_MK,600.0_MK,700.0_MK,800.0_MK,900.0_MK /)
        
      REAL(MK),DIMENSION(10)    :: job2_sigma = &
     &              (/1.0_MK,2.0_MK,1.5_MK,1.5_MK,1.0_MK,&
     &                1.0_MK,1.5_MK,1.5_MK,2.0_MK,2.0_MK /)
                    
      REAL(MK),DIMENSION(10)    :: job2_lambda = &
     &      (/2.0_MK*5.0_MK/32.0_MK,5.0_MK/32.0_MK,2.0_MK,1.0_MK,&
     &        2.0_MK*5.0_MK/100.0_MK,5.0_MK/100.0_MK,2.0_MK*10.0_MK,10.0_MK,&
     &        2.0_MK*5.0_MK/60.0_MK,5.0_MK/60.0_MK /)                 
 
 
      REAL(MK),DIMENSION(10)    :: job19_sigma = &
     &      (/0.10_MK,2.0_MK,1.5_MK,1.5_MK,1.0_MK,&
     &        1.0_MK,1.5_MK,1.5_MK,2.0_MK,2.0_MK /)
      REAL(MK),DIMENSION(10)    :: job19_lambda = &
     &      (/0.10_MK*5.0_MK/32.0_MK,5.0_MK/32.0_MK,2.0_MK,1.0_MK,&
     &        2.0_MK*5.0_MK/100.0_MK,5.0_MK/100.0_MK,2.0_MK*10.0_MK,10.0_MK,&
     &        2.0_MK*5.0_MK/60.0_MK,5.0_MK/60.0_MK /)
       
      REAL(MK),DIMENSION(10)    :: job21_lambda = &
     &      (/5.0_MK*5.0_MK/100.0_MK,5.0_MK/100.0_MK,5.0_MK,1.0_MK,5.0_MK,&
     &   1.0_MK,5.0_MK*10.0_MK,10.0_MK,5.0_MK*5.0_MK/200.0_MK,5.0_MK/200.0_MK/)            
       
      REAL(MK),DIMENSION(10)    :: job21_sigma = &
     &                  (/1.0_MK,1.0_MK,1.0_MK,1.0_MK,1.0_MK,&
     &                    2.0_MK,2.0_MK,2.0_MK,2.0_MK,2.0_MK /)            
 
      REAL(MK),DIMENSION(10)    :: job24_sigma = &
     &                  (/2.0_MK,2.0_MK,2.0_MK,2.0_MK,2.0_MK,&
     &                    2.0_MK,2.0_MK,2.0_MK,2.0_MK,2.0_MK /)
 
      REAL(MK),DIMENSION(10)    :: job24_lambda = &
     &      (/10.0_MK,5.0_MK/20.0_MK,1.0_MK,5.0_MK/32.0_MK,1.0_MK,&
     &   5_MK/100.0_MK,5.0_MK/50.0_MK,1.0_MK,5.0_MK/100.0_MK,5.0_MK/100.0_MK /)            


                            
      !-------------------------------------------------------------------------
      !  Variables for the Job Function (used in the hybrid functions)
      ! 
      !-------------------------------------------------------------------------
      ! Predefined constant
      REAL(MK)    :: job_C
      ! Estimated fmax
      REAL(MK),POINTER,DIMENSION(:)   :: job_fmax
      ! Shift global optimum for each basic function
      REAL(MK),POINTER,DIMENSION(:,:)   :: job_o
      ! Linear transformation matrix for each basic function
      REAL(MK),POINTER,DIMENSION(:,:,:)   :: job_M
    
      ! Working areas 
      REAL(MK),POINTER,DIMENSION(:) ::    job_w
      REAL(MK),POINTER,DIMENSION(:,:) ::  job_z
      REAL(MK),POINTER,DIMENSION(:,:) ::  job_zM
        
      INTEGER                         ::  job_num_dim
        
        
        
        
      CONTAINS
      !-------------------------------------------------------------------------
      ! benchmark
      !-------------------------------------------------------------------------
      !  Purpose      : runs benchmark on a given functionnr
      ! 
      !  Input       :  funcnr  (I) the number of the function to be benchmarked 
      !                 m       (I) m Dimension of the matrix (Rows)
      !                               = Dimension of Vectors
      !                 n       (I) n Dimension of the matrix (Columns)
      !                               = # of Vectors to process
      !  Input(optional): datafile (S)  TRIM(adjustl(path)) to the datafile used with the function
      !                   rotatefile(S) TRIM(adjustl(path)) to the rotatematrix used with the function
      !                   rotatesuffix(S) wich datasuffix (e.g. .txt) to use for the rotatefile
      !                   
      !
      !  Input/Output:  vars    (R) the matrix with the input values of size m*n
      !
      !  Output:        name    (S) string containing the name of the benchmark function
      !res(R) the vector with the results
      !
      !-------------------------------------------------------------------------
      
      
      !biases,shiftmatrix,rotatematrix
      SUBROUTINE benchmark(res,vars_in,m,n,boundarray_down,boundarray_up,&
     &               biasfile,datafile,rotatefile,rotatesuffix)
      !parameters

      INTEGER,INTENT(in)                            :: m
      INTEGER,INTENT(in)                            :: n
      REAL(MK),DIMENSION(n),INTENT(out)           :: res
      REAL(MK),DIMENSION(m,n),INTENT(IN)            :: vars_in
      
      CHARACTER(*),OPTIONAL,INTENT(IN)               :: biasfile 
      !REAL(MK),DIMENSION(:),POINTER                  :: shiftmatrix
      !REAL(MK),DIMENSION(:,:),POINTER                :: rotatematrix
      CHARACTER(*),OPTIONAL,INTENT(IN)               :: datafile
      CHARACTER(*),OPTIONAL,INTENT(IN)               :: rotatefile
      CHARACTER(*),OPTIONAL,INTENT(IN)               :: rotatesuffix
      REAL(MK),OPTIONAL,DIMENSION(m)      :: boundarray_up
      REAL(MK),OPTIONAL,DIMENSION(m)      :: boundarray_down
      !localy used variables    
      INTEGER                                        :: funcnr
      REAL(MK),DIMENSION(m,n)                        :: vars
      CHARACTER(len=200)                             :: usedfile
      CHARACTER(len=200)                             :: usedrotate
      CHARACTER(len=200)                             :: usedsuffix
      CHARACTER(len=5)                               :: nrstring
      INTEGER                                        :: alloc_status
      INTEGER                                        :: dealloc_status,i,j,k
      REAL(MK)                                       :: dummy      
      REAL(MK),DIMENSION(:,:),POINTER                :: tmp
      REAL(MK),DIMENSION(m)                          :: tmp_vec
      REAL(MK),DIMENSION(m)                          :: tmp_result

      LOGICAL :: t1,t2
      funcnr = prepared_function_nr
             
      vars = vars_in
      
      t1 = PRESENT(boundarray_up)
      t2 = PRESENT(boundarray_down)

      IF(.NOT. PRESENT(boundarray_up) .OR. .NOT. PRESENT(boundarray_down)) THEN
      !cause some fortran compiler would call with non existent optional arguments
      !we have to make 2 different versions of the call
      IF (funcnr .EQ. 1) THEN        
        CALL F01_shifted_sphere(res,vars,m,n) 
      ENDIF
      IF (funcnr .EQ. 2) THEN
        CALL F02_shifted_schwefel(res,vars,m,n)
      ENDIF
      IF (funcnr .EQ. 3) THEN
        CALL F03_shifted_rotated_high_cond_elliptic(res,vars,m,n)
      ENDIF
      IF (funcnr .EQ. 4) THEN
        CALL F04_shifted_schwefel_noise(res,vars,m,n)
      ENDIF
      IF (funcnr .EQ. 5) THEN
        CALL F05_schwefel_global_opt_bound(res,vars,m,n)
      ENDIF
      IF (funcnr .EQ. 6) THEN
        CALL F06_shifted_rosenbrock(res,vars,m,n)
      ENDIF
      IF (funcnr .EQ. 7) THEN
        CALL F07_shifted_rotated_griewank(res,vars,m,n)
      ENDIF    
      IF (funcnr .EQ. 8) THEN
        CALL F08_shifted_rotated_ackley_global_opt_bound(res,vars,m,n)
      ENDIF    
      IF (funcnr .EQ. 9) THEN
        CALL F09_shifted_rastrigin(res,vars,m,n)
      ENDIF    
      IF (funcnr .EQ. 10) THEN
        CALL F10_shifted_rotated_rastrigin(res,vars,m,n)
      ENDIF    
      IF (funcnr .EQ. 11) THEN
        CALL F11_shifted_rotated_weierstrass(res,vars,m,n)
      ENDIF
      IF (funcnr .EQ. 12) THEN
        CALL F12_schwefel(res,vars,m,n)
      ENDIF
      IF (funcnr .EQ. 13) THEN
        CALL F13_shifted_expanded_griewank_rosenbrock(res,vars,m,n)
      ENDIF
      IF (funcnr .EQ. 14) THEN
        CALL F14_shifted_rotated_expanded_scaffer(res,vars,m,n)
      ENDIF
      IF (funcnr .EQ. 15) THEN
        DO i = 1,n        
            CALL hybrid_composition(vars(:,i),dummy,m,15)        
            res(i) = dummy
        END DO
      ENDIF
      IF (funcnr .EQ. 16) THEN
        DO i = 1,n        
            CALL hybrid_composition(vars(:,i),dummy,m,16)        
            res(i) = dummy
        END DO
      ENDIF
      IF (funcnr .EQ. 17) THEN
        DO i = 1,n        
            CALL hybrid_composition(vars(:,i),dummy,m,17)        
            res(i) = dummy
        END DO
      ENDIF
      IF (funcnr .EQ. 18) THEN
        DO i = 1,n        
            CALL hybrid_composition(vars(:,i),dummy,m,18)        
            res(i) = dummy
        END DO
      ENDIF
      IF (funcnr .EQ. 19) THEN
        DO i = 1,n        
            CALL hybrid_composition(vars(:,i),dummy,m,19)        
            res(i) = dummy
        END DO
      ENDIF      
      IF (funcnr .EQ. 20) THEN
        DO i = 1,n        
            CALL hybrid_composition(vars(:,i),dummy,m,20)        
            res(i) = dummy
        END DO
      ENDIF            
      IF (funcnr .EQ. 21) THEN
        DO i = 1,n        
            CALL hybrid_composition(vars(:,i),dummy,m,21)        
            res(i) = dummy
        END DO
      ENDIF      
      IF (funcnr .EQ. 22) THEN
        

        DO i = 1,n        
            CALL hybrid_composition(vars(:,i),dummy,m,22)        
            res(i) = dummy
        END DO
      ENDIF      
      IF (funcnr .EQ. 23) THEN
        
        DO i = 1,n        
            CALL hybrid_composition(vars(:,i),dummy,m,23)        
            res(i) = dummy
        END DO
      ENDIF       
      IF (funcnr .EQ. 24) THEN
        
        DO i = 1,n        
            CALL hybrid_composition(vars(:,i),dummy,m,24)        
            res(i) = dummy
        END DO
      ENDIF
      
      IF (funcnr .EQ. 25) THEN
        
        DO i = 1,n        
            CALL hybrid_composition(vars(:,i),dummy,m,25)        
            res(i) = dummy
        END DO
      ENDIF            
      ELSE
      !give the name and call the function
      IF (funcnr .EQ. 1) THEN
        
        CALL F01_shifted_sphere(res,vars,m,n,boundarray_down,boundarray_up) 
      ENDIF
      IF (funcnr .EQ. 2) THEN
        CALL F02_shifted_schwefel(res,vars,m,n,boundarray_down,boundarray_up)
      ENDIF
      IF (funcnr .EQ. 3) THEN
        CALL F03_shifted_rotated_high_cond_elliptic(res,vars,m,n,&
     &                                           boundarray_down,boundarray_up)
      ENDIF
      IF (funcnr .EQ. 4) THEN
        CALL F04_shifted_schwefel_noise(res,vars,m,n,boundarray_down,&
     &                                                           boundarray_up)
      ENDIF
      IF (funcnr .EQ. 5) THEN
        CALL F05_schwefel_global_opt_bound(res,vars,m,n,boundarray_down,&
     &                                                           boundarray_up)
      ENDIF
      IF (funcnr .EQ. 6) THEN
        CALL F06_shifted_rosenbrock(res,vars,m,n,boundarray_down,boundarray_up)
      ENDIF
      IF (funcnr .EQ. 7) THEN
        CALL F07_shifted_rotated_griewank(res,vars,m,n,boundarray_down,&
     &                                                           boundarray_up)
      ENDIF    
      IF (funcnr .EQ. 8) THEN
        CALL F08_shifted_rotated_ackley_global_opt_bound(res,vars,m,n,&
     &                                           boundarray_down,boundarray_up)
      ENDIF    
      IF (funcnr .EQ. 9) THEN
        CALL F09_shifted_rastrigin(res,vars,m,n,boundarray_down,boundarray_up)
      ENDIF    
      IF (funcnr .EQ. 10) THEN
        CALL F10_shifted_rotated_rastrigin(res,vars,m,n,boundarray_down,&
     &                                                           boundarray_up)
      ENDIF    
      IF (funcnr .EQ. 11) THEN
        CALL F11_shifted_rotated_weierstrass(res,vars,m,n,boundarray_down,&
     &                                                           boundarray_up)
      ENDIF
      IF (funcnr .EQ. 12) THEN
        CALL F12_schwefel(res,vars,m,n,boundarray_down,boundarray_up)
      ENDIF
      IF (funcnr .EQ. 13) THEN
        CALL F13_shifted_expanded_griewank_rosenbrock(res,vars,m,n,&
     &                                           boundarray_down,boundarray_up)
      ENDIF
      IF (funcnr .EQ. 14) THEN
        CALL F14_shifted_rotated_expanded_scaffer(res,vars,m,n,boundarray_down&
     &                                                          ,boundarray_up)
      ENDIF
      IF (funcnr .EQ. 15) THEN
        DO i = 1,n        
            CALL hybrid_composition(vars(:,i),dummy,m,15,boundarray_down&
     &                                                          ,boundarray_up)        
            res(i) = dummy
        END DO
      ENDIF
      IF (funcnr .EQ. 16) THEN
        DO i = 1,n        
            CALL hybrid_composition(vars(:,i),dummy,m,16,boundarray_down&
     &                                                          ,boundarray_up)       
            res(i) = dummy
        END DO
      ENDIF
      IF (funcnr .EQ. 17) THEN
        DO i = 1,n        
            CALL hybrid_composition(vars(:,i),dummy,m,17,boundarray_down&
     &                                                          ,boundarray_up)        
            res(i) = dummy
        END DO
      ENDIF
      IF (funcnr .EQ. 18) THEN
        DO i = 1,n        
            CALL hybrid_composition(vars(:,i),dummy,m,18,boundarray_down&
     &                                                          ,boundarray_up)        
            res(i) = dummy
        END DO
      ENDIF
      IF (funcnr .EQ. 19) THEN
        DO i = 1,n        
            CALL hybrid_composition(vars(:,i),dummy,m,19,boundarray_down&
     &                                                          ,boundarray_up)        
            res(i) = dummy
        END DO
      ENDIF      
      IF (funcnr .EQ. 20) THEN
        DO i = 1,n        
            CALL hybrid_composition(vars(:,i),dummy,m,20,boundarray_down&
     &                                                          ,boundarray_up)        
            res(i) = dummy
        END DO
      ENDIF            
      IF (funcnr .EQ. 21) THEN
        DO i = 1,n        
            CALL hybrid_composition(vars(:,i),dummy,m,21,boundarray_down&
     &                                                          ,boundarray_up)       
            res(i) = dummy
        END DO
      ENDIF      
      IF (funcnr .EQ. 22) THEN
        
        
        DO i = 1,n        
            CALL hybrid_composition(vars(:,i),dummy,m,22,boundarray_down&
     &                                                          ,boundarray_up)
            res(i) = dummy
        END DO
      ENDIF      
      IF (funcnr .EQ. 23) THEN
        
        DO i = 1,n        
            CALL hybrid_composition(vars(:,i),dummy,m,23,boundarray_down&
     &                                                          ,boundarray_up)
            res(i) = dummy
        END DO
      ENDIF       
      IF (funcnr .EQ. 24) THEN
        
        DO i = 1,n        
            CALL hybrid_composition(vars(:,i),dummy,m,24,boundarray_down&
     &                                                          ,boundarray_up)
            res(i) = dummy
        END DO
      ENDIF
      
      IF (funcnr .EQ. 25) THEN
        
        DO i = 1,n        
            CALL hybrid_composition(vars(:,i),dummy,m,25,boundarray_down&
     &                                                          ,boundarray_up)
            res(i) = dummy
        END DO
      ENDIF            
      !end if - if the bounds are here
      ENDIF
      
      
#ifdef __CMA
      options%funcName = function_name
#endif

      
      !deallocate(shiftmatrix,stat=dealloc_status)
      !deallocate(rotatematrix,stat=dealloc_status)

      END SUBROUTINE benchmark
      
      !-------------------------------------------------------------------------
      ! prepare_benchmark
      !-------------------------------------------------------------------------
      !  Purpose      : runs benchmark on a given functionnr
      ! 
      !  Input       :  funcnr  (I) the number of the function to be benchmarked 
      !                 m       (I) m Dimension of the matrix (Rows)
      !                               = Dimension of Vectors
      !                 n       (I) n Dimension of the matrix (Columns)
      !                               = # of Vectors to process
      !  Input(optional): datafile (S)  TRIM(adjustl(path)) to the datafile used with the function
      !                   rotatefile(S) TRIM(adjustl(path)) to the rotatematrix used with the function
      !                   rotatesuffix(S) wich datasuffix (e.g. .txt) to use for the rotatefile
      !                   
      !
      !  Input/Output:  
      !
      !  Output:        name    (S) string containing the name of the benchmark function
      !
      !
      !-------------------------------------------------------------------------
      
      
      !biases,shiftmatrix,rotatematrix
      SUBROUTINE prepare_benchmark(funcnr,m,n,biasfile,datafile,rotatefile,&
     &                                                            rotatesuffix)
      !parameters
      !REAL(MK),DIMENSION(:),POINTER                  :: shiftmatrix
      INTEGER,INTENT(IN)                             :: funcnr
                      
      INTEGER,INTENT(in)                            :: m
      INTEGER,INTENT(in)                            :: n
      CHARACTER(*),OPTIONAL,INTENT(IN)               :: biasfile 
      CHARACTER(*),OPTIONAL,INTENT(IN)               :: datafile
      CHARACTER(*),OPTIONAL,INTENT(IN)               :: rotatefile
      CHARACTER(*),OPTIONAL,INTENT(IN)               :: rotatesuffix
          
      !localy used variables    
      CHARACTER(len=200)                             :: usedfile
      CHARACTER(len=200)                             :: usedrotate
      CHARACTER(len=200)                             :: usedsuffix
      CHARACTER(len=200)                             :: usedbias
      CHARACTER(len=5)                               :: nrstring
      INTEGER                                        :: alloc_status
      INTEGER                                        :: dealloc_status,i,j,k
      REAL(MK)                                       :: dummy      
      REAL(MK),DIMENSION(:,:),POINTER                :: tmp
      REAL(MK),DIMENSION(m)                          :: tmp_vec
      REAL(MK),DIMENSION(m)                          :: tmp_result   
      REAL(MK),DIMENSION(:,:),POINTER                :: f12_data
      
      
      !variables used for the hybrid functions
      REAL(MK),DIMENSION(m)                        :: job_testPoint
      REAL(MK),DIMENSION(m)                        :: job_testPointM
      
      !save in module var, so that benchmark nows which function
      prepared_function_nr = funcnr
      

      !Seed the random Number Generator (0 = system time)
#ifdef __CMA
      CALL ZBQLINI(options%seed)
      path = options%CECFolders // '/'
#else
      CALL ZBQLINI(1)
      path = ''
#endif
      
                                      


      allocate(shiftmatrix(m),stat=alloc_status)
      allocate(rotatematrix(m,m),stat=alloc_status)
      allocate(tmp(m+1,m),stat=alloc_status)
      
      IF (funcnr .EQ. 12) THEN
        allocate(schwefel_m_a(m,m),stat=alloc_status)
        allocate(schwefel_m_b(m,m),stat=alloc_status)
        allocate(schwefel_mA(m),stat=alloc_status)
      ENDIF   
      allocate(biases(25),stat=alloc_status)
      
      IF (.NOT. PRESENT(biasfile)) THEN
            usedbias = TRIM(adjustl(path)) // 'supportData/fbias_data.txt'
      ELSE
            usedbias = biasfile
      ENDIF
      
        
      
      CALL loadRowVectorFromFile(usedbias,25,biases)  
      IF (funcnr .EQ. 1) THEN
        function_name = 'Shifted Sphere Function'
        END IF
      IF (funcnr .EQ. 2) THEN
        function_name = 'Shifted Schwefels Problem 1.2'
        END IF
      IF (funcnr .EQ. 3) THEN
        function_name = 'Shifted Rotated High Conditioned Elliptic Function'
        END IF
      IF (funcnr .EQ. 4) THEN
        function_name = 'Shifted Schwefels Problem 1.2 with Noise in Fitness'
        END IF
      IF (funcnr .EQ. 5) THEN
        function_name = 'Schwefels Problem 2.6 with Global Optimum on Bounds'
        END IF
      IF (funcnr .EQ. 6) THEN
        function_name = 'Shifted Rosenbrocks Function'
        END IF
      IF (funcnr .EQ. 7) THEN
        function_name = 'Shifted Rotated Griewanks Function without Bounds'
        END IF    
      IF (funcnr .EQ. 8) THEN
        function_name = 'Shifted Rotated Ackleys Function with Global Optimum'
        function_name = function_name//'on Bounds'
        END IF    
      IF (funcnr .EQ. 9) THEN
        function_name = 'Shifted Rastrigins Function'
        END IF    
      IF (funcnr .EQ. 10) THEN
        function_name = 'Shifted Rotated Rastrigins Function'
        END IF    
      IF (funcnr .EQ. 11) THEN
        function_name = 'Shifted Rotated Weierstrass Function'
        END IF
      IF (funcnr .EQ. 12) THEN
        function_name = 'Schwefels Problem 2.13'
        END IF
      IF (funcnr .EQ. 13) THEN
        function_name = 'Shifted Expanded Griewanks plus Rosenbrocks Function'
        END IF
      IF (funcnr .EQ. 14) THEN
        function_name = 'Shifted Rotated Expanded Scaffers F6 Function'
        END IF
      IF (funcnr .EQ. 15) THEN
        function_name = 'Hybrid Composition Function 1'
      END IF
      IF (funcnr .EQ. 16) then
        function_name = 'Rotated Hybrid Composition Function 1'
      END IF
      IF (funcnr .EQ. 17) THEN
        function_name = 'Rotated Hybrid Composition Function 1 '
        function_name = function_name // 'with Noise in Fitness'
      END IF
      IF (funcnr .EQ. 18) THEN
        function_name = 'Rotated Hybrid Composition Function 2'
      END IF
      IF (funcnr .EQ. 19) THEN
        function_name = 'Rotated Hybrid Composition Function 2 with narrow'
        function_name = function_name // ' basin global optimum'
      END IF      
      IF (funcnr .EQ. 20) THEN
        function_name = 'Rotated Hybrid Composition Function 2 with Global'
        function_name = function_name // ' Optimimum on the Bounds'
      END IF            
      IF (funcnr .EQ. 21) THEN
        function_name = 'Rotated Hybrid Composition Function 3'
      END IF      
      IF (funcnr .EQ. 22) THEN
        function_name = 'Rotated Hybrid Composition Function 3 with High'
        function_name = function_name //' Condition Number Matrix'
      END IF      
      IF (funcnr .EQ. 23) THEN
        function_name = 'Non-Continuous Rotated Hybrid Composition Function 3'
      END IF       
      IF (funcnr .EQ. 24) THEN
        function_name = 'Rotated Hybrid Composition Function 4'
      END IF
      
      IF (funcnr .EQ. 25) THEN
        function_name = 'Rotated Hybrid Composition Function 4 without bounds'
      END IF            



      IF (.NOT. PRESENT(datafile)) THEN
        IF (funcnr .EQ. 1) THEN
            usedfile = TRIM(adjustl(path)) // 'supportData/sphere_func_data.txt'
        ENDIF
        IF (funcnr .EQ. 2) THEN
            usedfile = TRIM(adjustl(path)) // 'supportData/schwefel_102_data.txt'
        ENDIF
        IF (funcnr .EQ. 3) THEN
            usedfile = TRIM(adjustl(path)) // 'supportData/high_cond_elliptic_rot_data.txt'
        ENDIF       
        IF (funcnr .EQ. 4) THEN
            usedfile = TRIM(adjustl(path)) // 'supportData/schwefel_102_data.txt'
        ENDIF
        IF (funcnr .EQ. 5) THEN
            usedfile = TRIM(adjustl(path)) // 'supportData/schwefel_206_data.txt'
        ENDIF
        IF (funcnr .EQ. 6) THEN
            usedfile = TRIM(adjustl(path)) // 'supportData/rosenbrock_func_data.txt'
        ENDIF
        IF (funcnr .EQ. 7) THEN
            usedfile = TRIM(adjustl(path)) // 'supportData/griewank_func_data.txt'
        ENDIF
        IF (funcnr .EQ. 8) THEN
            usedfile = TRIM(adjustl(path)) // 'supportData/ackley_func_data.txt'
        ENDIF
        IF (funcnr .EQ. 9 .OR. funcnr .EQ. 10) THEN
            usedfile = TRIM(adjustl(path)) // 'supportData/rastrigin_func_data.txt'
        ENDIF
        IF (funcnr .EQ. 11) THEN
            usedfile = TRIM(adjustl(path)) // 'supportData/weierstrass_data.txt'
        ENDIF
        IF (funcnr .EQ. 12) THEN
            usedfile = TRIM(adjustl(path)) // 'supportData/schwefel_213_data.txt'
        ENDIF
        IF (funcnr .EQ. 13) THEN
            usedfile = TRIM(adjustl(path)) // 'supportData/EF8F2_func_data.txt'
        ENDIF
        IF (funcnr.EQ. 14) THEN
            usedfile = TRIM(adjustl(path)) // 'supportData/E_ScafferF6_func_data.txt'
        ENDIF
        IF (funcnr .EQ. 15 .OR. funcnr .EQ. 16 .OR. funcnr .EQ. 17) THEN
            usedfile = TRIM(adjustl(path)) // 'supportData/hybrid_func1_data.txt'
        ENDIF
        IF (funcnr .EQ. 18 .OR. funcnr .EQ. 19 .OR. funcnr .EQ. 20) THEN
            usedfile = TRIM(adjustl(path)) // 'supportData/hybrid_func2_data.txt'
        ENDIF
        IF (funcnr .EQ. 21 .OR. funcnr .EQ. 22 .OR. funcnr .EQ. 23) THEN
            usedfile = TRIM(adjustl(path)) // 'supportData/hybrid_func3_data.txt'
        ENDIF
        IF (funcnr .EQ. 25 .OR. funcnr .EQ. 24 ) THEN
            usedfile = TRIM(adjustl(path)) // 'supportData/hybrid_func4_data.txt'
        ENDIF
        
        
      ELSE
        usedfile = datafile
      END IF
      
      IF (.NOT. PRESENT(rotatesuffix)) THEN
        IF (funcnr .EQ. 3 .OR. funcnr .EQ. 7 .OR. funcnr .EQ. 8 .OR. funcnr &
     &      .EQ. 10 .OR. funcnr .EQ. 11 .OR. funcnr .EQ. 14 .OR. funcnr .EQ. &
     &         16 .OR. funcnr .EQ. 17 .OR. funcnr .EQ. 18 .OR. funcnr .EQ. 19 &
     &       .OR. funcnr .EQ. 20 .OR. funcnr .EQ. 21 .OR. funcnr .EQ. 22 .OR. &
     &            funcnr .EQ. 23 .OR. funcnr .EQ. 24 .OR. funcnr .EQ. 25) THEN
            usedsuffix = '.txt'
        ENDIF
      ELSE
        usedsuffix = rotatesuffix
      ENDIF   
      write(nrstring,'(I3)') m    
      IF (.NOT. PRESENT(rotatefile)) THEN
        IF (funcnr .EQ. 3) THEN
           usedrotate = TRIM(adjustl(path)) // 'supportData/elliptic_M_D'//&
     &           TRIM(adjustl(nrstring))//usedsuffix           
        ENDIF
        IF (funcnr .EQ. 7) THEN
            usedrotate = TRIM(adjustl(path)) // 'supportData/griewank_M_D'//&
     &            TRIM(adjustl(nrstring))//usedsuffix
        ENDIF
        IF (funcnr .EQ. 8) THEN
            usedrotate = TRIM(adjustl(path)) // 'supportData/ackley_M_D'//&
     &            TRIM(adjustl(nrstring))//usedsuffix
        ENDIF
        IF (funcnr .EQ. 10) THEN
            usedrotate = TRIM(adjustl(path)) // 'supportData/rastrigin_M_D'//&
     &            TRIM(adjustl(nrstring))//usedsuffix
        ENDIF
        IF (funcnr .EQ. 11) THEN
            usedrotate = TRIM(adjustl(path)) // 'supportData/weierstrass_M_D'//&
     &            TRIM(adjustl(nrstring))//usedsuffix
        ENDIF
        IF (funcnr .EQ. 14) THEN
            usedrotate = TRIM(adjustl(path)) // 'supportData/E_ScafferF6_M_D'//&
     &            TRIM(adjustl(nrstring))//usedsuffix
        ENDIF
        IF (funcnr .EQ. 16 .OR. funcnr .EQ. 17) THEN
           usedrotate = TRIM(adjustl(path)) // 'supportData/hybrid_func1_M_D'//&
     &           TRIM(adjustl(nrstring))//usedsuffix
        ENDIF
        IF (funcnr .EQ. 18 .OR. funcnr .EQ. 19 .OR. funcnr .EQ. 20) THEN
            usedrotate = TRIM(adjustl(path)) // 'supportData/hybrid_func2_M_D'//&
     &            TRIM(adjustl(nrstring))//usedsuffix
        ENDIF
        IF (funcnr .EQ. 21 .OR. funcnr .EQ. 23) THEN
            usedrotate = TRIM(adjustl(path)) // 'supportData/hybrid_func3_M_D'//&
     &            TRIM(adjustl(nrstring))//usedsuffix
        ENDIF
        IF (funcnr .EQ. 22) THEN
            usedrotate = TRIM(adjustl(path)) // 'supportData/hybrid_func3_HM_D'//&
     &            TRIM(adjustl(nrstring))//usedsuffix
        ENDIF    
        IF (funcnr .EQ. 24 .OR. funcnr .EQ. 25) THEN
            usedrotate = TRIM(adjustl(path)) // 'supportData/hybrid_func4_M_D'//&
     &            TRIM(adjustl(nrstring))//usedsuffix
        ENDIF
        
      ELSE
        usedrotate = rotatefile
      ENDIF
      IF (funcnr .EQ. 3 .OR. funcnr .EQ. 7 .OR. funcnr .EQ. 8 .OR. funcnr &
     &         .EQ. 10 .OR. funcnr .EQ. 11 .OR. funcnr .EQ. 14) THEN
         CALL loadMatrixFromFile(usedrotate,m,m,rotatematrix)
      END IF

      !shift the input 
      IF ((funcnr .LT. 5 .OR. funcnr .GT. 5) .AND. funcnr .NE. 12 .AND. &
     &          funcnr .LT. 15) THEN
        CALL loadRowVectorFromFile(usedfile,m,shiftmatrix)  
        IF (funcnr .EQ. 6 .OR. funcnr .EQ. 13) THEN
            DO i=1,m
                shiftmatrix(i) = shiftmatrix(i) - 1.0_MK;
            END DO   
        ENDIF
        IF (funcnr .EQ. 8) THEN
            DO i=1,m,2
                shiftmatrix(i) = -32.0_MK
            END DO      
        ENDIF    
      ELSE
        IF (funcnr .EQ. 5)THEN
            CALL loadMatrixFromFile(usedfile,m+1,m,tmp)
            DO i=1,m
                IF ((i) .GT. ceiling(m/4.0_MK)) THEN           
                    IF ((i) .LT. FLOOR((3.0_MK * m)/4.0_MK)) THEN
                        tmp_vec(i) = tmp(1,i)
                    ELSE
                        tmp_vec(i) = 100.0_MK
                    ENDIF                  
                ELSE
                    tmp_vec(i) = -100.0_MK
                ENDIF                                   
            END DO          
            DO i=1,m
                DO j=1,m
                    rotatematrix(j,i) = tmp(j+1,i)
                END DO
            END DO          
            CALL Ax(m,shiftmatrix,tmp_vec,rotatematrix)         
        ENDIF
        
        IF (funcnr .EQ. 12) THEN
            allocate(f12_data(201,m),stat=alloc_status)
            CALL loadMatrixFromFile(usedfile,201,m,f12_data)
            
            DO i = 1,m
                DO j = 1,m
                    schwefel_m_a(i,j) = f12_data(i,j)
                    schwefel_m_b(i,j) = f12_data(100+i,j)
                END DO
                shiftmatrix(i) = f12_data(201,i)
            END DO
            
            DO i = 1,m
                schwefel_mA(i) = 0.0d0
                DO j = 1,m
                 schwefel_mA(i)=schwefel_mA(i)+(schwefel_m_a(i,j)* &
     &               sin(shiftmatrix(j))+schwefel_m_b(i,j)*cos(shiftmatrix(j)))
                END DO
            END DO
            deallocate(f12_data,stat=dealloc_status)
        ENDIF
        
        
        
        
        !handle funcnr 15+ shift loading here
        
        IF (funcnr .GT. 14) THEN
           allocate(hybrid_shiftmatrix(job_num_func,m),stat=alloc_status)
           allocate(hybrid_identity(job_num_func,m,m),stat=alloc_status)
           
           allocate(job_o(job_num_func,m),stat=alloc_status)
           allocate(job_M(job_num_func,m,m),stat=alloc_status)
           
           !allocate(hybrid_testPoint(m),stat=alloc_status)
           !allocate(hybrid_testPointM(m),stat=alloc_status)
           allocate(job_fmax(job_num_func),stat=alloc_status)
           
           allocate(job_w(job_num_func),stat=alloc_status)
           allocate(job_z(job_num_func,m),stat=alloc_status)
           allocate(job_zM(job_num_func,m),stat=alloc_status)
           
           
           CALL loadMatrixFromFile(usedfile,job_num_func,m,job_o)
         
           IF (funcnr .EQ. 18 .OR. funcnr .EQ. 19 .OR. funcnr .EQ. 20) THEN
            DO i = 1,m
                job_o(10,i) = 0.0
            END DO
           ENDIF
           IF (funcnr .EQ. 20) THEN
            DO i = 2,m,2
                job_o(1,i) = 5.0
            END DO
           ENDIF
                
        
           
           
           
           IF (funcnr .EQ. 15) THEN
            ! Generate identiy matrices
            DO i = 1,job_num_func
                    DO j = 1,m
                        DO k = 1,m
                            job_M(i,j,k) = 0.0_MK
                        END DO
                    END DO
                    DO j = 1,m
                        job_M(i,j,j) = 1.0_MK
                    END DO
            END DO
           ELSE
                CALL loadNMatrixFromFile(usedrotate,job_num_func,m,m,job_M)
           END IF
           

           job_num_dim = m
           job_C = 2000.0_MK
      
           !Calculate/estimate the fmax for all the functions involed
           DO i = 1,job_num_func
                DO j=1,m
                    IF (funcnr .LT. 18) THEN
                        job_testPoint(j) = (5.0_MK / job_lambda(i))
                    ELSE
                        IF (funcnr .EQ. 18 .OR. funcnr .EQ. 20) THEN
                            job_testPoint(j) = (5.0_MK / job2_lambda(i))
                        ELSE
                            IF (funcnr .EQ. 19) THEN
                                job_testPoint(j) = (5.0_MK / job19_lambda(i))
                            ELSE
                             IF (funcnr .EQ. 21 .OR. funcnr .EQ. 22 .OR. &
     &                                 funcnr .EQ. 23) THEN
                                 job_testPoint(j) = (5.0_MK / job21_lambda(i))
                             ELSE
                                 job_testPoint(j) = (5.0_MK / job24_lambda(i))
                             ENDIF
                            END IF
                        END IF
                    END IF
                END DO
                
                
                !benchmark.rotate(m_testPointM,m_testPoint,m_M[i]);
                CALL xA(m,job_testPointM,job_testPoint,job_M(i,:,:))
              
                If(funcnr .EQ. 15 .OR. funcnr .EQ. 16 .OR. funcnr .EQ. 17) THEN
                    CALL Job15(i,m,job_testPointM,dummy)
                ELSE
                    IF (funcnr .LT. 21) THEN
                        CALL Job18(i,m,job_testPointM,dummy)
                    ELSE
                        IF (funcnr .GT. 23) THEN
                            CALL Job24(i,m,job_testPointM,dummy)
                        ELSE
                            CALL Job21(i,m,job_testPointM,dummy)
                        END IF
                    ENDIF
                END IF
                job_fmax(i) = abs(dummy)              
           END DO
           i = i
           
            
        ENDIF
        
        
      ENDIF
 
  
      
      !all deallocates
      deallocate(tmp,stat=dealloc_status)

      
      
      
      END SUBROUTINE prepare_benchmark
      
      
      !-------------------------------------------------------------------------
      ! close_benchmark
      !-------------------------------------------------------------------------
      !  Purpose      : just deallocates stuff allocated during init_benchmark
      !-------------------------------------------------------------------------
      
      
      
      SUBROUTINE close_benchmark()
      !parameters
      !REAL(MK),POINTER,DIMENSION(:)               :: biases
      !localy used variables
      INTEGER                                   :: alloc_status,dealloc_status  
      
      deallocate(biases,stat=dealloc_status)
      deallocate(shiftmatrix,stat=dealloc_status)
      deallocate(rotatematrix,stat=dealloc_status)
      deallocate(schwefel_m_a,stat=dealloc_status)
      deallocate(schwefel_m_b,stat=dealloc_status)
      deallocate(schwefel_mA,stat=dealloc_status)
      deallocate(hybrid_shiftmatrix,stat=alloc_status)
      deallocate(hybrid_identity,stat=alloc_status)
      deallocate(job_o,stat=alloc_status)
      deallocate(job_M,stat=alloc_status)           
      deallocate(job_fmax,stat=alloc_status)
      deallocate(job_w,stat=alloc_status)
      deallocate(job_z,stat=alloc_status)
      deallocate(job_zM,stat=alloc_status)
      
      END SUBROUTINE close_benchmark    
      
      
      !-------------------------------------------------------------------------
      ! loadRowVectorFromFile
      !-------------------------------------------------------------------------
      !  Purpose      : loads a row vector from a file (space seperated)
      ! 
      !  Input       :  filename (S)   TRIM(adjustl(path)) to the File from which to read a Row Vector
      !                 m        (I)   m Dimension of the vector
      !  Input/Output:  vector   (P,R) Pointer to the resulting vector
      !-------------------------------------------------------------------------
      
      SUBROUTINE loadRowVectorFromFile (filename,m,vector) 
      !Parameters 
      CHARACTER(*),INTENT(IN)                       :: filename
      INTEGER,INTENT(IN)                           :: m
      REAL(MK),POINTER,DIMENSION(:)                :: vector
      
      INTEGER :: io_error
      INTEGER :: n
      INTEGER :: i
      real    :: y 

      write(*,'(a)') 'File name is:'
      write(*,'(a)') TRIM(filename)
    
      open(unit=20,file=filename,status='old',action='read',&
     &           iostat=io_error) 

      IF ( io_error == 0) then
         read(20,*) vector
         !write(*,*) vector
      else
           IF (io_error > 0) THEN
           write(*,*) 'Error while opening File ', filename, ' - Error Nr:',& 
                       io_error,' occured in Function loadRowVectorFromFile' 
           STOP
           END IF
      end IF
      close(unit=20) 
      END SUBROUTINE loadRowVectorFromFile
      
      !-------------------------------------------------------------------------
      ! loadNMatrixFromFile
      !-------------------------------------------------------------------------
      !  Purpose      : loads a matrix from a file (space / linebreak seperated)
      !                 optional loads a column vector with results appending the matrix      
      ! 
      !  Input       : filename (S) TRIM(adjustl(path)) to the File from which to read a Row Vector
      !                 m       (I) m Dimension of the matrix
      !                 n       (I) n Dimension of the matrix
      !                 repeat  (I) input that tells how many matrices should be read out of a single file
      !                             default = 1
      !  Input/Output:  matrix(P,R) Pointer to the resulting matrix
      !                 res   (P,R) Pointer to a vector with corresponding results
      !                             OPTIONAL Parameter - only use if theres results also in the
      !                                                  matrixfile - dimensions are similar to matrix length m
      !-------------------------------------------------------------------------
      
      SUBROUTINE loadNMatrixFromFile (filename,repeat,n,m,matrix,res)
      !Parameters 
      CHARACTER(*),INTENT(IN)                                  :: filename
      INTEGER,INTENT(IN)                                      :: m
      INTEGER,INTENT(IN)                                      :: n
      INTEGER,INTENT(IN)                                      :: repeat
      REAL(MK),POINTER,DIMENSION(:,:,:)                       :: matrix
      REAL(MK),POINTER,DIMENSION(:),OPTIONAL                  :: res
      
      
      integer :: io_error
      
      integer :: i,j,k
      
    
      open(unit=20,file=filename,status='old',action='read',&
     &           iostat=io_error) 

      IF ( io_error == 0) then
        DO j = 1,repeat
         DO i = 1,n 
            read(20,*) matrix(j,i,:)
            !write(*,*) vector(1:m,n)
         END DO
        END DO
         IF (PRESENT(res)) THEN   
            DO i=1,n
                read(20,*) res(i)
                !write(*,*) res(i)
            END DO
         END IF
      else
           IF (io_error > 0) THEN          
               write(*,*) 'Error while opening File ', filename, ' - Error Nr:',& 
                           io_error,' occured in Function loadNMatrixFrom File'   
               STOP
           END IF
      end IF
      close(unit=20) 
      END SUBROUTINE loadNMatrixFromFile
      
      
      !-------------------------------------------------------------------------
      ! loadMatrixFromFile
      !-------------------------------------------------------------------------
      !  Purpose      : loads a matrix from a file (space / linebreak seperated)
      !                 optional loads a column vector with results appending the matrix      
      ! 
      !  Input       : filename (S) TRIM(adjustl(path)) to the File from which to read a Row Vector
      !                 m       (I) m Dimension of the matrix
      !                 n       (I) n Dimension of the matrix
      !                 
      !  Input/Output:  matrix(P,R) Pointer to the resulting matrix
      !                 res   (P,R) Pointer to a vector with corresponding results
      !                             OPTIONAL Parameter - only use if theres results also in the
      !                                                  matrixfile - dimensions are similar to matrix length m
      !-------------------------------------------------------------------------
      
      
      SUBROUTINE loadMatrixFromFile (filename,n,m,matrix,res)
      !Parameters 
      CHARACTER(*),INTENT(IN)                                 :: filename
      INTEGER,INTENT(IN)                                      :: m
      INTEGER,INTENT(IN)                                      :: n
      REAL(MK),POINTER,DIMENSION(:,:)                         :: matrix
      REAL(MK),POINTER,DIMENSION(:),OPTIONAL                  :: res
      
      

      integer :: io_error
      
      integer :: i,j,k
      
    
      open(unit=20,file=filename,status='old',action='read',&
     &           iostat=io_error) 

      IF ( io_error == 0) then
         DO i = 1,n 
            read(20,*) matrix(i,:)
            !write(*,*) vector(1:m,n)
         END DO
         
         IF (PRESENT(res)) THEN   
            DO i=1,n
                read(20,*) res(i)
                !write(*,*) res(i)
            END DO
         END IF
      else
           write(*,*) 'Error while opening File ', filename, ' - Error Nr:',& 
                       io_error,' occured' 
           STOP
      end IF
      close(unit=20) 
      END SUBROUTINE loadMatrixFromFile
          
      
  
      

      
      

      
      
      SUBROUTINE hybrid_composition(vars,res,m,hybrid_nr,lbounds,ubounds)
      INTEGER,INTENT(in)                               :: m
      REAL(MK),INTENT(out)                             :: res
      REAL(MK),DIMENSION(m),INTENT(INOUT)              :: vars
      INTEGER,INTENT(in)                               :: hybrid_nr
      REAL(MK),DIMENSION(m),OPTIONAL                   :: lbounds
      REAL(MK),DIMENSION(m),OPTIONAL                   :: ubounds  
      DOUBLE PRECISION                                 :: ZBQLNOR !the gaus gen

      
      !local used vars
      REAL(MK)                              :: wMax,sumSqr,wSum,w1mMaxPow
      INTEGER                               :: i,j,k
      REAL(MK)                              :: sumF,t_res
      REAL(MK),DIMENSION(m)                        :: ub,lb
      LOGICAL                      :: boundviolation
      
              
      res = 0.0_MK

           


       
       
       
      IF (.NOT. PRESENT(lbounds)) THEN
            lb = -5_MK
      ELSE
            lb = lbounds
      END IF
       
      IF (.NOT. PRESENT(ubounds)) THEN
        ub = 5_MK
      ELSE
        ub = ubounds
      END IF
       
      boundviolation = .FALSE.
      DO j=1,m
        IF (vars(j) .GT. ub(j) .OR. vars(j) .LT. lb(j)) THEN
            boundviolation = .TRUE.
        END IF            
      END DO
      
      IF (boundviolation ) THEN
        res = testf_NAN
      ELSE 
      
          
          IF (hybrid_nr .EQ. 23) THEN
           DO j = 1,m
               !myXround
               IF (ABS(vars(j) - job_o(1,j)) .GT. 0.5_MK) THEN
                      vars(j) = ANINT (vars(j) *2.0_MK) /2.0_MK
               END IF
           END DO
          
          
          END IF
          
          
          
          
          !get the raw weights
          wMax = - huge(wMax)
          DO i = 1,job_num_func
            sumSqr = 0.0_MK
            !Shift the Input
            DO j=1,m
               job_z(i,j) = vars(j) - job_o(i,j)
            END DO
            DO j=1,m
                sumSqr = sumSqr + (job_z(i,j) * job_z(i,j))
            END DO
            IF (hybrid_nr .EQ. 15 .OR. hybrid_nr .EQ. 16 .OR. &
     &                hybrid_nr .EQ. 17) THEN
                job_w(i) = exp(-1.0_MK * sumSqr / &
     &                          (2.0_MK * m * job_sigma(i)*job_sigma(i)))
            ELSE
                IF (hybrid_nr .EQ. 18 .OR. hybrid_nr .EQ. 20) THEN
                    job_w(i) = exp(-1.0_MK * sumSqr / (2.0_MK * m * &
     &                               job2_sigma(i)*job2_sigma(i)))
                ELSE
                    IF (hybrid_nr .EQ. 19) THEN
                        job_w(i) = exp(-1.0_MK * sumSqr / (2.0_MK * m * &
     &                                   job19_sigma(i)*job19_sigma(i)))
                    ELSE
                        IF (hybrid_nr .EQ. 21 .OR. hybrid_nr .EQ. 22 .OR. &
     &                            hybrid_nr .EQ. 23) THEN
                            job_w(i) = exp(-1.0_MK * sumSqr / (2.0_MK * m * &
     &                                       job21_sigma(i)*job21_sigma(i)))
                        ELSE
                            job_w(i) = exp(-1.0_MK * sumSqr / (2.0_MK * m * &
     &                                       job24_sigma(i)*job24_sigma(i)))
                        ENDIF
                    END IF
                END IF
            END IF
            IF (wMax < job_w(i)) THEN
                wMax = job_w(i)
            END IF  
          END DO
          
          !Modify the weights
          wSum = 0.0_MK
          w1mMaxPow = 1.0_MK - wMax**10
          DO i = 1,job_num_func
            IF (job_w(i) .NE. wMax) THEN
                job_w(i) = job_w(i)* w1mMaxPow
            END IF
            wSum = wSum + job_w(i)
          END DO
          
          
          !Normalize the weights
          DO i = 1,job_num_func
            job_w(i) = job_w(i) / wSum
          END DO    
          
          sumF = 0.0_MK
          DO i = 1,job_num_func
            DO j = 1,m
                IF (hybrid_nr .EQ. 15 .OR. hybrid_nr .EQ. 16 .OR. hybrid_nr &
     &                    .EQ. 17 ) THEN
                    job_z(i,j) = job_z(i,j) / job_lambda(i)
                ELSE
                    IF (hybrid_nr .EQ. 18 .OR. hybrid_nr .EQ. 20) THEN
                        job_z(i,j) = job_z(i,j) / job2_lambda(i)
                    ELSE
                        IF (hybrid_nr .EQ. 19) THEN
                            job_z(i,j) = job_z(i,j) / job19_lambda(i)
                        ELSE
                            IF ( hybrid_nr .EQ. 21 .OR. hybrid_nr .EQ. 22 .OR.&
     &                                 hybrid_nr .EQ. 23) THEN
                                job_z(i,j) = job_z(i,j) / job21_lambda(i)
                            ELSE
                                job_z(i,j) = job_z(i,j) / job24_lambda(i)
                            END IF
                        END IF
                    END IF
                END IF
            END DO
            !rotate 
            CALL xA(m,job_zM(i,:),job_z(i,:),job_M(i,:,:))
            
            !calling the basic functions
            IF (hybrid_nr .EQ. 15 .OR. hybrid_nr .EQ. 16 .OR. hybrid_nr &
     &                .EQ. 17) THEN
                CALL Job15(i,m,job_zM(i,:),t_res)
            ELSE
                IF (hybrid_nr .LT. 21) THEN
                    CALL Job18(i,m,job_zM(i,:),t_res)
                ELSE
                    IF (hybrid_nr .LT. 24) THEN
                        CALL Job21(i,m,job_zM(i,:),t_res)
                    ELSE
                        CALL Job24(i,m,job_zM(i,:),t_res)
                    ENDIF
                ENDIF
            END IF
            
            
            sumF = sumF + job_w(i) * (job_C*t_res/job_fmax(i) + job_biases(i) )
            
          END DO
          res = sumF + biases(hybrid_nr)
          
          
          IF (hybrid_nr .EQ. 17) THEN
            IF (noise ) THEN
                res = res * (1.0_MK + 0.2_MK * ABS(ZBQLNOR(0.0_MK,1.0_MK)))! * THE NOISE FUNCTION HERE!)
            END IF
          END IF
      END IF    
        
      
      
      END SUBROUTINE hybrid_composition
      
      
      
      SUBROUTINE Job24(funcnr,m,vars,res)
      !Parameters
      INTEGER,INTENT(IN)                    :: funcnr,m
      REAL(MK),INTENT(out)                   :: res
      REAL(MK),DIMENSION(m),INTENT(in)       :: vars
      
      !local used vars
      REAL(MK)                               :: sum,sum1,sum2,prod,e1,e2,t1,t2&
     &                                                ,temp1,temp2,temp3,a1
      INTEGER                                :: i,j,k
      !weierstrass vars
      INTEGER                                :: Kmax = 20
      REAL(MK)                               :: a = 0.5_MK
      REAL(MK)                               :: b = 3.0_MK
      REAL(MK)                               :: prevX,currX 
      !copy for non cont
      REAL(MK),DIMENSION(m)                  :: vars_t
      REAL(MK)                               :: ZBQLNOR !the gaus gen
      
      
      IF (funcnr .EQ. 1) THEN
         !weierstrass
         sum1 = 0.0_MK
         sum2 = 0.0_MK
         DO i = 1,m
            DO k = 0,Kmax
                sum1 = sum1 + a**k * cos((PI)*2.0_MK *b**k *(vars(i) + 0.5_MK))            
            END DO        
         END DO
         DO k = 0,Kmax
            sum2 = sum2 + a**k * cos((PI)*2.0_MK * b**k * (0.5_MK))
         END DO
         res = sum1 - sum2*m    

      ENDIF    
      IF (funcnr .EQ. 2) THEN
        !EScafferF6
        sum = 0.0_MK 
            DO i = 1,m-1
                temp1 = vars(i)**2 + vars(i+1)**2
                temp2 = sin(sqrt(temp1))
                temp3 = 1.0_MK + 0.001_MK * temp1
                sum = sum + (0.5_MK + ((temp2**2 - 0.5_MK)/(temp3**2)))
            END DO     
            temp1 = vars(m)**2 + vars(1)**2
            temp2 = sin(sqrt(temp1))
            temp3 = 1.0_MK + 0.001_MK * temp1
            sum = sum + (0.5_MK + ((temp2**2 - 0.5_MK)/(temp3**2)))
            res = sum
      ENDIF
      IF (funcnr .EQ. 3) THEN
      !F8F2
         sum = 0.0_MK
         DO i = 1,m-1
            !F2 part
            t1 = vars(i)**2 - vars(i+1)
            t2 = vars(i) - 1.0_MK
            t1 = ((100.0_MK * t1 * t1) + (t2*t2))
            !F8 part
            t2 = ((t1**2) / 4000.0_MK) - cos(t1) + 1.0_MK
            sum = sum + t2
         END DO
         
        !F2 part
        t1 = vars(m)**2 - vars(1)
        t2 = vars(m) - 1.0_MK
        t1 = ((100.0_MK * t1 * t1) + (t2*t2))
        !F8 part
        t2 = ((t1**2) / 4000.0_MK) - cos(t1) + 1.0_MK
        sum = sum + t2  
        res = sum
      ENDIF
      IF (funcnr .EQ. 4) THEN
      !ackley
        e1 = 0.0d0
        e2 = 0.0d0
        DO i = 1,m
            e1 = e1 + vars(i)**2
            e2 = e2 + cos(2.0d0*Pi*vars(i))
        END DO
        res = exp(1.0d0) + 20.0d0 - 20d0*exp(-0.2d0*sqrt(e1/real(m)))  
        res = res - exp(e2/real(m)) 
      ENDIF
      IF (funcnr .EQ. 5) THEN
        !rastrigin      
         sum = 10.0d0 * m
         DO i = 1,m
            sum=sum+vars(i)**2
            sum=sum-10.0d0*cos(2*Pi*vars(i))                      
         END DO
         res = sum
      
      END IF
      IF (funcnr .EQ. 6) THEN
        !griewank
         prod=1
         sum=0.0_MK
         DO i = 1,m
           sum= sum+(vars(i)**2)/4000
           prod=prod* (cos(vars(i)/(i**(1.0/2))))        
         END DO
         res = sum-prod+1 
      ENDIF
      IF (funcnr .EQ. 7) THEN
      !EScafferF6NonCont
       DO j = 1,m
           !myXround
           IF (ABS(vars(j)) .GT. 0.5_MK) THEN
                  vars_t(j) = ANINT (vars(j) *2.0_MK) /2.0_MK
           ELSE
                  vars_t(j) = vars(j)
           END IF
       END DO     
            sum = 0.0_MK 
            DO i = 1,m-1
                temp1 = vars_t(i)**2 + vars_t(i+1)**2
                temp2 = sin(sqrt(temp1))
                temp3 = 1.0_MK + 0.001_MK * temp1
                sum = sum + (0.5_MK + ((temp2**2 - 0.5_MK)/(temp3**2)))
            END DO     
            temp1 = vars_t(m)**2 + vars_t(1)**2
            temp2 = sin(sqrt(temp1))
            temp3 = 1.0_MK + 0.001_MK * temp1
            sum = sum + (0.5_MK + ((temp2**2 - 0.5_MK)/(temp3**2)))
            res = sum
      ENDIF 
      IF(funcnr .EQ. 8) THEN
      !rastriginNonCont
       DO j = 1,m
           !myXround
           IF (ABS(vars(j)) .GT. 0.5_MK) THEN
                  vars_t(j) = ANINT (vars(j) *2.0_MK) /2.0_MK
           ELSE
                  vars_t(j) = vars(j)
           END IF
       END DO 
      
         sum = 10.0d0 * m
         DO i = 1,m
            sum=sum+vars_t(i)**2
            sum=sum-10.0d0*cos(2*Pi*vars_t(i))                      
         END DO
         res = sum
      
      ENDIF
      IF (funcnr .EQ. 9) THEN
      !elliptic
      sum = 0.0_MK
      a1 = 1e6
      
      do i = 1,m
        sum = sum + a1**((i-1)/((m-1)*1.0_MK)) * vars(i)**2
      END DO
      res = sum
    
      ENDIF
      IF (funcnr .EQ. 10) THEN
      !sphere with noise
        sum = 0.0_MK
        DO i = 1,m
            sum = sum + vars(i)**2
        END DO
        IF (noise ) THEN
            sum = sum * (1.0_MK + 0.1 * ABS(ZBQLNOR(0.0_MK,1.0_MK)))
        END IF      
        res = sum    
      ENDIF
      
      res = res
      END SUBROUTINE Job24
      
      
      
      
      
      SUBROUTINE Job21(funcnr,m,vars,res)
      !Parameters
      INTEGER,INTENT(IN)                    :: funcnr,m
      REAL(MK),INTENT(out)                   :: res
      REAL(MK),DIMENSION(m),INTENT(in)       :: vars
      
      !local used vars
      REAL(MK)                               :: sum,sum1,sum2,prod,e1,e2,t1,t2&
     &                                                ,temp1,temp2,temp3
      INTEGER                                :: i,j,k
      !weierstrass vars
      INTEGER                                :: Kmax = 20
      REAL(MK)                               :: a = 0.5_MK
      REAL(MK)                               :: b = 3.0_MK
      
      
      
      IF (funcnr .LT. 3) THEN
        !EScafferF6
        sum = 0.0_MK 
            DO i = 1,m-1
                temp1 = vars(i)**2 + vars(i+1)**2
                temp2 = sin(sqrt(temp1))
                temp3 = 1.0_MK + 0.001_MK * temp1
                sum = sum + (0.5_MK + ((temp2**2 - 0.5_MK)/(temp3**2)))
            END DO     
            temp1 = vars(m)**2 + vars(1)**2
            temp2 = sin(sqrt(temp1))
            temp3 = 1.0_MK + 0.001_MK * temp1
            sum = sum + (0.5_MK + ((temp2**2 - 0.5_MK)/(temp3**2)))
            res = sum
      ENDIF
      IF (funcnr .EQ. 3 .OR. funcnr .EQ. 4) THEN
        !rastrigin      
         sum = 10.0d0 * m
         DO i = 1,m
            sum=sum+vars(i)**2
            sum=sum-10.0d0*cos(2*Pi*vars(i))                      
         END DO
         res = sum
      ENDIF
      IF (funcnr .EQ. 5 .OR. funcnr .EQ. 6) THEN
         !F8F2
         sum = 0.0_MK
         DO i = 1,m-1
            !F2 part
            t1 = vars(i)**2 - vars(i+1)
            t2 = vars(i) - 1.0_MK
            t1 = ((100.0_MK * t1 * t1) + (t2*t2))
            !F8 part
            t2 = ((t1**2) / 4000.0_MK) - cos(t1) + 1.0_MK
            sum = sum + t2
         END DO
         
        !F2 part
        t1 = vars(m)**2 - vars(1)
        t2 = vars(m) - 1.0_MK
        t1 = ((100.0_MK * t1 * t1) + (t2*t2))
        !F8 part
        t2 = ((t1**2) / 4000.0_MK) - cos(t1) + 1.0_MK
        sum = sum + t2  
        res = sum
      ENDIF
      IF (funcnr .EQ. 7 .OR. funcnr .EQ. 8) THEN
         !weierstrass
         sum1 = 0.0_MK
         sum2 = 0.0_MK
         DO i = 1,m
            DO k = 0,Kmax
                sum1 = sum1 + a**k * cos((PI)*2.0_MK *b**k *(vars(i) + 0.5_MK))            
            END DO        
         END DO
         DO k = 0,Kmax
            sum2 = sum2 + a**k * cos((PI)*2.0_MK * b**k * (0.5_MK))
         END DO
         res = sum1 - sum2*m    

      ENDIF
      IF (funcnr .EQ. 9 .OR. funcnr .EQ. 10) THEN
         !griewank
         prod=1
         sum=0.0_MK
         DO i = 1,m
           sum= sum+(vars(i)**2)/4000
           prod=prod* (cos(vars(i)/(i**(1.0/2))))        
         END DO
         res = sum-prod+1
      ENDIF
      
      
      END SUBROUTINE Job21
      
      
      
      SUBROUTINE Job18(funcnr,m,vars,res)
      !Parameters
      INTEGER,INTENT(IN)                    :: funcnr,m
      REAL(MK),INTENT(out)                   :: res
      REAL(MK),DIMENSION(m),INTENT(in)       :: vars
      
      !local used vars
      REAL(MK)                               :: sum,sum1,sum2,prod,e1,e2
      INTEGER                                :: i,j,k
      !weierstrass vars
      INTEGER                                :: Kmax = 20
      REAL(MK)                               :: a = 0.5_MK
      REAL(MK)                               :: b = 3.0_MK
      
      
      IF (funcnr .LT. 3) THEN
        !ackley
        e1 = 0.0d0
        e2 = 0.0d0
        DO i = 1,m
            e1 = e1 + vars(i)**2
            e2 = e2 + cos(2.0d0*Pi*vars(i))
        END DO
        res = exp(1.0d0) + 20.0d0 - 20d0*exp(-0.2d0*sqrt(e1/real(m)))  
        res = res - exp(e2/real(m))

      ENDIF
      IF (funcnr .EQ. 3 .OR. funcnr .EQ. 4) THEN
        !rastrigin      
         sum = 10.0d0 * m
         DO i = 1,m
            sum=sum+vars(i)**2
            sum=sum-10.0d0*cos(2*Pi*vars(i))                   
            res = sum
         END DO
      ENDIF
      IF (funcnr .EQ. 5 .OR. funcnr .EQ. 6) THEN
        !sphere
        sum = 0.0_MK
        DO i = 1,m
            sum = sum + vars(i)**2
        END DO
        res = sum    
       
      ENDIF
      IF (funcnr .EQ. 7 .OR. funcnr .EQ. 8) THEN
         !weierstrass
         sum1 = 0.0_MK
         sum2 = 0.0_MK
         DO i = 1,m
            DO k = 0,Kmax
                sum1 = sum1 + a**k * cos((PI)*2.0_MK * b**k *(vars(i)+0.5_MK))            
            END DO        
         END DO
         DO k = 0,Kmax
            sum2 = sum2 + a**k * cos((PI)*2.0_MK * b**k * (0.5_MK))
         END DO
         res = sum1 - sum2*m    

      ENDIF
      IF (funcnr .EQ. 9 .OR. funcnr .EQ. 10) THEN
         !griewank
         prod=1
         sum=0.0_MK
         DO i = 1,m
           sum= sum+(vars(i)**2)/4000
           prod=prod* (cos(vars(i)/(i**(1.0/2))))        
         END DO
         res = sum-prod+1
      ENDIF
      
      
      END SUBROUTINE Job18
      
      
      
      
      
      SUBROUTINE Job15(funcnr,m,vars,res)
      !Parameters
      INTEGER,INTENT(IN)                    :: funcnr,m
      REAL(MK),INTENT(out)                   :: res
      REAL(MK),DIMENSION(m),INTENT(in)       :: vars
      
      !local used vars
      REAL(MK)                               :: sum,sum1,sum2,prod,e1,e2
      INTEGER                                :: i,j,k
      !weierstrass vars
      INTEGER                                :: Kmax = 20
      REAL(MK)                               :: a = 0.5_MK
      REAL(MK)                               :: b = 3.0_MK
      
      
      IF (funcnr .LT. 3) THEN
        !rastrigin      
         sum = 10.0d0 * m
         DO i = 1,m
            sum=sum+vars(i)**2
            sum=sum-10.0d0*cos(2*Pi*vars(i))                   
            res = sum
         END DO
      ENDIF
      IF (funcnr .EQ. 3 .OR. funcnr .EQ. 4) THEN
        !weierstrass
         sum1 = 0.0_MK
         sum2 = 0.0_MK
         DO i = 1,m
            DO k = 0,Kmax
                sum1 = sum1 + a**k * cos((PI)*2.0_MK * b**k * (vars(i)+0.5_MK))            
            END DO        
         END DO
         DO k = 0,Kmax
            sum2 = sum2 + a**k * cos((PI)*2.0_MK * b**k * (0.5_MK))
         END DO
         res = sum1 - sum2*m    
      ENDIF
      IF (funcnr .EQ. 5 .OR. funcnr .EQ. 6) THEN
        !griewank
         prod=1
         sum=0.0_MK
         DO i = 1,m
           sum= sum+(vars(i)**2)/4000
           prod=prod* (cos(vars(i)/(i**(1.0/2))))        
         END DO
         res = sum-prod+1
      ENDIF
      IF (funcnr .EQ. 7 .OR. funcnr .EQ. 8) THEN
        !ackley
        e1 = 0.0d0
        e2 = 0.0d0
        DO i = 1,m
            e1 = e1 + vars(i)**2
            e2 = e2 + cos(2.0d0*Pi*vars(i))
        END DO
        res = exp(1.0d0) + 20.0d0 - 20d0*exp(-0.2d0*sqrt(e1/real(m)))  
        res = res - exp(e2/real(m))
      ENDIF
      IF (funcnr .EQ. 9 .OR. funcnr .EQ. 10) THEN
        !sphere
        sum = 0.0_MK
        DO i = 1,m
            sum = sum + vars(i)**2
        END DO
        res = sum
      ENDIF
      
      
      END SUBROUTINE Job15
      
      !-------------------------------------------------------------------------
      ! xA
      !
      
      SUBROUTINE xA (m,res,row,matrix)
       !Parameters
       INTEGER,INTENT(in)                            :: m
       REAL(MK),DIMENSION(m),INTENT(out)           :: res
       REAL(MK),DIMENSION(m,m),INTENT(in)            :: matrix
       REAL(MK),DIMENSION(m),INTENT(in)            :: row
       
       INTEGER                                        :: i,j
       
       DO i = 1,m
        res(i) = 0
        DO j= 1,m
            res(i) = res(i) + (row(j) * matrix(j,i))
        END DO
       END DO
      
      END SUBROUTINE xA
      

      !-------------------------------------------------------------------------
      ! Ax
      !
      
      SUBROUTINE Ax (m,res,row,matrix)
       !Parameters
       INTEGER,INTENT(in)                            :: m
       REAL(MK),DIMENSION(m),INTENT(out)           :: res
       REAL(MK),DIMENSION(m,m),INTENT(in)            :: matrix
       REAL(MK),DIMENSION(m),INTENT(in)            :: row
       
       INTEGER                                        :: i,j
       
       DO i = 1,m
        res(i) = 0
        DO j= 1,m
            res(i) = res(i) + (row(j) * matrix(i,j))
        END DO
       END DO
      
      END SUBROUTINE Ax
      
      !-------------------------------------------------------------------------
      ! F01_shifted_sphere
      !-------------------------------------------------------------------------
      ! References   :     De Jong,K.: An Analysis of the
      !                    Behaviour of a Class of Genetic Adaptive Systems. PhD
      !                    thesis,University of Michigan (1975)
      !
      !  Input        : vars    (R) a matrix containing an array of vectors
      !                             the function should process
      !                             vectors are saved per column at the
      !                             moment cause for a big problem size           
      !                             the population is smaller then the
      !                             dimension
      !                 m       (I) m Dimension of the matrix (Rows)
      !                               = Dimension of Vectors
      !                 n       (I) n Dimension of the matrix (Columns)
      !                               = # of Vectors to process
      !                 bias    (R) a bias wich is added to the result              
      !                 lbounds (R) lower bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result
      !                 ubounds (R) upper bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result      
      !
      !
      !  Output       : res     (R) a row vector with the result of each
      !                             column vector of the vars input,has
      !                             size n
      !-------------------------------------------------------------------------
      SUBROUTINE F01_shifted_sphere(res,vars,m,n,lbounds,ubounds) 
       
      !Parameters
      INTEGER,INTENT(in)                           :: m
      INTEGER,INTENT(in)                           :: n
      REAL(MK),DIMENSION(n),INTENT(out)      :: res
      REAL(MK),DIMENSION(m,n),INTENT(inout)        :: vars
          
      REAL(MK),DIMENSION(m),OPTIONAL              :: lbounds
      REAL(MK),DIMENSION(m),OPTIONAL              :: ubounds
           
      ! Locally used variables
      REAL(MK),DIMENSION(m)                        :: ub,lb
      INTEGER               :: var_size
      INTEGER                                       :: i  
      INTEGER                                       :: j
      LOGICAL,DIMENSION(n)                  :: boundviolation
       
       
       
      IF (.NOT. PRESENT(lbounds)) THEN
        lb = -100_MK
      ELSE
        lb = lbounds
      END IF
       
      IF (.NOT. PRESENT(ubounds)) THEN
        ub = 100_MK
      ELSE
        ub = ubounds
      END IF
       
      boundviolation = .FALSE.
       
       
      !Shift the Input
      DO j=1,n
        DO i=1,m
            IF (vars(i,j) .GT. ub(i) .OR. vars(i,j) .LT. lb(i)) THEN
                boundviolation(j) = .TRUE.   
            ELSE
                 vars(i,j) = vars(i,j) - shiftmatrix(i)
            END IF
        END DO
      END DO
       
      DO j = 1,n
               
        res(j) = 0.0_MK
        IF (boundviolation(j) ) THEN
            res(j) = testf_NAN
        ELSE
            DO i = 1,m      
                res(j) = res(j) + vars(i,j) * vars(i,j)
            END DO
            res(j) = res(j) + biases(1)  
        END IF
       END DO
      
        
       
      END SUBROUTINE F01_shifted_sphere
      
      !-------------------------------------------------------------------------
      ! Schwefel's Double Sum Function (f2)
      !-------------------------------------------------------------------------
      !  References   :     H. P. Schwefel.
      !                    Evolution and Optimum Seeking.
      !                    John Wiley & Sons,1995.
      !
      !  Input        : vars    (R) a matrix containing an array of vectors
      !                             the function should process
      !                             vectors are saved per column at the
      !                             moment cause for a big problem size           
      !                             the population is smaller then the
      !                             dimension
      !                 m       (I) m Dimension of the matrix
      !                               = Dimension of Vectors
      !                 n       (I) n Dimension of the matrix
      !                               = # of Vectors to process 
      !                 bias    (R) a bias wich is added to the result                           
      !                 lbounds (R) lower bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result
      !                 ubounds (R) upper bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result      
      !
      !
      !  Output       : res     (R) a row vector with the result of each
      !                             column vector of the vars input,has
      !                             size n      
      !-------------------------------------------------------------------------
      SUBROUTINE F02_shifted_schwefel (res,vars,m,n,lbounds,ubounds)   
      
      INTEGER,INTENT(in)                       :: m
      INTEGER,INTENT(in)                       :: n
      REAL(MK),DIMENSION(n),INTENT(out)     :: res
      REAL(MK),DIMENSION(m,n),INTENT(inout)    :: vars
      REAL(MK),DIMENSION(m),OPTIONAL          :: lbounds
      REAL(MK),DIMENSION(m),OPTIONAL          :: ubounds
           
      ! Locally used variables
      REAL(MK),DIMENSION(m)                    :: ub,lb
      INTEGER                         :: h
      INTEGER                              :: i
      INTEGER                              :: j
      LOGICAL,DIMENSION(n)                   :: boundviolation
       
      IF (.NOT. PRESENT(lbounds)) THEN
       lb = -100_MK
      ELSE
       lb = lbounds
      END IF
       
      IF (.NOT. PRESENT(ubounds)) THEN
       ub = 100_MK
      ELSE
       ub = ubounds
      END IF
      
      boundviolation = .FALSE.
            
       
      !Shift the Input
      DO j=1,n
        DO i=1,m
            IF (vars(i,j) .GT. ub(i) .OR. vars(i,j) .LT. lb(i)) THEN
                boundviolation(j) = .TRUE.   
            ELSE
                 vars(i,j) = vars(i,j) - shiftmatrix(i)
            END IF
        END DO
      END DO
 
      res = 0.0_MK
      DO j = 1,n
      
        IF (boundviolation(j) ) THEN
            res(j) = testf_NAN
        ELSE
            DO i = 1,m
                res(j) = res(j) + sum(vars(1:i,j))**2
            END DO
            res(j) = res(j) + biases(2)
        END IF
      END DO    
      
       RETURN
      END SUBROUTINE F02_shifted_schwefel

      !-------------------------------------------------------------------------
      ! F03_shifted_rotated_high_cond_elliptic
      !-------------------------------------------------------------------------
      ! References   :     
      !
      !  Input        : vars    (R) a matrix containing an array of vectors
      !                             the function should process
      !                             vectors are saved per column at the
      !                             moment cause for a big problem size           
      !                             the population is smaller then the
      !                             dimension
      !                 m       (I) m Dimension of the matrix (Rows)
      !                               = Dimension of Vectors
      !                 n       (I) n Dimension of the matrix (Columns)
      !                               = # of Vectors to process
      !                 bias    (R) a bias wich is added to the result              
      !                 lbounds (R) lower bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result
      !                 ubounds (R) upper bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result      
      !
      !
      !  Output       : res     (R) a row vector with the result of each
      !                             column vector of the vars input,has
      !                             size n
      !-------------------------------------------------------------------------
      SUBROUTINE F03_shifted_rotated_high_cond_elliptic(res,vars,m,n,&
     &                                                        lbounds,ubounds) 
       
      !Parameters
      INTEGER,INTENT(in)                            :: m
      INTEGER,INTENT(in)                            :: n
      REAL(MK),DIMENSION(n),INTENT(out)           :: res
      REAL(MK),DIMENSION(m,n),INTENT(inout)         :: vars
       
       
      REAL(MK),DIMENSION(m),OPTIONAL               :: lbounds
      REAL(MK),DIMENSION(m),OPTIONAL               :: ubounds
       
      ! Locally used variables
      REAL(MK),DIMENSION(m)                         :: ub,lb
      INTEGER             :: var_size
      INTEGER                                        :: i  
      INTEGER                                        :: j,k
      REAL(MK)                                       :: constant
      REAL(MK),DIMENSION(m)                          :: tmp_result
      LOGICAL,DIMENSION(n)               :: boundviolation
      
       
      IF (.NOT. PRESENT(lbounds)) THEN
        lb = -100_MK
      ELSE
        lb = lbounds
      END IF
       
      IF (.NOT. PRESENT(ubounds)) THEN
        ub = 100_MK
      ELSE
        ub = ubounds
      END IF
      
      boundviolation = .FALSE.
      
       
      !Shift the Input
      DO j=1,n
        DO i=1,m
            IF (vars(i,j) .GT. ub(i) .OR. vars(i,j) .LT. lb(i)) THEN
                boundviolation(j) = .TRUE.   
            ELSE
                vars(i,j) = vars(i,j) - shiftmatrix(i)
            END IF
        END DO
      END DO
      
      !rotate the input
        
      DO i=1,n
        IF (.NOT. boundviolation(i)) THEN
            DO j=1,m
                tmp_result(j) = 0.0_MK
                DO k=1,m
                    tmp_result(j) = tmp_result(j)+ vars(k,i) *rotatematrix(k,j)
                END DO
            END DO
            vars(:,i) = tmp_result
        END IF          
        !vars(:,i) = matmul(tmp_vec,tmp)
      END DO

      constant = 1.0e6_MK**(1.0_MK/(m-1.0_MK))
              
      DO j = 1,n          
        res(j) = 0.0_MK
        IF (boundviolation(j) ) THEN
            res(j) = testf_NAN
        ELSE
            DO i = 1,m      
                res(j) = res(j) + constant**(i-1) * vars(i,j) * vars(i,j)
            END DO
            res(j) = res(j) + biases(3)
        END IF
      END DO
       
       
      END SUBROUTINE F03_shifted_rotated_high_cond_elliptic
      
      !-------------------------------------------------------------------------
      ! F04_shifted_schwefel_noise
      !-------------------------------------------------------------------------
      ! References   :     
      !
      !  Input        : vars    (R) a matrix containing an array of vectors
      !                             the function should process
      !                             vectors are saved per column at the
      !                             moment cause for a big problem size           
      !                             the population is smaller then the
      !                             dimension
      !                 m       (I) m Dimension of the matrix (Rows)
      !                               = Dimension of Vectors
      !                 n       (I) n Dimension of the matrix (Columns)
      !                               = # of Vectors to process
      !                 bias    (R) a bias wich is added to the result              
      !                 lbounds (R) lower bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result
      !                 ubounds (R) upper bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result      
      !
      !
      !  Output       : res     (R) a row vector with the result of each
      !                             column vector of the vars input,has
      !                             size n
      !-------------------------------------------------------------------------

      SUBROUTINE F04_shifted_schwefel_noise(res,vars,m,n,lbounds,ubounds) 
      !Parameters
      INTEGER,INTENT(in)                            :: m
      INTEGER,INTENT(in)                            :: n
      REAL(MK),DIMENSION(n),INTENT(out)           :: res
      REAL(MK),DIMENSION(m,n),INTENT(inout)         :: vars
       
       
      REAL(MK),DIMENSION(m),OPTIONAL               :: lbounds
      REAL(MK),DIMENSION(m),OPTIONAL               :: ubounds
       
      ! Locally used variables
      INTEGER                                        :: i
      REAL(MK)                                       :: ZBQLNOR !the gaus gen
      
      CALL F02_shifted_schwefel(res,vars,m,n,lbounds,ubounds)
      
      
      ! add noise function here
      DO i = 1,n
        IF (res(i) .NE. testf_NAN) THEN
            res(i) = (res(i) - biases(2))
            IF (noise ) THEN
                res(i) = res(i)*(1.0_MK + 0.4_MK *ABS(ZBQLNOR(0.0_MK,1.0_MK)))! * THE NOISE FUNCTION HERE!)
            END IF
            res(i) = res(i) + biases(4)
        END IF
      END DO
            
      END SUBROUTINE F04_shifted_schwefel_noise
      !-------------------------------------------------------------------------
      ! F05_schwefel_global_opt_bound
      !-------------------------------------------------------------------------
      ! References   :     
      !  Input        : vars    (R) a matrix containing an array of vectors
      !                             the function should process
      !                             vectors are saved per column at the
      !                             moment cause for a big problem size           
      !                             the population is smaller then the
      !                             dimension
      !                 m       (I) m Dimension of the matrix (Rows)
      !                               = Dimension of Vectors
      !                 n       (I) n Dimension of the matrix (Columns)
      !                               = # of Vectors to process
      !                 bias    (R) a bias wich is added to the result              
      !                 lbounds (R) lower bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result
      !                 ubounds (R) upper bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result      
      !
      !
      !  Output       : res     (R) a row vector with the result of each
      !                             column vector of the vars input,has
      !                             size n
      !-------------------------------------------------------------------------
      SUBROUTINE F05_schwefel_global_opt_bound(res,vars,m,n,lbounds,ubounds) 
      
      INTEGER,INTENT(in)                           :: m
      INTEGER,INTENT(in)                           :: n
      REAL(MK),DIMENSION(n),INTENT(out)            :: res
      REAL(MK),DIMENSION(m,n),INTENT(in)           :: vars
      REAL(MK),DIMENSION(m),OPTIONAL              :: lbounds
      REAL(MK),DIMENSION(m),OPTIONAL              :: ubounds
       
      ! Locally used variables
      REAL(MK),DIMENSION(m)                        :: ub,lb
      INTEGER                             :: h
      INTEGER                                  :: i
      INTEGER                                  :: j
      LOGICAL,DIMENSION(n)                       :: boundviolation
      REAL(MK),DIMENSION(m)                        :: tmparray
      REAL(MK),DIMENSION(m)                        :: tmparray2
      
      IF (.NOT. PRESENT(lbounds)) THEN
        lb = -100_MK
      ELSE
        lb = lbounds
      END IF
      
      IF (.NOT. PRESENT(ubounds)) THEN
        ub = 100_MK
      ELSE
        ub = ubounds
      END IF
      
      boundviolation = .FALSE.
          
      DO j=1,n
        DO i=1,m
            IF (vars(i,j) .GT. ub(i) .OR. vars(i,j) .LT. lb(i)) THEN
                boundviolation(j) = .TRUE.   
            END IF
        END DO
      END DO
      
      
      res = 0.0_MK
      DO i = 1,n
        IF (boundviolation(i) ) THEN
            res(i) = testf_NAN
        ELSE
            CALL Ax(m,tmparray,vars(:,i),rotatematrix)      
            !tmparray = matmul(vars(:,i),rotatematrix)
            DO j = 1,m
                tmparray2(j) = shiftmatrix(j)  
                tmparray2(j) = tmparray(j) - shiftmatrix(j)
                tmparray2(j) = abs(tmparray2(j))     
            END DO
            res(i) = maxval(tmparray2) + biases(5)
        END IF
      END DO    
      
       RETURN
              
      END SUBROUTINE F05_schwefel_global_opt_bound
      
      
      
      !-------------------------------------------------------------------------
      ! F06_shifted_rosenbrock
      !-------------------------------------------------------------------------
      ! References   :     
      !
      !  Input        : vars    (R) a matrix containing an array of vectors
      !                             the function should process
      !                             vectors are saved per column at the
      !                             moment cause for a big problem size           
      !                             the population is smaller then the
      !                             dimension
      !                 m       (I) m Dimension of the matrix (Rows)
      !                               = Dimension of Vectors
      !                 n       (I) n Dimension of the matrix (Columns)
      !                               = # of Vectors to process
      !                 bias    (R) a bias wich is added to the result              
      !                 lbounds (R) lower bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result
      !                 ubounds (R) upper bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result      
      !
      !
      !  Output       : res     (R) a row vector with the result of each
      !                             column vector of the vars input,has
      !                             size n
      !-------------------------------------------------------------------------
      SUBROUTINE F06_shifted_rosenbrock(res,vars,m,n,lbounds,ubounds) 
      
      INTEGER,INTENT(in)                            :: m
      INTEGER,INTENT(in)                            :: n
      REAL(MK),DIMENSION(n),INTENT(out)             :: res
      REAL(MK),DIMENSION(m,n),INTENT(inout)         :: vars
      REAL(MK),DIMENSION(m),OPTIONAL               :: lbounds
      REAL(MK),DIMENSION(m),OPTIONAL               :: ubounds
       
      ! Locally used variables
      REAL(MK),DIMENSION(m)                         :: ub,lb
      INTEGER                                        :: i  
      INTEGER                                        :: j
      LOGICAL,DIMENSION(n)                           :: boundviolation
       

      IF (.NOT. PRESENT(lbounds)) THEN
        lb = -100.0_MK
      ELSE
        lb = lbounds
      END IF
       
      IF (.NOT. PRESENT(ubounds)) THEN
        ub = 100.0_MK
      ELSE
        ub = ubounds
      END IF
      
      boundviolation = .FALSE.
      
       
      !Shift the Input
      DO j=1,n
        DO i=1,m
            IF (vars(i,j) .GT. ub(i) .OR. vars(i,j) .LT. lb(i)) THEN
                boundviolation(j) = .TRUE.   
            ELSE
                vars(i,j) = vars(i,j) - shiftmatrix(i)
            END IF
        END DO
      END DO
          
      res = 0.0_MK
      DO j = 1,n 
        IF (boundviolation(j) ) THEN
            res(j) = testf_NAN
        ELSE 
            DO i = 1,m - 1
            res(j) = res(j)+100.0_MK*(vars(i+1,j)-vars(i,j)**2)**2
            res(j) = res(j) + (vars(i,j) - 1)**2
            END DO
            res(j) = res(j) + biases(6)
        END IF
      END DO
      RETURN

      END SUBROUTINE F06_shifted_rosenbrock
      
      !-------------------------------------------------------------------------
      ! F07_shifted_rotated_griewank
      !-------------------------------------------------------------------------
      ! References   :     Griewangk,A.O.: Generalized Descent of Global Optimization. 
      !                    Journal of Optimization Theory and Applications,34: 11.39 (1981)
      !
      !  Input        : vars    (R) a matrix containing an array of vectors
      !                             the function should process
      !                             vectors are saved per column at the
      !                             moment cause for a big problem size           
      !                             the population is smaller then the
      !                             dimension
      !                 m       (I) m Dimension of the matrix
      !                               = Dimension of Vectors
      !                 n       (I) n Dimension of the matrix
      !                               = # of Vectors to process 
      !                 bias    (R) a bias wich is added to the result                                 
      !                 lbounds (R) lower bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result
      !                 ubounds (R) upper bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result      
      !
      !
      !  Output       : res     (R) a row vector with the result of each
      !                             column vector of the vars input,has
      !                             size n      
      !-------------------------------------------------------------------------
      SUBROUTINE F07_shifted_rotated_griewank (res,vars,m,n,lbounds,ubounds)   
      
      INTEGER,INTENT(in)                            :: m
      INTEGER,INTENT(in)                            :: n
      REAL(MK),DIMENSION(n),INTENT(out)             :: res
      REAL(MK),DIMENSION(m,n),INTENT(inout)         :: vars
           
      REAL(MK),DIMENSION(m),OPTIONAL               :: lbounds
      REAL(MK),DIMENSION(m),OPTIONAL               :: ubounds
       
      ! Locally used variables
      REAL(MK),DIMENSION(m)                         :: ub,lb
      INTEGER                                   :: i  
      REAL(MK)                                  :: sum
      REAL(MK)                                  :: prod
      INTEGER                                     :: j,k 
      REAL(MK),DIMENSION(m)                       :: tmp_result       
      LOGICAL,DIMENSION(n)                        :: boundviolation


      IF (.NOT. PRESENT(lbounds)) THEN
        lb = -HUGE(sum)
      ELSE
        lb = lbounds
      END IF
       
      IF (.NOT. PRESENT(ubounds)) THEN
        ub = HUGE(sum)
      ELSE
        ub = ubounds
      END IF 
      
      boundviolation = .FALSE.
      
       
      !Shift the Input
      DO j=1,n
        DO i=1,m
            IF (vars(i,j) .GT. ub(i) .OR. vars(i,j) .LT. lb(i)) THEN
                boundviolation(j) = .TRUE.   
            ELSE
                vars(i,j) = vars(i,j) - shiftmatrix(i)
            END IF
        END DO
      END DO
       
      !rotate the input
      DO i=1,n
        IF (.NOT. boundviolation(i)) THEN
            DO j=1,m
                tmp_result(j) = 0.0_MK
                DO k=1,m
                    tmp_result(j) = tmp_result(j)+ vars(k,i) *rotatematrix(k,j)
                END DO
            END DO
            vars(:,i) = tmp_result
        END IF
        !vars(:,i) = matmul(tmp_vec,tmp)
      END DO           
  
      DO j = 1,n
        prod=1
        sum=0
        IF (boundviolation(j) ) THEN
            res(j) = testf_NAN
        ELSE
            DO i = 1,m
                sum= sum+(vars(i,j)**2)/4000
                prod=prod* (cos(vars(i,j)/(i**(1.0/2))))
            END DO
            res(j) = sum-prod+1+biases(7)
        END IF
      END DO
             
      RETURN
      END SUBROUTINE F07_shifted_rotated_griewank
      
      !-------------------------------------------------------------------------
      ! F08_shifted_rotated_ackley_global_opt_bound
      !-------------------------------------------------------------------------
      ! References   :      Ackley,D.: An Empirical Study of Bit Vector Function Optimization.
      !                     Genetic Algorithms and Simulated Annealing,(1987) 170-215
      !
      !  Input        : vars    (R) a matrix containing an array of vectors
      !                             the function should process
      !                             vectors are saved per column at the
      !                             moment cause for a big problem size           
      !                             the population is smaller then the
      !                             dimension
      !                 m       (I) m Dimension of the matrix
      !                               = Dimension of Vectors
      !                 n       (I) n Dimension of the matrix
      !                               = # of Vectors to process              
      !                 lbounds (R) lower bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result
      !                 ubounds (R) upper bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result      
      !
      !
      !  Output       : res     (R) a row vector with the result of each
      !                             column vector of the vars input,has
      !                             size n      
      !-------------------------------------------------------------------------
      SUBROUTINE F08_shifted_rotated_ackley_global_opt_bound (res,vars,m,n,&
     &                                                         lbounds,ubounds)   
       
      INTEGER,INTENT(in)                      :: m
      INTEGER,INTENT(in)                      :: n
      REAL(MK),DIMENSION(n),INTENT(out)       :: res
      REAL(MK),DIMENSION(m,n),INTENT(inout)   :: vars
      REAL(MK),DIMENSION(m),OPTIONAL         :: lbounds
      REAL(MK),DIMENSION(m),OPTIONAL         :: ubounds
       
      ! Locally used variables
      REAL(MK),DIMENSION(m)                   :: ub,lb
      INTEGER                                  :: i  
      INTEGER                                  :: j,k 
      REAL(MK)                                 :: e1
      REAL(MK)                                 :: e2
        
      LOGICAL,DIMENSION(n)                     :: boundviolation
      REAL(MK),DIMENSION(m)                    :: tmp_result

      IF (.NOT. PRESENT(lbounds)) THEN
        lb = -32.000_MK
      ELSE
        lb = lbounds
      END IF

      IF (.NOT. PRESENT(ubounds)) THEN
        ub = 32.000_MK
      ELSE
        ub = ubounds
      END IF  

      boundviolation = .FALSE.

      !Shift the Input
      DO j=1,n
        DO i=1,m
            IF (vars(i,j) .GT. ub(i) .OR. vars(i,j) .LT. lb(i)) THEN
                boundviolation(j) = .TRUE.   
            ELSE
                vars(i,j) = vars(i,j) - shiftmatrix(i)
            END IF
        END DO
      END DO
       
       
      !rotate the input
        
      DO i=1,n
        IF (.NOT. boundviolation(i) ) THEN
            DO j=1,m
                tmp_result(j) = 0.0d0
                DO k=1,m
                    tmp_result(j) = tmp_result(j) + vars(k,i)*rotatematrix(k,j)
                END DO
            END DO
            vars(:,i) = tmp_result
        END IF
            
            !vars(:,i) = matmul(tmp_vec,tmp)
      END DO          
  
      DO j = 1,n
        e1 = 0.0_MK
        e2 = 0.0_MK
        
        IF (boundviolation(j) ) THEN
            res(j) = testf_NAN
        ELSE
            DO i = 1,m
                e1 = e1 + vars(i,j)**2
                e2 = e2 + cos(2.0_MK*Pi*vars(i,j))
            END DO
            res(j) = exp(1.0_MK) + 20.0_MK-20_MK*exp(-0.2_MK*sqrt(e1/real(m)))  
            res(j) = res(j) - exp(e2/real(m))
            res(j) = res(j) + biases(8)
        END IF                
      END DO

       RETURN
      END  SUBROUTINE F08_shifted_rotated_ackley_global_opt_bound
      
      !-------------------------------------------------------------------------
      ! F09_shifted_rastrigin
      !-------------------------------------------------------------------------
      ! References   :     Rastrigin,L. A.: Extremal Control Systems. In 
      !                    Theoretical Foundations of Engineering Cybernetics Series,
      !                    Moscow,Nauka,Russian (1974)
      !
      !  Input        : vars    (R) a matrix containing an array of vectors
      !                             the function should process
      !                             vectors are saved per column at the
      !                             moment cause for a big problem size           
      !                             the population is smaller then the
      !                             dimension
      !                 m       (I) m Dimension of the matrix
      !                               = Dimension of Vectors
      !                 n       (I) n Dimension of the matrix
      !                               = # of Vectors to process              
      !                 lbounds (R) lower bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result
      !                 ubounds (R) upper bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result      
      !
      !
      !  Output       : res     (R) a row vector with the result of each
      !                             column vector of the vars input,has
      !                             size n      
      !-------------------------------------------------------------------------
      SUBROUTINE F09_shifted_rastrigin (res,vars,m,n,lbounds,ubounds)   
       
      INTEGER,INTENT(in)                            :: m
      INTEGER,INTENT(in)                            :: n
      REAL(MK),DIMENSION(n),INTENT(out)             :: res
      REAL(MK),DIMENSION(m,n),INTENT(inout)         :: vars
      REAL(MK),DIMENSION(m),OPTIONAL               :: lbounds
      REAL(MK),DIMENSION(m),OPTIONAL               :: ubounds
       
      ! Locally used variables
      REAL(MK),DIMENSION(m)                         :: ub,lb
      INTEGER                                :: i  
      REAL(MK)                                  :: sum
      LOGICAL,DIMENSION(n)                          :: boundviolation
      INTEGER                                 :: j 

      IF (.NOT. PRESENT(lbounds)) THEN
        lb = -5.00_MK
      ELSE
        lb = lbounds
      END IF
       
      IF (.NOT. PRESENT(ubounds)) THEN
        ub = 5.00_MK
      ELSE
        ub = ubounds
      END IF  
       
      boundviolation = .FALSE.
     
      !Shift the Input
      DO j=1,n
        DO i=1,m
            IF (vars(i,j) .GT. ub(i) .OR. vars(i,j) .LT. lb(i)) THEN
                boundviolation(j) = .TRUE.        
            ELSE
                vars(i,j) = vars(i,j) - shiftmatrix(i)
            END IF
        END DO
      END DO   
    
      DO j = 1,n
         
        sum=0
        IF(boundviolation(j) ) THEN        
            res(j) = testf_NAN
        ELSE      
            sum = 10.0d0 * m
            DO i = 1,m
                sum=sum+vars(i,j)**2
                sum=sum-10.0d0*cos(2*Pi*vars(i,j))
            END DO
            res(j) = sum + biases(9)
         END IF   
      END DO
      
       RETURN  
      END SUBROUTINE F09_shifted_rastrigin
      
      !-------------------------------------------------------------------------
      ! F09_shifted_rastrigin_grad
      !-------------------------------------------------------------------------
      ! Remarks:           This Subroutine calculates the analytical Gradient of the
      !                    F09 CEC2005 Function. Used for testing reasons with BFGS.
      !
      ! References   :     Rastrigin,L. A.: Extremal Control Systems. In 
      !                    Theoretical Foundations of Engineering Cybernetics Series,
      !                    Moscow,Nauka,Russian (1974)
      !
      !  Input        : vars    (R) a matrix containing an array of vectors
      !                             the function should process
      !                             vectors are saved per column at the
      !                             moment cause for a big problem size           
      !                             the population is smaller then the
      !                             dimension
      !                 m       (I) m Dimension of the matrix
      !                               = Dimension of Vectors
      !                 n       (I) n Dimension of the matrix
      !                               = # of Vectors to process              
      !                 lbounds (R) lower bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result
      !                 ubounds (R) upper bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result      
      !
      !
      !  Output       : res     (R) a row vector with the result of each
      !                             column vector of the vars input,has
      !                             size n      
      !-------------------------------------------------------------------------
      SUBROUTINE F09_grad (grad,vars1,m,n,lbounds,ubounds)   
       
      INTEGER,INTENT(in)                            :: m
      INTEGER,INTENT(in)                            :: n
      REAL(MK),DIMENSION(m),INTENT(out)             :: grad
      REAL(MK),DIMENSION(m),INTENT(in)              :: vars1
      REAL(MK),DIMENSION(m),OPTIONAL               :: lbounds
      REAL(MK),DIMENSION(m),OPTIONAL               :: ubounds
       
      ! Locally used variables
      REAL(MK),DIMENSION(m)                       :: vars
      REAL(MK),DIMENSION(m)                         :: ub,lb
      INTEGER                                 :: i  
      REAL(MK)                                :: sum
      LOGICAL                                 :: boundviolation
      INTEGER                                 :: j 

      IF (.NOT. PRESENT(lbounds)) THEN
        lb = -5.00_MK
      ELSE
        lb = lbounds
      END IF
       
      IF (.NOT. PRESENT(ubounds)) THEN
        ub = 5.00_MK
      ELSE
        ub = ubounds
      END IF  
       
      boundviolation = .FALSE.
      vars = vars1
      !Shift the Input
      
        DO i=1,m
            IF (vars(i) .GT. ub(i) .OR. vars(i) .LT. lb(i)) THEN
                !Bboundviolation = .TRUE.        
            ELSE
                vars(i) = vars(i) - shiftmatrix(i)
            END IF
        END DO
      
    
        sum=0
        IF(boundviolation ) THEN        
            grad = testf_NAN
        ELSE            
            DO i = 1,m
                grad(i) = vars(i)*2+20*Pi*sin(2*Pi*vars(i))
            END DO
         END IF   

      
       RETURN  
      END SUBROUTINE F09_grad
      
      !-------------------------------------------------------------------------
      ! F10_shifted_rotated_rastrigin
      !-------------------------------------------------------------------------
      ! References   :     Rastrigin,L. A.: Extremal Control Systems. In 
      !                    Theoretical Foundations of Engineering Cybernetics Series,
      !                    Moscow,Nauka,Russian (1974)
      !
      !  Input        : vars    (R) a matrix containing an array of vectors
      !                             the function should process
      !                             vectors are saved per column at the
      !                             moment cause for a big problem size           
      !                             the population is smaller then the
      !                             dimension
      !                 m       (I) m Dimension of the matrix
      !                               = Dimension of Vectors
      !                 n       (I) n Dimension of the matrix
      !                               = # of Vectors to process              
      !                 lbounds (R) lower bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result
      !                 ubounds (R) upper bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result      
      !
      !
      !  Output       : res     (R) a row vector with the result of each
      !                             column vector of the vars input,has
      !                             size n      
      !-------------------------------------------------------------------------
      SUBROUTINE F10_shifted_rotated_rastrigin (res,vars,m,n,lbounds,ubounds)   
       
      INTEGER,INTENT(in)                            :: m
      INTEGER,INTENT(in)                            :: n
      REAL(MK),DIMENSION(n),INTENT(out)             :: res
      REAL(MK),DIMENSION(m,n),INTENT(inout)         :: vars
      REAL(MK),DIMENSION(m),OPTIONAL               :: lbounds
      REAL(MK),DIMENSION(m),OPTIONAL               :: ubounds
       
      ! Locally used variables
      REAL(MK),DIMENSION(m)                         :: ub,lb
      INTEGER                                   :: i  
      REAL(MK)                                  :: sum    
      LOGICAL,DIMENSION(n)                          :: boundviolation
      INTEGER                                   :: j,k 
      REAL(MK),DIMENSION(m)                     :: tmp_result

      IF (.NOT. PRESENT(lbounds)) THEN
        lb = -5.00_MK
      ELSE
        lb = lbounds
      END IF
       
      IF (.NOT. PRESENT(ubounds)) THEN
        ub = 5.00_MK
      ELSE
        ub = ubounds
      END IF  

      boundviolation = .FALSE.

      !Shift the Input
      DO j=1,n
        DO i=1,m
            IF (vars(i,j) .GT. ub(i) .OR. vars(i,j) .LT. lb(i)) THEN
                boundviolation(j) = .TRUE.   
            ELSE
                vars(i,j) = vars(i,j) - shiftmatrix(i)
            END IF
        END DO
      END DO        

       !rotate the input
        
      DO i=1,n
        IF (.NOT. boundviolation(i)) THEN
            DO j=1,m
                tmp_result(j) = 0.0_MK
                    DO k=1,m
                        tmp_result(j)=tmp_result(j)+vars(k,i)*rotatematrix(k,j)
                    END DO
            END DO
            vars(:,i) = tmp_result
        END IF
            
            !vars(:,i) = matmul(tmp_vec,tmp)
      END DO           
    
      DO j = 1,n
         
         sum=0
         IF(boundviolation(j) ) THEN
            res(j) = testf_NAN
         ELSE
            sum = 10.0d0 * m
            DO i = 1,m
                sum=sum+vars(i,j)**2
                sum=sum-10.0d0*cos(2*Pi*vars(i,j))
            END DO
            res(j) = sum + biases(9)
         END IF
      END DO
      
      END SUBROUTINE F10_shifted_rotated_rastrigin
      
      
      !-------------------------------------------------------------------------
      ! F11_Shifted_Rotated_Weierstrass
      !-------------------------------------------------------------------------
      ! References   :     2DO
      !
      !  Input        : vars    (R) a matrix containing an array of vectors
      !                             the function should process
      !                             vectors are saved per column at the
      !                             moment cause for a big problem size           
      !                             the population is smaller then the
      !                             dimension
      !                 m       (I) m Dimension of the matrix
      !                               = Dimension of Vectors
      !                 n       (I) n Dimension of the matrix
      !                               = # of Vectors to process              
      !                 lbounds (R) lower bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result
      !                 ubounds (R) upper bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result      
      !
      !
      !  Output       : res     (R) a row vector with the result of each
      !                             column vector of the vars input,has
      !                             size n      
      !-------------------------------------------------------------------------
      SUBROUTINE F11_shifted_rotated_weierstrass (res,vars,m,n,lbounds,ubounds)   
      
      INTEGER,INTENT(in)                            :: m
      INTEGER,INTENT(in)                            :: n
      REAL(MK),DIMENSION(n),INTENT(out)             :: res
      REAL(MK),DIMENSION(m,n),INTENT(inout)         :: vars
      REAL(MK),DIMENSION(m),OPTIONAL               :: lbounds
      REAL(MK),DIMENSION(m),OPTIONAL               :: ubounds
       
      ! Locally used variables
      REAL(MK),DIMENSION(m)                         :: ub,lb
      REAL(MK)                                  :: sum1,sum2
      INTEGER                                   :: Kmax = 20
      REAL(MK)                                  :: a = 0.5_MK
      REAL(MK)                                  :: b = 3.0_MK
      INTEGER                                   :: i  
      LOGICAL,DIMENSION(n)                          :: boundviolation
      INTEGER                                   :: j,k 
      REAL(MK),DIMENSION(m)                     :: tmp_result

      IF (.NOT. PRESENT(lbounds)) THEN
        lb = -0.5_MK
      ELSE
        lb = lbounds
      END IF
       
      IF (.NOT. PRESENT(ubounds)) THEN
        ub = 0.5_MK
      ELSE
        ub = ubounds
      END IF  

      boundviolation = .FALSE.


      ! Shift the Input
      DO j=1,n
        DO i=1,m
            IF (vars(i,j) .GT. ub(i) .OR. vars(i,j) .LT. lb(i)) THEN
                boundviolation(j) = .TRUE.   
            ELSE
                vars(i,j) = vars(i,j) - shiftmatrix(i)
            END IF
        END DO
      END DO         

       !rotate the input
        
      DO i=1,n
        IF (.NOT. boundviolation(i)) THEN
            DO j=1,m
                tmp_result(j) = 0.0_MK
                    DO k=1,m
                        tmp_result(j)=tmp_result(j)+vars(k,i)*rotatematrix(k,j)
                    END DO
            END DO
            vars(:,i) = tmp_result
        END IF
        
        !vars(:,i) = matmul(tmp_vec,tmp)
      END DO           
    
      DO j = 1,n
          
        sum1 = 0.0_MK
        sum2 = 0.0_MK
         
        IF (boundviolation(j) ) THEN
            res(j) = testf_NAN
        ELSE
            DO i = 1,m
                DO k = 0,Kmax
                    sum1 = sum1+a**k*cos((PI)*2.0_MK*b**k*(vars(i,j)+0.5_MK))
                END DO        
            END DO
         
            DO k = 0,Kmax
                sum2 = sum2 + a**k * cos((PI)*2.0_MK * b**k * (0.5_MK))
            END DO
            res(j) = sum1 - sum2*m + biases(11)
        END IF   
      END DO
      
      END SUBROUTINE F11_shifted_rotated_weierstrass
      
      
      
      !-------------------------------------------------------------------------
      ! F12_schwefel
      !-------------------------------------------------------------------------
      ! References   :     2DO
      !
      !  Input        : vars    (R) a matrix containing an array of vectors
      !                             the function should process
      !                             vectors are saved per column at the
      !                             moment cause for a big problem size           
      !                             the population is smaller then the
      !                             dimension
      !                 m       (I) m Dimension of the matrix
      !                               = Dimension of Vectors
      !                 n       (I) n Dimension of the matrix
      !                               = # of Vectors to process              
      !                 lbounds (R) lower bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result
      !                 ubounds (R) upper bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result      
      !
      !
      !  Output       : res     (R) a row vector with the result of each
      !                             column vector of the vars input,has
      !                             size n      
      !-------------------------------------------------------------------------
      SUBROUTINE F12_schwefel (res,vars,m,n,lbounds,ubounds)   
       
       
      INTEGER,INTENT(in)                            :: m
      INTEGER,INTENT(in)                            :: n
      REAL(MK),DIMENSION(n),INTENT(out)             :: res
      REAL(MK),DIMENSION(m,n),INTENT(inout)         :: vars
      REAL(MK),DIMENSION(m),OPTIONAL               :: lbounds
      REAL(MK),DIMENSION(m),OPTIONAL               :: ubounds
       
      ! Locally used variables
      REAL(MK),DIMENSION(m)                         :: ub,lb
      INTEGER                                        :: i  
      LOGICAL,DIMENSION(n)                          :: boundviolation
      INTEGER                                        :: j,k 
      REAL(MK),DIMENSION(m)                          :: tmp_result
      REAL(MK)                                       :: sum,tmp
      REAL(MK),DIMENSION(m)                          :: mB

      IF (.NOT. PRESENT(lbounds)) THEN
        lb = -PI
      ELSE
        lb = lbounds
      END IF
       
      IF (.NOT. PRESENT(ubounds)) THEN
        ub = PI
      ELSE
        ub = ubounds
      END IF  

      boundviolation = .FALSE.

      DO j = 1,n
         
         sum = 0.0_MK
           
         DO i = 1,m
         
            IF (vars(i,j) .GT. ub(i) .OR. vars(i,j) .LT. lb(i)) THEN
                boundviolation(j) = .TRUE. 
                EXIT   
            ELSE
                   
                mB(i) = 0.0d0
                DO k = 1,m
                    mB(i)=mB(i)+(schwefel_m_a(i,k)*sin(vars(k,j))+ &
     &                          schwefel_m_b(i,k) * cos(vars(k,j))) 
                END DO 
                tmp = schwefel_mA(i) - mB(i)
                sum = sum + (tmp * tmp)      
            END IF
         END DO
         
         IF (boundviolation(j) ) THEN
             res(j) = testf_NAN
         ELSE
             res(j) = sum + biases(12)
         END IF   
      END DO
      
      END SUBROUTINE F12_schwefel
    
    
      !-------------------------------------------------------------------------
      ! F13_shifted_expanded_griewank_rosenbrock
      !-------------------------------------------------------------------------
      ! References   :     2DO
      !
      !  Input        : vars    (R) a matrix containing an array of vectors
      !                             the function should process
      !                             vectors are saved per column at the
      !                             moment cause for a big problem size           
      !                             the population is smaller then the
      !                             dimension
      !                 m       (I) m Dimension of the matrix
      !                               = Dimension of Vectors
      !                 n       (I) n Dimension of the matrix
      !                               = # of Vectors to process              
      !                 lbounds (R) lower bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result
      !                 ubounds (R) upper bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result      
      !
      !
      !  Output       : res     (R) a row vector with the result of each
      !                             column vector of the vars input,has
      !                             size n      
      !-------------------------------------------------------------------------
      SUBROUTINE F13_shifted_expanded_griewank_rosenbrock (res,vars,m,n,&
     &                                                         lbounds,ubounds)   
       
      INTEGER,INTENT(in)                            :: m
      INTEGER,INTENT(in)                            :: n
      REAL(MK),DIMENSION(n),INTENT(out)             :: res
      REAL(MK),DIMENSION(m,n),INTENT(inout)         :: vars
      REAL(MK),DIMENSION(m),OPTIONAL               :: lbounds
      REAL(MK),DIMENSION(m),OPTIONAL               :: ubounds
       
      ! Locally used variables
      REAL(MK),DIMENSION(m)                         :: ub,lb
      INTEGER                                  :: i  
      LOGICAL,DIMENSION(n)                          :: boundviolation
      INTEGER                                  :: j,k 
      REAL(MK),DIMENSION(m)                    :: tmp_result
      REAL(MK)                                 :: sum,tmp
      REAL(MK),DIMENSION(m)                    :: mB
      REAL(MK)                                 :: t1,t2

      IF (.NOT. PRESENT(lbounds)) THEN
        lb = -3.0_MK
      ELSE
        lb = lbounds
      END IF
       
      IF (.NOT. PRESENT(ubounds)) THEN
        ub = 1.0_MK
      ELSE
        ub = ubounds
      END IF  

      boundviolation = .FALSE.


      ! Shift the Input
      DO j=1,n
        DO i=1,m
            IF (vars(i,j) .GT. ub(i) .OR. vars(i,j) .LT. lb(i)) THEN
                boundviolation(j) = .TRUE.   
            ELSE
                vars(i,j) = vars(i,j) - shiftmatrix(i)
            END IF
        END DO
      END DO      

      
    
      DO j = 1,n
        sum = 0.0_MK
              
        IF (boundviolation(j) ) THEN
            res(j) = testf_NAN
        ELSE
            
            DO i = 1,m-1
                !F2 part
                t1 = vars(i,j)**2 - vars(i+1,j)
                t2 = vars(i,j) - 1.0_MK
                t1 = ((100.0_MK * t1 * t1) + (t2*t2))
                !F8 part
                t2 = ((t1**2) / 4000.0_MK) - cos(t1) + 1.0_MK
                sum = sum + t2
            END DO
                     
            !F2 part
            t1 = vars(m,j)**2 - vars(1,j)
            t2 = vars(m,j) - 1.0_MK
            t1 = ((100.0_MK * t1 * t1) + (t2*t2))
            !F8 part
            t2 = ((t1**2) / 4000.0_MK) - cos(t1) + 1.0_MK
            sum = sum + t2
            res(j) = sum + biases(13)
         
         END IF   
      END DO
       
      END SUBROUTINE F13_shifted_expanded_griewank_rosenbrock  
      
      
      !-------------------------------------------------------------------------
      ! F14_shifted_rotated_expanded_scaffer
      !-------------------------------------------------------------------------
      ! References   :     2DO
      !
      !  Input        : vars    (R) a matrix containing an array of vectors
      !                             the function should process
      !                             vectors are saved per column at the
      !                             moment cause for a big problem size           
      !                             the population is smaller then the
      !                             dimension
      !                 m       (I) m Dimension of the matrix
      !                               = Dimension of Vectors
      !                 n       (I) n Dimension of the matrix
      !                               = # of Vectors to process              
      !                 lbounds (R) lower bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result
      !                 ubounds (R) upper bound of the area to process,
      !                             if the vector lies outside of this
      !                             bound it will return a penalty value
      !                             as result      
      !
      !
      !  Output       : res     (R) a row vector with the result of each
      !                             column vector of the vars input,has
      !                             size n      
      !-------------------------------------------------------------------------
      SUBROUTINE F14_shifted_rotated_expanded_scaffer (res,vars,m,n,&
     &                                                       lbounds,ubounds)   
       
      INTEGER,INTENT(in)                            :: m
      INTEGER,INTENT(in)                            :: n
      REAL(MK),DIMENSION(n),INTENT(out)             :: res
      REAL(MK),DIMENSION(m,n),INTENT(inout)         :: vars
      REAL(MK),DIMENSION(m),OPTIONAL               :: lbounds
      REAL(MK),DIMENSION(m),OPTIONAL               :: ubounds
       
      ! Locally used variables
      REAL(MK),DIMENSION(m)                         :: ub,lb
      INTEGER                                   :: i  
      LOGICAL,DIMENSION(n)                          :: boundviolation
      INTEGER                                   :: j,k 
      REAL(MK),DIMENSION(m)                     :: tmp_result
      REAL(MK)                                  :: sum,tmp
      REAL(MK),DIMENSION(m)                     :: mB
      REAL(MK)                                  :: temp1,temp2,temp3

      IF (.NOT. PRESENT(lbounds)) THEN
        lb = -100.0_MK
      ELSE
        lb = lbounds
      END IF
       
      IF (.NOT. PRESENT(ubounds)) THEN
        ub = 100.0_MK
      ELSE
        ub = ubounds
      END IF  

      boundviolation = .FALSE.


      ! Shift the Input
      DO j=1,n
        DO i=1,m
            IF (vars(i,j) .GT. ub(i) .OR. vars(i,j) .LT. lb(i)) THEN
                boundviolation(j) = .TRUE.   
            ELSE
                vars(i,j) = vars(i,j) - shiftmatrix(i)
            END IF
        END DO
      END DO    
      
      !rotate the input
        
      DO i=1,n
        IF (.NOT. boundviolation(i)) THEN
            DO j=1,m
                tmp_result(j) = 0.0_MK
                DO k=1,m
                    tmp_result(j) = tmp_result(j)+ vars(k,i) *rotatematrix(k,j)
                END DO
            END DO
            
            vars(:,i) = tmp_result
            !vars(:,i) = matmul(tmp_vec,tmp)
        END IF
      END DO
    
      DO j = 1,n
        sum = 0.0_MK 
          
        IF (boundviolation(j) ) THEN
             res(j) = testf_NAN
        ELSE
         
            DO i = 1,m-1
                temp1 = vars(i,j)**2 + vars(i+1,j)**2
                temp2 = sin(sqrt(temp1))
                temp3 = 1.0_MK + 0.001_MK * temp1
                sum = sum + (0.5_MK + ((temp2**2 - 0.5_MK)/(temp3**2)))
            END DO
         
            temp1 = vars(m,j)**2 + vars(1,j)**2
            temp2 = sin(sqrt(temp1))
            temp3 = 1.0_MK + 0.001_MK * temp1
            sum = sum + (0.5_MK + ((temp2**2 - 0.5_MK)/(temp3**2)))
            
            res(j) = sum + biases(14)
        END IF   
      
      END DO
       
      END SUBROUTINE F14_shifted_rotated_expanded_scaffer
      
      END MODULE CEC2005
     
