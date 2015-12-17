        !-------------------------------------------------------------------------
        !  Module       :            testCEC2005
        !-------------------------------------------------------------------------
        !
        !  Purpose      : Module to check correctness of the CEC2005 test suit implementation
        !
        !  Remarks      : 
        !
        !  References   : 
        !
        !  Revisions    :
        !-------------------------------------------------------------------------
        
        MODULE testCEC2005
        USE CEC2005
        IMPLICIT NONE
        
        
        REAL(MK),DIMENSION(25)    :: lower_boundary = &
     &            (/-100.0_MK, -100.0_MK, -100.0_MK, -100.0_MK, -100.0_MK, &
     &              -100.0_MK, 0.0_MK, -32.0_MK, -5.0_MK, -5.0_MK,&
     &              -0.5_MK, &
     &              -3.14159265358979323846264338327950288419716939937510d0,&
     &               -3.0_MK, -100.0_MK, -5.0_MK,&
     &              -5.0_MK, -5.0_MK, -5.0_MK, -5.0_MK, -5.0_MK,&
     &              -5.0_MK, -5.0_MK, -5.0_MK, -5.0_MK, 2.0_MK /)
        REAL(MK),DIMENSION(25)    :: upper_boundary = &
     &            (/100.0_MK, 100.0_MK, 100.0_MK, 100.0_MK, 100.0_MK, &
     &              100.0_MK, 600.0_MK, 32.0_MK, 5.0_MK, 5.0_MK,&
     &              0.5_MK,&
     &              3.14159265358979323846264338327950288419716939937510d0,&
     &              1.0_MK, 100.0_MK, 5.0_MK,&
     &              5.0_MK, 5.0_MK, 5.0_MK, 5.0_MK, 5.0_MK,&
     &              5.0_MK, 5.0_MK, 5.0_MK, 5.0_MK, 5.0_MK /)                    
        
        
        
        CONTAINS
        
      !----------------------------------------------------------------------
      ! test_minima
      !-------------------------------------------------------------------------
      !  Purpose      : tests the benchmark on a given functionnr (against given results)
      ! 
      !  Input       :  nr  (I) the number of the function to be benchmarked - 0 = all
      ! 
      !  References:
      !-------------------------------------------------------------------------  
      
      recursive SUBROUTINE test_minima(nr)
      !parameters
      INTEGER, INTENT(IN)                           :: nr
      
      !local used vars
      REAL(MK),POINTER,DIMENSION(:,:)   :: res_matrix
      INTEGER                           :: alloc_status, dealloc_status,i
      REAL(MK),DIMENSION(1)             :: res
      
      IF (nr .EQ. 0) THEN
        DO i=1, 25
            CALL test_minima(i)
        END DO
      ELSE        
          allocate(res_matrix(25,50),stat=alloc_status)                  
          CALL loadMatrixFromFile ('supportData/global_optima.txt',25,50,&
     &    res_matrix)          
          CALL prepare_benchmark(nr,50,1)   !first initialize the benchmark  
          CALL benchmark(res,res_matrix(nr,1:50),50,1)
          
          write (*,*) '------------------------------ Test for global optimum'  
          write (*,"(I2)") nr
          write (*,*) TRIM(adjustl(function_name))
          
          IF (nr .EQ. 8 .OR. nr .EQ. 5 .OR. nr .EQ. 20) THEN
              write (*,*) 'skipped global Optima test (outside of Bounds)'
          ELSE
              IF (res(1) .EQ. biases(nr)) THEN            
                  write (*,*) 'Global Optimum correct!'            
              ELSE
                  write(*,*) 'Global Optimum differs by ', biases(nr)-res(1)            
              END IF        
          END IF
          
          
          CALL close_benchmark()                             
          deallocate(res_matrix,stat=dealloc_status)
      END IF
      END SUBROUTINE test_minima
        
        
        
      !-------------------------------------------------------------------------
      ! test_boundaries
      !-------------------------------------------------------------------------
      !  Purpose      : tests the benchmark on a given functionnr (against given results)
      ! 
      !  Input       :  nr  (I) the number of the function to be benchmarked - 0 = all
      ! 
      !  References:
      !-------------------------------------------------------------------------  
      
      recursive SUBROUTINE test_boundaries(nr,m)
      !parameters
      INTEGER, INTENT(IN)                           :: nr
      INTEGER, INTENT(IN)                            :: m
      
      !locals
      REAL(MK)                                      :: dummy      
      INTEGER                                       :: i,j,k ,l    
      REAL(MK), DIMENSION(m)          :: res
      REAL(MK), DIMENSION(m,m)        :: vars
      INTEGER                         :: alloc_status, dealloc_status  
      
     
      IF (nr .EQ. 0) THEN
        DO i=1, 25
            CALL test_boundaries(i,m)
        END DO
      ELSE
      
           
            !check if lower boundary pass works on every dimension
            
            DO k=1,m
                vars(k,1) = NEAREST(lower_boundary(nr),1.0_MK)                            
            END DO
            
            
            CALL prepare_benchmark(nr,m,1)   !first initialize the benchmark    
            CALL benchmark(res,vars,m,1)
            CALL close_benchmark()          
            
            IF (nr .EQ. 25 .OR. nr .EQ. 7) THEN
            write (*,*) '------------------------------ Boundary Test'  
            write (*,"(I2)") nr
            write (*,*) TRIM(adjustl(function_name))
            write (*,*)'Skipped the test cause this function is without bounds'
            RETURN
            END IF   
      
          
            write (*,*) '------------------------------ Boundary Test'  
            write (*,"(I2)") nr
            write (*,*) TRIM(adjustl(function_name))

       
            
           
            
            IF (res(1) .EQ. testf_NAN) THEN
                write (*,*) ' Lower Boundary pass failed!' 
            ELSE     
                write (*,*) 'Lower Boundary pass works ok!' 
            ENDIF
            
            
            
            !check if lower boundary violation works on every dimension
            DO j = 1, m
                DO k=1,m
                    vars(k,j) = NEAREST(lower_boundary(nr),1.0_MK)                          
                END DO
                l = NEAREST(3.0_MK,1.0_MK)
                vars(j,j) = NEAREST(lower_boundary(nr),-1.0_MK)
            END DO
            CALL prepare_benchmark(nr,m,m)   !first initialize the benchmark    
            CALL benchmark(res,vars,m,m)
            CALL close_benchmark()
            
            k=0
            DO j=1,m
                IF (res(j) .NE. testf_NAN) THEN
                    write (*,*) j,' Dim. Lower Boundary violation failed!' 
                    k = 1
                ENDIF
            END DO
            IF (k .EQ. 0) THEN
                    write (*,*) 'Lower Boundary violoation works ok!' 
            ENDIF
                    
            
            !check if upper boundary pass works on every dimension
            
            DO k=1,m
                vars(k,1) = NEAREST(upper_boundary(nr),-1.0_MK)
            END DO
            
            CALL prepare_benchmark(nr,m,1)   !first initialize the benchmark    
            CALL benchmark(res,vars,m,1)
            CALL close_benchmark()          
            IF (res(1) .EQ. testf_NAN) THEN
                write (*,*) 'Upper Boundary pass failed!' 
            ELSE     
                write (*,*) 'Upper Boundary pass works ok!' 
            ENDIF
            
            
            
            !check if Upper Boundary violation works on every dimension
            DO j = 1, m
                DO k=1,m
                    vars(k,j) = NEAREST(upper_boundary(nr),-1.0_MK)
                END DO                          
                vars(j,j) = NEAREST(upper_boundary(nr),1.0_MK)
            END DO
            CALL prepare_benchmark(nr,m,m)   !first initialize the benchmark    
            CALL benchmark(res,vars,m,m)
            CALL close_benchmark()
            
            k=0
            DO j=1,m
                IF (res(j) .NE. testf_NAN) THEN
                    write (*,*) j,' Dim. Upper Boundary violation failed!' 
                    k = 1
                ENDIF
            END DO
            IF (k .EQ. 0) THEN
                    write (*,*) 'Upper Boundary violoation works ok!' 
            ENDIF
            
            
            
            
      
      END IF
      END SUBROUTINE test_boundaries
      



        
        
      !-------------------------------------------------------------------------
      ! test_benchmark
      !-------------------------------------------------------------------------
      !  Purpose      : tests the benchmark on a given functionnr (against given results)
      !                 Boundaries set to +/- Infinity.
      !                 all functions that include noise might give bad results - for that
      !                 reason compile without noise to get proper results
      ! 
      !  Input       :  nr  (I) the number of the function to be benchmarked - 0 = all
      ! 
      !  References:
      !-------------------------------------------------------------------------
      
      
      
      recursive SUBROUTINE test_benchmark(nr)
      !parameters
      INTEGER, INTENT(IN)                            :: nr
!     REAL(MK),POINTER,DIMENSION(:)                  :: biases
!   REAL(MK),DIMENSION(:),POINTER                  :: shiftmatrix
!   REAL(MK),DIMENSION(:,:),POINTER                :: rotatematrix
      !localy used variables
      INTEGER                                       :: i
      CHARACTER(len=200)                            :: testfile
      CHARACTER(len=4)                              :: nrstring,istring,&
     &                                                 numstring,dimstring
      INTEGER                                       :: num_test_points = 10
      INTEGER                                       :: test_dimension =  50
      
      REAL(MK),POINTER,DIMENSION(:,:)               :: transpose_vars
      REAL(MK),POINTER,DIMENSION(:,:)               :: test_vars
      REAL(MK),POINTER,DIMENSION(:)                 :: test_res
      REAL(MK),POINTER,DIMENSION(:)                 :: calculated_res
      REAL(MK),DIMENSION(50)                        :: boundarray_up
      REAL(MK),DIMENSION(50)                        :: boundarray_down
      REAL(MK)                                      :: dummy
      REAL(MK)                                      :: tstart,tend,t1,t2
      
      INTEGER                                       :: alloc_status,&
     &                                                 dealloc_status    
      dummy = HUGE(dummy)
      
      boundarray_up = dummy
      boundarray_down = -dummy
      
      
      
      IF (nr .EQ. 0) THEN
        DO i=1, 25
            CALL test_benchmark(i)
        END DO
      ELSE
      

      
      write(nrstring,'(I3)') nr
      testfile = 'testData/test_data_func'//TRIM(adjustl(nrstring))//'.txt';


      CALL prepare_benchmark(nr,test_dimension,num_test_points)         
      allocate(test_res(num_test_points),stat=alloc_status)
      allocate(calculated_res(num_test_points),stat=alloc_status)
      allocate(test_vars(test_dimension,num_test_points),stat=alloc_status)
      allocate(transpose_vars(num_test_points,test_dimension),stat=alloc_status)
          
          
      CALL loadMatrixFromFile (testfile,num_test_points,test_dimension,&
     &transpose_vars,test_res)
      
        CALL CPU_TIME(tstart)  
        DO i = 1, 100   !to get a time that is meaningful
          test_vars = transpose(transpose_vars)
          CALL benchmark(calculated_res,test_vars,test_dimension,&
     &    num_test_points,boundarray_down,boundarray_up)
      END DO
      CALL CPU_TIME(tend) 
      t1 = tend - tstart
      write(istring,'(I3)') i-1
      write(numstring,'(I3)') num_test_points
      write(dimstring,'(I3)') test_dimension     
     
      write (*,*) ''                  
      write (*,*) '------------------------------ Test against known results'  
      write (*,"(I2)") nr
      write (*,*) TRIM(adjustl(function_name))
      write (*,*) 'calling ',TRIM(adjustl(istring))&
     &,' times with ',TRIM(adjustl(numstring)),&
     &' Test points ',TRIM(adjustl(dimstring)),&
     &' dimensional took ',t1,' seconds'
          
      DO i=1, num_test_points
      
          IF (test_res(i)-calculated_res(i) .EQ. 0.0_MK) THEN
            write (*,*) i,' Testpoints OK!' 
          ELSE
            write (*,*) i,' Testpoint has a difference of:'
            write (*,*) i,' absolut: ',abs((test_res(i)-calculated_res(i)))
            write (*,*) i,' ratio: ',abs((test_res(i)-calculated_res(i)))/&
     &      abs(calculated_res(i)), '%'
              
              
          END IF
      
            !write (*,*) calculated_res(i)
      END DO
      
      CALL close_benchmark()
      deallocate(test_res,stat=alloc_status)
      deallocate(calculated_res,stat=alloc_status)
      deallocate(test_vars,stat=alloc_status)
      deallocate(transpose_vars,stat=alloc_status)
      
      END IF         
      END SUBROUTINE test_benchmark
      
      END MODULE testCEC2005
