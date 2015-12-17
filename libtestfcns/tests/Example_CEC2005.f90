      !****************************************************************************
      !
      !  PROGRAM: Example_CEC2005
      !
      !  PURPOSE:  Entry point for the console application.
      !            Demonstrates various ways how to use the CEC2005 Test Suit
      !
      ! Georg Ofenbeck
      ! MOSAIC Group, ETH Zurich, Switzerland
      !****************************************************************************

      program Example_CEC2005
    
      !USE testf_module
      USE CEC2005
      USE testCEC2005
      implicit none

      INTEGER, PARAMETER                :: n = 1
      INTEGER, PARAMETER                :: m = 50
      CHARACTER(len=200)                :: name
      REAL(MK), DIMENSION(n)            :: res
        REAL(MK), DIMENSION(m,n)          :: vars
        REAL(MK)                          :: dummy
        REAL(MK)                          :: ZBQLNOR


    
      ! First a few functions that are implemented in testCEC2005.F90 (not needed in practice!!)
      !--------------------------------------------------------------------
    
    
      !To test the implemented Functions vs the Results provided by Matlab:
      !Note: turn of Noise (Noise Parameter in CEC2005.F90) to get usefull results for the Functions
      !      containing Noise
    
      CALL test_benchmark(0) !This function calls the Benchmark Function given
      !by the number handed over as parameter. If 0 is given
      !to the function it will sequentially call all 25 Functions and 
      !compare the calculated results to the given Results.
      CALL test_boundaries(0,50) !This function tests if the boundary handling is implemented
      !correctly. The First Parameter is again the function nr - the 
      !second set in how many dimensions to test it.
    
      CALL test_minima(0)   !Test_minimas compares the Minimas given also as supporting Material
      !against the Shifts for each function (Functions are 0 leveld and by that
      ! it should be the same value)
    
    
      ! How to use the CEC2005 Suit in Practice !
      !------------------------------------------------------------------------   
    
    
    
      CALL prepare_benchmark(12,m,n)   !first initialize the benchmark
      ! Parameter1: Integer that tells which Function Nr to initialize
      ! Parameter2: Dimension of a Test Point
      ! Parameter3: # of Test points to be calculated at once 
                                    
      !optional parameters are datafile,rotatefile,rotatesuffix which give the option to
      !   use different Supporting Files for the Functions
    
         
      CALL benchmark(res,vars,m,n) !Then call the benchmark function as often as needed
      ! Parameter1: Integer that tell wich function Nr to use
      ! Parameter2: Vector of Dimension n that will recive the Result(s)
      ! Parameter3: Matrix of Dimension m*n wich carries n m-dimensional Input Variables
      !       !!!! the Input Variables will get modified during the call (INOUT Parameter)!!!!
      ! Parameter4: the dimension of the Input Variables
      ! Parameter5: the number of input Variables and accordingly also the number of Results
                                    
      !optional: 
      !Parameter6: m-dim Array with the Lower Boundaries for each Dimension
      !Parameter7: m-dim Array with the Upper Boundaries for each Dimension
      !Parameter8: alternative BiasFile
      !Parameter9: alternative DataFile (Shifts)
      !Parameter10: alternative Rotationfile (without suffix)
      !Parmaeter11: alternative Suffix for the Rotationfile (default .txt)
    
    
      CALL close_benchmark()          !After all work with the Benchmark Function is done call close_benchmark to
      !release all allocated Memory
    
    
      end program Example_CEC2005

