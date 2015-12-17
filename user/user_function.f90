      !-------------------------------------------------------------------------
      !  Routine       :            user_function
      !-------------------------------------------------------------------------
      !
      !  Purpose      :  
      !
      !  Remark       : at the moment pcmalib will only call single values at once
      !                 meaning it will send 'vars' of size m x 1 and expect a vector of size
      !                 n=1 back. This might change in the future and therefore it is recommended
      !                 to write the function in a way that it can handle multiple input vectors in form
      !                 of a matrix and return a vector of results
      !
      !  Input       :  
      !                 m       (I) m Dimension of the matrix (Rows)
      !                               = Dimension of Vectors
      !                 n       (I) n Dimension of the matrix (Columns)
      !                               = # of Vectors to process
      !  Input(optional): lbounds / ubounds - m dimensional array with boundaries given to CMA-ES
      !                   
      !
      !  Input/Output:  vars    (R) the matrix with the input values of size m*n
      !
      !  Output:        res     (R) the vector with the results of size n
      !
      !-------------------------------------------------------------------------


      SUBROUTINE user_function(res,vars,m,n,lbounds,ubounds) 
      USE cmaes_param_mod
      !-------------------------------------------------------------------------
      !  Parameters
      !-------------------------------------------------------------------------
      REAL(MK),DIMENSION(n),INTENT(out)            :: res
      REAL(MK),DIMENSION(m,n),INTENT(in)           :: vars
      INTEGER,INTENT(in)                           :: m
      INTEGER,INTENT(in)                           :: n
      REAL(MK),DIMENSION(m),OPTIONAL              :: lbounds
      REAL(MK),DIMENSION(m),OPTIONAL              :: ubounds
      
      
      
      !your code here!!!!!
      WRITE(*,*) 'no code provided in the user_function yet'
      STOP
      
      
      
      END SUBROUTINE user_function
