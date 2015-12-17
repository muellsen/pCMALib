 	  !-------------------------------------------------------------------------
      !  Subroutine   :          tool_realtxt2arr
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Read data (REAL numbers) from .txt file to an array
      !
      !  Input		  : filename	(C) Name of file to read from
      !					N			(I)	Number of rows of output array
      !					M			(I) Number of columns of output array
	  !
	  !  Output		  : arr(:,:)	(R) Matrix containing values from file
	  !         
      !  Remarks      : txt file must be in ascii format and each value in its
      !					own line. The Routine reads the file line by line and
      !					creates an Array (column-wise) with specified 
      !					dimensions.
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
	  SUBROUTINE tool_realtxt2arr(filename,N,M,arr)
	  USE cmaes_param_mod	  
	  IMPLICIT NONE
	  CHARACTER(len=30),INTENT(in)			:: filename
	  INTEGER,INTENT(in)		 			:: N
	  INTEGER,INTENT(in)					:: M
	  REAL(MK),DIMENSION(N,M),INTENT(out) 	:: arr
	  
	  INTEGER 								:: i,j,k
	  REAL(MK),DIMENSION(N*M) 				:: vec
	  INTEGER 								:: ioerror
	  
	  !WRITE(*,*) filename
	  OPEN(unit=50,file=trim(filename),status='old',action='read',iostat=ioerror)	  
	  IF(ioerror .EQ. 0) THEN
	  	DO i = 1, N*M
	  		READ(50,*) vec(i)
	  	END DO
	  ELSE
	  	WRITE(*,*) 'Error reading from file.',filename
	  END IF
	  CLOSE(unit=50)

	  k = 1
	  DO j=1,M
	  	DO i=1,N
	  		arr(i,j) = vec(k)
	  		k = k+1
	  	END DO
	  END DO

	  RETURN
	  END SUBROUTINE tool_realtxt2arr