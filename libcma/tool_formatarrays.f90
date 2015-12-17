      !-------------------------------------------------------------------------
      !  Subroutine      :                    tool_formatarrays 
      !-------------------------------------------------------------------------
      ! 
      !  Purpose      : sets format descriptor to enable a nice looking output
      !                    of matrices and vectors
      !
      !  Remarks      : values are displayed in exponential form 
      !                    with 3 decimal places in a field of width 12
      !
      !  Input          : D        (I) Number of columns (matrix) or vector length
      !
      !  Input/Output : form    (C) format descriptor
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
      SUBROUTINE tool_formatarrays(form,D)
      IMPLICIT NONE
      INTEGER,INTENT(in)    :: D
      CHARACTER(len=10) :: dimension
      CHARACTER(len=20),INTENT(inout) :: form
      form = '(ES12.3))'
      WRITE(dimension,'(I10)') D ! Convert D to String
      form = '(3X,' // trim(adjustl(dimension)) // form ! Set format descriptor
      RETURN
      END SUBROUTINE tool_formatarrays
