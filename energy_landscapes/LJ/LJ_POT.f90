 	  !-------------------------------------------------------------------------
      !  SUBROUTINE   :                    LJ_POT
      !-------------------------------------------------------------------------
      !
      !  Purpose      : calculates LJ Potential 
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


       SUBROUTINE LJ_POT(res,vars,m,n,lbounds,ubounds) 
	   USE cmaes_param_mod
	   !USE deriv_class      ! make the module accessible
	   !Parameters
	   REAL(MK), DIMENSION(n),  INTENT(out)           :: res
	   REAL(MK), DIMENSION(m,n),INTENT(in)            :: vars
	   INTEGER, INTENT(in)                            :: m ! (3*LJ_N-6) dim vector
	   INTEGER, INTENT(in)                            :: n
	   
	   REAL(MK), OPTIONAL                             :: lbounds
	   REAL(MK), OPTIONAL                             :: ubounds
	   	   
       ! Locally used variables
       REAL(MK)										  :: dist
       REAL(MK), DIMENSION((m+6)/3,3)				  :: part_mat
       REAL(MK), DIMENSION(1,3)				  		  :: part_mean
       REAL(MK)                                       :: ub,lb
       REAL(MK)							              :: LJsigma,LJeps
       INTEGER                                        :: i,j,k,LJ_N
       
       IF (.NOT. PRESENT(lbounds)) THEN
       	lb = -3_MK
       ELSE
        lb = lbounds	
       END IF
       
       IF (.NOT. PRESENT(ubounds)) THEN
        ub = 3_MK
       ELSE
        ub = ubounds
       END IF
       
       !-------------------------------------------------------------------------
       !  Lennard Jones Parameters
       !------------------------------------------------------------------------- 
       
       ! Lennard-Jones Potential
       !------------------------------------------------------------------------- 
       LJsigma=1.0_MK
       LJeps=1.0_MK
       
       
       ! number of particles 
       LJ_N = (m+6)/3
       
       
       res(:) = 0.0_MK
       
       ! loop over all particle configurations
       ! Sequential computation of n solutions
       DO i=1,n
       	 
       	 ! place first particle at the origin
      	 part_mat(:,:)=0.0_MK
       	 ! place second particle at the x-axis
       	 part_mat(2,1)=vars(1,i)       
         part_mat(2,2)=0.0_MK       
         part_mat(2,3)=0.0_MK
         ! place third particle in the x-y plane
       	 part_mat(3,1)=vars(2,i)       
         part_mat(3,2)=vars(3,i)       
         part_mat(3,3)=0.0_MK
                
         ! build particle position matrix
         DO j=4,m-2,3
         	k=(j+8)/3
         	part_mat(k,1)=vars(j,i)
         	part_mat(k,2)=vars(j+1,i)
         	part_mat(k,3)=vars(j+2,i)
         	
         END DO    		
         
         ! initialize distance
         dist=0.0_MK
         
         ! compute distance matrix and the Lennard-Jones energy
         DO k=1,LJ_N-1		
        	DO l=k+1,LJ_N
        		dist=(sum((part_mat(k,:)-part_mat(l,:))**2))**(0.5)  
        		res(i) = res(i) + 4.0_MK*LJeps*((LJsigma/dist)**12 - (LJsigma/dist)**6) 
        	END DO  		
         END DO
         
       END DO
       	   
	  END SUBROUTINE LJ_POT
