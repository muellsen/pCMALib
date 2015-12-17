 	  !-------------------------------------------------------------------------
      !  SUBROUTINE   :                    LJ_POT_COMP
      !-------------------------------------------------------------------------
      !
      !  Purpose      : calculates LJ Potential with applied compression term
      !                
      !  Remarks      : Minima for LJ 38 are fcc like with slight variation
      !                 dependent on compression factor mu.
      !                 Minimal values for LJ 38 for mu=0:0.25:5
      !                 -173.928426589861 -151.101188107311	-128.429136646439
      !                 -105.900692687039 -83.5059543852499	-61.2363482726032
      !                 -39.0843723641102 -17.0434022950991	 4.89245748612146
      !                  26.7284904120414  48.4694603891765	 70.1196863760991
      !                  91.6831031448869  113.163311320461	 134.563618980143
      !                  155.887076549895  177.136506231594	 198.314527184472
      !                  219.423576974386  240.465930108432	 261.443714169179
      !
      !  References   : J.P.K Doye, Phys.Rev. E 62, 2000
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  pCMALib: a parallel fortran 90 library for the evolution strategy with
      !           covariance matrix adaptation
      !  Christian L. Mueller, Benedikt Baumgartner, Georg Ofenbeck
      !  MOSAIC group, ETH Zurich, Switzerland
      !-------------------------------------------------------------------------


       SUBROUTINE LJ_POT_COMP(res,vars,m,n,lbounds,ubounds) 
	   USE cmaes_param_mod
	   USE cmaes_opts_mod
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
       REAL(MK)							              :: LJsigma,LJeps,mu_comp
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
       
       ! Compression weight
	   mu_comp = options%LJ_comp*LJeps
       
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
        
         part_mean(1,1) = 1./LJ_N*sum(part_mat(:,1))
         part_mean(1,2) = 1./LJ_N*sum(part_mat(:,2))
         part_mean(1,3) = 1./LJ_N*sum(part_mat(:,3))
         	
         DO k=1,LJ_N
         	res(i) = res(i) + mu_comp * sum((part_mat(k,:)-part_mean(1,:))**2) / LJsigma**2 
         END DO
                
        END DO
       	   
	  END SUBROUTINE LJ_POT_COMP
