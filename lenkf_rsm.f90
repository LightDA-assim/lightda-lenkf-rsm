! Copyright (c) 2004-2018 Lars Nerger
!
! This file is part of PDAF.
!
! PDAF is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License
! as published by the Free Software Foundation, either version
! 3 of the License, or (at your option) any later version.
!
! PDAF is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with PDAF.  If not, see <http://www.gnu.org/licenses/>.
!
!$Id: PDAF-D_lenkf_analysis_rsm.F90 29 2018-03-09 20:06:46Z lnerger $
!BOP
!
! !ROUTINE: PDAF_lenkf_analysis_rsm --- Perform LEnKF analysis step
!
! !INTERFACE:
SUBROUTINE lenkf_analysis_rsm(step,ind_p,dim_p, dim_obs_p, dim_obs, dim_ens, rank_ana, &
     state_p, ens_p, resid, U_add_obs_err, U_localize, forget, flag)


! !DESCRIPTION:
! Analysis step of ensemble Kalman filter with 
! representer-type formulation.  In this version 
! HP is explicitly computed.  This variant is 
! optimal if the number of observations is 
! smaller than or equal to half of the ensemble 
! size.
! The final ensemble update uses a block
! formulation to reduce memory requirements.
!
! Variant for domain decomposition.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2003-10 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
  ! (Defines BLAS/LAPACK routines and MPI_REALTYPE)

  USE iso_c_binding

  IMPLICIT NONE

  ! !ARGUMENTS:
  INTEGER(c_int32_t), INTENT(in) :: step      ! Iteration number
  INTEGER(c_int32_t), INTENT(in) :: ind_p      ! Integer index of PE
  INTEGER(c_int32_t), INTENT(in)  :: dim_p     ! PE-local dimension of model state
  INTEGER(c_int32_t), INTENT(in) :: dim_obs_p  ! PE-local dimension of observation vector
  INTEGER(c_int32_t), INTENT(in) :: dim_obs    ! Global dimension of observation vector
  INTEGER(c_int32_t), INTENT(in)  :: dim_ens   ! Size of state ensemble
  INTEGER(c_int32_t), INTENT(in)  :: rank_ana  ! Rank to be considered for inversion of HPH
  REAL(c_double), INTENT(inout)  :: state_p(dim_p)        ! PE-local ensemble mean state
  REAL(c_double), INTENT(inout)  :: ens_p(dim_p, dim_ens) ! PE-local state ensemble
  REAL(c_double), INTENT(in)     :: resid(dim_obs, dim_ens) ! Global array of residuals
  REAL(c_double), INTENT(in)     :: forget    ! Forgetting factor
  INTEGER(c_int32_t), INTENT(inout) :: flag    ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_add_obs_err, &  ! Add observation error covariance matrix
       U_localize               ! Apply localization to HP and HPH^T

! !CALLING SEQUENCE:
! Called by: PDAF_enkf_update
! Calls: U_add_obs_err
! Calls: U_localize
! Calls: gemmTYPE (BLAS; dgemm or sgemm dependent on precision)
! Calls: gesvTYPE (LAPACK; dgesv or sgesv dependent on precision)
! Calls: syevxTYPE (LAPACK; dsyevx or ssyevx dependent on precision)
!EOP

! *** local variables ***
  INTEGER :: i, j, member  ! counters
  REAL(c_double) :: invdim_ens       ! inverse of ensemble size
  REAL(c_double) :: invdim_ensm1     ! inverse of ensemble size minus 1
  REAL(c_double) :: sqrtinvforget    ! square root of inverse forgetting factor
  INTEGER, SAVE :: allocflag = 0       ! Flag for first-time allocation
  INTEGER, SAVE :: allocflag_b = 0     ! Flag for first-time allocation
  REAL(c_double), ALLOCATABLE :: HP_p(:,:)       ! Temporary matrix for analysis
  REAL(c_double), ALLOCATABLE :: HPH(:,:)        ! Temporary matrix for analysis
  REAL(c_double), ALLOCATABLE :: XminMean_p(:,:) ! Temporary matrix for analysis
  REAL(c_double), ALLOCATABLE :: m_state_p(:)    ! PE-local observed state
  REAL(c_double), ALLOCATABLE :: HXmean_p(:)     ! Temporary vector for analysis
  INTEGER, ALLOCATABLE :: ipiv(:)      ! vector of pivot indices
  INTEGER :: sgesv_info                ! output flag of SGESV

  ! *** Variables for variant using pseudo inverse with eigendecompositon
  REAL(c_double), ALLOCATABLE :: eval(:)         ! vector of eigenvalues
  REAL(c_double), ALLOCATABLE :: rwork(:)        ! workarray for eigenproblem
  REAL(c_double), ALLOCATABLE :: evec(:,:)       ! matrix of eigenvectors
  REAL(c_double), ALLOCATABLE :: evec_temp(:,:)  ! matrix of eigenvectors
  REAL(c_double), ALLOCATABLE :: repres(:,:)     ! matrix of representer vectors
  INTEGER :: syev_info     ! output flag of eigenproblem routine
  REAL(c_double)    :: VL, VU         ! temporary variables for SYEVX (never really used)
  INTEGER :: Ilower, Iupper ! variables defining the interval of eigenvalues
  REAL(c_double)    :: abstol         ! error tolerance for eigenvalue problem
  INTEGER :: nEOF           ! number of EOFs as computed by SYEVX
  INTEGER, ALLOCATABLE :: iwork(:)     ! workarray for SYEVX
  INTEGER, ALLOCATABLE :: ifail(:)     ! workarray for SYEVX
  REAL(c_double),EXTERNAL :: DLAMCH   ! function to specify tolerance of SYEVX
  REAL(c_double)    :: eval_inv       ! inverse of an eigenvalue


! **********************
! *** INITIALIZATION ***
! **********************

  ! init numbers
  invdim_ens = 1.0 / REAL(dim_ens)
  invdim_ensm1 = 1.0 / (REAL(dim_ens - 1))
  sqrtinvforget = SQRT(1.0 / forget)


! **********************************
! *** Compute representer vector ***
! ***                            ***
! *** We compute the ensemble of ***
! *** representer vectors b by   ***
! *** solving                    ***
! ***        T                   ***
! *** (H P H  + R) b  = y - H x  ***
! **********************************

  ! *******************************************************
  ! *** Compute mean forecasted state                   ***
  ! *** for normal EnKF using ensemble mean as forecast ***
  ! *******************************************************

  state_p = 0.0
  DO member = 1, dim_ens
     DO i = 1, dim_p
        state_p(i) = state_p(i) + invdim_ens * ens_p(i, member)
     END DO
  END DO

  ! **********************************************
  ! *** We directly compute the matrices       ***
  ! ***                                T       ***
  ! ***   HP = H P     and  HPH = H P H        ***
  ! *** as ensemble means by just projecting   ***
  ! *** the state ensemble onto observation    ***
  ! *** space. The covariance matrix is not    ***
  ! *** explicitly computed.                   ***
  ! **********************************************
  ALLOCATE(XminMean_p(dim_p, dim_ens))
  

  ! ***                             T ***
  ! *** get HP = H P and HPH = H P H  ***
  ! *** as ensemble means             ***
  ENSa: DO member = 1, dim_ens

     ! spread out state ensemble according to forgetting factor
     IF (forget /= 1.0) THEN
        ens_p(:, member) = state_p(:) &
             + (ens_p(:, member) - state_p(:)) * sqrtinvforget
     END IF

     ! initialize XMINMEAN
     XminMean_p(:, member) = ens_p(:, member) - state_p(:)

  END DO ENSa

  ! compute mean of ensemble projected on obseration space
  ALLOCATE(HXmean_p(dim_obs_p))
    
  ! Finish computation of HP and HPH
  ALLOCATE(HP_p(dim_obs, dim_p))
  ALLOCATE(HPH(dim_obs, dim_obs))
  
  CALL gemmTYPE('n', 't', dim_obs, dim_p, dim_ens, &
       invdim_ensm1, resid, dim_obs, XminMean_p, dim_p, &
       0.0, HP_p, dim_obs)

  CALL gemmTYPE('n', 't', dim_obs, dim_obs, dim_ens, &
       invdim_ensm1, resid, dim_obs, resid, dim_obs, &
       0.0, HPH, dim_obs)

  DEALLOCATE(XminMean_p)

  ! Apply localization
  CALL U_localize(step, ind_p, dim_p, dim_obs, HP_p, HPH)

  ! *** Add observation error covariance ***
  ! ***       HPH^T = (HPH + R)          ***
  CALL U_add_obs_err(step, ind_p, dim_obs, HPH)

  whichupdate: IF (rank_ana > 0) THEN
! **************************************************
! *** Update using pseudo inverse of HPH         ***
! *** by performing incomplete eigendecompostion ***
! *** and using Moore-Penrose inverse of this    ***
! *** matrix                                     ***
! **************************************************

     ! *** Initialization ***
     ALLOCATE(repres(dim_obs, dim_ens))
     ALLOCATE(eval(rank_ana))
     ALLOCATE(evec(dim_obs, rank_ana))
     ALLOCATE(evec_temp(dim_obs, rank_ana))
     ALLOCATE(rwork(8 * dim_obs))
     ALLOCATE(iwork(5 * dim_obs))
     ALLOCATE(ifail(dim_obs))
    
     ! **************************************
     ! *** compute pseudo inverse of HPH  ***
     ! *** using Moore-Penrose inverse    ***
     ! *** of rank reduced matrix         ***
     ! **************************************

     Iupper = dim_obs
     Ilower = dim_obs - rank_ana + 1
     abstol = 2 * DLAMCH('S')

     ! *** Decompose HPH = eigenvec ev eigenvec^T by   ***
     ! *** computing the RANK_ANA largest eigenvalues  ***
     ! *** and the corresponding eigenvectors          ***
     ! *** We use the LAPACK routine SYEVX             ***
     CALL syevxTYPE('v', 'i', 'u', dim_obs, HPH, &
          dim_obs, VL, VU, Ilower, Iupper, &
          abstol, nEOF, eval, evec, dim_obs, &
          rwork, 8 * dim_obs, iwork, ifail, syev_info)

     ! check if eigendecomposition was successful
     EVPok: IF (syev_info == 0) THEN
        ! Eigendecomposition OK, continue

        ! *** store V ***
        evec_temp = evec

        ! *** compute  V diag(ev^(-1)) ***
        DO j = 1, rank_ana
           eval_inv = 1.0 / eval(j)
           DO i = 1, dim_obs
              evec(i, j) = eval_inv * evec(i, j)
           END DO
        END DO
      
        ! *** compute HPH^(-1) ~ V evinv V^T ***
        ! *** HPH^(-1) is stored in HPH      ***
        CALL gemmTYPE('n', 't', dim_obs, dim_obs, rank_ana, &
             1.0, evec, dim_obs, evec_temp, dim_obs, &
             0.0, HPH, dim_obs)

        DEALLOCATE(eval, evec, evec_temp, rwork, iwork, ifail)
        
        ! ****************************************
        ! *** Compute ensemble of representer  ***
        ! *** vectors b as the product         ***
        ! ***           b = invHPH d           ***
        ! ****************************************

        CALL gemmTYPE('n', 'n', dim_obs, dim_ens, dim_obs, &
             1.0, HPH, dim_obs, resid, dim_obs, &
             0.0, repres, dim_obs)

        ! **********************************
        ! *** Update model state ensemble ***
        ! ***          a   f              ***
        ! ***         x = x + K d         ***
        ! ***********************************


        CALL gemmTYPE('t', 'n', dim_p, dim_ens, dim_obs, &
             1.0, HP_p, dim_obs, repres, dim_obs, &
             1.0, ens_p, dim_p)
 
     END IF EVPok

     ! *** Clean up ***
     DEALLOCATE(repres)


  ELSE whichupdate
! *******************************************
! *** Update using matrix HPH directly to ***      
! *** compute representer amplitudes b by ***
! *** solving HPH b = d for b.            ***
! *******************************************

     ! ****************************************
     ! *** Compute ensemble of representer  ***
     ! *** vectors b by solving             ***
     ! ***              HPH b = d           ***
     ! *** We use the LAPACK routine GESV   ***
     ! ****************************************
     ALLOCATE(ipiv(dim_obs))

     CALL gesvTYPE(dim_obs, dim_ens, HPH, dim_obs, ipiv, &
          resid, dim_obs, sgesv_info)
    

     ! *** check if solve was successful
     update: IF (sgesv_info /= 0) THEN
        WRITE (*, '(/a, 3x, a/)') 'PDAF', '!!! Problem in solve for Kalman gain !!!'
        flag = 2
     ELSE

        ! ***********************************
        ! *** Update model state ensemble ***
        ! ***    a   f         f    T     ***
        ! ***   x = x + K d = x + HP b    ***
        ! ***********************************

        CALL gemmTYPE('t', 'n', dim_p, dim_ens, dim_obs, &
             1.0, HP_p, dim_obs, resid, dim_obs, &
             1.0, ens_p, dim_p)

     END IF update
      
     DEALLOCATE(ipiv)
     
  END IF whichupdate


! ********************
! *** Finishing up ***
! ********************

  DEALLOCATE(HP_p)
  DEALLOCATE(HPH)

  IF (allocflag == 0) allocflag = 1

END SUBROUTINE lenkf_analysis_rsm
