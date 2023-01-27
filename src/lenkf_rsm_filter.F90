module mod_lenkf_rsm_filter

  use exceptions, ONLY: error_container, throw, new_exception
  use mod_assimilation_filter, ONLY: assimilation_filter
  use lenkf_rsm, ONLY: lenkf_analysis_rsm
  use, intrinsic::iso_fortran_env, ONLY: real64
  use mod_base_assimilation_manager, ONLY: base_assimilation_manager
  use util, ONLY: str

  implicit none

  type, extends(assimilation_filter)::lenkf_rsm_filter
    real(kind=8)::forget = 1.0
    real(real64), allocatable::innovations(:, :)
  contains
    procedure::assimilate
  end type lenkf_rsm_filter

contains
  function get_innovations( &
    observations, predictions, obs_errors, status) result(innovations)

    !! Compute innovations for a given subset of the model domain.

    use random

    implicit none

    ! Arguments
    real(kind=8), intent(in)::observations(:)
        !! Observation values for the subset
    real(kind=8), intent(in)::obs_errors(:)
        !! Observation errors for the subset
    real(kind=8), intent(in)::predictions(:, :)
        !! Predictions for the subset
    real(kind=8), allocatable::innovations(:, :)
        !! Innovations for the subset
    type(error_container), intent(out), optional::status
        !! Error status

    real(kind=8), allocatable::obs_perturbations(:, :)
        !! Obs perturbations

    integer::imember, iobs, obs_count, n_ensemble
    real(kind=8) mean_perturbation, perturbation_std

    obs_count = size(observations)
    n_ensemble = size(predictions, 2)

    if (size(observations) /= obs_count) then
      call throw(status, new_exception( &
                 'Observations array has wrong length. Expected '// &
                 str(obs_count)//', got '//str(size(observations))//'.', &
                 'get_innovations'))
      return
    end if

    if (size(predictions, 1) /= obs_count .or. &
        size(predictions, 2) /= n_ensemble) then
      call throw(status, new_exception( &
                 'Predictions array has wrong shape. Expected ('// &
                 str(obs_count)//','//str(n_ensemble)//'),&
                 & got ('// &
                 str(size(predictions, 1))//',' &
                 //str(size(predictions, 2))//').', &
                 'get_innovations'))
      return
    end if

    allocate (innovations(obs_count, n_ensemble))
    allocate (obs_perturbations(obs_count, n_ensemble))

    ! Compute obs perturbations
    do imember = 1, n_ensemble
      do iobs = 1, obs_count
        obs_perturbations(iobs, imember) = random_normal()*obs_errors(iobs)
      end do
    end do

    ! Subtract off the mean from obs perturbations (so the new mean perturbation
    ! is zero)
    do iobs = 1, obs_count
      mean_perturbation = sum(obs_perturbations(iobs, :))/n_ensemble
      obs_perturbations(iobs, :) = obs_perturbations(iobs, :) - mean_perturbation
    end do

    ! Scale the perturbations so their standard deviation equals the
    ! observation error
    do iobs = 1, obs_count
      perturbation_std = &
        sqrt(sum(obs_perturbations(iobs, :)**2)/(n_ensemble - 1))
      obs_perturbations(iobs, :) = &
        obs_perturbations(iobs, :)*obs_errors(iobs)/perturbation_std

    end do

    do imember = 1, n_ensemble
      do iobs = 1, obs_count
        innovations(iobs, imember) = observations(iobs) &
                                     - predictions(iobs, imember) &
                                     + obs_perturbations(iobs, imember)
      end do
    end do

  end function get_innovations

  subroutine assimilate( &
    this, ibatch, dim_p, dim_obs, dim_ens, &
    ens_p, predictions, observations, obs_errors, &
    mgr, status)

    class(lenkf_rsm_filter) :: this
    integer, intent(in)::ibatch
    integer, intent(in)::dim_p
    integer, intent(in)::dim_obs
    integer, intent(in)::dim_ens
    real(kind=8), intent(inout)::ens_p(dim_p, dim_ens)
    real(kind=8), intent(in)::predictions(dim_obs, dim_ens)
    real(kind=8), intent(in)::observations(dim_obs)
    real(kind=8), intent(in)::obs_errors(dim_obs)
    real(kind=8), allocatable::innovations(:, :)

    class(base_assimilation_manager)::mgr
    type(error_container), intent(out), optional :: status

    real(kind=8)::state_p(dim_p)
    integer::flag

    innovations = get_innovations(observations, predictions, obs_errors)
    this%innovations = innovations

    call lenkf_analysis_rsm( &
      ibatch, dim_p, dim_obs, dim_obs, &
      dim_ens, int(0), state_p, ens_p, predictions, innovations, &
      this%forget, flag, mgr)

    if (flag == 2) then
      call throw(status, new_exception( &
           'lenkf_analysis_rsm encountered a problem in solve for &
           &Kalman gain', 'assimilate'))
      return
    else if (flag /= 0) then
      call throw(status, new_exception( &
                 'Error in lenkf_analysis_rsm', 'assimilate'))
      return
    end if

  end subroutine assimilate

end module mod_lenkf_rsm_filter
