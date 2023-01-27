module single_obs_test

  use mod_base_assimilation_manager, ONLY: base_assimilation_manager
  use exceptions, ONLY: throw, new_exception, error_container

  implicit none

  type, extends(base_assimilation_manager)::simple_assimilation_manager
   contains
    procedure::localize
    procedure::add_obs_err
  end type simple_assimilation_manager

contains

  SUBROUTINE localize(this, ibatch, dim_p, dim_obs, HP_p, HPH, status)
    ! Apply localization to HP and HPH^T
    USE iso_c_binding
    class(simple_assimilation_manager), target::this
    INTEGER(c_int32_t), INTENT(in), value :: ibatch, dim_p, dim_obs
    REAL(c_double), INTENT(inout) :: HP_p(dim_obs, dim_p)
    REAL(c_double), INTENT(inout) :: HPH(dim_obs, dim_obs)
    type(error_container), intent(out), optional :: status

  end SUBROUTINE localize

  SUBROUTINE add_obs_err(this, ibatch, dim_obs, HPH, status)
    ! Add observation error covariance matrix
    USE iso_c_binding
    class(simple_assimilation_manager), target::this
    INTEGER(c_int32_t), INTENT(in), value :: ibatch, dim_obs
    REAL(c_double), INTENT(inout) :: HPH(dim_obs, dim_obs)
    type(error_container), intent(out), optional :: status

    integer::iobs

    do iobs = 1, dim_obs
      HPH(iobs, iobs) = HPH(iobs, iobs) + 0.1**2
    end do


  end SUBROUTINE add_obs_err
    
  subroutine test_single_obs()

    use random
    use mod_lenkf_rsm_filter, only: lenkf_rsm_filter
    use mod_assimilation_manager, ONLY: &
         assimilation_manager
    use iso_fortran_env, only: real64, error_unit
    use util, ONLY: str
    use localization, ONLY: base_localizer
    use system_mpi
  use random_observations, ONLY: &
    random_observation_set, new_random_observation_set

    integer, parameter::ibatch = 1
    integer, parameter::n_obs = 1
    integer, parameter:: n_ens = 15
    integer, parameter::batch_size = 1
    integer, parameter::istep = 1
    integer, parameter::report_interval=1
    real(real64)::observations(n_obs)
    real(real64)::obs_errors(n_obs)
    real(real64)::predictions(n_obs, n_ens)
    real(real64)::state(n_obs, n_ens)
    real(real64)::prior_state(n_obs, n_ens)
    real(real64)::posterior_state(n_obs, n_ens)
    type(lenkf_rsm_filter)::filter
    type(simple_assimilation_manager)::mgr
    real(real64)::prior_mean, posterior_mean

    integer::imember

    real(real64)::prior_innovations(n_obs, n_ens)
    real(real64)::prior_variance, expected_increment(n_obs, n_ens), actual_increment(n_obs, n_ens), mean_perturbation
    real(real64)::obs_perturbations(n_obs,n_ens)

    observations = 1
    obs_errors = 0.1
    
    do imember = 1, n_ens

       predictions(1,imember) = random_normal() + 2
       prior_innovations(:,imember) = observations - predictions(:,imember)

    end do

    prior_mean = sum(predictions)/n_ens

    ! Extremely simple test, model state is identical to predictions
    state = predictions

    prior_state = state

    call filter%assimilate(ibatch, n_obs, n_obs, n_ens, state, predictions, observations, obs_errors, mgr)

    posterior_state = state

    posterior_mean = sum(state)/n_ens

    if (posterior_mean <= min(prior_mean, observations(1)) .or. &
         posterior_mean >= max(prior_mean, observations(1))) then

       write (error_unit, *) "Posterior mean "// &
          str(real(posterior_mean), '(F0.2)')// &
          " out of range (should be between "// &
          str(real(min(prior_mean, observations(1))), '(F0.2)')//" and "// &
          str(real(max(prior_mean, observations(1))), '(F0.2)')//")."
        error stop
    end if

    prior_variance = sum((predictions - prior_mean)**2)/(n_ens-1)
    obs_perturbations = filter%innovations - prior_innovations

    do imember = 1, n_ens
       expected_increment(1,imember) = &
         prior_variance/(prior_variance + obs_errors(1)**2)*filter%innovations(1,imember)
    end do

    actual_increment = posterior_state - prior_state

    if (any(abs((actual_increment - expected_increment)/expected_increment) > 1e-15)) then
       write (error_unit, *) "Wrong analysis increment. Expected",&
            expected_increment,"; got ", &
            actual_increment,"(error", actual_increment - expected_increment,")"
       error stop
    end if

  end subroutine test_single_obs

end module single_obs_test

program run_single_obs_test

  use single_obs_test, only: test_single_obs
  use system_mpi

  implicit none

  integer::ierr

  call mpi_init(ierr)

  call test_single_obs()

  call mpi_finalize(ierr)

end program run_single_obs_test
