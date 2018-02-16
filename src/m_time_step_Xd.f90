module m_time_step_$Dd
  use m_a$D_all
  use m_particle_core
  use m_globals_$Dd

  implicit none
  public

  ! Current time step
  real(dp) :: ST_dt

  ! Maximum allowed time step
  real(dp), protected :: ST_dt_max = 1.0e-10_dp

  ! Minimum allowed time step
  real(dp), protected :: ST_dt_min = 1.0e-14_dp

  ! Time between writing output
  real(dp), protected :: ST_dt_output = 1.0e-10_dp

contains

  subroutine time_step_init(cfg)
    use m_config
    type(CFG_t), intent(inout) :: cfg

    call CFG_add_get(cfg, "dt_output", ST_dt_output, &
         "The timestep for writing output (s)")
    call CFG_add_get(cfg, "dt_max", ST_dt_max, &
         "The maximum timestep (s)")
    call CFG_add_get(cfg, "dt_min", ST_dt_min, &
         "The minimum timestep (s)")

  end subroutine time_step_init

  !> This routine can be used to estimate the error in the electric field. It
  !> should first be called with store_samples = .true., then a later call with
  !> store_samples = .false. sets the error in fld_err.
  subroutine PM_fld_error(tree, pc, rng, n_samples, fld_err, store_samples)
    use m_random
    type(a$D_t), intent(in)    :: tree
    class(PC_t), intent(in)    :: pc
    type(RNG_t), intent(inout) :: rng
    integer, intent(in)        :: n_samples
    real(dp), intent(out)      :: fld_err
    logical, intent(in)        :: store_samples

    integer                     :: i, n
    real(dp)                    :: fld($D), this_err, avg_norm
    real(dp), allocatable, save :: pos_samples(:, :)
    real(dp), allocatable, save :: fld_samples(:, :)
    logical, save               :: have_samples = .false.

    if (store_samples) then
       ! Store samples of the electric field

       if (pc%n_part < 1) then
          ! Not enough particles to store samples
          have_samples = .false.
       else
          have_samples = .true.
          if (allocated(fld_samples)) deallocate(fld_samples)
          if (allocated(pos_samples)) deallocate(pos_samples)
          allocate(fld_samples($D, n_samples))
          allocate(pos_samples($D, n_samples))

          do i = 1, n_samples
             ! Randomly select a particle
             n = floor(rng%unif_01() * pc%n_part) + 1
             pos_samples(:,i) = pc%particles(n)%x(1:$D)
#if $D == 2
             fld_samples(:,i) = a$D_interp1(tree, pos_samples(:,i), &
                  [i_Ex, i_Ey], $D)
#elif $D == 3
             fld_samples(:,i) = a$D_interp1(tree, pos_samples(:,i), &
                  [i_Ex, i_Ey, i_Ez], $D)
#endif
          end do
       end if

    else if (.not. store_samples) then
       ! Compute the error
       fld_err = 0

       if (have_samples) then
          avg_norm = norm2(fld_samples) / sqrt(1.0_dp * n_samples)

          do i = 1, n_samples
#if $D == 2
             fld = a$D_interp1(tree, pos_samples(:,i), &
                  [i_Ex, i_Ey], $D)
#elif $D == 3
             fld = a$D_interp1(tree, pos_samples(:,i), &
                  [i_Ex, i_Ey, i_Ez], $D)
#endif
             this_err = norm2(fld - fld_samples(:,i)) / &
                  max(norm2(fld), norm2(fld_samples(:,i)), avg_norm)
             fld_err  = max(fld_err, this_err)
          end do
       else
          error stop "PM_fld_error: no samples available"
       end if
    end if

  end subroutine PM_fld_error

  function PM_get_max_dt(pc, rng, n_samples, cfl_num) result(dt_max)
    use m_random
    use m_mrgrnk
    class(PC_t), intent(in)    :: pc
    type(RNG_t), intent(inout) :: rng
    integer, intent(in)        :: n_samples
    real(dp), intent(in)       :: cfl_num
    real(dp)                   :: dt_max
    real(dp)                   :: vel_est, min_dr
    real(dp), allocatable      :: velocities(:)
    integer, allocatable       :: ix_list(:)
    integer                    :: n, ix

    allocate(velocities(n_samples))
    allocate(ix_list(n_samples))

    if (pc%n_part > 0) then
       ! Estimate maximum velocity of particles
       do n = 1, n_samples
          ix = floor(rng%unif_01() * pc%n_part) + 1
          velocities(n) = norm2(pc%particles(ix)%v)
       end do

       call mrgrnk(velocities, ix_list)
       velocities = velocities(ix_list)

       vel_est = velocities(nint(n_samples * 0.9_dp))

       ! Get smallest grid delta
       min_dr = a$D_min_dr(tree)
       dt_max = cfl_num * min_dr / max(vel_est, 1e-10_dp)
    else
       dt_max = huge(1.0_dp)
    end if

  end function PM_get_max_dt

  !> Given the old stepsize 'old_dt', the error 'err', the maximum allowed error
  !> 'max_err', and the relative change in number of particles, return a new
  !> stepsize.
  real(dp) function get_new_dt(old_dt, err, max_err)
    real(dp), intent(in) :: old_dt
    real(dp), intent(in) :: err, max_err

    ! Adjust the timesteps to obtain an error close to the desired one. The
    ! interpolation error should be ~ interp_errFrac * max_fld_err
    if (err < 0.5D0 * max_err) then
       get_new_dt = min(2 * old_dt, &
            old_dt * (max_err/(err+epsilon(err)))**0.1D0)
    else if (err > max_err) then
       get_new_dt = old_dt * (max_err/err)
    else
       get_new_dt = old_dt
    end if

    get_new_dt = min(max(get_new_dt, ST_dt_min), ST_dt_max)

  end function get_new_dt

end module m_time_step_$Dd
