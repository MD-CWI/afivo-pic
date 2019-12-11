module m_particles
  use m_particle_core
  use m_af_all
  use m_globals

  implicit none
  public

  ! Magnetic field vector
  real(dp), protected :: B_vec(3) = [0.0_dp, 0.0_dp, 0.0_dp]

  real(dp) :: min_merge_increase = 1.25_dp

  real(dp), parameter :: array_incr_fac = 1.25_dp

  real(dp), protected :: steps_per_period = 30.0_dp

  real(dp) :: bfield_delay = 0.0_dp
  real(dp) :: bfield_onset_duration = 0.0_dp

contains

  subroutine init_particle(cfg, pc)
    use m_units_constants
    use m_gas
    use m_cross_sec
    use m_config
    use m_domain
    use m_field
    type(CFG_t), intent(inout) :: cfg
    type(PC_t), intent(inout) :: pc

    integer                        :: nn, tbl_size, max_num_part
    integer                        :: n_gas_comp, n_gas_frac
    real(dp)                       :: temperature, max_ev
    character(len=200)             :: cs_file
    character(len=20), allocatable :: gas_names(:)
    real(dp), allocatable          :: gas_fracs(:)
    type(CS_t), allocatable        :: cross_secs(:)

    ! Gas parameters
    call CFG_add(cfg, "gas%temperature", 300.0_dp, &
         "The gas temperature (Kelvin)")
    call CFG_add(cfg, "gas%components", (/"N2"/), &
         "The names of the gases used in the simulation", .true.)
    call CFG_add(cfg, "gas%file", "input/cs_example.txt", &
         "The file in which to find cross section data")
    call CFG_add(cfg, "gas%fractions", (/1.0_dp /), &
         & "The partial pressure of the gases (as if they were ideal gases)", .true.)

    ! Particle model related parameters
    call CFG_add_get(cfg, "particle%per_cell", particle_per_cell, &
         "Desired number of particles per cell")
    call CFG_add_get(cfg, "particle%min_weight", particle_min_weight, &
         "Minimum weight for simulation particles")
    call CFG_add_get(cfg, "particle%max_weight", particle_max_weight, &
         "Maximum weight for simulation particles")

    call CFG_add(cfg, "particle%max_number", 50*1000*1000, &
         "Maximum number of particles")
    call CFG_add_get(cfg, "particle%min_merge_increase", min_merge_increase, &
         "Minimum increase in particle count before merging")
    call CFG_add(cfg, "particle%lkptbl_size", 1000, &
         "The size of the lookup table for the collision rates")
    call CFG_add(cfg, "particle%max_energy_ev", 500.0_dp, &
         "The maximum energy in eV for particles in the simulation")
    call CFG_add_get(cfg, "boris_steps_per_period", steps_per_period, &
         "The minimum number of time steps per gyration")

    call CFG_get_size(cfg, "gas%components", n_gas_comp)
    call CFG_get_size(cfg, "gas%fractions", n_gas_frac)
    if (n_gas_comp /= n_gas_frac) &
         error stop "gas%components and gas%fractions have unequal size"
    allocate(gas_names(n_gas_comp))
    allocate(gas_fracs(n_gas_comp))

    call CFG_get(cfg, "gas%components", gas_names)
    call CFG_get(cfg, "gas%fractions", gas_fracs)
    call CFG_get(cfg, "gas%temperature", temperature)

    call CFG_get(cfg, "gas%file", cs_file)
    call CFG_get(cfg, "particle%max_energy_ev", max_ev)
    call CFG_get(cfg, "particle%max_number", max_num_part)
    call CFG_add_get(cfg, "magnetic_field", B_vec, &
         "Magnetic field vector (T)")
    call CFG_add_get(cfg, "magnetic_field_delay", bfield_delay, &
         "Time (s) before applying the magnetic field")
    call CFG_add_get(cfg, "magnetic_field_onset_duration", &
         bfield_onset_duration, &
         "Time (s) during which B goes from zero to final value")

    ! Initialize gas and electric field module
    call GAS_initialize(gas_names, gas_fracs, ST_gas_pressure, temperature)

    do nn = 1, n_gas_comp
       call CS_add_from_file(trim(cs_file), &
            trim(gas_names(nn)), gas_fracs(nn) * &
            GAS_number_dens, max_ev, cross_secs)
    end do

    call CS_write_summary(cross_secs, &
         trim(ST_output_dir) // "/" // trim(ST_simulation_name) // "_cs_summary.txt")

    call CFG_get(cfg, "particle%lkptbl_size", tbl_size)

    call pc%initialize(UC_elec_mass, max_num_part, &
         rng_seed=get_random_seed())
    call pc%use_cross_secs(max_ev, tbl_size, cross_secs)

    pc%accel_function => get_accel
    pc%outside_check => outside_check

    if (any(abs(B_vec) > 0)) then
       print *, "Using Boris mover"
       pc%B_vec = B_vec
       pc%particle_mover => PC_boris_advance
       pc%after_mover => PC_after_dummy
       pc%dt_max = 2 * UC_pi / (steps_per_period * &
            abs(UC_elec_q_over_m) * norm2(B_vec))
       write(*, "(A,E10.2)") " Max time step Boris (s): ", pc%dt_max
    else
       print *, "Using Verlet mover"
    end if

    where (pc%colls(:)%type == CS_ionize_t .or. &
         pc%colls(:)%type == CS_attach_t)
       pc%coll_is_event(:) = .true.
    end where

  end subroutine init_particle

  subroutine update_bfield(time)
    use m_units_constants
    real(dp), intent(in) :: time

    if (time < Bfield_delay) then
       pc%B_vec = 0.0_dp
    else if (time < Bfield_delay + Bfield_onset_duration) then
       pc%B_vec = B_vec * (time - Bfield_delay) / Bfield_onset_duration
    else
       pc%B_vec = B_vec
    end if

    pc%dt_max = 2 * UC_pi / (steps_per_period * &
         abs(UC_elec_q_over_m) * norm2(pc%B_vec))

  end subroutine update_bfield

  subroutine adapt_weights(tree, pc)
    type(af_t), intent(in)   :: tree
    type(PC_t), intent(inout) :: pc
    integer                   :: id, n_part_id
    integer, allocatable      :: id_ipart(:)

    call sort_by_id(tree, pc, id_ipart)
    ! print *, "before: ", pc%get_num_sim_part(), pc%get_num_real_part()
    !$omp parallel do private(id, n_part_id) schedule(dynamic)
    do id = 1, tree%highest_id
       n_part_id = id_ipart(id+1) - id_ipart(id)
       if (n_part_id > 0) then
          call pc%merge_and_split_range(id_ipart(id), id_ipart(id+1)-1, &
               (/.true., .true., .false./), 1e-12_dp, &
               .true., get_desired_weight, 1.0e10_dp, &
               PC_merge_part_rxv, PC_split_part)
       end if
    end do
    !$omp end parallel do

    call pc%clean_up()
    ! print *, "after:  ", pc%get_num_sim_part(), pc%get_num_real_part()
  end subroutine adapt_weights

  subroutine sort_by_id(tree, pc, id_ipart)
    type(af_t), intent(in)              :: tree
    type(PC_t), intent(inout)           :: pc
    integer, intent(inout), allocatable :: id_ipart(:)

    integer                      :: n, id, new_ix
    integer, allocatable         :: id_count(:)
    integer, allocatable         :: id_ix(:)
    type(PC_part_t), allocatable, save :: p_copy(:)

    if (allocated(id_ipart)) deallocate(id_ipart)
    allocate(id_ipart(tree%highest_id+1))
    allocate(id_ix(tree%highest_id+1))
    allocate(id_count(tree%highest_id))

    if (.not. allocated(p_copy)) then
       allocate(p_copy(nint(pc%n_part * array_incr_fac)))
    else if (size(p_copy) < pc%n_part) then
       deallocate(p_copy)
       allocate(p_copy(nint(pc%n_part * array_incr_fac)))
    end if

    id_count(:) = 0

    do n = 1, pc%n_part
       id           = pc%particles(n)%id
       id_count(id) = id_count(id) + 1
       p_copy(n)    = pc%particles(n)
    end do

    id_ix(1) = 1
    do id = 2, tree%highest_id+1
       id_ix(id) = id_ix(id-1) + id_count(id-1)
    end do

    id_ipart(:) = id_ix(:)

    do n = 1, pc%n_part
       id                   = p_copy(n)%id
       new_ix               = id_ix(id)
       id_ix(id)            = id_ix(id) + 1
       pc%particles(new_ix) = p_copy(n)
    end do
  end subroutine sort_by_id

  subroutine particles_to_density_and_events(tree, pc, init_cond)
    use m_cross_sec
    type(af_t), intent(inout)        :: tree
    type(PC_t), intent(inout)        :: pc
    logical, intent(in)              :: init_cond

    integer                     :: n, i, n_part
    real(dp), allocatable, save :: coords(:, :)
    real(dp), allocatable, save :: weights(:)
    integer, allocatable, save  :: id_guess(:)

    n_part = pc%get_num_sim_part()
    n = max(n_part, pc%n_events)

    if (.not. allocated(weights)) then
       n = nint(n * array_incr_fac)
       allocate(coords(NDIM, n))
       allocate(weights(n))
       allocate(id_guess(n))
    else if (size(weights) < n) then
       n = nint(n * array_incr_fac)
       deallocate(coords)
       deallocate(weights)
       deallocate(id_guess)
       allocate(coords(NDIM, n))
       allocate(weights(n))
       allocate(id_guess(n))
    end if

    !$omp parallel do
    do n = 1, n_part
       coords(:, n) = pc%particles(n)%x(1:NDIM)
       weights(n) = pc%particles(n)%w
       id_guess(n) = pc%particles(n)%id
    end do
    !$omp end parallel do

    if (init_cond) then
       call af_tree_clear_cc(tree, i_electron)
       call af_particles_to_grid(tree, i_electron, coords(:, 1:n_part), &
            weights(1:n_part), n_part, 1, id_guess(1:n_part))
       call af_particles_to_grid(tree, i_pos_ion, coords(:, 1:n_part), &
            weights(1:n_part), n_part, 1, id_guess(1:n_part))
    else
       call af_tree_clear_cc(tree, i_electron)
       call af_particles_to_grid(tree, i_electron, coords(:, 1:n_part), &
            weights(1:n_part), n_part, 1, id_guess(1:n_part))
    end if

    pc%particles(1:n_part)%id = id_guess(1:n_part)

    i = 0
    do n = 1, pc%n_events
       if (pc%event_list(n)%ctype == CS_ionize_t) then
          i = i + 1
          coords(:, i) = pc%event_list(n)%part%x(1:NDIM)
          weights(i) = pc%event_list(n)%part%w
          id_guess(i) = pc%event_list(n)%part%id
       else if (pc%event_list(n)%ctype == CS_attach_t) then
          i = i + 1
          coords(:, i) = pc%event_list(n)%part%x(1:NDIM)
          weights(i) = -pc%event_list(n)%part%w
          id_guess(i) = pc%event_list(n)%part%id
       end if
    end do

    if (i > 0) then
       call af_particles_to_grid(tree, i_pos_ion, coords(:, 1:i), &
            weights(1:i), i, 1, id_guess(1:i))
    end if

    pc%n_events = 0

  end subroutine particles_to_density_and_events

  function get_accel_pos(x) result(accel)
    use m_units_constants
    real(dp), intent(in) :: x(NDIM)
    real(dp)             :: accel(3)
    logical              :: success

#if NDIM == 2
    accel(1:NDIM) = af_interp1(tree, x, [i_Ex, i_Ey], success)
    accel(3) = 0.0_dp
#elif NDIM == 3
    accel(1:NDIM) = af_interp1(tree, x, [i_Ex, i_Ey, i_Ez], success)
#endif
    accel(:) = accel(:) * UC_elec_q_over_m
  end function get_accel_pos

  function get_accel(my_part) result(accel)
    use m_units_constants
    type(PC_part_t), intent(inout) :: my_part
    real(dp)                       :: accel(3)
    accel = get_accel_pos(my_part%x(1:NDIM))
  end function get_accel

  function get_desired_weight(my_part) result(weight)
    type(PC_part_t), intent(in) :: my_part
    real(dp)                    :: weight, n_elec
    type(af_loc_t)              :: loc
    integer                     :: id, ix(NDIM)

    loc = af_get_loc(tree, my_part%x(1:NDIM), my_part%id)
    id = loc%id
    ix = loc%ix

#if NDIM == 2
    n_elec = tree%boxes(id)%cc(ix(1), ix(2), i_electron) * &
         product(tree%boxes(id)%dr)
#elif NDIM == 3
    n_elec = tree%boxes(id)%cc(ix(1), ix(2), ix(3), i_electron) * &
         product(tree%boxes(id)%dr)
#endif

    weight = n_elec / particle_per_cell
    weight = max(particle_min_weight, min(particle_max_weight, weight))
  end function get_desired_weight

end module m_particles
