module m_particles
  use m_particle_core
  use m_af_all
  use m_globals
  use m_photons

  implicit none
  public

  real(dp) :: min_merge_increase = 1.25_dp

  !> Maximum number of iterations between calling adapt_weights
  integer :: iterations_between_merge_split = 1000

  real(dp), parameter :: array_incr_fac = 1.25_dp

  real(dp), protected :: steps_per_period = 30.0_dp

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
    call CFG_add(cfg, "gas%components", ["N2"], &
         "The names of the gases used in the simulation", .true.)
    call CFG_add(cfg, "gas%file", "input/cs_example.txt", &
         "The file in which to find cross section data")
    call CFG_add(cfg, "gas%fractions", [1.0_dp ], &
         & "The partial pressure of the gases (as if they were ideal gases)", .true.)

    ! Particle model related parameters
    call CFG_add_get(cfg, "particle%per_cell", particle_per_cell, &
         "Desired number of particles per cell")
    call CFG_add_get(cfg, "particle%min_weight", particle_min_weight, &
         "Minimum weight for simulation particles")
    call CFG_add_get(cfg, "particle%max_weight", particle_max_weight, &
         "Maximum weight for simulation particles")

    call CFG_add(cfg, "particle%max_number", 15*1000*1000, &
         "Maximum number of particles")
    call CFG_add_get(cfg, "particle%min_merge_increase", min_merge_increase, &
         "Minimum increase in particle count before merging")

    if (GL_cylindrical) then
       ! Adapt weights more frequently to prevent fluctuations near the axis
       iterations_between_merge_split = 10
    end if
    call CFG_add_get(cfg, "particle%iterations_between_merge_split", &
         iterations_between_merge_split, &
         "Maximum number of iterations between calling adapt_weights")

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

    ! Initialize gas and electric field module
    call GAS_initialize(gas_names, gas_fracs, GL_gas_pressure, temperature)

    do nn = 1, n_gas_comp
       call CS_add_from_file(trim(cs_file), &
            trim(gas_names(nn)), gas_fracs(nn) * &
            GAS_number_dens, max_ev, cross_secs)
    end do

    call CS_write_summary(cross_secs, &
         trim(GL_output_dir) // "/" // trim(GL_simulation_name) // "_cs_summary.txt")

    call CFG_get(cfg, "particle%lkptbl_size", tbl_size)

    pc%particle_mover => PC_verlet_advance
    pc%after_mover => PC_verlet_correct_accel

    call pc%initialize(UC_elec_mass, max_num_part)
    call pc%use_cross_secs(max_ev, tbl_size, cross_secs)

    pc%accel_function => get_accel
    pc%outside_check => outside_check

    where (pc%colls(:)%type == CS_ionize_t .or. &
         pc%colls(:)%type == CS_attach_t)
       pc%coll_is_event(:) = .true.
    end where

  end subroutine init_particle

  !> Adjust the weights of the particles
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
               [.true., .true., .false.], 1e-12_dp, &
               .true., get_desired_weight, 1.0e10_dp, &
               PC_merge_part_rxv, PC_split_part)
       end if
    end do
    !$omp end parallel do

    call pc%clean_up()
    ! print *, "after:  ", pc%get_num_sim_part(), pc%get_num_real_part()
  end subroutine adapt_weights

  !> Sort the particles by their id
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
    use m_particle_core
    use m_domain, only: outside_check
    type(af_t), intent(inout)        :: tree
    type(PC_t), intent(inout)        :: pc
    logical, intent(in)              :: init_cond

    integer                     :: n, i, n_part, n_photons
    ! type(PC_part_t)             :: new_part
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
       coords(:, n) = get_coordinates(pc%particles(n))
       weights(n) = pc%particles(n)%w
       id_guess(n) = pc%particles(n)%id
    end do
    !$omp end parallel do

    if (init_cond) then
       call af_tree_clear_cc(tree, i_electron)
       call af_particles_to_grid(tree, i_electron, coords(:, 1:n_part), &
            weights(1:n_part), n_part, interpolation_order_to_density, &
            id_guess(1:n_part), iv_tmp=i_tmp_dens)
       call af_particles_to_grid(tree, i_pos_ion, coords(:, 1:n_part), &
            weights(1:n_part), n_part, interpolation_order_to_density, &
            id_guess(1:n_part), iv_tmp=i_tmp_dens)
    else
       call af_tree_clear_cc(tree, i_electron)
       call af_particles_to_grid(tree, i_electron, coords(:, 1:n_part), &
            weights(1:n_part), n_part, interpolation_order_to_density, &
            id_guess(1:n_part), iv_tmp=i_tmp_dens)
    end if

    pc%particles(1:n_part)%id = id_guess(1:n_part)

    i = 0
    do n = 1, pc%n_events
       if (pc%event_list(n)%ctype == CS_ionize_t) then
          i = i + 1
          coords(:, i) = get_coordinates(pc%event_list(n)%part)
          weights(i) = pc%event_list(n)%part%w
          id_guess(i) = pc%event_list(n)%part%id
       else if (pc%event_list(n)%ctype == CS_attach_t) then
          i = i + 1
          coords(:, i) = get_coordinates(pc%event_list(n)%part)
          weights(i) = -pc%event_list(n)%part%w
          id_guess(i) = pc%event_list(n)%part%id

       else if (pc%event_list(n)%ctype == PC_particle_went_out .and. &
            pc%event_list(n)%cIx == inside_dielectric) then
          ! Now we map the particle to surface charge
          call particle_to_surface_charge(tree, pc%event_list(n)%part, &
               i_surf_elec)
       end if
    end do

    if (i > 0) then ! only for the events that created an ion
       call af_particles_to_grid(tree, i_pos_ion, coords(:, 1:i), &
            weights(1:i), i, interpolation_order_to_density, &
            id_guess(1:i), iv_tmp=i_tmp_dens)
    end if

    if (photoi_enabled) then
      call photoionization(tree, pc, coords, weights, n_photons)
      call af_particles_to_grid(tree, i_pos_ion, coords(:, 1:n_photons), &
           weights(1:n_photons), n_photons, interpolation_order_to_density, &
           id_guess(1:n_photons), iv_tmp=i_tmp_dens)
    end if
    if (photoe_enabled) then
      call photoemission(tree, pc)
    end if

    pc%n_events = 0

    if (GL_use_dielectric) then
       ! Map densities inside the first dielectric layers to surface charge.
       ! Electrons can still move away from the dielectric.
       call dielectric_inside_layer_to_surface(tree, diel, i_electron, &
            i_surf_elec_close, 1.0_dp, clear_cc=.true., clear_surf=.true.)
       call dielectric_inside_layer_to_surface(tree, diel, i_pos_ion, &
            i_surf_pos_ion, 1.0_dp, clear_cc=.true., clear_surf=.false.)
       ! Sum densities together
       call dielectric_set_weighted_sum(diel, i_surf_sum_dens, &
            [i_surf_elec, i_surf_elec_close, i_surf_pos_ion], &
            [-1.0_dp, -1.0_dp, 1.0_dp])
    end if

  end subroutine particles_to_density_and_events

  subroutine particle_to_surface_charge(tree, my_part, i_surf)
    ! Input: a particle that is found in the dielectric (after timestep, and is
    ! thus flagged for removal) This particle will be mapped to the surface
    ! charge of the corresponding cell
    use m_units_constants

    type(af_t), intent(in)      :: tree
    type(PC_part_t), intent(in) :: my_part
    integer, intent(in)         :: i_surf !< Surface variable
    integer                     :: ix_surf, ix_cell(NDIM-1)

    call dielectric_get_surface_cell(tree, diel, get_coordinates(my_part), &
         ix_surf, ix_cell)
    ! Update the charge in the surface cell
#if NDIM == 2
    diel%surfaces(ix_surf)%sd(ix_cell(1), i_surf) = &
         diel%surfaces(ix_surf)%sd(ix_cell(1), i_surf) + &
         my_part%w / diel%surfaces(ix_surf)%dr(1)
#elif NDIM == 3
!FLAG
    diel%surfaces(ix_surf)%sd(ix_cell(1), ix_cell(2), i_surf) = &
         diel%surfaces(ix_surf)%sd(ix_cell(1), ix_cell(2), i_surf) + &
         my_part%w / product(diel%surfaces(ix_surf)%dr)
#endif
  end subroutine particle_to_surface_charge

  function get_accel(my_part) result(accel)
    use m_units_constants
    type(PC_part_t), intent(inout) :: my_part
    real(dp)                       :: accel(3), coord(NDIM)
    logical                        :: success

    ! Set acceleration for extra dimensions to zero
    accel(NDIM+1:) = 0.0_dp
    coord = get_coordinates(my_part)

    ! Interpolation of face-centered fields
    accel(1:NDIM) = af_interp1_fc(tree, coord, ifc_E, &
         success, id_guess=my_part%id) * UC_elec_q_over_m

    if (GL_cylindrical) then
       ! Convert back to xyz coordinates
       accel([1, 3]) = accel(1) * my_part%x([1, 3]) / coord(1)
    end if

  end function get_accel

  function get_desired_weight(my_part) result(weight)
    type(PC_part_t), intent(in) :: my_part
    real(dp)                    :: weight, n_elec, coord(NDIM)
    type(af_loc_t)              :: loc
    integer                     :: id, ix(NDIM)

    coord = get_coordinates(my_part)
    loc = af_get_loc(tree, coord, my_part%id)
    id = loc%id
    ix = loc%ix

#if NDIM == 2
    if (GL_cylindrical) then
       n_elec = tree%boxes(id)%cc(ix(1), ix(2), i_electron) * &
            af_cyl_volume_cc(tree%boxes(id), ix(1))
    else
       n_elec = tree%boxes(id)%cc(ix(1), ix(2), i_electron) * &
            product(tree%boxes(id)%dr)
    end if
#elif NDIM == 3
    n_elec = tree%boxes(id)%cc(ix(1), ix(2), ix(3), i_electron) * &
         product(tree%boxes(id)%dr)
#endif

    weight = n_elec / particle_per_cell
    weight = max(particle_min_weight, min(particle_max_weight, weight))
  end function get_desired_weight

  function write_EEDF_as_curve(pc) result(curve_dat)
    !> Make a histogram of electron energies and save it pass a curve-object (can be added to Silo-file)
    type(PC_t), intent(in)  :: pc
    integer                 :: i, num_bins
    real(dp), allocatable   :: bins(:), bin_values(:)
    real(dp)                :: n_part, max_en, en_step = 0.75
    real(dp), allocatable   :: curve_dat(:, :, :)

    n_part = pc%get_num_real_part()
    max_en = get_max_energy(pc)
    num_bins = ceiling(max_en / en_step) ! Generate bins of en_step (eV) each up until max energy
    if (num_bins < 1) num_bins = 1 ! At least one bin

    allocate(bins(num_bins))
    allocate(bin_values(num_bins))
    allocate(curve_dat(1, 2, num_bins))

    do i = 1, num_bins
      bins(i) = en_step * (i-1)
    end do

    call pc%histogram(calc_elec_energy, is_alive, [0.0_dp], bins, bin_values)
    ! Convert histogram to density and save as curve-object
    curve_dat(1, 1, :) = bins
    curve_dat(1, 2, :) = bin_values/(n_part * en_step) + 1e-9 ! Add regularization parameter to prevent errors when converting to semilogy plots
  end function write_EEDF_as_curve

  function calc_elec_energy(part) result(energy)
    use m_units_constants
    type(PC_part_t), intent(in) :: part
    real(dp)       :: energy

    energy = PC_v_to_en(part%v, UC_elec_mass) / UC_elec_volt
  end function calc_elec_energy

  logical function is_alive(part, real_args) result(alive)
    type(PC_part_t), intent(in) :: part
    real(dp), intent(in)        :: real_args(:) ! This basically does nothing but it seems non-optional

    alive = (part%w > 0.0_dp)
  end function is_alive

  function get_max_energy(pc) result(max_en)
    use m_units_constants
    type(PC_t), intent(in)  :: pc
    real(dp)                :: max_en
    integer                 :: ll
    max_en = 0.0_dp
    do ll = 1, pc%n_part
      if (pc%particles(ll)%w > 0.0_dp) then
        max_en = max(max_en, PC_v_to_en(pc%particles(ll)%v, UC_elec_mass))
      end if
    end do
    max_en = max_en / UC_elec_volt
  end function get_max_energy

end module m_particles
