!> Module for ion movement as simple particles with only a drift velocity.
module m_ions
  use m_particle_core
  use m_particles
  use m_dielectric
  use m_lookup_table

  implicit none
  public

  !> parameter used for merging ion-particles
  real(dp) :: ion_min_merge_increase = 1.25_dp
  !> Secondary electron coefficient
  real(dp), parameter :: se_coefficient = 1.0e-1_dp
  !> Index of Ar^+ mobility
  integer, parameter :: i_td_mu = 1
  ! Table with electric field vs mobility
  type(LT_t), protected :: td_tbl_mu
contains

  subroutine init_pc_ion(cfg, pc_ions)
    !> initialize a PC instance for (Argon) ions.
    type(CFG_t), intent(inout)  :: cfg
    type(PC_t), intent(inout)   :: pc_ions

    integer   :: max_num_ions

    ! TODO Different settings for the ions?
    ! Particle model related parameters for ions
    ! call CFG_add_get(cfg, "ions%min_weight", ions_min_weight, &
    !      "Minimum weight for simulation particles")
    ! call CFG_add_get(cfg, "ions%max_weight", ions_max_weight, &
    !      "Maximum weight for simulation particles")

    call CFG_add(cfg, "ions%max_number", 15*1000*1000, &
         "Maximum number of particles")
    call CFG_add_get(cfg, "ions%min_merge_increase", ion_min_merge_increase, &
         "Minimum increase in particle count before merging")
    call CFG_get(cfg, "ions%max_number", max_num_ions)

    call load_ion_mobility_data(cfg)

    ! set ion-particles as a tracer and pass drift_velocity
    pc_ions%tracer_velocity => drift_velocity
    pc_ions%particle_mover => PC_tracer_advance_midpoint

    call pc_ions%initialize(UC_Ar_mass, max_num_ions)
    pc_ions%outside_check => outside_check

  end subroutine init_pc_ion

  subroutine particles_and_ions_to_density_and_events(tree, pc, pc_ions, init_cond)
    use m_cross_sec
    use m_random
    use m_particle_core
    use m_domain, only: outside_check
    type(af_t), intent(inout)        :: tree
    type(PC_t), intent(inout)        :: pc
    type(PC_t), intent(inout)        :: pc_ions
    logical, intent(in)              :: init_cond

    integer                     :: n, i, n_photons
    integer                     :: ix_surf, ix_cell

    real(dp), allocatable, save :: coords(:, :)
    real(dp), allocatable, save :: weights(:)
    integer, allocatable, save  :: id_guess(:)

    n = pc%n_events
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
          weights(i) = -pc%event_list(n)%part%w !TODO negative weights are not allowed for tracer particles!!!
          id_guess(i) = pc%event_list(n)%part%id

       else if (pc%event_list(n)%ctype == PC_particle_went_out .and. &
            pc%event_list(n)%cIx == inside_dielectric) then
          ! Now we map the particle to surface charge
          call particle_to_surface_charge(tree, pc%event_list(n)%part, &
               i_surf_elec)
       end if
    end do

    if (i > 0) then ! only for the events that created an ion
      call generate_ions_as_particles(pc_ions, i_pos_ion, coords(:, 1:i), &
          weights(1:i), i, interpolation_order_to_density, &
          id_guess(1:i))
    end if

    if (photoi_enabled) then
      call photoionization(tree, pc, coords, weights, n_photons)
      call generate_ions_as_particles(pc_ions, i_pos_ion, coords(:, 1:n_photons), &
          weights(1:n_photons), n_photons, interpolation_order_to_density)!id_guess(1:n_photons)
    end if

    if (photoe_enabled) then
      call photoemission(tree, pc)
    end if

    pc%n_events = 0

    !Do the events (i.e. secondary emission) for outside ions
    do n = 1, pc_ions%n_events
      if (pc_ions%event_list(n)%ctype == PC_particle_went_out .and. &
          pc_ions%event_list(n)%cIx == inside_dielectric) then
        call particle_to_surface_charge(tree, pc_ions%event_list(n)%part, &
              i_surf_pos_ion)
        call secondary_electron_emission(tree, diel, pc, pc_ions%event_list(n)%part)
      end if
    end do

    pc_ions%n_events = 0

    ! update densities (as electrons and ions might have been created or destroyed)
    call update_density(tree, pc_ions, i_pos_ion)
    call update_density(tree, pc, i_electron)

    if (GL_use_dielectric) then
       ! Map densities inside the first dielectric layers to surface charge.
       ! Particles can still move away from the dielectric.
       call dielectric_inside_layer_to_surface(tree, diel, i_electron, &
            i_surf_elec_close, 1.0_dp, clear_cc=.true., clear_surf=.true.)
       call dielectric_inside_layer_to_surface(tree, diel, i_pos_ion, &
            i_surf_pos_ion_close, 1.0_dp, clear_cc=.true., clear_surf=.true.)
       ! Sum densities together
       call dielectric_set_weighted_sum(diel, i_surf_sum_dens, &
            [i_surf_elec, i_surf_elec_close, i_surf_pos_ion,i_surf_pos_ion_close], &
            [-1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp])
    end if
  end subroutine particles_and_ions_to_density_and_events

  subroutine update_density(tree, pc, iv)
    !> Update the ion density (tree-variable)
    type(af_t), intent(inout) ::  tree
    type(pc_t), intent(inout) ::  pc
    integer, intent(in)       :: iv
    integer                   :: n_part, n
    real(dp), allocatable, save :: coords(:, :)
    real(dp), allocatable, save :: weights(:)
    integer, allocatable, save  :: id_guess(:)

    n_part = pc%get_num_sim_part()
    n = nint(n_part * array_incr_fac)
    if (.not. allocated(weights)) then
      allocate(coords(NDIM, n))
      allocate(weights(n))
      allocate(id_guess(n))
    else if (size(weights) < n_part) then
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

    call af_tree_clear_cc(tree, iv)
    call af_particles_to_grid(tree, iv, coords(:, 1:n_part), &
        weights(1:n_part), n_part, interpolation_order_to_density, &
        id_guess(1:n_part))
  end subroutine update_density

  subroutine generate_ions_as_particles(pc_ions, i_pos_ion, coords, weights, n_new_ions, &
       order, id_guess)
    !> Generation of ions as tracer particles
    type(pc_t), intent(inout)        :: pc_ions
    integer, intent(in)              :: i_pos_ion            !< Variable to store density
    integer, intent(in)              :: n_new_ions           !< The number of particles
    real(dp), intent(in)             :: coords(NDIM, n_new_ions) !< The particle coordinates
    real(dp), intent(in)             :: weights(n_new_ions)    !< Weights for the particles
    integer, intent(in)              :: order                   !< Order of interpolation
    !> Guess for box id containing particle, set to 0 where no guess is available
    integer, intent(inout), optional :: id_guess(n_new_ions)

    integer :: i
    type(PC_part_t)         :: new_part
    do i = 1, n_new_ions !TODO do it with buffer!  or in parallel idk
      new_part%x(1:NDIM) = coords(1:NDIM, i)
      new_part%v = drift_velocity(new_part)
      new_part%w = weights(i)
      if (present(id_guess)) then
        new_part%id = id_guess(i)
      end if

      call pc_ions%add_part(new_part)
    end do
  end subroutine generate_ions_as_particles

  subroutine secondary_electron_emission(tree, diel, pc, ion)
    type(af_t), intent(inout)         :: tree
    type(dielectric_t), intent(inout) :: diel
    type(pc_t), intent(inout)         :: pc
    type(PC_part_t)                   :: ion, new_electron

    type(af_loc_t)  :: loc
    integer         :: ix_surf, nb, id

    ! Find location
    loc = af_get_loc(tree, ion%x(1:NDIM))
    id = loc%id
    ! Check if id is valid
    if (id == -1) error stop "Coordinate out of domain"
    ! Check for a surface (only from the inside)
    ix_surf = diel%box_id_in_to_surface_ix(id)
    if (ix_surf == -1) error stop "No surface found for secondary electron emission"
    ! Get the direction of the surface
    nb = diel%surfaces(ix_surf)%direction

    ! The new electron will be created at the cell-center of the closest
    ! gas-phase surface cell (i.e. one of the ghost cells)
    select case (nb)
#if NDIM == 2
    case (af_neighb_lowx)
      new_electron%x(1:NDIM) = af_r_cc(tree%boxes(id), loc%ix + [1, 0])
    case (af_neighb_highx)
      new_electron%x(1:NDIM) = af_r_cc(tree%boxes(id), loc%ix + [-1, 0])
    case (af_neighb_lowy)
      new_electron%x(1:NDIM) = af_r_cc(tree%boxes(id), loc%ix + [0, 1])
    case (af_neighb_highy)
      new_electron%x(1:NDIM) = af_r_cc(tree%boxes(id), loc%ix + [0, -1])
#elif NDIM == 3
    case default
    error stop
#endif
  end select
    new_electron%v(:) = 0
    new_electron%a(:) = pc%accel_function(new_electron)
    new_electron%w    = se_coefficient * ion%w

    call pc%add_part(new_electron) !TODO do it with buffer! or in parallel idk

  end subroutine secondary_electron_emission

  function drift_velocity(ion)
    !> Calculate the drift velocity for a particle treated as a tracer
    type(PC_part_t), intent(inout) :: ion
    real(dp)        :: drift_velocity(3)
    real(dp)        :: E(NDIM)
    real(dp)        :: mu
    logical         :: success
    ! Get mu from the lookuptable
    E  = af_interp1_fc(tree, ion%x(1:NDIM), ifc_E, success, id_guess=ion%id)
    if (.not. success) error stop "Can not find E for the calculation of ion drift velocity."
    mu = LT_get_col(td_tbl_mu, i_td_mu, norm2(E))

    drift_velocity = 0
    drift_velocity(1:NDIM) = mu * E

    ! ! > test debug
    ! if (norm2(drift_velocity) >= 1000) &
    ! write(*, "(A20,E12.4)") "the ion velocity is", norm2(drift_velocity
    ! ! < end test debug_flag
  end function drift_velocity

  subroutine load_ion_mobility_data(cfg)
    ! Initialize the ion mobility data and write to a lookup table from a datafile
    ! Assumptions on data file: x-axis in V and y-axis in m^2/(V s).
    use m_transport_data
    use m_config
    use m_gas
    type(CFG_t), intent(inout) :: cfg

    character(len=GL_slen)     :: td_file_ions = "../../input/transport_data_ar0.txt"
    integer                    :: table_size_ions           = 500
    real(dp)                   :: max_electric_fld = 3.5e7_dp
    real(dp), allocatable      :: x_data(:), y_data(:)
    character(len=GL_slen)     :: data_name

    call CFG_add_get(cfg, "ions%transport_data_file", td_file_ions, &
         "Input file with ion transport data")
    call CFG_add_get(cfg, "ions%lookup_table_size", table_size_ions, &
         "The ion transport data table size")
    call CFG_add_get(cfg, "ions%lookup_table_max_efield", max_electric_fld, &
         "The maximum reduced electric field for the ion mobility")

    ! Create a lookup table for the model coefficients
    td_tbl_mu = LT_create(0.0_dp, max_electric_fld, table_size_ions, 1)

    ! Fill table with data
    data_name = "efield[V/m]_vs_ion_mu[m2/Vs]"
    call CFG_add_get(cfg, "td_mu_name", data_name, &
         "The name of the ion mobility coefficient.")
    call TD_get_from_file(td_file_ions, GL_gas_name, &
         trim(data_name), x_data, y_data)

    ! Convert mobility in 1 atm to the low pressure mobility with propotional relation
    call LT_set_col(td_tbl_mu, i_td_mu, x_data, y_data * (1.0_dp/GL_gas_pressure))
  end subroutine load_ion_mobility_data

end module m_ions
