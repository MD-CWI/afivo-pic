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
    integer                        :: particle_bytes
    real(dp)                       :: temperature, max_ev
    character(len=200)             :: cs_file
    character(len=20), allocatable :: gas_names(:)
    real(dp), allocatable          :: gas_fracs(:)
    type(CS_t), allocatable        :: cross_secs(:)
    integer                        :: ii
    type(PC_part_t)                :: dummy_part

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

    ! Initialize gas and electric field module
    call GAS_initialize(gas_names, gas_fracs, GL_gas_pressure, temperature)

    do nn = 1, n_gas_comp
       call CS_add_from_file(trim(cs_file), &
            trim(gas_names(nn)), gas_fracs(nn) * &
            GAS_number_dens, max_ev, cross_secs)
    end do

    call CS_write_summary(cross_secs, &
         trim(GL_output_dir) // "/" // trim(GL_simulation_name) // "_cs_summary.txt")

    call CS_create_ledger(cross_secs, &
         trim(GL_output_dir) // "/" // trim(GL_simulation_name) // "_cs_ledger.txt")

    call CFG_get(cfg, "particle%lkptbl_size", tbl_size)

    if (GL_cylindrical) then
       pc%particle_mover => PC_verlet_cyl_advance
       pc%after_mover => PC_verlet_cyl_correct_accel
    else
       pc%particle_mover => PC_verlet_advance
       pc%after_mover => PC_verlet_correct_accel
    end if

    ! How many bytes are required per particle (accounting for overhead, for
    ! example when they are mapped to a grid)
    particle_bytes = (storage_size(dummy_part) + 2 * storage_size(0)) / 8
    max_num_part = nint(GL_memory_particles_GB * 2.0_dp**30 / particle_bytes)
    write(*, "(A,I12)")   " Max number of particles: ", max_num_part

    call pc%initialize(UC_elec_mass, max_num_part)
    call pc%use_cross_secs(max_ev, tbl_size, cross_secs)

    pc%accel_function => get_accel
    pc%outside_check => outside_check

    where (pc%colls(:)%type == CS_ionize_t .or. &
         pc%colls(:)%type == CS_attach_t)
         pc%coll_is_event(:) = .true.
    end where

    if (GL_track_CAS) then
      do ii = 1, num_cIx_to_track
        pc%coll_is_event(GL_cIx_to_track(ii)) = .true.
      end do
    end if

  end subroutine init_particle

  !> Set particle tags to cell they are in. Assumes that the 'id' of each
  !> particle is correctly set, and that there are no 'dead' particles
  subroutine set_particle_tags(tree, pc, n_velocity_bins)
    type(af_t), intent(in)    :: tree
    type(PC_t), intent(inout) :: pc
    integer, intent(out)      :: n_velocity_bins
    integer                   :: n, nc, ix(NDIM), mapping(NDIM)
    integer                   :: max_ix_tag, v_to_int
    real(dp)                  :: inv_max_vel

    nc = tree%n_cell
    inv_max_vel = 1 / pc%max_vel

    do n = 1, NDIM
       mapping(n) = nc**(n-1)
    end do

    ! Maximal tag due to cell index
    max_ix_tag = sum((nc-1) * mapping + 1)

    ! Determine number of bins for velocity
    n_velocity_bins = huge(pc%particles(1)%tag)/(max_ix_tag + 1)

    !$omp parallel do private(ix, v_to_int)
    do n = 1, pc%n_part
       ix = af_cc_ix(tree%boxes(pc%particles(n)%id), &
            get_coordinates(pc%particles(n)))
       ! Convert cell index to number
       pc%particles(n)%tag = sum((ix-1) * mapping) + 1

       ! Then the velocity is mapped to an integer
       v_to_int = floor((n_velocity_bins - 1) * &
            min(1.0_dp, norm2(pc%particles(n)%v) * inv_max_vel))

       ! Combine the two pieces into one number
       pc%particles(n)%tag = pc%particles(n)%tag * n_velocity_bins + &
            v_to_int
    end do
    !$omp end parallel do
  end subroutine set_particle_tags

  !> Adjust the weights of the particles
  subroutine adapt_weights(tree, pc, t_sort, t_rest)
    use omp_lib
    type(af_t), intent(in)    :: tree
    type(PC_t), intent(inout) :: pc
    real(dp), intent(out)     :: t_sort, t_rest
    real(dp)                  :: t0, t1, t2
    integer                   :: n, n_threads, thread_id
    integer, allocatable      :: ix_thread(:)
    integer                   :: cell_tag, cell_i0, cell_i1
    integer                   :: i, j, i_buffer, n_part_prev, N_vb
    real(dp)                  :: v_new, w_min, w_max, desired_weight
    real(dp)                  :: w_new, remainder
    integer, parameter        :: buffer_size = 1024
    type(PC_part_t)           :: pbuffer(buffer_size)
    type(prng_t)              :: prng

    ! print *, "before: ", pc%get_num_sim_part(), pc%get_num_real_part(), &
         ! pc%get_mean_energy() / UC_elec_volt

    ! Set tags that contain information about the cell index and velocity
    call set_particle_tags(tree, pc, N_vb)
    t0 = omp_get_wtime()
    call pc%sort_in_place_by_id_tag()
    t1 = omp_get_wtime()

    n_threads = af_get_max_threads()
    allocate(ix_thread(0:n_threads))

    ! Determine which threads work on which particles
    ix_thread(0) = 1
    ix_thread(n_threads) = pc%n_part + 1
    do n = 1, n_threads-1
       ix_thread(n) = nint(pc%n_part * real(n, dp)/n_threads)
    end do

    ! Correct so that the boundaries between threads occur at tag boundaries
    do n = 1, n_threads-1
       do while (pc%particles(ix_thread(n))%tag/N_vb == &
            pc%particles(ix_thread(n)-1)%tag/N_vb .and. &
            ix_thread(n) > ix_thread(n-1))
          ix_thread(n) = ix_thread(n) - 1
       end do
    end do

    call prng%init_parallel(omp_get_max_threads(), GL_rng)

    !$omp parallel private(thread_id, cell_tag, cell_i0, cell_i1, &
    !$omp i, j, i_buffer, v_new, w_min, w_max, desired_weight, &
    !$omp pbuffer, w_new, remainder)
    i_buffer = 0
    thread_id = omp_get_thread_num()

    cell_i0 = ix_thread(thread_id)
    do while (cell_i0 < ix_thread(thread_id+1))
       ! Find indices of particles with current tag
       cell_tag = pc%particles(cell_i0)%tag/N_vb
       do cell_i1 = cell_i0+1, ix_thread(thread_id+1) - 1
          if (pc%particles(cell_i1)%tag/N_vb /= cell_tag) exit
       end do
       ! Went one index too far, so subtract one
       cell_i1 = cell_i1 - 1

       desired_weight = get_desired_weight(pc%particles(cell_i0))
       w_min = desired_weight / 1.5_dp
       w_max = desired_weight * 1.5_dp

       ! Merge particles
       i = cell_i0
       do while (i < cell_i1)
          if (pc%particles(i)%w > 0 .and. pc%particles(i)%w < w_min) then
             ! Look for candidate to merge with
             do j = i + 1, cell_i1
                if (pc%particles(j)%w > 0 .and. pc%particles(j)%w < w_min) exit
             end do

             if (j == cell_i1 + 1) exit ! No candidate found

             ! Merge particles
             associate (pa => pc%particles(i), pb => pc%particles(j))
               if (prng%rngs(thread_id+1)%unif_01() > pa%w/(pa%w + pb%w)) then
                  ! Keep particle b
                  pb%w = pa%w + pb%w
                  pa%w = PC_dead_weight
               else
                  ! Keep particle a
                  pa%w = pa%w + pb%w
                  pb%w = PC_dead_weight
               end if
             end associate

             ! Jump to next possible particle that can be merged
             i = j + 1
          else
             ! Proceed with next particle
             i = i + 1
          end if
       end do

       ! Split particles
       do i = cell_i0, cell_i1
          if (pc%particles(i)%w > 0 .and. pc%particles(i)%w > w_max) then
             ! Determine new weights that are multiples of particle_min_weight,
             ! there is a correction below for odd multiples
             w_new = floor(0.5_dp * pc%particles(i)%w/particle_min_weight) * &
                  particle_min_weight
             remainder = pc%particles(i)%w - 2 * w_new
             pc%particles(i)%w = w_new

             ! Empty buffer if necessary
             if (i_buffer == buffer_size) then
                !$omp critical
                j = pc%n_part
                pc%n_part = pc%n_part + i_buffer
                !$omp end critical
                call pc%check_space(j+i_buffer)
                pc%particles(j+1:j+i_buffer) = pbuffer(1:i_buffer)
                i_buffer = 0
             end if

             i_buffer = i_buffer + 1
             pbuffer(i_buffer) = pc%particles(i)

             ! Correct for original weight that could not evenly be divided in two
             if (remainder > 0) then
                if (prng%rngs(thread_id+1)%unif_01() > 0.5_dp) then
                   pc%particles(i)%w = pc%particles(i)%w + remainder
                else
                   pbuffer(i_buffer)%w = pbuffer(i_buffer)%w + remainder
                end if
             end if
          end if
       end do

       ! Go to next cell
       cell_i0 = cell_i1 + 1
    end do

    ! Empty buffer at end
    if (i_buffer > 0) then
       !$omp critical
       j = pc%n_part
       pc%n_part = pc%n_part + i_buffer
       !$omp end critical
       call pc%check_space(j+i_buffer)
       pc%particles(j+1:j+i_buffer) = pbuffer(1:i_buffer)
       i_buffer = 0
    end if
    !$omp end parallel

    ! Clean up dead particles at the end
    i = 1
    do while (i <= pc%n_part)
       if (pc%particles(i)%w <= PC_dead_weight) then
          n_part_prev = pc%n_part
          ! This is overridden if a replacement is found
          pc%n_part = min(pc%n_part, i-1)

          ! Find the last "alive" particle in the list
          do j = n_part_prev, i+1, -1
             if (pc%particles(j)%w > PC_dead_weight) then
                pc%particles(i) = pc%particles(j)
                pc%n_part       = j-1
                exit
             end if
          end do
       end if
       i = i + 1
    end do

    call prng%update_seed(GL_rng)

    t2 = omp_get_wtime()
    t_sort = t1 - t0
    t_rest = t2 - t1

    ! print *, "after: ", pc%get_num_sim_part(), pc%get_num_real_part(), &
         ! pc%get_mean_energy() / UC_elec_volt
  end subroutine adapt_weights

  subroutine particles_to_density_and_events(tree, pc, init_cond)
    use m_cross_sec
    use m_particle_core
    use m_domain, only: outside_check
    type(af_t), intent(inout)        :: tree
    type(PC_t), intent(inout)        :: pc
    logical, intent(in)              :: init_cond

    integer                     :: n, i, j, ii, n_part, n_photons, n_part_before

    real(dp), allocatable, save :: coords(:, :), coords_CAS(:, :)
    real(dp), allocatable, save :: weights(:), weights_CAS(:), mask_var(:)
    integer, allocatable, save  :: id_guess(:), id_guess_CAS(:), mask(:)
    real(dp), allocatable, save :: energy_lost(:)

    n_part = pc%get_num_sim_part()
    n = pc%n_events

    if (.not. allocated(weights)) then
       n = nint(n * array_incr_fac)
       allocate(coords(NDIM, n), weights(n), id_guess(n), energy_lost(n))
       if (GL_track_CAS) then
         allocate(coords_CAS(NDIM, n), weights_CAS(n), id_guess_CAS(n), mask(n), mask_var(n))
       end if
    else if (size(weights) < n) then
       n = nint(n * array_incr_fac)
       deallocate(coords, weights, id_guess, energy_lost)
       allocate(coords(NDIM, n), weights(n), id_guess(n), energy_lost(n))
       if (GL_track_CAS) then
         deallocate(coords_CAS, mask_var, mask, id_guess_CAS, weights_CAS)
         allocate(coords_CAS(NDIM, n), weights_CAS(n), id_guess_CAS(n), mask(n), mask_var(n))
       end if
    end if

    if (init_cond) then
       call af_tree_clear_cc(tree, i_electron)
       call af_particles_to_grid(tree, i_electron, n_part, &
            get_id, get_rw, interpolation_order_to_density, &
            iv_tmp=i_tmp_dens)
       call af_particles_to_grid(tree, i_pos_ion, n_part, &
            get_id, get_rw, interpolation_order_to_density, &
            iv_tmp=i_tmp_dens)
    else
      call af_tree_clear_cc(tree, i_electron)
      call af_particles_to_grid(tree, i_electron, n_part, &
           get_id, get_rw, interpolation_order_to_density, &
           iv_tmp=i_tmp_dens)
      ! FLAG This routine needs to be rewritten in new syntax
      ! call af_particles_to_grid(tree, i_energy_dep, coords(:, 1:n_part), &
      !     weights(1:n_part)*energy_lost(1:n_part), n_part, 1, &
      !     id_guess(1:n_part))
      call af_particles_to_grid(tree, i_energy_dep, n_part, &
          get_id, get_r_energydeposit, interpolation_order_to_density, &
          iv_tmp=i_tmp_dens)
      ! END FLAG
    end if

    !$omp parallel do
    do n = 1, n_part
       ! coords(:, n)   = pc%particles(n)%x(1:NDIM)
       ! weights(n)     = pc%particles(n)%w
       ! id_guess(n)    = pc%particles(n)%id
       ! energy_lost(n) = pc%particles(n)%en_loss
       pc%particles(n)%en_loss = 0.0_dp ! Reset the energy_loss of all particles
    end do
    !$omp end parallel do

    i = 0 ! Counter for ionizing events
    j = 0 ! Counter for CAS tracker

    do n = 1, pc%n_events
       if (pc%event_list(n)%ctype == CS_ionize_t) then
          i = i + 1
          coords(:, i) = get_coordinates(pc%event_list(n)%part)
          weights(i) = pc%event_list(n)%part%w
          id_guess(i) = pc%event_list(n)%part%id
          energy_lost(i) = 0.0_dp ! Energy loss for this collision is already stored in the particles
       else if (pc%event_list(n)%ctype == CS_attach_t) then
          i = i + 1
          coords(:, i) = get_coordinates(pc%event_list(n)%part)
          weights(i) = -pc%event_list(n)%part%w
          id_guess(i) = pc%event_list(n)%part%id
          energy_lost(i) = PC_v_to_en(pc%event_list(n)%part%v, UC_elec_mass) ! All of the particles energy is deposited in the gas

       else if (pc%event_list(n)%ctype == PC_particle_went_out .and. &
            pc%event_list(n)%cIx == inside_dielectric) then
          ! Now we map the particle to surface charge
          call particle_to_surface_charge(tree, pc%event_list(n)%part, &
               i_surf_elec)
       end if

       if (GL_track_CAS) then
         if(any(pc%event_list(n)%cIx == GL_cIx_to_track)) then
           j = j + 1
           coords_CAS(:, j) = pc%event_list(n)%part%x(1:NDIM)
           weights_CAS(j)   = pc%event_list(n)%part%w
           id_guess_CAS(j)  = pc%event_list(n)%part%id
           mask(j)          = pc%event_list(n)%cIx
         end if
       end if
    end do

    if (i > 0) then ! only for the events that created an ion
       call af_particles_to_grid(tree, i_pos_ion, i, &
           get_event_id, get_event_rw, interpolation_order_to_density, &
           iv_tmp=i_tmp_dens)
      ! FLAG rewrite in terms of new syntax
       ! call af_particles_to_grid(tree, i_energy_dep, coords(:, 1:i), &
       !     abs(weights(1:i)*energy_lost(1:i)), i, 1, &
       !     id_guess(1:i))
       call af_particles_to_grid(tree, i_energy_dep, i, &
           get_event_id, get_event_r_energydeposit, interpolation_order_to_density, &
           iv_tmp=i_tmp_dens)
    end if

    if (j > 0) then
      do ii = 1, num_cIx_to_track
        mask_var = 0.0_dp
        where (mask == GL_cIx_to_track(ii)) mask_var = 1.0_dp
        call af_particles_to_grid(tree, i_tracked_cIx(ii), j, &
              get_event_id_CAS, get_event_rw_CAS, interpolation_order_to_density, &
              iv_tmp=i_tmp_dens)
      end do
      ! END FLAG
    end if

    if (photoi_enabled) then
       ! Photo-electrons are added to the end of the particle list
       n_part_before = pc%n_part
       call photoionization(tree, pc, n_photons)
       call af_particles_to_grid(tree, i_pos_ion, n_photons, &
            get_id, get_rw, interpolation_order_to_density, &
            iv_tmp=i_tmp_dens, offset_particles=n_part_before)
    end if

    if (photoe_enabled) then
       call photoemission(tree, pc)
    end if

    pc%n_events = 0

    if (GL_use_dielectric) then
       ! Map densities inside the first dielectric layers to surface charge.
       ! Electrons can still move away from the dielectric.
       call surface_inside_layer_to_surface(tree, diel, i_electron, &
            i_surf_elec_close, 1.0_dp, clear_cc=.true., clear_surf=.true.)
       call surface_inside_layer_to_surface(tree, diel, i_pos_ion, &
            i_surf_pos_ion, 1.0_dp, clear_cc=.true., clear_surf=.false.)
       ! Sum densities together
       call surface_set_weighted_sum(diel, i_surf_sum_dens, &
            [i_surf_elec, i_surf_elec_close, i_surf_pos_ion], &
            [-1.0_dp, -1.0_dp, 1.0_dp])
    end if

  contains

    subroutine get_event_id(n, id)
      integer, intent(in)  :: n
      integer, intent(out) :: id

      id = af_get_id_at(tree, coords(:, n), guess=id_guess(n))
    end subroutine get_event_id

    subroutine get_event_id_CAS(n, id)
      integer, intent(in)  :: n
      integer, intent(out) :: id

      id = af_get_id_at(tree, coords_CAS(:, n), guess=id_guess_CAS(n))
    end subroutine get_event_id_CAS

    !> Get particle position and weight
    subroutine get_event_rw(n, r, w)
      integer, intent(in)   :: n
      real(dp), intent(out) :: r(NDIM)
      real(dp), intent(out) :: w

      r = coords(:, n)
      w = weights(n)
    end subroutine get_event_rw

    !> Get energy deposited to the gas
    subroutine get_event_r_energydeposit(n, r, w)
      integer, intent(in)   :: n
      real(dp), intent(out) :: r(NDIM)
      real(dp), intent(out) :: w

      r = coords(:, n)
      w = abs(weights(n) * energy_lost(n))
    end subroutine get_event_r_energydeposit

    !> Get position and reactive species produced
    subroutine get_event_rw_CAS(n, r, w)
      integer, intent(in)   :: n
      real(dp), intent(out) :: r(NDIM)
      real(dp), intent(out) :: w

      r = coords_CAS(:, n)
      w = mask_var(n) * weights_CAS(n)
    end subroutine get_event_rw_CAS

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

    call surface_get_surface_cell(tree, diel, get_coordinates(my_part), &
         ix_surf, ix_cell)
    ! Update the charge in the surface cell
#if NDIM == 2
    diel%surfaces(ix_surf)%sd(ix_cell(1), i_surf) = &
         diel%surfaces(ix_surf)%sd(ix_cell(1), i_surf) + &
         my_part%w / diel%surfaces(ix_surf)%dr(1)
#elif NDIM == 3
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
    curve_dat(1, 2, :) = bin_values/(n_part * en_step)
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

  !> Get particle id
  subroutine get_id(n, id)
    integer, intent(in)  :: n
    integer, intent(out) :: id

    id = af_get_id_at(tree, get_coordinates(pc%particles(n)), &
         guess=pc%particles(n)%id)
  end subroutine get_id

  !> Get particle position and weight
  subroutine get_rw(n, r, w)
    integer, intent(in)   :: n
    real(dp), intent(out) :: r(NDIM)
    real(dp), intent(out) :: w

    r = get_coordinates(pc%particles(n))
    w = pc%particles(n)%w
  end subroutine get_rw

  !> Get energy deposited to the gas
  subroutine get_r_energydeposit(n, r, w)
    integer, intent(in)   :: n
    real(dp), intent(out) :: r(NDIM)
    real(dp), intent(out) :: w

    r = get_coordinates(pc%particles(n))
    w = pc%particles(n)%en_loss * pc%particles(n)%w
  end subroutine get_r_energydeposit

  !> Get particle position and energy
  subroutine get_r_energy(n, r, w)
    integer, intent(in)   :: n
    real(dp), intent(out) :: r(NDIM)
    real(dp), intent(out) :: w

    r = get_coordinates(pc%particles(n))
    w = pc%particles(n)%w * PC_v_to_en(pc%particles(n)%v, UC_elec_mass) / &
            UC_elec_volt
  end subroutine get_r_energy

  !> Get particle position and unit density (for particle-per-cell computation)
  subroutine get_r_unit_dens(n, r, w)
    integer, intent(in)   :: n
    real(dp), intent(out) :: r(NDIM)
    real(dp), intent(out) :: w

    r = get_coordinates(pc%particles(n))
    w = 1.0_dp
  end subroutine get_r_unit_dens

end module m_particles
