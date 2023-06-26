#include "../afivo/src/cpp_macros.h"
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
    type(PC_part_t)                :: dummy_part
    real(dp)                       :: boris_dt_factor = 0.1_dp

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

    call CFG_add_get(cfg, "particle%boris_dt_factor", boris_dt_factor, &
         "With B-field: limit time step to boris_dt_factor / cyclotron freq.")

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

    call CFG_get(cfg, "particle%lkptbl_size", tbl_size)

    if (magnetic_field_used) then
       if (GL_cylindrical) error stop "Not implemented"

       pc%B_vec = magnetic_field
       pc%particle_mover => PC_boris_advance

       pc%dt_max = boris_dt_factor * 2 * UC_pi / &
            (norm2(magnetic_field) * abs(UC_elec_q_over_m))
       write(*, "(A,E12.2)") " Boris max. time step:    ", pc%dt_max
    else
       if (GL_cylindrical) then
          pc%particle_mover => PC_verlet_cyl_advance
          pc%after_mover => PC_verlet_cyl_correct_accel
       else
          pc%particle_mover => PC_verlet_advance
          pc%after_mover => PC_verlet_correct_accel
       end if
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
    real(dp)                  :: inv_particle_min_weight, w_relative
    integer, parameter        :: buffer_size = 1024
    type(PC_part_t)           :: pbuffer(buffer_size)
    type(prng_t)              :: prng

    t_sort = 0
    t_rest = 0

    ! Don't attempt to adapt weights without particles
    if (pc%n_part == 0) return

    ! print *, "before: ", pc%get_num_sim_part(), pc%get_num_real_part(), &
         ! pc%get_mean_energy() / UC_elec_volt
    inv_particle_min_weight = 1/particle_min_weight

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
       ix_thread(n) = max(ix_thread(n-1), nint(pc%n_part * real(n, dp)/n_threads))
    end do

    ! Correct so that the boundaries between threads occur at tag boundaries
    do n = 1, n_threads-1
       ! Check if the thread has particles assigned to it
       do while (ix_thread(n) > ix_thread(n-1))
          if (pc%particles(ix_thread(n))%tag/N_vb == &
               pc%particles(ix_thread(n)-1)%tag/N_vb) then
             ! Adjust boundary
             ix_thread(n) = ix_thread(n) - 1
          else
             exit
          end if
       end do
    end do

    call prng%init_parallel(omp_get_max_threads(), GL_rng)

    !$omp parallel private(thread_id, cell_tag, cell_i0, cell_i1, &
    !$omp i, j, i_buffer, v_new, w_min, w_max, desired_weight, &
    !$omp pbuffer, w_new, w_relative, remainder)
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

       ! Now pc%particles(cell_i0:cell_i1) are in the same grid cell

       desired_weight = get_desired_weight(pc%particles(cell_i0))
       w_min = desired_weight / 1.5_dp
       w_max = desired_weight * 1.5_dp

       ! Merge particles
       i = cell_i0
       do while (i < cell_i1)
          if (pc%particles(i)%w >= 0 .and. pc%particles(i)%w < w_min) then
             ! Look for candidate to merge with
             do j = i + 1, cell_i1
                if (pc%particles(j)%w >= 0 .and. pc%particles(j)%w < w_min) exit
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
          w_relative = pc%particles(i)%w * inv_particle_min_weight
          if (w_relative >= 2 .and. pc%particles(i)%w > w_max) then
             ! Determine new weights that are multiples of particle_min_weight,
             ! there is a correction below for odd multiples
             w_new = floor(0.5_dp * w_relative) * particle_min_weight
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
    use m_domain
    type(af_t), intent(inout)        :: tree
    type(PC_t), intent(inout)        :: pc
    logical, intent(in)              :: init_cond

    integer                     :: n, i, n_part, n_photons
    integer                     :: n_part_before
    real(dp), allocatable, save :: coords(:, :)
    real(dp), allocatable, save :: weights(:)
    integer, allocatable, save  :: id_guess(:)

    n_part = pc%get_num_sim_part()
    n = pc%n_events

    if (.not. allocated(weights)) then
       n = nint(n * array_incr_fac)
       allocate(coords(NDIM, n), weights(n), id_guess(n))
    else if (size(weights) < n) then
       n = nint(n * array_incr_fac)
       deallocate(coords, weights, id_guess)
       allocate(coords(NDIM, n), weights(n), id_guess(n))
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
    end if

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
       call af_particles_to_grid(tree, i_pos_ion, i, &
            get_event_id, get_event_rw, interpolation_order_to_density, &
            iv_tmp=i_tmp_dens)
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

    !> Get particle position and weight
    subroutine get_event_rw(n, r, w)
      integer, intent(in)   :: n
      real(dp), intent(out) :: r(NDIM)
      real(dp), intent(out) :: w

      r = coords(:, n)
      w = weights(n)
    end subroutine get_event_rw

  end subroutine particles_to_density_and_events

  !> Add background ionization due to a continues volumetric source
  subroutine particles_background_ionization(tree, pc, time_elapsed)
    use m_particle_core
    use m_domain
    type(af_t), intent(inout) :: tree
    type(PC_t), intent(inout) :: pc
    real(dp), intent(in)      :: time_elapsed
    integer                   :: n, n_new, n_part_before, n_ionizations
    real(dp)                  :: num_electrons, x(NDIM)
    type(PC_part_t)           :: new_part

    if (GL_background_ionization_rate > 0) then
       num_electrons = GL_background_ionization_rate * domain_volume * &
            time_elapsed / particle_min_weight

       ! Draw Poisson random number
       n_ionizations = GL_rng%poisson(num_electrons)

       new_part%a = 0
       new_part%v = 0
       new_part%x = 0
       new_part%w = particle_min_weight
       n_part_before = pc%n_part

       do n = 1, n_ionizations
          ! Get location (TODO: cylindrical sampling)

          if (NDIM == 2 .and. GL_cylindrical) then
             ! r, z
             x(1:2) = [sqrt(GL_rng%unif_01()), GL_rng%unif_01()]
          else
             x = [DTIMES(GL_rng%unif_01())]
          end if

          new_part%x(1:NDIM) = x * domain_len
          new_part%id = af_get_id_at(tree, new_part%x(1:NDIM))

          if (outside_check(new_part) == 0) then
             call pc%add_part(new_part)
          end if
       end do

       ! Add positive ions in the end
       n_new = pc%n_part - n_part_before
       if (n_new > 0) then
          call af_particles_to_grid(tree, i_pos_ion, n_new, &
               get_id, get_rw, interpolation_order_to_density, &
               iv_tmp=i_tmp_dens, offset_particles=n_part_before)
       end if
    end if
  end subroutine particles_background_ionization

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

    if (n_part > 0) then
       call pc%histogram(calc_elec_energy, is_alive, [0.0_dp], bins, bin_values)
       ! Convert histogram to density and save as curve-object
       curve_dat(1, 1, :) = bins
       curve_dat(1, 2, :) = bin_values/(n_part * en_step) + 1e-9 ! Add regularization parameter to prevent errors when converting to semilogy plots
    else
       curve_dat(1, 1, :) = bins
       curve_dat(1, 2, :) = 0.0_dp
    end if
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
