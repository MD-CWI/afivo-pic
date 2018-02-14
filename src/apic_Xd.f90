#include "afivo/src/cpp_macros_$Dd.h"
!> Program to perform $Dd discharge simulations in Cartesian and cylindrical coordinates
program apic_$Dd

  use m_a$D_all
  use m_globals
  use m_field_$Dd
  use m_init_cond_$Dd
  use m_particle_core
  use m_photoi_$Dd
  use m_domain_$Dd
  use m_refine_$Dd
  use m_time_step_$Dd

  implicit none

  integer, parameter     :: int8 = selected_int_kind(18)
  integer, parameter     :: ndim = $D
  integer(int8)          :: t_start, t_current, count_rate
  real(dp)               :: dt, fld_err
  real(dp)               :: wc_time, inv_count_rate, time_last_print
  integer                :: it, id, n_part_id
  character(len=ST_slen) :: fname
  logical                :: write_out
  type(CFG_t)            :: cfg  ! The configuration for the simulation
  type(a$D_t)            :: tree ! This contains the full grid information
  type(mg$D_t)           :: mg   ! Multigrid option struct
  type(PC_t)             :: pc
  type(PC_events_t)      :: events
  type(ref_info_t)       :: ref_info

  integer :: output_cnt = 0 ! Number of output files written

  call CFG_update_from_arguments(cfg)
  call domain_init(cfg)
  call refine_init(cfg, ndim)
  call time_step_init(cfg)
  call ST_initialize(cfg, ndim)

  call field_initialize(cfg, mg)
  call init_cond_initialize(cfg, ndim)

  call init_particle(cfg, pc)

  call pi_initialize(cfg)

  fname = trim(ST_output_dir) // "/" // trim(ST_simulation_name) // "_out.cfg"
  call CFG_write(cfg, trim(fname))

  ! Initialize the tree (which contains all the mesh information)
  call init_tree(tree)

  ! Set the multigrid options. First define the variables to use
  mg%i_phi = i_phi
  mg%i_tmp = i_Ex
  mg%i_rhs = i_rhs

  ! This automatically handles cylindrical symmetry
  mg%box_op => mg$D_auto_op
  mg%box_gsrb => mg$D_auto_gsrb
  mg%box_corr => mg$D_auto_corr

  ! This routine always needs to be called when using multigrid
  call mg$D_init_mg(mg)

  call a$D_set_cc_methods(tree, i_electron, a$D_bc_neumann_zero)
  call a$D_set_cc_methods(tree, i_pos_ion, a$D_bc_neumann_zero)
  call a$D_set_cc_methods(tree, i_phi, mg%sides_bc, mg%sides_rb)
  call a$D_set_cc_methods(tree, i_ppc, a$D_bc_neumann_zero)

  output_cnt      = 0         ! Number of output files written
  ST_time         = 0         ! Simulation time (all times are in s)

  ! Set up the initial conditions
  call init_cond_particles(tree, pc)

  do
     call a$D_tree_clear_cc(tree, i_pos_ion)
     call a$D_loop_box(tree, init_cond_set_box)
     call particles_to_density(tree, pc, events, .true.)
     call field_compute(tree, mg, .false.)
     call a$D_adjust_refinement(tree, refine_routine, ref_info, &
          refine_buffer_width, .true.)
     if (ref_info%n_add == 0) exit
  end do

  call pc%set_accel()

  print *, "Number of threads", af_get_max_threads()
  call a$D_print_info(tree)

  ! Start from small time step
  ST_dt   = ST_dt_min

  ! Initial wall clock time
  call system_clock(t_start, count_rate)
  inv_count_rate = 1.0_dp / count_rate
  time_last_print = -1e10_dp

  do it = 1, huge(1)-1
     if (ST_time >= ST_end_time) exit

     call system_clock(t_current)
     wc_time = (t_current - t_start) * inv_count_rate

     ! Every ST_print_status_interval, print some info about progress
     if (wc_time - time_last_print > ST_print_status_sec) then
        call print_status()
        time_last_print = wc_time
     end if

     ! Every ST_dt_output, write output
     if (output_cnt * ST_dt_output <= ST_time + ST_dt) then
        write_out  = .true.
        dt         = output_cnt * ST_dt_output - ST_time
        output_cnt = output_cnt + 1
     else
        write_out = .false.
        dt        = ST_dt
     end if

     call PM_fld_error(tree, pc, ST_rng, 1000, fld_err, .true.)

     call pc%advance_openmp(dt, events)

     ST_time = ST_time + dt
     call particles_to_density(tree, pc, events, .false.)

     ! Compute field with new density
     call field_compute(tree, mg, .true.)

     call PM_fld_error(tree, pc, ST_rng, 1000, fld_err, .false.)
     print *, "fld err", fld_err, ST_dt
     ST_dt = get_new_dt(ST_dt, fld_err, 10.0e-2_dp)

     call PC_verlet_correct_accel(pc, dt)

     if (modulo(it, 10) == 0) then
        call adapt_weights(tree, pc)
     end if

     if (write_out) then
        call prepare_output_variables()

        write(fname, "(A,I6.6)") trim(ST_simulation_name) // "_", output_cnt
        call a$D_write_silo(tree, fname, output_cnt, ST_time, &
             vars_for_output, dir=ST_output_dir)
        call print_info()
     end if

     if (mod(it, refine_per_steps) == 0) then
        call a$D_adjust_refinement(tree, refine_routine, ref_info, &
             refine_buffer_width, .true.)

        if (ref_info%n_add + ref_info%n_rm > 0) then
           call adapt_weights(tree, pc)
           ! Compute the field on the new mesh
           call particles_to_density(tree, pc, events, .false.)
           call field_compute(tree, mg, .true.)
        end if
        ! print *, "after"
        ! call print_info()
     end if
  end do

contains

  subroutine adapt_weights(tree, pc)
    type(a$D_t), intent(in)    :: tree
    type(PC_t), intent(inout) :: pc
    integer, allocatable      :: id_ipart(:)

    call sort_by_id(tree, pc, id_ipart)
    print *, "before: ", pc%get_num_sim_part()
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
    print *, "after:  ", pc%get_num_sim_part()

    call pc%clean_up()
  end subroutine adapt_weights

  subroutine sort_by_id(tree, pc, id_ipart)
    type(a$D_t), intent(in)              :: tree
    type(PC_t), intent(inout)           :: pc
    integer, intent(inout), allocatable :: id_ipart(:)

    integer                      :: n, id, new_ix
    integer, allocatable         :: id_count(:)
    integer, allocatable         :: id_ix(:)
    type(PC_part_t), allocatable :: p_copy(:)

    if (allocated(id_ipart)) deallocate(id_ipart)
    allocate(id_ipart(tree%highest_id+1))
    allocate(id_ix(tree%highest_id+1))
    allocate(id_count(tree%highest_id))
    allocate(p_copy(pc%n_part))

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

  !> Initialize the AMR tree
  subroutine init_tree(tree)
    type(a$D_t), intent(inout) :: tree

    ! Variables used below to initialize tree
    real(dp)                  :: dr
    integer                   :: id
    integer                   :: ix_list($D, 1) ! Spatial indices of initial boxes

    dr = domain_len / box_size

    ! Initialize tree
    if (ST_cylindrical) then
       call a$D_init(tree, box_size, n_var_cell, n_var_face, dr, &
            coarsen_to=2, coord=af_cyl, cc_names=ST_cc_names)
    else
       call a$D_init(tree, box_size, n_var_cell, n_var_face, dr, &
            coarsen_to=2, cc_names=ST_cc_names)
    end if

    ! Set up geometry
    id             = 1          ! One box ...
    ix_list(:, id) = 1          ! With index 1,1 ...

    ! Create the base mesh
    call a$D_set_base(tree, 1, ix_list)

  end subroutine init_tree

  subroutine init_particle(cfg, pc)
    use m_units_constants
    use m_gas
    use m_cross_sec
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
    call CFG_add(cfg, "particle%max_number", 50*1000*1000, &
         "Maximum number of particles")
    call CFG_add(cfg, "particle%lkptbl_size", 1000, &
         "The size of the lookup table for the collision rates")
    call CFG_add(cfg, "particle%max_energy_ev", 500.0_dp, &
         "The maximum energy in eV for particles in the simulation")
    call CFG_add(cfg, "particle%mover", "verlet", &
         "Which particle mover to use. Options: analytic, verlet, boris")
    call CFG_add(cfg, "boris_dt_factor", 0.1_dp, &
         "The maximum time step in terms of the cylotron frequency")

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
    call GAS_initialize(gas_names, gas_fracs, ST_gas_pressure, temperature)

    do nn = 1, n_gas_comp
       call CS_add_from_file(trim(cs_file), &
            trim(gas_names(nn)), gas_fracs(nn) * &
            GAS_number_dens, max_ev, cross_secs)
    end do

    call CS_write_summary(cross_secs, &
         trim(ST_output_dir) // "/" // trim(ST_simulation_name) // "_cs_summary.txt")

    call CFG_get(cfg, "particle%lkptbl_size", tbl_size)

    call pc%initialize(UC_elec_mass, cross_secs, &
         tbl_size, max_ev, max_num_part, get_random_seed())

    pc%accel_function => get_accel
    pc%outside_check => outside_check

    where (pc%colls(:)%type == CS_ionize_t .or. &
         pc%colls(:)%type == CS_attach_t)
       pc%coll_is_event(:) = .true.
    end where

  end subroutine init_particle

  function get_accel_pos(x) result(accel)
    use m_units_constants
    real(dp), intent(in) :: x($D)
    real(dp)             :: accel(3)

#if $D == 2
    accel(1:$D) = a$D_interp1(tree, x, [i_Ex, i_Ey], $D)
    accel(1:$D) = accel(1:$D) * UC_elec_q_over_m
    accel(3) = 0.0_dp
#elif $D == 3
    accel(1:$D) = a$D_interp1(tree, x, [i_Ex, i_Ey, i_Ez], $D)
    accel(1:$D) = accel(1:$D) * UC_elec_q_over_m
#endif
  end function get_accel_pos

  function get_accel(my_part) result(accel)
    use m_units_constants
    type(PC_part_t), intent(in) :: my_part
    real(dp)                    :: accel(3)

#if $D == 2
    accel(1:$D) = a$D_interp1(tree, my_part%x(1:$D), [i_Ex, i_Ey], $D, my_part%id)
    accel(3) = 0.0_dp
#elif $D == 3
    accel(1:$D) = a$D_interp1(tree, my_part%x(1:$D), [i_Ex, i_Ey, i_Ez], &
         $D, my_part%id)
#endif

    accel(:) = accel(:) * UC_elec_q_over_m

  end function get_accel

  subroutine particles_to_density(tree, pc, events, init_cond)
    use m_cross_sec
    use m_photoi_$Dd
    type(a$D_t), intent(inout)        :: tree
    type(PC_t), intent(inout)        :: pc
    type(PC_events_t), intent(inout) :: events
    logical, intent(in)              :: init_cond

    integer               :: n, i, n_part, n_events, n_photons, id
    real(dp)              :: x(3), v(3), a(3)
    real(dp), allocatable :: coords(:, :)
    real(dp), allocatable :: weights(:)
    integer, allocatable  :: id_guess(:)

    n_part = pc%get_num_sim_part()
    n_events = events%n_stored
    n = max(n_part, n_events)

    allocate(coords($D, n))
    allocate(weights(n))
    allocate(id_guess(n))

    !$omp parallel do
    do n = 1, n_part
       coords(:, n) = pc%particles(n)%x(1:$D)
       weights(n) = pc%particles(n)%w
       id_guess(n) = pc%particles(n)%id
    end do
    !$omp end parallel do

    if (init_cond) then
       call a$D_tree_clear_cc(tree, i_electron)
       call a$D_particles_to_grid(tree, i_electron, coords(:, 1:n_part), &
            weights(1:n_part), n_part, 1, id_guess(1:n_part))
       call a$D_particles_to_grid(tree, i_pos_ion, coords(:, 1:n_part), &
            weights(1:n_part), n_part, 1, id_guess(1:n_part))
    else
       call a$D_tree_clear_cc(tree, i_electron)
       call a$D_particles_to_grid(tree, i_electron, coords(:, 1:n_part), &
            weights(1:n_part), n_part, 1, id_guess(1:n_part))
    end if

    pc%particles(1:n_part)%id = id_guess(1:n_part)

    i = 0
    do n = 1, n_events
       if (events%list(n)%ctype == CS_ionize_t) then
          i = i + 1
          coords(:, i) = events%list(n)%part%x(1:$D)
          weights(i) = events%list(n)%part%w
          id_guess(i) = events%list(n)%part%id
       else if (events%list(n)%ctype == CS_attach_t) then
          i = i + 1
          coords(:, i) = events%list(n)%part%x(1:$D)
          weights(i) = -events%list(n)%part%w
          id_guess(i) = events%list(n)%part%id
       end if
    end do

    if (i > 0) then
       call a$D_particles_to_grid(tree, i_pos_ion, coords(:, 1:i), &
            weights(1:i), i, 1, id_guess(1:i))
    end if

    if (photoi_enabled) then
       call get_photoionization(events, coords, weights, n_photons)
       call a$D_particles_to_grid(tree, i_pos_ion, coords(:, 1:n_photons), &
            weights(1:n_photons), n_photons, 1)
       print *, "n_photons", n_photons, n_events
       do n = 1, n_photons
          x(1:$D) = coords(:, n)
#if $D == 2
          x(3)   = 0
#endif
          v      = 0
          a      = get_accel_pos(x(1:$D))
          id     = a$D_get_id_at(tree, x(1:$D))

          call pc%create_part(x, v, a, weights(n), 0.0_dp, id=id)
       end do
    end if

    events%n_stored = 0

  end subroutine particles_to_density


  function get_desired_weight(my_part) result(weight)
    type(PC_part_t), intent(in) :: my_part
    real(dp)                    :: weight, n_elec
    type(a$D_loc_t)              :: loc
    integer                     :: id, ix($D)

    loc = a$D_get_loc(tree, my_part%x(1:$D), my_part%id)
    id = loc%id
    ix = loc%ix

#if $D == 2
    n_elec = tree%boxes(id)%cc(ix(1), ix(2), i_electron) * &
         tree%boxes(id)%dr**2
#elif $D == 3
    n_elec = tree%boxes(id)%cc(ix(1), ix(2), ix(3), i_electron) * &
         tree%boxes(id)%dr**2
#endif

    weight = n_elec / particle_per_cell
    weight = max(particle_min_weight, min(particle_max_weight, weight))
  end function get_desired_weight

  subroutine print_status()
    write(*, "(F7.2,A,I0,A,E10.3,A,E10.3,A,E10.3,A,E10.3,A,E10.3)") &
         100 * ST_time / ST_end_time, "% it: ", it, &
         " t:", ST_time, " dt:", ST_dt, " wc:", wc_time, &
         " ncell:", real(a$D_num_cells_used(tree), dp), &
         " npart:", real(pc%get_num_sim_part(), dp)
  end subroutine print_status

  subroutine print_info()
    use m_units_constants
    real(dp) :: max_fld, max_elec, max_pion
    real(dp) :: sum_elec, sum_pos_ion
    real(dp) :: mean_en, n_elec
    integer  :: n_part
    call a$D_tree_max_cc(tree, i_E, max_fld)
    call a$D_tree_max_cc(tree, i_electron, max_elec)
    call a$D_tree_max_cc(tree, i_pos_ion, max_pion)
    call a$D_tree_sum_cc(tree, i_electron, sum_elec)
    call a$D_tree_sum_cc(tree, i_pos_ion, sum_pos_ion)
    mean_en = pc%get_mean_energy()
    n_part  = pc%get_num_sim_part()
    n_elec  = pc%get_num_real_part()

    print *, "max field", max_fld
    print *, "max elec/pion", max_elec, max_pion
    print *, "sum elec/pion", sum_elec, sum_pos_ion
    print *, "mean energy", mean_en / UC_elec_volt
    print *, "n_part, n_elec", n_part, n_elec
    print *, "mean weight", n_elec/n_part
  end subroutine print_info

  subroutine prepare_output_variables()
    integer :: n, n_part
    real(dp), allocatable :: coords(:, :)
    real(dp), allocatable :: weights(:)
    integer, allocatable  :: id_guess(:)

    n_part = pc%get_num_sim_part()
    allocate(coords($D, n_part))
    allocate(weights(n_part))
    allocate(id_guess(n_part))

    !$omp parallel do
    do n = 1, n_part
       coords(:, n) = pc%particles(n)%x(1:$D)
       weights(n) = 1.0_dp
       id_guess(n) = pc%particles(n)%id
    end do
    !$omp end parallel do

    call a$D_tree_clear_cc(tree, i_ppc)
    call a$D_particles_to_grid(tree, i_ppc, coords(:, 1:n_part), &
         weights(1:n_part), n_part, 0, id_guess(1:n_part), .false.)

    ! Fill ghost cells before writing output
    call a$D_gc_tree(tree, i_electron)
    call a$D_gc_tree(tree, i_pos_ion)

  end subroutine prepare_output_variables

end program apic_$Dd
