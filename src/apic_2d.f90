!> Program to perform 2d discharge simulations in Cartesian and cylindrical coordinates
program apic_2d

  use m_a2_all
  use m_streamer
  use m_field_2d
  use m_init_cond_2d
  use m_particle_core

  implicit none

  integer, parameter     :: int8 = selected_int_kind(18)
  integer(int8)          :: t_start, t_current, count_rate
  real(dp)               :: dt, sum_elec, sum_pos_ion
  real(dp)               :: wc_time, inv_count_rate, time_last_print
  integer                :: it, id, n_part_id
  character(len=ST_slen) :: fname
  logical                :: write_out
  type(CFG_t)            :: cfg  ! The configuration for the simulation
  type(a2_t)             :: tree ! This contains the full grid information
  type(mg2_t)            :: mg   ! Multigrid option struct
  type(PC_t)             :: pc
  type(PC_events_t)      :: events
  type(ref_info_t)       :: ref_info

  real(dp), parameter :: PM_part_per_cell = 100
  real(dp), parameter :: PM_max_weight    = 1e20_dp

  integer :: output_cnt = 0 ! Number of output files written

  call CFG_update_from_arguments(cfg)
  call ST_initialize(cfg, 2)

  call field_initialize(cfg, mg)
  call init_cond_initialize(cfg, 2)

  call init_particle(cfg, pc)

  fname = trim(ST_output_dir) // "/" // trim(ST_simulation_name) // "_out.cfg"
  call CFG_write(cfg, trim(fname))

  ! Initialize the tree (which contains all the mesh information)
  call init_tree(tree)

  ! Set the multigrid options. First define the variables to use
  mg%i_phi = i_phi
  mg%i_tmp = i_Ex
  mg%i_rhs = i_rhs

  ! This automatically handles cylindrical symmetry
  mg%box_op => mg2_auto_op
  mg%box_gsrb => mg2_auto_gsrb
  mg%box_corr => mg2_auto_corr

  ! This routine always needs to be called when using multigrid
  call mg2_init_mg(mg)

  call a2_set_cc_methods(tree, i_electron, a2_bc_neumann_zero)
  call a2_set_cc_methods(tree, i_pos_ion, a2_bc_neumann_zero)
  call a2_set_cc_methods(tree, i_phi, mg%sides_bc, mg%sides_rb)

  output_cnt      = 0         ! Number of output files written
  ST_time         = 0         ! Simulation time (all times are in s)

  ! Set up the initial conditions
  call init_cond_particles(tree, pc)

  do
     call a2_tree_clear_cc(tree, i_pos_ion)
     call a2_loop_box(tree, init_cond_set_box)
     call particles_to_density(tree, pc, events, .true.)
     call field_compute(tree, mg, .false.)
     call a2_adjust_refinement(tree, refine_routine, ref_info, &
          ST_refine_buffer_width, .true.)
     if (ref_info%n_add == 0) exit
  end do

  call pc%set_accel()

  print *, "Number of threads", af_get_max_threads()
  call a2_print_info(tree)

  ! Start from small time step
  ST_dt   = ST_dt_min

  ! Initial wall clock time
  call system_clock(t_start, count_rate)
  inv_count_rate = 1.0_dp / count_rate
  time_last_print = -1e10_dp

  do it = 1, huge(1)-1
     if (ST_time >= ST_end_time) exit
     ! call a2_tree_sum_cc(tree, i_electron, sum_elec)
     ! call a2_tree_sum_cc(tree, i_pos_ion, sum_pos_ion)

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

     call pc%advance_openmp(dt, events)

     ST_time = ST_time + dt
     call particles_to_density(tree, pc, events, .false.)

     ! Compute field with new density
     call field_compute(tree, mg, .true.)

     call PC_verlet_correct_accel(pc, dt)

     if (modulo(it, 10) == 0) then
        call adapt_weights(tree, pc)
     end if

     if (write_out) then
        ! Fill ghost cells before writing output
        call a2_gc_tree(tree, i_electron, a2_gc_interp_lim, a2_bc_neumann_zero)
        call a2_gc_tree(tree, i_pos_ion, a2_gc_interp_lim, a2_bc_neumann_zero)

        write(fname, "(A,I6.6)") trim(ST_simulation_name) // "_", output_cnt
        call a2_write_silo(tree, fname, output_cnt, ST_time, &
             vars_for_output, dir=ST_output_dir)
     end if

     if (mod(it, ST_refine_per_steps) == 0) then
        call a2_adjust_refinement(tree, refine_routine, ref_info, &
             ST_refine_buffer_width, .true.)

        if (ref_info%n_add > 0) then
           ! Compute the field on the new mesh
           call particles_to_density(tree, pc, events, .false.)
           call field_compute(tree, mg, .true.)
        end if
     end if
  end do

contains

  subroutine adapt_weights(tree, pc)
    type(a2_t), intent(in)    :: tree
    type(PC_t), intent(inout) :: pc
    integer, allocatable      :: id_ipart(:)

    call sort_by_id(tree, pc, id_ipart)

    ! print *, "befor", pc%get_num_sim_part(), pc%get_num_real_part()
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
    ! print *, "after", pc%get_num_sim_part(), pc%get_num_real_part()
    ! print *, "after", pc%get_num_real_part()/pc%get_num_sim_part()
  end subroutine adapt_weights

  subroutine sort_by_id(tree, pc, id_ipart)
    type(a2_t), intent(in)              :: tree
    type(PC_t), intent(inout)           :: pc
    integer, intent(inout), allocatable :: id_ipart(:)

    integer                      :: n, n_part, id, new_ix
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
    type(a2_t), intent(inout) :: tree

    ! Variables used below to initialize tree
    real(dp)                  :: dr
    integer                   :: id
    integer                   :: ix_list(2, 1) ! Spatial indices of initial boxes

    dr = ST_domain_len / ST_box_size

    ! Initialize tree
    if (ST_cylindrical) then
       call a2_init(tree, ST_box_size, n_var_cell, n_var_face, dr, &
            coarsen_to=2, coord=af_cyl, &
            cc_names=ST_cc_names)
    else
       call a2_init(tree, ST_box_size, n_var_cell, n_var_face, dr, &
            coarsen_to=2, cc_names=ST_cc_names)
    end if

    ! Set up geometry
    id             = 1          ! One box ...
    ix_list(:, id) = 1          ! With index 1,1 ...

    ! Create the base mesh
    call a2_set_base(tree, 1, ix_list)

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

    where (pc%colls(:)%type == CS_ionize_t)
       pc%coll_is_event(:) = .true.
    end where

  end subroutine init_particle

  function get_accel(my_part) result(accel)
    use m_units_constants
    type(PC_part_t), intent(in) :: my_part
    real(dp)                    :: accel(3)

    accel(1:2) = a2_interp1(tree, my_part%x(1:2), [i_Ex, i_Ey], 2)
    accel(1:2) = accel(1:2) * UC_elec_q_over_m
    accel(3) = 0.0_dp
  end function get_accel

  subroutine particles_to_density(tree, pc, events, init_cond)
    use m_cross_sec
    type(a2_t), intent(inout)        :: tree
    type(PC_t), intent(inout)        :: pc
    type(PC_events_t), intent(inout) :: events
    logical, intent(in)              :: init_cond

    integer               :: n, i, n_part, n_events
    real(dp), allocatable :: coords(:, :)
    real(dp), allocatable :: weights(:)
    integer, allocatable  :: id_guess(:)

    n_part = pc%get_num_sim_part()
    n_events = events%n_stored
    n = max(n_part, n_events)

    allocate(coords(2, n))
    allocate(weights(n))
    allocate(id_guess(n))

    !$omp parallel do
    do n = 1, n_part
       coords(:, n) = pc%particles(n)%x(1:2)
       weights(n) = pc%particles(n)%w
       id_guess(n) = pc%particles(n)%id
    end do
    !$omp end parallel do

    if (init_cond) then
       call a2_tree_clear_cc(tree, i_electron)
       call a2_particles_to_grid(tree, i_electron, coords(:, 1:n_part), &
            weights(1:n_part), n_part, 1, id_guess(1:n_part))
       call a2_particles_to_grid(tree, i_pos_ion, coords(:, 1:n_part), &
            weights(1:n_part), n_part, 1, id_guess(1:n_part))
    else
       call a2_tree_clear_cc(tree, i_electron)
       call a2_particles_to_grid(tree, i_electron, coords(:, 1:n_part), &
            weights(1:n_part), n_part, 1, id_guess(1:n_part))
    end if

    pc%particles(1:n_part)%id = id_guess(1:n_part)

    i = 0
    do n = 1, n_events
       if (events%list(n)%ctype == CS_ionize_t) then
          i = i + 1
          coords(:, i) = events%list(n)%part%x(1:2)
          weights(i) = events%list(n)%part%w
          id_guess(i) = events%list(n)%part%id
       end if
    end do

    if (i > 0) then
       call a2_particles_to_grid(tree, i_pos_ion, coords(:, 1:i), &
            weights(1:i), i, 1, id_guess(1:i))
    end if

    events%n_stored = 0

  end subroutine particles_to_density

  ! This routine sets the cell refinement flags for box
  subroutine refine_routine(box, cell_flags)
    use m_geometry
    use m_init_cond_2d
    type(box2_t), intent(in) :: box
    ! Refinement flags for the cells of the box
    integer, intent(out)     :: &
         cell_flags(box%n_cell, box%n_cell)
    integer                  :: i, j, n, nc
    real(dp)                 :: dx, dx2, fld, alpha, adx, cphi, elec_dens
    real(dp)                 :: dist, rmin(2), rmax(2)

    nc = box%n_cell
    dx = box%dr
    dx2 = dx**2

    do j = 1, nc; do i = 1, nc
       fld   = box%cc(i, j, i_E)
       alpha = LT_get_col(ST_td_tbl, i_td_alpha, fld)
       adx   = box%dr * alpha

       ! The refinement is also based on the intensity of the source term.
       ! Here we estimate the curvature of phi (given by dx**2 *
       ! Laplacian(phi))
       cphi = dx2 * abs(box%cc(i, j, i_rhs))

       elec_dens = box%cc(i, j, i_electron)

       if (adx > ST_refine_adx .and. elec_dens > ST_refine_elec_dens) then
          cell_flags(i, j) = af_do_ref
       else if (adx < 0.125_dp * ST_refine_adx .or. &
            elec_dens < 0.125_dp * ST_refine_elec_dens) then
          cell_flags(i, j) = af_rm_ref
       else
          cell_flags(i, j) = af_keep_ref
       end if
    end do; end do

    ! Check fixed refinements
    rmin = box%r_min
    rmax = box%r_min + box%dr * box%n_cell

    do n = 1, size(ST_refine_regions_dr)
       if (ST_time <= ST_refine_regions_tstop(n) .and. &
            dx > ST_refine_regions_dr(n) .and. all(&
            rmax >= ST_refine_regions_rmin(:, n) .and. &
            rmin <= ST_refine_regions_rmax(:, n))) then
          ! Mark just the center cell to prevent refining neighbors
          cell_flags(nc/2, nc/2) = af_do_ref
       end if
    end do

    ! Make sure we don't have or get a too fine or too coarse grid
    if (dx > ST_refine_max_dx) then
       cell_flags = af_do_ref
    else if (dx < 2 * ST_refine_min_dx) then
       where(cell_flags == af_do_ref) cell_flags = af_keep_ref
    end if

  end subroutine refine_routine

  function get_desired_weight(my_part) result(weight)
    type(PC_part_t), intent(in) :: my_part
    real(dp)                    :: weight, n_elec
    type(a2_loc_t)              :: loc
    integer                     :: id, ix(2)

    loc = a2_get_loc(tree, my_part%x(1:2), my_part%id)
    id = loc%id
    ix = loc%ix

    n_elec = tree%boxes(id)%cc(ix(1), ix(2), i_electron) * tree%boxes(id)%dr**2
    weight = n_elec / PM_part_per_cell
    weight = max(1.0_dp, min(PM_max_weight, weight))
    ! print *, n_elec, weight, my_part%w
  end function get_desired_weight

  real(dp) function id_as_real(my_part)
    type(PC_part_t), intent(in) :: my_part
    id_as_real = my_part%id
  end function id_as_real

  subroutine print_status()
    write(*, "(F7.2,A,I0,A,E10.3,A,E10.3,A,E10.3,A,E10.3,A,E10.3)") &
         100 * ST_time / ST_end_time, "% it: ", it, &
         " t:", ST_time, " dt:", ST_dt, " wc:", wc_time, &
         " ncell:", real(a2_num_cells_used(tree), dp), &
         " npart:", real(pc%get_num_sim_part(), dp)
  end subroutine print_status

end program apic_2d
