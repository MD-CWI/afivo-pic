#include "../afivo/src/cpp_macros.h"
!> Program to perform discharge simulations with particle-in-cell method
program apic

  use omp_lib
  use m_af_all
  use m_globals
  use m_field
  use m_particle_core
  use m_domain
  use m_refine
  use m_time_step
  use m_particles
  use m_user
  use m_user_methods
  use m_cross_sec
  use m_init_cond
  use m_output

  implicit none

  integer, parameter     :: int8 = selected_int_kind(18)
  integer, parameter     :: ndim = NDIM
  integer(int8)          :: t_start, t_current, count_rate
  real(dp)               :: dt
  real(dp)               :: wc_time, inv_count_rate
  real(dp)               :: time_last_print, time_last_generate
  integer                :: it, it_last_merge
  integer                :: n_part, n_prev_merge, n_samples
  integer                :: lvl, i, id
  integer, allocatable   :: ref_links(:, :)
  character(len=GL_slen) :: fname
  logical                :: write_out
  type(ref_info_t)       :: ref_info
  real(dp)               :: n_elec, n_elec_prev, max_elec_dens
  real(dp)               :: dt_cfl, dt_growth, dt_drt

  real(dp) :: t0, t1, t_sort, t_rest
  real(dp) :: wtime_start
  real(dp) :: wtime_run = 0.0_dp
  real(dp) :: wtime_advance = 0.0_dp
  real(dp) :: wtime_mergesplit = 0.0_dp
  real(dp) :: wtime_sort = 0.0_dp
  real(dp) :: wtime_todensity = 0.0_dp
  real(dp) :: wtime_field = 0.0_dp
  real(dp) :: wtime_io = 0.0_dp
  real(dp) :: wtime_accel = 0.0_dp
  real(dp) :: wtime_amr = 0.0_dp
  integer :: output_cnt = 0 ! Number of output files written

  ! Read command line arguments and configuration files
  call CFG_update_from_arguments(cfg)

  ! Initialize the user's code
  call user_initialize(cfg)

  ! Initialize other modules
  call domain_init(cfg)
  call time_step_init(cfg)
  call GL_initialize(cfg, ndim)
  call check_path_writable(trim(GL_output_dir), trim(GL_simulation_name))
  call field_initialize(cfg, mg)
  call init_particle(cfg, pc)
  call refine_init(cfg, ndim)
  call photons_initialize(cfg)
  call init_cond_initialize(cfg)

  ! Write configuration to output
  fname = trim(GL_output_dir) // "/" // trim(GL_simulation_name) // "_out.cfg"
  call CFG_write(cfg, trim(fname))

  ! Set the multigrid options. First define the variables to use
  mg%i_phi = i_phi
  mg%i_tmp = i_residual
  mg%i_rhs = i_rhs
  if (GL_use_dielectric) tree%mg_i_eps = i_eps

  if (GL_use_electrode) then
     if (any(field_rod_r0 <= -1.0_dp)) &
          error stop "field_rod_r0 not set correctly"
     if (any(field_rod_r1 <= -1.0_dp)) &
          error stop "field_rod_r1 not set correctly"
     if (field_rod_radius <= 0) &
          error stop "field_rod_radius not set correctly"

     call af_set_cc_methods(tree, i_lsf, af_bc_neumann_zero, &
          af_gc_prolong_copy, af_prolong_zeroth, &
          funcval=field_set_lsf_box)

     mg%lsf => field_get_lsf
     tree%mg_i_lsf = i_lsf

     mg%lsf_dist => mg_lsf_dist_gss
     mg%lsf_length_scale = field_rod_radius
  end if

  ! Initialize the tree (which contains all the mesh information)
  call init_tree(tree)

  ! Set the dielectric permittivity before initializing multigrid
  if (GL_use_dielectric) then
     if (.not. associated(user_set_dielectric_eps)) &
          error stop "user_set_dielectric_eps not defined"
     call af_loop_box(tree, user_set_dielectric_eps)
  end if

  ! This routine always needs to be called when using multigrid
  call mg_init(tree, mg)

  call af_set_cc_methods(tree, i_electron, af_bc_neumann_zero)
  call af_set_cc_methods(tree, i_pos_ion, af_bc_neumann_zero, &
       prolong=af_prolong_limit)
  call af_set_cc_methods(tree, i_E, af_bc_neumann_zero)

  if (GL_use_dielectric) then
     ! Initialize dielectric surfaces at the third refinement level
     call af_refine_up_to_lvl(tree, 3)
     call surface_initialize(tree, i_eps, diel, n_surf_vars)
  end if

  output_cnt         = 0 ! Number of output files written
  GL_time            = 0 ! Simulation time (all times are in s)
  time_last_generate = GL_time

  ! Set up the initial conditions
  if (associated(user_initial_particles)) &
       call user_initial_particles(pc)

  ! Perform additional refinement
  do i = 1, 100
     call af_tree_clear_cc(tree, i_pos_ion)
     call particles_to_density_and_events(tree, pc, .true.)

     if (associated(user_initial_ion_density)) &
          call af_loop_box(tree, user_initial_ion_density)

     call field_compute(tree, mg, .false.)

     if (GL_use_dielectric) then
        ! Make sure there are no refinement jumps across the dielectric
        call surface_get_refinement_links(diel, ref_links)
        call af_adjust_refinement(tree, refine_routine, ref_info, &
             refine_buffer_width, ref_links)
        call surface_update_after_refinement(tree, diel, ref_info)
     else
        call af_adjust_refinement(tree, refine_routine, ref_info, &
             refine_buffer_width)
     end if
     if (ref_info%n_add == 0) exit
  end do

  ! Add particles from the initial condition. Note that this cannot be done in
  ! parallel at the moment
  if (any(init_conds%seed_density > 0) .or. &
       init_conds%background_density > 0) then
     do lvl = 1, tree%highest_lvl
        do i = 1, size(tree%lvls(lvl)%leaves)
           id = tree%lvls(lvl)%leaves(i)
           call init_cond_set_box(tree%boxes(id))
        end do
     end do

     call af_tree_clear_cc(tree, i_pos_ion)
     call particles_to_density_and_events(tree, pc, .true.)
  end if

  ! Set an initial surface charge
  if (associated(user_set_surface_charge)) then
     call surface_set_values(tree, diel, i_surf_pos_ion, user_set_surface_charge)
  end if

  call pc%set_accel()

  write(*, "(A,I12)") " Number of threads:       ", af_get_max_threads()
  call af_print_info(tree)

  ! Start from small time step
  GL_dt   = GL_dt_min

  ! Initial wall clock time
  call system_clock(t_start, count_rate)
  inv_count_rate = 1.0_dp / count_rate
  time_last_print = -1e10_dp
  it_last_merge = 0
  n_prev_merge = pc%get_num_sim_part()
  n_elec_prev = pc%get_num_real_part()

  wtime_start = omp_get_wtime()

  ! Start of time integration
  do it = 1, huge(1)-1
     if (GL_time >= GL_end_time) exit

     if (associated(user_generate_particles)) then
        call user_generate_particles(pc, GL_time, GL_time - time_last_generate)
        time_last_generate = GL_time
     end if

     if (pc%get_num_sim_part() == 0) then
        print *, "No particles, end of simulation"
        exit
     end if

     call system_clock(t_current)
     wc_time = (t_current - t_start) * inv_count_rate

     ! Every GL_print_status_interval, print some info about progress
     if (wc_time - time_last_print > GL_print_status_sec) then
        call print_status()
        time_last_print = wc_time
     end if

     ! Every GL_dt_output, write output
     if (output_cnt * GL_dt_output <= GL_time + GL_dt) then
        write_out  = .true.
        dt         = output_cnt * GL_dt_output - GL_time
        output_cnt = output_cnt + 1
     else
        write_out = .false.
        dt        = GL_dt
     end if

     t0 = omp_get_wtime()
     call pc%advance_openmp(dt)
     if (.not. magnetic_field_used) call pc%after_mover(dt)
     t1 = omp_get_wtime()
     wtime_advance = wtime_advance + (t1 - t0)

     GL_time = GL_time + dt

     call particles_to_density_and_events(tree, pc, .false.)
     t0 = omp_get_wtime()
     wtime_todensity = wtime_todensity + (t0 - t1)

     n_part = pc%get_num_sim_part()
     if (n_part > n_prev_merge * min_merge_increase .or. &
          it - it_last_merge >= iterations_between_merge_split) then
        call adapt_weights(tree, pc, t_sort, t_rest)

        n_prev_merge     = pc%get_num_sim_part()
        it_last_merge    = it
        wtime_sort       = wtime_sort + t_sort
        wtime_mergesplit = wtime_mergesplit + t_rest
     end if

     ! Compute field with new density
     t0 = omp_get_wtime()
     call field_compute(tree, mg, .true.)
     t1 = omp_get_wtime()
     wtime_field = wtime_field + (t1 - t0)

     call pc%set_accel()
     t0 = omp_get_wtime()
     wtime_accel = wtime_accel + (t0 - t1)

     n_samples = min(n_part, 1000)
     call af_tree_max_cc(tree, i_electron, max_elec_dens)
     n_elec      = pc%get_num_real_part()
     dt_cfl      = PM_get_max_dt(pc, GL_rng, n_samples, cfl_particles)
     dt_drt      = dielectric_relaxation_time(max_elec_dens)
     dt_growth   = get_new_dt(GL_dt, abs(1-n_elec/n_elec_prev), 20.0e-2_dp)
     GL_dt       = min(dt_cfl, dt_growth, dt_drt)
     n_elec_prev = n_elec

     if (write_out) then
        t0 = omp_get_wtime()
        call set_output_variables()

        write(fname, "(A,I6.6)") trim(GL_output_dir) // "/" // &
             trim(GL_simulation_name) // "_", output_cnt
        call af_write_silo(tree, fname, output_cnt, GL_time, &
             add_curve_names = ["EEDF"], &
             add_curve_dat = write_EEDF_as_curve(pc))
        call print_info()
        call CS_write_ledger(pc%coll_ledger, &
        trim(GL_output_dir) // "/" // trim(GL_simulation_name) // "_cs_ledger.txt", &
        GL_time)

        ! output the log file
        write(fname, "(A,I6.6)") trim(GL_output_dir) // "/" // trim(GL_simulation_name) // "_log.txt"
        if (associated(user_write_log)) then
          call user_write_log(tree, fname, output_cnt)
        else
          call output_log(tree, fname, output_cnt, wc_time)
        end if

        if (GL_write_to_dat) then
          if (GL_write_to_dat_interval(1) .le. GL_time .and. GL_time .le. GL_write_to_dat_interval(2)) then
            call af_write_tree(tree, trim(GL_output_dir) // "/" // trim(fname), write_sim_data)
          end if
        end if
        t1 = omp_get_wtime()
        wtime_io = wtime_io + (t1 - t0)
     end if

     if (mod(it, refine_per_steps) == 0) then
        t0 = omp_get_wtime()
        if (GL_use_dielectric) then
           ! Make sure there are no refinement jumps across the dielectric
           call surface_get_refinement_links(diel, ref_links)
           call af_adjust_refinement(tree, refine_routine, ref_info, &
                refine_buffer_width, ref_links)
           call surface_update_after_refinement(tree, diel, ref_info)
        else
           call af_adjust_refinement(tree, refine_routine, ref_info, &
                refine_buffer_width)
        end if
        t1 = omp_get_wtime()
        wtime_amr = wtime_amr + (t1 - t0)

        if (ref_info%n_add + ref_info%n_rm > 0) then
           ! Compute the field on the new mesh
           call particles_to_density_and_events(tree, pc, .false.)
           t0 = omp_get_wtime()
           wtime_todensity = wtime_todensity + (t0 - t1)
           call field_compute(tree, mg, .true.)
           t1 = omp_get_wtime()
           wtime_field = wtime_field + (t1 - t0)
           call adapt_weights(tree, pc, t_sort, t_rest)

           n_prev_merge     = pc%get_num_sim_part()
           it_last_merge    = it
           wtime_sort       = wtime_sort + t_sort
           wtime_mergesplit = wtime_mergesplit + t_rest
        end if
     end if
  end do

contains

  !> Initialize the AMR tree
  subroutine init_tree(tree)
    type(af_t), intent(inout) :: tree
    integer                   :: coord

    if (GL_cylindrical) then
       coord = af_cyl
    else
       coord = af_xyz
    end if

    ! Initialize tree
    call af_init(tree, box_size, domain_len, coarse_grid_size, &
         coord=coord, mem_limit_gb=GL_memory_afivo_GB)
  end subroutine init_tree

  subroutine print_status()
    write(*, "(F7.2,A,I0,A,E9.2,A,E9.2,A,E9.2,A,E9.2,A,E9.2)") &
         100 * GL_time / GL_end_time, "% it: ", it, &
         " t:", GL_time, " dt:", GL_dt, " wc:", wc_time, &
         " ncell:", real(af_num_cells_used(tree), dp), &
         " npart:", real(pc%get_num_sim_part(), dp)
  end subroutine print_status

  subroutine write_sim_data(my_unit)
    integer, intent(in) :: my_unit
    real(dp) :: time, global_time, photoi_prev_time, global_dt

    time = GL_time
    global_time = GL_time
    photoi_prev_time = GL_time
    global_dt = GL_dt

    write(my_unit) output_cnt
    write(my_unit) time
    write(my_unit) global_time
    write(my_unit) photoi_prev_time
    write(my_unit) global_dt
  end subroutine write_sim_data

  subroutine print_info()
    use m_units_constants
    real(dp) :: max_fld, max_elec, max_pion
    real(dp) :: sum_elec, sum_pos_ion
    real(dp) :: mean_en, n_elec, n_part
    real(dp) :: surf_int

    call af_tree_max_cc(tree, i_E, max_fld)
    call af_tree_max_cc(tree, i_electron, max_elec)
    call af_tree_max_cc(tree, i_pos_ion, max_pion)
    call af_tree_sum_cc(tree, i_electron, sum_elec)
    call af_tree_sum_cc(tree, i_pos_ion, sum_pos_ion)
    mean_en = pc%get_mean_energy()
    n_part  = pc%get_num_sim_part()
    n_elec  = pc%get_num_real_part()

    write(*, "(A20,E12.4)") "dt", GL_dt
    write(*, "(A20,E12.4)") "max field", max_fld
    write(*, "(A20,2E12.4)") "max elec/pion", max_elec, max_pion
    write(*, "(A20,2E12.4)") "sum elec/pion", sum_elec, sum_pos_ion
    call surface_get_integral(diel, i_surf_sum_dens, surf_int)
    write(*, "(A20,E12.4)") "Net charge", sum_pos_ion - sum_elec + surf_int
    write(*, "(A20,E12.4)") "mean energy", mean_en / UC_elec_volt
    write(*, "(A20,2E12.4)") "n_part, n_elec", n_part, n_elec
    write(*, "(A20,E12.4)") "mean weight", n_elec/n_part

    wtime_run = omp_get_wtime() - wtime_start
    write(*, "(8A10)") "advance", "to_grid", "weights", "sort", &
         "field", "accel", "amr", "io"
    write(*, "(8F10.2)") 1e2 * [wtime_advance, wtime_todensity, &
         wtime_mergesplit, wtime_sort, wtime_field, wtime_accel, &
         wtime_amr, wtime_io] / wtime_run
  end subroutine print_info

  subroutine set_output_variables()
    use m_units_constants
    integer :: n_part

    n_part = pc%get_num_sim_part()

    call af_tree_clear_cc(tree, i_ppc)
    ! Don't divide by cell volume (last .false. argument)
    call af_particles_to_grid(tree, i_ppc, n_part, get_id, get_r_unit_dens, &
         0, density=.false., fill_gc=.false.)

    call af_tree_clear_cc(tree, i_energy)
    call af_particles_to_grid(tree, i_energy, n_part, get_id, &
         get_r_energy, 1, fill_gc=.false., iv_tmp=i_tmp_dens)
    call af_tree_apply(tree, i_energy, i_electron, '/', 1e-10_dp)

    if (GL_use_dielectric) then
       call clear_cc_inside_dielectric(tree, i_energy)
    end if

    ! Fill ghost cells before writing output
    call af_gc_tree(tree, [i_electron, i_pos_ion])

  end subroutine set_output_variables

  subroutine check_path_writable(pathname, filename)
    character(len=*), intent(in) :: pathname, filename
    integer                      :: my_unit, iostate

    open(newunit=my_unit, file=trim(pathname)//"/"//trim(filename)// &
         "_DUMMY", iostat=iostate)
    if (iostate /= 0) then
       print *, "Output directory: " // trim(pathname)
       error stop "Directory not writable (does it exist?)"
    else
       close(my_unit, status='delete')
    end if
  end subroutine check_path_writable

  subroutine clear_cc_inside_dielectric(tree, iv)
    type(af_t), intent(inout) :: tree
    integer, intent(in)        :: iv !< Variable to clear
    integer                    :: lvl, i, id

    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
            if (tree%boxes(id)%cc(1, 1, i_eps) > 1) then
              call af_box_clear_cc(tree%boxes(id), iv)
            end if
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine clear_cc_inside_dielectric

end program apic
