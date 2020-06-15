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
  use m_dielectric
  use m_ions

  implicit none

  integer, parameter     :: int8 = selected_int_kind(18)
  integer, parameter     :: ndim = NDIM
  integer(int8)          :: t_start, t_current, count_rate
  real(dp)               :: dt, dt_ions
  real(dp)               :: wc_time, inv_count_rate
  real(dp)               :: time_last_print, time_last_generate
  integer                :: it
  integer                :: n_part, n_prev_merge, n_ion_prev_merge, n_samples
  integer, allocatable   :: ref_links(:, :)
  character(len=GL_slen) :: fname
  logical                :: write_out, update_ions
  type(ref_info_t)       :: ref_info
  real(dp)               :: n_elec, n_elec_prev, max_elec_dens, max_ion_dens
  real(dp)               :: dt_cfl, dt_growth, dt_drt
  real(dp)               :: dt_ions_cfl, dt_ions_drt
  real(dp)               :: time_ions = 0.0_dp
  ! real(dp) :: max_surf_density,min_surf_density, max_potential


  real(dp) :: t0
  real(dp) :: wtime_start
  real(dp) :: wtime_run = 0.0_dp
  real(dp) :: wtime_advance = 0.0_dp
  integer :: output_cnt = 0 ! Number of output files written

  wtime_start = omp_get_wtime()

  ! Read command line arguments and configuration files
  call CFG_update_from_arguments(cfg)

  ! Initialize the user's code
  call user_initialize(cfg)

  ! Initialize other modules
  call domain_init(cfg)
  call GL_initialize(cfg, ndim)
  call refine_init(cfg, ndim)
  call time_step_init(cfg)
  call check_path_writable(trim(GL_output_dir))
  call field_initialize(cfg, mg)
  call init_particle(cfg, pc)
  call init_pc_ion(cfg, pc_ions)
  call photons_initialize(cfg)

  ! Write configuration to output
  fname = trim(GL_output_dir) // "/" // trim(GL_simulation_name) // "_out.cfg"
  call CFG_write(cfg, trim(fname))

  ! Initialize the tree (which contains all the mesh information)
  call init_tree(tree)

  ! Set the multigrid options. First define the variables to use
  mg%i_phi = i_phi
  mg%i_tmp = i_Ex
  mg%i_rhs = i_rhs
  if (GL_use_dielectric) mg%i_eps = i_eps

  ! This automatically handles cylindrical symmetry
  mg%box_op => mg_auto_op
  mg%box_gsrb => mg_auto_gsrb
  mg%box_corr => mg_auto_corr
  mg%box_stencil => mg_auto_stencil

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
  call af_set_cc_methods(tree, i_Ex, af_bc_neumann_zero)
  call af_set_cc_methods(tree, i_Ey, af_bc_neumann_zero)
#if NDIM == 3
  call af_set_cc_methods(tree, i_Ez, af_bc_neumann_zero)
#endif

if (GL_use_dielectric) then
   ! Initialize dielectric surfaces at the third refinement level
   call af_refine_up_to_lvl(tree, 3)
   call dielectric_initialize(tree, i_eps, diel, n_surf_vars)

   ! Initialize the dielectric surface charge
   if (associated(user_set_dielectric_charge)) &
      call dielectric_set_values(tree, diel, i_surf_pos_ion, user_set_dielectric_charge)

end if

  output_cnt         = 0 ! Number of output files written
  GL_time            = 0 ! Simulation time (all times are in s)
  time_ions          = 0
  time_last_generate = GL_time

  ! Set up the initial conditions
  if (.not. associated(user_initial_particles) .and. .not. associated(user_initial_particles_and_ions)) &
    error stop "user routine for particles (or ions) not defined"

  if (associated(user_initial_particles)) then
    call user_initial_particles(pc)
  end if
  if (associated(user_initial_particles_and_ions)) then
    call user_initial_particles_and_ions(pc, pc_ions)
  end if

  ! Perform additional refinement
  do
     call af_tree_clear_cc(tree, i_pos_ion)
     call particles_and_ions_to_density_and_events(tree, pc, pc_ions, .false., .true.)
     call field_compute(tree, mg, .false.)

     if (GL_use_dielectric) then
        ! Make sure there are no refinement jumps across the dielectric
        call dielectric_get_refinement_links(diel, ref_links)
        call af_adjust_refinement(tree, refine_routine, ref_info, &
             refine_buffer_width, ref_links)
        call dielectric_update_after_refinement(tree, diel, ref_info)
     else
        call af_adjust_refinement(tree, refine_routine, ref_info, &
             refine_buffer_width)
     end if
     if (ref_info%n_add == 0) exit
  end do


  call pc%set_accel()

  print *, "Number of threads", af_get_max_threads()
  call af_print_info(tree)

  ! Start from small time step
  GL_dt   = GL_dt_min
  dt_ions = 2*GL_dt_min

  ! Initial wall clock time
  call system_clock(t_start, count_rate)
  inv_count_rate = 1.0_dp / count_rate
  time_last_print = -1e10_dp
  n_prev_merge = pc%get_num_sim_part()
  n_elec_prev = pc%get_num_real_part()
  n_ion_prev_merge = pc_ions%get_num_sim_part()

  ! Start of time integration
  do it = 1, huge(1)-1
     if (GL_time >= GL_end_time) exit

     if (associated(user_generate_particles)) then
        call user_generate_particles(pc, GL_time, GL_time - time_last_generate)
        ! time_last_generate = GL_time
     end if

     if (pc%get_num_sim_part() == 0 .and. pc_ions%get_num_sim_part() == 0) then
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
        write_out   = .true.
        update_ions = .true.
        dt         = output_cnt * GL_dt_output - GL_time
        dt_ions    = output_cnt * GL_dt_output - time_ions
        output_cnt = output_cnt + 1
        if (abs((GL_time + dt) - (time_ions + dt_ions)) > 1.0e-20_dp) then
          print *, "Ion and electron times are not correctly synchronized within allowed tolerance"
          print *, "The difference equals: ", abs((GL_time + dt) - (time_ions + dt_ions))
          error stop
        end if
     else
        write_out = .false.
        dt        = GL_dt
        if (time_ions + dt_ions < GL_time + GL_dt) then
          update_ions = .true.
        else
          update_ions = .false.
        end if
     end if

     t0 = omp_get_wtime()
     call pc%advance_openmp(dt)
     call pc%after_mover(dt)
     if (update_ions) then
       call pc_ions%advance_openmp(dt_ions)
     end if
     wtime_advance = wtime_advance + omp_get_wtime() - t0

     GL_time = GL_time + dt
     if (update_ions) time_ions = time_ions + dt_ions

     call particles_and_ions_to_density_and_events(tree, pc, pc_ions, .false., update_ions)

     n_part = pc_ions%get_num_sim_part()
     ! if (n_part > n_ion_prev_merge * ion_min_merge_increase .and. update_ions) then
     if (update_ions) then ! Only when ions have moved do we need to reconsider weights
        call adapt_weights_ions(tree, pc_ions)
        n_ion_prev_merge = pc_ions%get_num_sim_part()
     end if

     n_part = pc%get_num_sim_part()
     if (n_part > n_prev_merge * min_merge_increase) then
        call adapt_weights(tree, pc)
        n_prev_merge = pc%get_num_sim_part()
     end if

     ! Compute field with new density
     call field_compute(tree, mg, .true.)
     call pc%set_accel()

     ! Time step for the ions
     n_samples = min(pc_ions%get_num_sim_part(), 1000)
     call af_tree_max_cc(tree, i_pos_ion, max_ion_dens)
     dt_ions_cfl = PM_get_max_dt(pc_ions, GL_rng, n_samples, 0.5_dp)
     dt_ions_drt = dielectric_relaxation_time(max_ion_dens, 1.5e-4_dp / GL_gas_pressure) ! max ion mobility
     dt_ions     = min(dt_ions_cfl, dt_ions_drt, GL_dt_output)

     ! Time step for the electrons
     n_samples = min(pc%get_num_sim_part(), 1000)
     call af_tree_max_cc(tree, i_electron, max_elec_dens)
     n_elec      = pc%get_num_real_part()
     if (n_elec > 0.0_dp) then
       dt_cfl      = PM_get_max_dt(pc, GL_rng, n_samples, cfl_particles)
       dt_drt      = dielectric_relaxation_time(max_elec_dens, 0.75_dp) ! max elec mobility
       dt_growth   = 1.0_dp!get_new_dt(GL_dt, abs(1-n_elec/n_elec_prev), 20.0e-2_dp)
       GL_dt       = min(dt_cfl, dt_growth, dt_drt)
     else
       GL_dt = dt_ions !if no electrons are present, relax dt to the ion time scale
     end if
     n_elec_prev = n_elec

     if (write_out) then
        call set_output_variables()

        write(fname, "(A,I6.6)") trim(GL_simulation_name) // "_", output_cnt
        call af_write_silo(tree, fname, output_cnt, GL_time, &
             dir=GL_output_dir, add_curve_names = ["EEDF"], &
             add_curve_dat = write_EEDF_as_curve(pc))
        call print_info()
     end if

     if (mod(it, refine_per_steps) == 0) then
        t0 = omp_get_wtime()
        if (GL_use_dielectric) then
           ! Make sure there are no refinement jumps across the dielectric
           call dielectric_get_refinement_links(diel, ref_links)
           call af_adjust_refinement(tree, refine_routine, ref_info, &
                refine_buffer_width, ref_links)
           call dielectric_update_after_refinement(tree, diel, ref_info)
        else
           call af_adjust_refinement(tree, refine_routine, ref_info, &
                refine_buffer_width)
        end if

        if (ref_info%n_add + ref_info%n_rm > 0) then
           ! Compute the field on the new mesh
           call particles_and_ions_to_density_and_events(tree, pc, pc_ions, .false., .false.)
           call field_compute(tree, mg, .true.)
           call adapt_weights(tree, pc)
           ! if (update_ions) call adapt_weights_ions(tree, pc_ions)
        end if
     end if
  end do

contains

  !> Initialize the AMR tree
  subroutine init_tree(tree)
    type(af_t), intent(inout) :: tree

    ! Initialize tree
    if (GL_cylindrical) then
       call af_init(tree, box_size, domain_len, &
            coarse_grid_size, coord=af_cyl)
    else
       call af_init(tree, box_size, domain_len, &
            coarse_grid_size)
    end if

  end subroutine init_tree

  subroutine print_status()
    write(*, "(F7.2,A,I0,A,E9.2,A,E9.2,A,E9.2,A,E9.2,A,E9.2)") &
         100 * GL_time / GL_end_time, "% it: ", it, &
         " t:", GL_time, " dt:", GL_dt, " wc:", wc_time, &
         " ncell:", real(af_num_cells_used(tree), dp), &
         " npart:", real(pc%get_num_sim_part(), dp)
  end subroutine print_status

  subroutine print_info()
    use m_units_constants
    real(dp) :: max_fld, max_elec, max_pion
    real(dp) :: sum_elec, sum_pos_ion
    real(dp) :: mean_en, n_elec, n_part
    real(dp) :: n_part_ions, n_ions
    real(dp) :: surf_int

    call af_tree_max_cc(tree, i_E, max_fld)
    call af_tree_max_cc(tree, i_electron, max_elec)
    call af_tree_max_cc(tree, i_pos_ion, max_pion)
    call af_tree_sum_cc(tree, i_electron, sum_elec)
    call af_tree_sum_cc(tree, i_pos_ion, sum_pos_ion)
    mean_en = pc%get_mean_energy()
    n_part  = pc%get_num_sim_part()
    n_elec  = pc%get_num_real_part()
    n_part_ions = pc_ions%get_num_sim_part()
    n_ions  = pc_ions%get_num_real_part()

    write(*, "(A20,E12.4)") "dt", GL_dt
    write(*, "(A20,E12.4)") "dt_cfl", dt_cfl
    write(*, "(A20,E12.4)") "dt_drt", dt_drt
    write(*, "(A20,E12.4)") "dt_growth", dt_growth
    write(*, "(A20,E12.4)") "dt_ions", dt_ions
    write(*, "(A20,E12.4)") "max field", max_fld
    write(*, "(A20,2E12.4)") "max elec/pion", max_elec, max_pion
    write(*, "(A20,2E12.4)") "sum elec/pion", sum_elec, sum_pos_ion
    call dielectric_get_integral(diel, i_surf_sum_dens, surf_int)
    write(*, "(A20,E12.4)") "Net charge", sum_pos_ion - sum_elec + surf_int
    write(*, "(A20,E12.4)") "mean energy", mean_en / UC_elec_volt
    write(*, "(A20,2E12.4)") "n_part_elec, n_elec", n_part, n_elec
    write(*, "(A20,2E12.4)") "n_part_ion, n_ion", n_part_ions, n_ions
    write(*, "(A20,E12.4)") "mean weight elec", n_elec/n_part
    write(*, "(A20,E12.4)") "mean weight ions", n_ions/n_part_ions

    wtime_run = omp_get_wtime() - wtime_start
    write(*, "(A20,F8.2,A)") "cost of advance", &
         1e2 * wtime_advance / wtime_run, "%"
  end subroutine print_info

  subroutine set_output_variables()
    use m_units_constants
    integer :: n, n_part
    real(dp), allocatable :: coords(:, :)
    real(dp), allocatable :: weights(:)
    real(dp), allocatable :: energy(:)
    integer, allocatable  :: id_guess(:)

    n_part = pc%get_num_sim_part()
    allocate(coords(NDIM, n_part))
    allocate(weights(n_part))
    allocate(energy(n_part))
    allocate(id_guess(n_part))

    !$omp parallel do
    do n = 1, n_part
       coords(:, n) = pc%particles(n)%x(1:NDIM)
       weights(n) = 1.0_dp
       energy(n) = pc%particles(n)%w * &
            PC_v_to_en(pc%particles(n)%v, UC_elec_mass) / &
            UC_elec_volt
       id_guess(n) = pc%particles(n)%id
    end do
    !$omp end parallel do

    call af_tree_clear_cc(tree, i_ppc)
    ! Don't divide by cell volume (last .false. argument)
    call af_particles_to_grid(tree, i_ppc, coords(:, 1:n_part), &
         weights(1:n_part), n_part, 0, id_guess(1:n_part), &
         density=.false., fill_gc=.false.)

    call af_tree_clear_cc(tree, i_energy)
    call af_particles_to_grid(tree, i_energy, coords(:, 1:n_part), &
         energy(1:n_part), n_part, 1, id_guess(1:n_part), &
         fill_gc=.false.)
    call af_tree_apply(tree, i_energy, i_electron, '/', 1e-10_dp)

    ! Fill ghost cells before writing output
    call af_gc_tree(tree, [i_electron, i_pos_ion])

  end subroutine set_output_variables

  subroutine check_path_writable(pathname)
    character(len=*), intent(in) :: pathname
    integer                      :: my_unit, iostate

    open(newunit=my_unit, file=trim(pathname)//"/DUMMY", iostat=iostate)
    if (iostate /= 0) then
       print *, "Output directory: " // trim(pathname)
       error stop "Directory not writable (does it exist?)"
    else
       close(my_unit, status='delete')
    end if
  end subroutine check_path_writable

end program apic
