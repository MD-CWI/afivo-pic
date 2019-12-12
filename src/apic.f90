#include "../afivo/src/cpp_macros.h"
!> Program to perform discharge simulations in Cartesian and cylindrical coordinates
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

  implicit none

  integer, parameter     :: int8 = selected_int_kind(18)
  integer, parameter     :: ndim = NDIM
  integer(int8)          :: t_start, t_current, count_rate
  real(dp)               :: dt
  real(dp)               :: wc_time, inv_count_rate, time_last_print
  integer                :: it
  integer                :: n_part, n_prev_merge, n_samples
  character(len=GL_slen) :: fname
  logical                :: write_out
  type(ref_info_t)       :: ref_info
  real(dp)               :: n_elec, n_elec_prev, max_elec_dens
  real(dp)               :: dt_cfl, dt_growth, dt_drt

  real(dp) :: t0
  real(dp) :: wtime_start
  real(dp) :: wtime_run = 0.0_dp
  real(dp) :: wtime_advance = 0.0_dp
  real(dp) :: wtime_density = 0.0_dp
  real(dp) :: wtime_merge = 0.0_dp
  real(dp) :: wtime_field = 0.0_dp
  real(dp) :: wtime_io = 0.0_dp
  real(dp) :: wtime_amr = 0.0_dp
  integer :: output_cnt = 0 ! Number of output files written

  wtime_start = omp_get_wtime()

  call CFG_update_from_arguments(cfg)
  call user_initialize(cfg)
  call domain_init(cfg)
  call refine_init(cfg, ndim)
  call time_step_init(cfg)
  call GL_initialize(cfg, ndim)
  call field_initialize(cfg, mg)
  call init_particle(cfg, pc)
  call check_path_writable(trim(GL_output_dir))

  fname = trim(GL_output_dir) // "/" // trim(GL_simulation_name) // "_out.cfg"
  call CFG_write(cfg, trim(fname))

  ! Initialize the tree (which contains all the mesh information)
  call init_tree(tree)

  ! Set the multigrid options. First define the variables to use
  mg%i_phi = i_phi
  mg%i_tmp = i_Ex
  mg%i_rhs = i_rhs

  ! This automatically handles cylindrical symmetry
  mg%box_op => mg_auto_op
  mg%box_gsrb => mg_auto_gsrb
  mg%box_corr => mg_auto_corr

  ! This routine always needs to be called when using multigrid
  call mg_init(tree, mg)

  call af_set_cc_methods(tree, i_electron, af_bc_neumann_zero)
  call af_set_cc_methods(tree, i_pos_ion, af_bc_neumann_zero)
  call af_set_cc_methods(tree, i_phi, mg%sides_bc, mg%sides_rb)

  output_cnt      = 0         ! Number of output files written
  GL_time         = 0         ! Simulation time (all times are in s)

  ! Set up the initial conditions
  if (.not. associated(user_initial_particles)) &
       error stop "user_initial_particles not defined"
  call user_initial_particles(pc)

  do
     call af_tree_clear_cc(tree, i_pos_ion)
     call particles_to_density_and_events(tree, pc, .true.)
     call field_compute(tree, mg, .false.)
     call af_adjust_refinement(tree, refine_routine, ref_info, &
          refine_buffer_width)
     if (ref_info%n_add == 0) exit
  end do

  call pc%set_accel()

  print *, "Number of threads", af_get_max_threads()
  call af_print_info(tree)

  ! Start from small time step
  GL_dt   = GL_dt_min

  ! Initial wall clock time
  call system_clock(t_start, count_rate)
  inv_count_rate = 1.0_dp / count_rate
  time_last_print = -1e10_dp
  n_prev_merge = pc%get_num_sim_part()
  n_elec_prev = pc%get_num_real_part()

  do it = 1, huge(1)-1
     if (GL_time >= GL_end_time) exit

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
     call pc%after_mover(dt)
     wtime_advance = wtime_advance + omp_get_wtime() - t0

     GL_time = GL_time + dt

     t0 = omp_get_wtime()
     call particles_to_density_and_events(tree, pc, .false.)
     wtime_density = wtime_density + omp_get_wtime() - t0

     n_part = pc%get_num_sim_part()
     if (n_part > n_prev_merge * min_merge_increase) then
        t0 = omp_get_wtime()
        call adapt_weights(tree, pc)
        wtime_merge = wtime_merge + omp_get_wtime() - t0
        n_prev_merge = pc%get_num_sim_part()
     end if

     ! Compute field with new density
     t0 = omp_get_wtime()
     call field_compute(tree, mg, .true.)
     wtime_field = wtime_field + omp_get_wtime() - t0

     call pc%set_accel()

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

        write(fname, "(A,I6.6)") trim(GL_simulation_name) // "_", output_cnt
        call af_write_silo(tree, fname, output_cnt, GL_time, &
             dir=GL_output_dir)
        call print_info()
        wtime_io = wtime_io + omp_get_wtime() - t0
     end if

     if (mod(it, refine_per_steps) == 0) then
        t0 = omp_get_wtime()
        call af_adjust_refinement(tree, refine_routine, ref_info, &
             refine_buffer_width)

        if (ref_info%n_add + ref_info%n_rm > 0) then
           ! Compute the field on the new mesh
           call particles_to_density_and_events(tree, pc, .false.)
           call field_compute(tree, mg, .true.)
           call adapt_weights(tree, pc)
        end if
     end if
     wtime_amr = wtime_amr + omp_get_wtime() - t0
  end do

contains

  !> Initialize the AMR tree
  subroutine init_tree(tree)
    type(af_t), intent(inout) :: tree
    integer                   :: coarse_grid(NDIM)

    coarse_grid = box_size

    ! Initialize tree
    if (GL_cylindrical) then
       call af_init(tree, box_size, domain_len, coarse_grid, coord=af_cyl)
    else
       call af_init(tree, box_size, domain_len, coarse_grid)
    end if

  end subroutine init_tree

  subroutine print_status()
    write(*, "(F7.2,A,I0,A,E10.3,A,E10.3,A,E10.3,A,E10.3,A,E10.3)") &
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
    write(*, "(A20,E12.4)") "mean energy", mean_en / UC_elec_volt
    write(*, "(A20,2E12.4)") "n_part, n_elec", n_part, n_elec
    write(*, "(A20,E12.4)") "mean weight", n_elec/n_part

    wtime_run = omp_get_wtime() - wtime_start
    write(*, "(A20,F10.4)") "advance", wtime_advance / wtime_run
    write(*, "(A20,F10.4)") "field", wtime_field / wtime_run
    write(*, "(A20,F10.4)") "merge", wtime_merge / wtime_run
    write(*, "(A20,F10.4)") "density", wtime_density / wtime_run
    write(*, "(A20,F10.4)") "output", wtime_io / wtime_run
    write(*, "(A20,F10.4)") "AMR", wtime_amr / wtime_run
    write(*, "(A20,F10.4)") "other", 1 - (wtime_advance + wtime_field + &
         wtime_merge + wtime_density + wtime_io + wtime_amr) / wtime_run
  end subroutine print_info

  subroutine set_output_variables()
    use m_units_constants
    integer :: n, n_part
    real(dp), allocatable :: coords(:, :)
    real(dp), allocatable :: weights(:)
    real(dp), allocatable :: energy(:)
    integer, allocatable  :: id_guess(:)
    ! integer, allocatable :: ionize_ix(:)

    n_part = pc%get_num_sim_part()
    allocate(coords(NDIM, n_part))
    allocate(weights(n_part))
    allocate(energy(n_part))
    allocate(id_guess(n_part))

    ! call pc%get_colls_of_type(CS_ionize_t, ionize_ix)

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
    call af_gc_tree(tree, i_electron)
    call af_gc_tree(tree, i_pos_ion)

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
