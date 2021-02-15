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

  implicit none

  integer, parameter     :: int8 = selected_int_kind(18)
  integer, parameter     :: ndim = NDIM
  integer(int8)          :: t_start, t_current, count_rate
  real(dp)               :: dt
  real(dp)               :: wc_time, inv_count_rate
  real(dp)               :: time_last_print, time_last_generate
  integer                :: it
  integer                :: n_part, n_prev_merge, n_samples
  integer, allocatable   :: ref_links(:, :)
  character(len=GL_slen) :: fname
  logical                :: write_out
  type(ref_info_t)       :: ref_info
  real(dp)               :: n_elec, n_elec_prev, max_elec_dens
  real(dp)               :: dt_cfl, dt_growth, dt_drt

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
  call refine_init(cfg, ndim)
  call time_step_init(cfg)
  call GL_initialize(cfg, ndim)
  call check_path_writable(trim(GL_output_dir))
  call field_initialize(cfg, mg)
  call init_particle(cfg, pc)
  call photons_initialize(cfg)

  ! Write configuration to output
  fname = trim(GL_output_dir) // "/" // trim(GL_simulation_name) // "_out.cfg"
  call CFG_write(cfg, trim(fname))

  ! Set the multigrid options. First define the variables to use
  mg%i_phi = i_phi
  mg%i_tmp = i_residual
  mg%i_rhs = i_rhs
  if (GL_use_dielectric) mg%i_eps = i_eps

  if (GL_use_electrode) then
     if (any(field_rod_r0 <= -1.0_dp)) &
          error stop "field_rod_r0 not set correctly"
     if (any(field_rod_r1 <= -1.0_dp)) &
          error stop "field_rod_r1 not set correctly"
     if (field_rod_radius <= 0) &
          error stop "field_rod_radius not set correctly"
     ! print *, "AT YOUR SERVICE"
     call af_set_cc_methods(tree, i_lsf, af_bc_neumann_zero, &
          af_gc_prolong_copy, af_prolong_zeroth, funcval=field_rod_lsf)
     mg%i_lsf = i_lsf
  end if

  ! This automatically handles cylindrical symmetry
  mg%box_op => mg_auto_op
  mg%box_gsrb => mg_auto_gsrb
  mg%box_corr => mg_auto_corr
  mg%box_stencil => mg_auto_stencil

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
     call dielectric_initialize(tree, i_eps, diel, n_surf_vars)
  end if

  output_cnt         = 0 ! Number of output files written
  GL_time            = 0 ! Simulation time (all times are in s)
  time_last_generate = GL_time

  ! Set up the initial conditions
  if (.not. associated(user_initial_particles)) &
       error stop "user_initial_particles not defined"
  call user_initial_particles(pc)

  ! Perform additional refinement
  do
     call af_tree_clear_cc(tree, i_pos_ion)
     call particles_to_density_and_events(tree, pc, .true.)

     if (associated(user_initial_ion_density)) &
          call af_loop_box(tree, user_initial_ion_density)

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

  ! Set an initial surface charge
  if (associated(user_set_surface_charge)) then
     call dielectric_set_values(tree, diel, i_surf_pos_ion, user_set_surface_charge)
  end if

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
     call pc%after_mover(dt)
     wtime_advance = wtime_advance + omp_get_wtime() - t0

     GL_time = GL_time + dt

     call particles_to_density_and_events(tree, pc, .false.)

     n_part = pc%get_num_sim_part()
     if (n_part > n_prev_merge * min_merge_increase) then
        call adapt_weights(tree, pc)
        n_prev_merge = pc%get_num_sim_part()
     end if

     ! Compute field with new density
     call field_compute(tree, mg, .true.)

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
        call set_output_variables()

        write(fname, "(A,I6.6)") trim(GL_simulation_name) // "_", output_cnt
        call af_write_silo(tree, fname, output_cnt, GL_time, &
             dir=GL_output_dir, add_curve_names = ["EEDF"], &
             add_curve_dat = write_EEDF_as_curve(pc))
        call print_info()
        call CS_write_ledger(pc%coll_ledger, &
        trim(GL_output_dir) // "/" // trim(GL_simulation_name) // "_cs_ledger.txt", &
        GL_time)
        write(fname, "(A,I6.6)") trim(GL_output_dir) // "_log.txt"
        
        call output_log(tree, trim(GL_output_dir) // "/" // trim(GL_simulation_name) // &
        &"_log.txt", output_cnt, wc_time)
        call write_EEDF_txt(pc, trim(GL_output_dir) // "/" // trim(GL_simulation_name) // &
        &"_EEDF.txt", output_cnt)
        
        if (GL_write_to_dat) then
          if (GL_write_to_dat_interval(1) .le. GL_time .and. GL_time .le. GL_write_to_dat_interval(2)) then
            call af_write_tree(tree, trim(GL_output_dir) // "/" // trim(fname), write_sim_data)
          end if
        end if

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
           call particles_to_density_and_events(tree, pc, .false.)
           call field_compute(tree, mg, .true.)
           call adapt_weights(tree, pc)
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
    call dielectric_get_integral(diel, i_surf_sum_dens, surf_int)
    write(*, "(A20,E12.4)") "Net charge", sum_pos_ion - sum_elec + surf_int
    write(*, "(A20,E12.4)") "mean energy", mean_en / UC_elec_volt
    write(*, "(A20,2E12.4)") "n_part, n_elec", n_part, n_elec
    write(*, "(A20,E12.4)") "mean weight", n_elec/n_part

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

  subroutine output_log(tree, filename, out_cnt, wc_time)
    !use m_field, only: field_voltage
    !use m_user_methods
    !use m_chemistry
    !use m_analysis
    use m_gas
    use m_domain
    type(af_t), intent(in)       :: tree
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: out_cnt !< Output number
    real(dp), intent(in)         :: wc_time !< Wallclock time
    character(len=50), save      :: fmt
    integer                      :: my_unit
    real(dp)                     :: velocity, dt
    real(dp), save               :: prev_pos(NDIM) = 0
    real(dp)                     :: sum_elec, sum_pos_ion
    real(dp)                     :: max_elec, max_field, max_Er, min_Er
    real(dp)                     :: sum_elem_charge
    !real(dp)                     :: elecdens_threshold, max_field_tip
    !real(dp)                     :: r0(NDIM), r1(NDIM)
    type(af_loc_t)               :: loc_elec, loc_field, loc_Er
    integer                      :: n_reals
    !character(len=name_len)      :: var_names(user_max_log_vars)
    !real(dp)                     :: var_values(user_max_log_vars)
    logical, save                :: first_time     = .true.
    !real(dp) :: density_threshold = 1e18_dp

    my_unit= 888
    call af_tree_sum_cc(tree, i_electron, sum_elec)
    call af_tree_sum_cc(tree, i_pos_ion, sum_pos_ion)
    call af_tree_max_cc(tree, i_electron, max_elec, loc_elec)
    call af_tree_max_cc(tree, i_E, max_field, loc_field)
    call af_tree_max_fc(tree, 1, ifc_E, max_Er, loc_Er)
    call af_tree_min_fc(tree, 1, ifc_E, min_Er)

    ! Scale threshold with gas number density
    !elecdens_threshold = density_threshold * &
    !     (GAS_number_dens/2.414e25_dp)**2
    !call analysis_zmin_zmax_threshold(tree, i_electron, elecdens_threshold, &
    !     [domain_len(NDIM), 0.0_dp], ne_zminmax)

    ! Assume streamer tip is located at the farthest plasma density (above a
    ! threshold) away from the electrode boundaries
    !r0(:) = 0.0_dp
    !r1 = r0 + domain_len

    !if (ne_zminmax(1)-  r0 < &
    !      r0  + domain_len(NDIM) - ne_zminmax(2)) then
    !   r0(NDIM) = ne_zminmax(2) - 0.02_dp * domain_len(NDIM)
    !   r1(NDIM) = ne_zminmax(2) + 0.02_dp * domain_len(NDIM)
    !else
    !   r0(NDIM) = ne_zminmax(1) - 0.02_dp * domain_len(NDIM)
    !   r1(NDIM) = ne_zminmax(1) + 0.02_dp * domain_len(NDIM)
    !end if

    !call analysis_max_var_region(tree, i_E, r0, r1, &
    !     max_field_tip, loc_tip)

    sum_elem_charge = 0
    call af_tree_sum_cc(tree, i_rhs, sum_elem_charge)

    dt = GL_dt

    if (first_time) then
       first_time = .false.

       open(newunit=my_unit, file=trim(filename), action="write")
#if NDIM == 1
       write(my_unit, "(A)", advance="no") "it time dt v sum(n_e) sum(n_i) &
            &sum(charge) max(E) x max(n_e) x wc_time n_cells min(dx) &
            &highest(lvl)"
#elif NDIM == 2
       write(my_unit, "(A)", advance="no") "it time dt v sum(n_e) sum(n_i) &
            &sum(charge) max(E) x y max(n_e) x y max(E_r) x y min(E_r) &
            &wc_time n_cells min(dx) highest(lvl)"
#elif NDIM == 3
       write(my_unit, "(A)", advance="no") "it time dt v sum(n_e) sum(n_i) &
            &sum(charge) max(E) x y z max(n_e) x y z &
            &wc_time n_cells min(dx) highest(lvl)"
#endif
       
       write(my_unit, *) ""
       close(my_unit)

       ! Start with velocity zero
       prev_pos = af_r_loc(tree, loc_field)
    end if

#if NDIM == 1
    n_reals = 10
#elif NDIM == 2
    n_reals = 17
#elif NDIM == 3
    n_reals = 16
#endif

    write(fmt, "(A,I0,A)") "(I6,", n_reals, "E16.8,I12,1E16.8)"
    
    velocity = norm2(af_r_loc(tree, loc_field) - prev_pos) / GL_dt_output
    prev_pos = af_r_loc(tree, loc_field)

    open(newunit=my_unit, file=trim(filename), action="write", &
         position="append")
#if NDIM == 1
    write(my_unit, fmt) out_cnt, GL_time, dt, velocity, sum_elec, &
         sum_pos_ion, sum_elem_charge, &
         max_field, af_r_loc(tree, loc_field), max_elec, &
         af_r_loc(tree, loc_elec), &
         wc_time, af_num_cells_used(tree), &
         af_min_dr(tree),tree%highest_lvl
         !var_values(1:n_user_vars)
#elif NDIM == 2
    write(my_unit, fmt) out_cnt, GL_time, dt, velocity, sum_elec, &
         sum_pos_ion, sum_elem_charge, &
         max_field, af_r_loc(tree, loc_field), max_elec, &
         af_r_loc(tree, loc_elec), max_Er, af_r_loc(tree, loc_Er), min_Er, &
         wc_time, af_num_cells_used(tree), af_min_dr(tree),tree%highest_lvl
         !var_values(1:n_user_vars)
#elif NDIM == 3
    write(my_unit, fmt) out_cnt, GL_time, dt, velocity, sum_elec, &
         sum_pos_ion, sum_elem_charge, &
         max_field, af_r_loc(tree, loc_field), max_elec, &
         af_r_loc(tree, loc_elec),&
         wc_time, af_num_cells_used(tree), &
         af_min_dr(tree),tree%highest_lvl
         !var_values(1:n_user_vars)
#endif
    close(my_unit)

  end subroutine output_log

  subroutine write_EEDF_txt(pc, filename, out_cnt)
    type(PC_t), intent(in)  :: pc
    integer                      :: my_unit  
    !character(len=50), save      :: fmt
    character(len=*), intent(in) :: filename
    integer                      :: n_reals, num_bins
    integer, intent(in)          :: out_cnt !< Output number
    real(dp)                     :: n_part, max_en, en_step = 0.75
    logical, save                :: first_time     = .true.
    real(dp), allocatable        :: EEDF_dat(:, :, :)
    
    my_unit= 888
    
    n_part = pc%get_num_real_part()
    max_en = get_max_energy(pc)
    num_bins = ceiling(max_en / en_step) ! Generate bins of en_step (eV) each up until max energy
    if (num_bins < 1) num_bins = 1 ! At least one bin


    allocate(EEDF_dat(1, 2, num_bins))
    
    
    
    EEDF_dat = write_EEDF_as_curve(pc)
    if (first_time) then
       first_time = .false.

       open(newunit=my_unit, file=trim(filename), action="write")
        
       write(my_unit, "(A)", advance="no") "it time electron_energy distribution "
       
       write(my_unit, *) ""
       close(my_unit)
    end if
    
    n_reals = 1
    
    open(newunit=my_unit, file=trim(filename), action="write", &
         position="append")
         
    write(my_unit, *) out_cnt, GL_time,EEDF_dat(1, 1, :), EEDF_dat(1, 2, :)
    close(my_unit)
  
  end subroutine write_EEDF_txt

    
end program apic
