#include "../afivo/src/cpp_macros.h"
!> Program to perform discharge simulations with particle-in-cell method
program test_photoemission

  ! use omp_lib
  ! use m_af_all
  ! use m_globals
  ! use m_domain
  ! use m_refine
  ! use m_user
  ! use m_user_methods

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

  ! Write configuration to output
  fname = trim(GL_output_dir) // "/" // trim(GL_simulation_name) // "_out.cfg"
  call CFG_write(cfg, trim(fname))

! integer, parameter     :: ndim = NDIM
!
!   ! Read command line arguments and configuration files
!   call CFG_update_from_arguments(cfg)
!
!   ! Initialize the user's code
!   call user_initialize(cfg)
!
!   ! Initialize other modules
!   call domain_init(cfg)
!   call refine_init(cfg, ndim)
!   call GL_initialize(cfg, ndim)

  ! Initialize the tree (which contains all the mesh information)
  call init_tree(tree)

  ! Set the dielectric permittivity before initializing multigrid
  if (GL_use_dielectric) then
     if (.not. associated(user_set_dielectric_eps)) &
          error stop "user_set_dielectric_eps not defined"
     call af_loop_box(tree, user_set_dielectric_eps)
  end if

  if (GL_use_dielectric) then
     ! Initialize dielectric surfaces at the third refinement level
     call af_refine_up_to_lvl(tree, 3)
     call dielectric_initialize(tree, i_eps, diel, 1)
  end if

call random_photoemission_event()

contains

  subroutine random_photoemission_event()
    use m_random
    use m_dielectric, only: bisect_line
    real(dp) :: start(2), final(2)
    real(dp) :: x_start(2), x_stop(2)
    logical  :: on_surface

    ! Generate two random numbers between [0.25, 0.75]
    start(1) = GL_rng%unif_01() * 0.5 + 0.25
    start(2) = GL_rng%unif_01() * 0.5 + 0.25
    x_start = start * domain_len

    ! Generate two random numbers between [-2, 2]
    final(1) = GL_rng%unif_01() * 4 - 2
    final(2) = GL_rng%unif_01() * 4 - 2
    x_stop = final * domain_len

    print *, "=========================="
    print *, " "
    print *, "start coordinates: ", start
    print *, "stop coordinates: ", final
    print *, " "

    call bisect_line(tree, x_start, x_stop, on_surface, i_eps)

    print *, "Calculated x_start: ", x_start / domain_len
    print *, "Calculated x_stop: ", x_stop / domain_len
    print *, "on_surface = ", on_surface
    print *, " "
    print *, "=========================="

  end subroutine



end program
