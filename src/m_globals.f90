!> This module contains several pre-defined variables like:
!! * Indices of cell-centered variables
!! * Names of the cell-centered variables
!! * Indices of face-centered variables
!! * Indices of transport data
module m_globals
  use m_af_all
  use m_particle_core
  use m_config
  use m_random

  implicit none
  public

  type(CFG_t) :: cfg  ! The configuration for the simulation
  type(af_t)  :: tree ! This contains the full grid information
  type(mg_t)  :: mg   ! Multigrid option struct
  type(PC_t)  :: pc

  ! Default length of strings
  integer, parameter :: GL_slen = 200

  ! Named integer for outside checks
  integer, parameter :: outside_domain = 1
  integer, parameter :: inside_dielectric = 2

  ! Which interpolation order is used
  integer, protected :: interpolation_order = 1

  ! ** Indices of cell-centered variables **
  integer, protected :: n_var_cell = 0  ! Number of variables
  integer, protected :: i_electron = -1 ! Electron density
  integer, protected :: i_pos_ion  = -1 ! Positive ion density
  integer, protected :: i_phi      = -1 ! Electrical potential
  integer, protected :: i_Ex       = -1 ! Electric field (x)
  integer, protected :: i_Ey       = -1 ! Electric field (y)
  integer, protected :: i_Ez       = -1 ! Electric field (z)
  integer, protected :: i_E        = -1 ! Electric field (x)
  integer, protected :: i_rhs      = -1 ! Source term Poisson
  integer, protected :: i_ppc      = -1 ! Particles per cell
  integer, protected :: i_energy   = -1 ! Energy density
  integer, protected :: i_eps      = -1 ! Dielectric permittivity
  integer, parameter :: name_len   = 12

  ! Whether cylindrical coordinates are used
  logical, protected :: GL_cylindrical = .false.

  ! Whether a dielectric is used
  logical, protected :: GL_use_dielectric = .false.

  ! Random number generator
  type(rng_t) :: GL_rng

  ! Name of the simulations
  character(len=GL_slen), protected :: GL_simulation_name = "sim"

  ! Output directory
  character(len=GL_slen), protected :: GL_output_dir = "output"

  ! Time between writing output
  real(dp), protected :: GL_dt_output = 1.0e-10_dp

  ! Print status every this many seconds
  real(dp), protected :: GL_print_status_sec = 60.0_dp

  ! Current time
  real(dp)  :: GL_time

  ! End time of the simulation
  real(dp), protected :: GL_end_time = 10e-9_dp

  ! Pressure of the gas in bar
  real(dp), protected :: GL_gas_pressure = 1.0_dp

  ! Name of the gas mixture
  character(len=GL_slen) :: GL_gas_name = "AIR"

  ! Fraction of O2
  real(dp), protected :: GL_gas_frac_O2 = 0.2_dp

  ! Position of initial electron seed
  real(dp) :: GL_init_seed_pos(3)   = [0.5_dp, 0.5_dp, 0.9_dp]
  real(dp) :: GL_init_seed_sigma    = 1.0e-3_dp
  integer  :: GL_init_num_particles = 10000

  real(dp) :: particle_min_weight = 1.0_dp
  real(dp) :: particle_max_weight = 1.0e20_dp
  real(dp) :: particle_per_cell = 100.0_dp

contains

  !> Create the configuration file with default values
  subroutine GL_initialize(cfg, ndim)
    use iso_fortran_env, only: int64
    use m_config
    use omp_lib
    type(CFG_t), intent(inout) :: cfg  !< The configuration for the simulation
    integer, intent(in)        :: ndim !< Number of dimensions

    call af_add_cc_variable(tree, "electron", .true., ix=i_electron)
    call af_add_cc_variable(tree, "pos_ion", .true., ix=i_pos_ion)
    call af_add_cc_variable(tree, "phi", .true., ix=i_phi)
    call af_add_cc_variable(tree, "Ex", .true., ix=i_Ex)
    call af_add_cc_variable(tree, "Ey", .true., ix=i_Ey)
    if (ndim > 2) then
       call af_add_cc_variable(tree, "Ez", .true., ix=i_Ez)
    end if
    call af_add_cc_variable(tree, "E", .true., ix=i_E)
    call af_add_cc_variable(tree, "rhs", .true., ix=i_rhs)
    call af_add_cc_variable(tree, "ppc", .true., ix=i_ppc)
    call af_add_cc_variable(tree, "energy", .true., ix=i_energy)

    call CFG_add_get(cfg, "cylindrical", GL_cylindrical, &
         "Whether cylindrical coordinates are used (only in 2D)")
    call CFG_add_get(cfg, "end_time", GL_end_time, &
         "The desired endtime (s) of the simulation")
    call CFG_add_get(cfg, "simulation_name", GL_simulation_name, &
         "The name of the simulation")
    call CFG_add_get(cfg, "output_dir", GL_output_dir, &
         "Directory where the output should be written")
    call CFG_add_get(cfg, "dt_output", GL_dt_output, &
         "The timestep for writing output (s)")
    call CFG_add_get(cfg, "print_status_sec", GL_print_status_sec, &
         "Print status every this many seconds")
    call CFG_add_get(cfg, "gas%pressure", GL_gas_pressure, &
         "The gas pressure (bar), used for photoionization")
    call CFG_add_get(cfg, "gas%name", GL_gas_name, &
         "The name of the gas mixture used")
    call CFG_add_get(cfg, "gas%frac_O2", GL_gas_frac_O2, &
         "Fraction of O2, used for photoionization")
    call CFG_add_get(cfg, "use_dielectric", GL_use_dielectric, &
         "Whether a dielectric is used")

    if (GL_use_dielectric) then
       interpolation_order = 0
       call af_add_cc_variable(tree, "eps", ix=i_eps)
       call af_set_cc_methods(tree, i_eps, af_bc_neumann_zero, &
            af_gc_prolong_copy, af_prolong_zeroth)
    end if

    call GL_rng%set_random_seed()

  end subroutine GL_initialize

end module m_globals
