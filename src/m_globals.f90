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

  type(CFG_t)            :: cfg  ! The configuration for the simulation
  type(af_t)            :: tree ! This contains the full grid information
  type(mg_t)           :: mg   ! Multigrid option struct
  type(PC_t)             :: pc

  ! Default length of strings
  integer, parameter :: ST_slen = 200

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
  integer, parameter :: name_len = 12

  ! Names of the cell-centered variables
  character(len=name_len), allocatable :: ST_cc_names(:)

  ! Indices of variables to be included in output
  integer, allocatable :: vars_for_output(:)

  ! ** Indices of face-centered variables **
  integer, protected :: n_var_face   = 0 ! Number of variables

  ! Whether cylindrical coordinates are used
  logical :: ST_cylindrical = .false.

  ! Random number generator
  type(rng_t) :: ST_rng

  ! Parallel random number generator
  type(prng_t) :: ST_prng

  ! Name of the simulations
  character(len=ST_slen), protected :: ST_simulation_name = "sim"

  ! Output directory
  character(len=ST_slen), protected :: ST_output_dir = "output"

  ! Print status every this many seconds
  real(dp), protected :: ST_print_status_sec = 60.0_dp

  ! Current time
  real(dp)  :: ST_time

  ! End time of the simulation
  real(dp), protected :: ST_end_time = 10e-9_dp

  ! Pressure of the gas in bar
  real(dp), protected :: ST_gas_pressure = 1.0_dp

  ! Name of the gas mixture
  character(len=ST_slen) :: ST_gas_name = "AIR"

  ! Fraction of O2
  real(dp), protected :: ST_gas_frac_O2 = 0.2_dp

  ! Position of initial electron seed
  real(dp) :: ST_init_seed_pos(3)   = [0.5_dp, 0.5_dp, 0.9_dp]
  real(dp) :: ST_init_seed_sigma    = 1.0e-3_dp
  integer  :: ST_init_num_particles = 10000

  real(dp) :: particle_min_weight = 1.0_dp
  real(dp) :: particle_max_weight = 1.0e20_dp
  real(dp) :: particle_per_cell = 100.0_dp

contains

  integer function ST_add_cc_variable(name, include_in_output)
    character(len=*), intent(in) :: name
    logical, intent(in)          :: include_in_output
    integer                      :: i, n

    ST_cc_names = [character(len=name_len) :: &
         (ST_cc_names(i), i=1,n_var_cell), name]

    if (include_in_output) then
       if (allocated(vars_for_output)) then
          n = size(vars_for_output)
       else
          n = 0
       end if
       vars_for_output = [(vars_for_output(i), i=1,n), n_var_cell+1]
    end if

    ST_add_cc_variable = n_var_cell + 1
    n_var_cell         = n_var_cell + 1
  end function ST_add_cc_variable

  integer function ST_cc_var_index(name)
    character(len=*), intent(in) :: name
    integer :: n

    ST_cc_var_index = -1
    do n = 1, size(ST_cc_names)
       if (ST_cc_names(n) == name) then
          ST_cc_var_index = n
          exit
       end if
    end do

  end function ST_cc_var_index

  integer function ST_add_fc_variable()
    ST_add_fc_variable = n_var_face + 1
    n_var_face         = n_var_face + 1
  end function ST_add_fc_variable

  !> Create the configuration file with default values
  subroutine ST_initialize(cfg, ndim)
    use iso_fortran_env, only: int64
    use m_config
    use omp_lib
    use m_afivo_types
    type(CFG_t), intent(inout) :: cfg  !< The configuration for the simulation
    integer, intent(in)        :: ndim !< Number of dimensions
    integer                    :: n_threads
    integer                    :: rng_int4_seed(4) = &
         [8123, 91234, 12399, 293434]
    integer(int64)             :: rng_int8_seed(2)

    i_electron = ST_add_cc_variable("electron", .true.)
    i_pos_ion = ST_add_cc_variable("pos_ion", .true.)
    i_phi = ST_add_cc_variable("phi", .true.)
    i_Ex = ST_add_cc_variable("Ex", .true.)
    i_Ey = ST_add_cc_variable("Ey", .true.)
    if (ndim > 2) then
       i_Ez = ST_add_cc_variable("Ez", .true.)
    end if
    i_E = ST_add_cc_variable("E", .true.)
    i_rhs = ST_add_cc_variable("rhs", .true.)
    i_ppc = ST_add_cc_variable("ppc", .true.)
    i_energy = ST_add_cc_variable("energy", .true.)

    call CFG_add_get(cfg, "cylindrical", ST_cylindrical, &
         "Whether cylindrical coordinates are used (only in 2D)")
    call CFG_add_get(cfg, "end_time", ST_end_time, &
         "The desired endtime (s) of the simulation")
    call CFG_add_get(cfg, "simulation_name", ST_simulation_name, &
         "The name of the simulation")
    call CFG_add_get(cfg, "output_dir", ST_output_dir, &
         "Directory where the output should be written")
    call CFG_add_get(cfg, "print_status_sec", ST_print_status_sec, &
         "Print status every this many seconds")
    call CFG_add_get(cfg, "gas%pressure", ST_gas_pressure, &
         "The gas pressure (bar), used for photoionization")
    call CFG_add_get(cfg, "gas%name", ST_gas_name, &
         "The name of the gas mixture used")
    call CFG_add_get(cfg, "gas%frac_O2", ST_gas_frac_O2, &
         "Fraction of O2, used for photoionization")

    call CFG_add_get(cfg, "rng_seed", rng_int4_seed, &
         "Seed for random numbers. If all zero, generate from clock.")

    if (all(rng_int4_seed == 0)) then
       rng_int4_seed = get_random_seed()
       print *, "RNG seed: ", rng_int4_seed
    end if

    rng_int8_seed = transfer(rng_int4_seed, rng_int8_seed)
    call ST_rng%set_seed(rng_int8_seed)
    n_threads = af_get_max_threads()
    call ST_prng%init_parallel(n_threads, ST_rng)

  end subroutine ST_initialize

  !> Get a random seed based on the current time
  function get_random_seed() result(seed)
    integer :: seed(4)
    integer :: time, i

    call system_clock(time)
    do i = 1, 4
       seed(i) = ishftc(time, i*8)
    end do
  end function get_random_seed

end module m_globals
