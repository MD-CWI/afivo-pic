!> This module contains several pre-defined variables like:
!! * Indices of cell-centered variables
!! * Names of the cell-centered variables
!! * Indices of face-centered variables
!! * Indices of transport data
module m_streamer

  use m_config
  use m_random
  use m_lookup_table

  implicit none
  public

  integer, parameter, private :: dp = kind(0.0d0)

  ! Default length of strings
  integer, parameter :: ST_slen = 200

  ! ** Indices of cell-centered variables **
  integer, protected :: n_var_cell = 0  ! Number of variables
  integer, protected :: i_electron = -1 ! Electron density
  integer, protected :: i_pos_ion  = -1 ! Positive ion density
  integer, protected :: i_phi      = -1 ! Electrical potential
  integer, protected :: i_Ex       = -1 ! Electric field (x)
  integer, protected :: i_Ey       = -1 ! Electric field (x)
  integer, protected :: i_rhs      = -1 ! Source term Poisson

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

  ! The refinement buffer width in cells (around flagged cells)
  integer, protected :: ST_refine_buffer_width = 2

  ! The number of steps after which the mesh is updated
  integer, protected :: ST_refine_per_steps = 2

  ! The grid spacing will always be larger than this value
  real(dp), protected :: ST_refine_min_dx = 1.0e-6_dp

  ! The grid spacing will always be smaller than this value
  real(dp), protected :: ST_refine_max_dx = 1.0e-3_dp

  ! Refine if alpha*dx is larger than this value
  real(dp), protected :: ST_refine_adx = 1.0_dp

  ! For refinement, use alpha(f * E)/f, where f is this factor
  real(dp), protected :: ST_refine_adx_fac = 1.0_dp

  ! Only derefine if grid spacing if smaller than this value
  real(dp), protected :: ST_derefine_dx = 1e-4_dp

  ! Refine around initial conditions up to this time
  real(dp), protected :: ST_refine_init_time = 10e-9_dp

  ! Refine until dx is smaller than this factor times the seed width
  real(dp), protected :: ST_refine_init_fac = 0.25_dp

  ! Refine a region up to this grid spacing
  real(dp), protected, allocatable :: ST_refine_regions_dr(:)

  ! Refine regions up to this simulation time
  real(dp), protected, allocatable :: ST_refine_regions_tstop(:)

  ! Minimum coordinate of the refinement regions
  real(dp), protected, allocatable :: ST_refine_regions_rmin(:,:)

  ! Maximum coordinate of the refinement regions
  real(dp), protected, allocatable :: ST_refine_regions_rmax(:,:)

  ! Current time step
  real(dp) :: ST_dt

  ! Maximum allowed time step
  real(dp), protected :: ST_dt_max = 1.0e-11_dp

  ! Minimum allowed time step
  real(dp), protected :: ST_dt_min = 1.0e-12_dp

  ! Time between writing output
  real(dp), protected :: ST_dt_output = 1.0e-10_dp

  ! Current time
  real(dp)  :: ST_time

  ! End time of the simulation
  real(dp), protected :: ST_end_time = 10e-9_dp

  ! The size of the boxes that we use to construct our mesh
  integer, protected :: ST_box_size = 8

  ! The length of the (square) domain
  real(dp), protected :: ST_domain_len = 4e-3_dp

  ! Pressure of the gas in bar
  real(dp), protected :: ST_gas_pressure = 1.0_dp

  ! Fraction of O2
  real(dp), protected :: ST_gas_frac_O2 = 0.2_dp

  ! Number of V-cycles to perform per time step
  integer, protected :: ST_multigrid_num_vcycles = 2

  real(dp) :: ST_init_seed_pos(3)   = [0.5_dp, 0.5_dp, 0.5_dp]
  real(dp) :: ST_init_seed_sigma    = 1.0e-5_dp
  integer  :: ST_init_num_particles = 10000

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
    integer                    :: n, n_threads
    real(dp)                   :: vec(ndim)
    real(dp), allocatable      :: dbuffer(:)
    integer                    :: rng_int4_seed(4) = &
         [8123, 91234, 12399, 293434]
    integer(int64)             :: rng_int8_seed(2)

    i_electron = ST_add_cc_variable("electron", .true.)
    i_pos_ion = ST_add_cc_variable("pos_ion", .true.)
    i_phi = ST_add_cc_variable("phi", .true.)
    i_Ex = ST_add_cc_variable("Ex", .true.)
    i_Ey = ST_add_cc_variable("Ey", .true.)
    i_rhs = ST_add_cc_variable("rhs", .true.)

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
    call CFG_add_get(cfg, "box_size", ST_box_size, &
         "The number of grid cells per coordinate in a box")
    call CFG_add_get(cfg, "domain_len", ST_domain_len, &
         "The length of the domain (m)")
    call CFG_add_get(cfg, "gas_pressure", ST_gas_pressure, &
         "The gas pressure (bar), used for photoionization")
    call CFG_add_get(cfg, "gas_frac_O2", ST_gas_frac_O2, &
         "Fraction of O2, used for photoionization")

    call CFG_add_get(cfg, "multigrid_num_vcycles", ST_multigrid_num_vcycles, &
         "Number of V-cycles to perform per time step")

    call CFG_add_get(cfg, "dt_output", ST_dt_output, &
         "The timestep for writing output (s)")
    call CFG_add_get(cfg, "dt_max", ST_dt_max, &
         "The maximum timestep (s)")
    call CFG_add_get(cfg, "dt_min", ST_dt_min, &
         "The minimum timestep (s)")

    call CFG_add_get(cfg, "refine_buffer_width", ST_refine_buffer_width, &
         "The refinement buffer width in cells (around flagged cells)")
    call CFG_add_get(cfg, "refine_per_steps", ST_refine_per_steps, &
         "The number of steps after which the mesh is updated")
    call CFG_add_get(cfg, "refine_min_dx", ST_refine_min_dx, &
         "The grid spacing will always be larger than this value")
    call CFG_add_get(cfg, "refine_max_dx", ST_refine_max_dx, &
         "The grid spacing will always be smaller than this value")

    if (ST_refine_min_dx > ST_refine_max_dx) &
         error stop "Cannot have refine_min_dx < refine_max_dx"

    call CFG_add_get(cfg, "refine_adx", ST_refine_adx, &
         "Refine if alpha*dx is larger than this value")
    call CFG_add_get(cfg, "refine_adx_fac", ST_refine_adx_fac, &
         "For refinement, use alpha(f * E)/f, where f is this factor")
    call CFG_add_get(cfg, "derefine_dx", ST_derefine_dx, &
         "Only derefine if grid spacing if smaller than this value")
    call CFG_add_get(cfg, "refine_init_time", ST_refine_init_time, &
         "Refine around initial conditions up to this time")
    call CFG_add_get(cfg, "refine_init_fac", ST_refine_init_fac, &
         "Refine until dx is smaller than this factor times the seed width")

    call CFG_add(cfg, "refine_regions_dr", [1.0e99_dp], &
         "Refine regions up to this grid spacing", .true.)
    call CFG_add(cfg, "refine_regions_tstop", [-1.0e99_dp], &
         "Refine regions up to this simulation time", .true.)
    vec = 0.0_dp
    call CFG_add(cfg, "refine_regions_rmin", vec, &
         "Minimum coordinate of the refinement regions", .true.)
    call CFG_add(cfg, "refine_regions_rmax", vec, &
         "Maximum coordinate of the refinement regions", .true.)

    call CFG_get_size(cfg, "refine_regions_dr", n)
    allocate(ST_refine_regions_dr(n))
    allocate(ST_refine_regions_tstop(n))
    allocate(ST_refine_regions_rmin(ndim, n))
    allocate(ST_refine_regions_rmax(ndim, n))
    allocate(dbuffer(ndim * n))

    call CFG_get(cfg, "refine_regions_dr", ST_refine_regions_dr)
    call CFG_get(cfg, "refine_regions_tstop", ST_refine_regions_tstop)
    call CFG_get(cfg, "refine_regions_rmin", dbuffer)
    ST_refine_regions_rmin = reshape(dbuffer, [ndim, n])
    call CFG_get(cfg, "refine_regions_rmax", dbuffer)
    ST_refine_regions_rmax = reshape(dbuffer, [ndim, n])

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

end module m_streamer
