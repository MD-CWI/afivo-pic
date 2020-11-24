!> This module contains several pre-defined variables like:
!! * Indices of cell-centered variables
!! * Names of the cell-centered variables
!! * Indices of face-centered variables
!! * Indices of transport data
module m_globals
  use m_af_all
  use m_particle_core
  use m_units_constants
  use m_config
  use m_random

  implicit none
  public

  type(CFG_t) :: cfg  ! The configuration for the simulation
  type(af_t)  :: tree ! This contains the full grid information
  type(mg_t)  :: mg   ! Multigrid option struct
  type(PC_t)  :: pc
  type(dielectric_t) :: diel ! To store dielectric surface

  ! Default length of strings
  integer, parameter :: GL_slen = 200

  ! Named integer for outside checks
  integer, parameter :: outside_domain = 1
  integer, parameter :: inside_dielectric = 2

  ! Interpolation order for the electric field
  integer, protected :: interpolation_order_field = 1

  ! Interpolation order for mapping particles to density
  integer, protected :: interpolation_order_to_density = 1

  ! ** Indices of cell-centered variables **
  integer, protected :: n_var_cell = 0  ! Number of variables
  integer, protected :: i_electron = -1 ! Electron density
  integer, protected :: i_pos_ion  = -1 ! Positive ion density
  integer, protected :: i_phi      = -1 ! Electrical potential
  integer, protected :: i_Ex       = -1 ! Electric field (x)
  integer, protected :: i_Ey       = -1 ! Electric field (y)
  integer, protected :: i_Ez       = -1 ! Electric field (z)
  integer, protected :: i_E_all(NDIM) = -1 ! Electric field components
  integer, protected :: ifc_E = -1 ! Face-centered electric field
  integer, protected :: i_E        = -1 ! norm(E) (cell-centered)
  integer, protected :: i_E_v2 = -1     ! norm(E) (face-centered)
  integer, protected :: i_rhs      = -1 ! Source term Poisson
  integer, protected :: i_ppc      = -1 ! Particles per cell
  integer, protected :: i_energy   = -1 ! Energy density
  integer, protected :: i_eps      = -1 ! Dielectric permittivity
  integer, protected :: i_P_dep    = -1 ! Density of deposited power
  integer, protected, allocatable :: i_tracked_cIx(:) ! CAS that will be tracked
  integer, protected :: i_lsf      = -1 ! level set function (for electrode)

  integer, parameter :: name_len   = 12

  ! Index of surface charge on dielectric
  integer, parameter :: n_surf_vars = 4
  integer, parameter :: i_surf_sum_dens = 1
  integer, parameter :: i_surf_elec_close = 2
  integer, parameter :: i_surf_elec = 3
  integer, parameter :: i_surf_pos_ion = 4

  ! Whether cylindrical coordinates are used
  logical, protected :: GL_cylindrical = .false.

  ! Whether a dielectric is used
  logical, protected :: GL_use_dielectric = .false.

  ! Whether a dielectric is used
  logical, protected :: GL_use_electrode = .false.

  ! Stop multigrid when residual is smaller than this factor times max(|rhs|)
  real(dp), public, protected :: GL_multigrid_max_rel_residual = 1e-4_dp

  !> Boundary condition for the plasma species
  procedure(af_subr_bc), public, protected, pointer :: &
       bc_species => null()

  ! Whether a the output is also written to a .dat file
  logical, protected :: GL_write_to_dat = .false.

  ! Wheter additional active species will be tracked
  logical, protected :: GL_track_CAS = .false.

  ! List of collision indices corresponding to species to track
  integer, allocatable, protected :: GL_cIx_to_track(:)

  character(len=GL_slen), allocatable, protected :: GL_cIx_labels(:)

  ! Total number of species to track
  integer :: num_cIx_to_track

  ! Interval for writing to .dat file
  real(dp), protected :: GL_write_to_dat_interval(2)


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

  integer  :: ii
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

#if NDIM == 2
    i_E_all = [i_Ex, i_Ey]
#elif NDIM == 3
    call af_add_cc_variable(tree, "Ez", .true., ix=i_Ez)
    i_E_all = [i_Ex, i_Ey, i_Ez]
#endif

    call af_add_fc_variable(tree, "E_fc", ix=ifc_E)

    call af_add_cc_variable(tree, "E", .false., ix=i_E)
    call af_add_cc_variable(tree, "E_v2", .true., ix=i_E_v2)
    call af_add_cc_variable(tree, "rhs", .true., ix=i_rhs)
    call af_add_cc_variable(tree, "ppc", .false., ix=i_ppc)
    call af_add_cc_variable(tree, "energy", .false., ix=i_energy)
    call af_add_cc_variable(tree, "power_deposition", .true., ix=i_P_dep)


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
    call CFG_add_get(cfg, "track_CAS", GL_track_CAS, &
          "Whether additional grid-data is used to track specific collision types")
    call CFG_add_get(cfg, "write_to_dat", GL_write_to_dat, &
          "Wheter the output is also written to a .dat file")
    call CFG_add_get(cfg, "write_to_dat_interval", GL_write_to_dat_interval, &
          "The time interval that is written to a .dat file")
    call CFG_add_get(cfg, "use_electrode", GL_use_electrode, &
             "Whether to include an electrode")


    if (GL_use_dielectric) then
       interpolation_order_field = 1
       interpolation_order_to_density = 1

       call af_add_cc_variable(tree, "eps", ix=i_eps)
       call af_set_cc_methods(tree, i_eps, af_bc_neumann_zero, &
            af_gc_prolong_copy, af_prolong_zeroth)
    end if

    if (GL_track_CAS) then
      call CFG_add(cfg, "cIx_to_track", [-1], &
            "List of collision indices that are tracked", .true.)

      call CFG_get_size(cfg, "cIx_to_track", num_cIx_to_track)

      allocate(GL_cIx_to_track(num_cIx_to_track))
      allocate(i_tracked_cIx(num_cIx_to_track))
      allocate(GL_cIx_labels(num_cIx_to_track))

      i_tracked_cIx = -1

      call CFG_get(cfg, "cIx_to_track", GL_cIx_to_track)
      call CFG_add_get(cfg, "cIx_labels", GL_cIx_labels, &
           "Output labels for tracked species")

      do ii = 1, num_cIx_to_track
        call af_add_cc_variable(tree, trim(GL_cIx_labels(ii)), .true., ix=i_tracked_cIx(ii))
        call af_set_cc_methods(tree, i_tracked_cIx(ii), af_bc_neumann_zero, &
            prolong=af_prolong_limit)
      end do
    end if

    if (GL_use_electrode) then
      call af_add_cc_variable(tree, "lsf", .true., ix=i_lsf)
      bc_species => bc_species_dirichlet_zero ! TODO why this instead of neumann zero?
    end if

    call GL_rng%set_random_seed()

  end subroutine GL_initialize

    !> Impose a Dirichlet zero boundary condition for plasma species in the last
    !> dimension, which is supposed to have the electrodes. We use Neumann
    !> conditions in the other dimensions. Note that this version avoids
    !> extrapolation (in contrast to the regular Dirichlet b.c.), which is more
    !> suitable for conserved species densities.
    subroutine bc_species_dirichlet_zero(box, nb, iv, coords, bc_val, bc_type)
      type(box_t), intent(in) :: box
      integer, intent(in)     :: nb
      integer, intent(in)     :: iv
      real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
      real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
      integer, intent(out)    :: bc_type

      if (af_neighb_dim(nb) == NDIM) then
         bc_type = af_bc_dirichlet_copy
         bc_val  = 0.0_dp
      else
         bc_type = af_bc_neumann
         bc_val  = 0.0_dp
      end if
    end subroutine bc_species_dirichlet_zero

end module m_globals
