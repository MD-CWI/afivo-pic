!> This module contains several pre-defined variables like:
!! * Indices of cell-centered variables
!! * Names of the cell-centered variables
!! * Indices of face-centered variables
!! * Indices of transport data
#include "../afivo/src/cpp_macros.h"
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
  type(surfaces_t) :: diel ! To store dielectric surface

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
  integer, protected :: ifc_E = -1 ! Face-centered electric field
  integer, protected :: i_E = -1     ! norm(E) (face-centered)
  integer, protected :: i_rhs      = -1 ! Source term Poisson
  integer, protected :: i_residual = -1 ! Multigrid residual
  integer, protected :: i_ppc      = -1 ! Particles per cell
  integer, protected :: i_energy   = -1 ! Energy density
  integer, protected :: i_eps      = -1 ! Dielectric permittivity
  integer, protected :: i_lsf      = -1 ! level set function (for electrode)
  integer, protected :: i_tmp_dens = -1 ! used in cylindrical coordinate system
  integer, parameter :: name_len   = 12

  ! Index of surface charge on dielectric
  integer, parameter :: n_surf_vars = 4
  integer, parameter :: i_surf_sum_dens = 1
  integer, parameter :: i_surf_elec_close = 2
  integer, parameter :: i_surf_elec = 3
  integer, parameter :: i_surf_pos_ion = 4

  !> How many gigabytes of memory to use for afivo
  real(dp), protected :: GL_memory_afivo_GB = 2.0_dp

  !> How many gigabytes of memory to use for particles
  real(dp), protected :: GL_memory_particles_GB = 4.0_dp

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

  ! Whether the mesh is written to a binary file
  logical, protected :: GL_write_mesh_binary = .false.

  ! Whether the particles are written to a binary file
  logical, protected :: GL_write_particles_binary = .false.

  ! Per how many regular outputs binary output is written
  integer, protected :: GL_binary_per_outputs = 1

  ! Random number generator
  type(rng_t) :: GL_rng

  ! Name of the simulations
  character(len=GL_slen), protected :: GL_simulation_name = "sim"

  ! Output directory
  character(len=GL_slen), protected :: GL_output_dir = "output"

  ! Whether to write Silo output
  logical, protected :: GL_write_silo = .true.

  ! Whether to write VTK unstructured output
  logical, protected :: GL_write_vtk = .false.

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

  ! Background ionization rate (m^-3 s^-1)
  real(dp), protected :: GL_background_ionization_rate = 0.0_dp

  ! Position of initial electron seed
  real(dp) :: GL_init_seed_pos(3)   = [0.5_dp, 0.5_dp, 0.9_dp]
  real(dp) :: GL_init_seed_sigma    = 1.0e-3_dp
  integer  :: GL_init_num_particles = 10000

  ! The length of the (square) domain
  ! TODO: cannot directly use domain_len in m_domain
  real(dp), protected, public :: GL_domain_len(NDIM) = 4e-3_dp

  real(dp) :: particle_min_weight = 1.0_dp
  real(dp) :: particle_max_weight = 1.0e20_dp
  real(dp) :: particle_per_cell = 100.0_dp

  ! Density threshold for detecting plasma regions
  real(dp) :: density_threshold = 1e18_dp
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
    call af_add_cc_variable(tree, "phi", .false., ix=i_phi)
    call af_add_fc_variable(tree, "E_fc", ix=ifc_E)

    call af_add_cc_variable(tree, "E", .true., ix=i_E)
    call af_add_cc_variable(tree, "rhs", .true., ix=i_rhs)
    call af_add_cc_variable(tree, "residual", .false., ix=i_residual)
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
    call CFG_add_get(cfg, "write_silo", GL_write_silo, &
         "Whether to write Silo output")
    call CFG_add_get(cfg, "write_vtk", GL_write_vtk, &
         "Whether to write VTK unstructured output")
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
    call CFG_add_get(cfg, "write_mesh_binary", GL_write_mesh_binary, &
          "Wheter the mesh is written to a binary file")
    call CFG_add_get(cfg, "write_particles_binary", GL_write_particles_binary, &
          "Wheter the particles are written to a binary file")
    call CFG_add_get(cfg, "binary_per_outputs", GL_binary_per_outputs, &
          "Per how many regular outputs binary output is written")
    call CFG_add_get(cfg, "use_electrode", GL_use_electrode, &
         "Whether to include an electrode")
    call CFG_add_get(cfg, "background_ionization_rate", &
         GL_background_ionization_rate, &
         "Background production rate (m^-3 s^-1) of electrons and positive ions")


    call CFG_add_get(cfg, "memory_afivo_GB", GL_memory_afivo_GB, &
         "How many gigabytes of memory to use for afivo")
    call CFG_add_get(cfg, "memory_particles_GB", GL_memory_particles_GB, &
         "How many gigabytes of memory to use for particles")

    write(*, "(A,F12.2)") " memory_particles_GB:     ", GL_memory_particles_GB
    write(*, "(A,F12.2)") " memory_afivo_GB:         ", GL_memory_afivo_GB


    if (GL_cylindrical) then
       call af_add_cc_variable(tree, "tmp_dens", .false., ix=i_tmp_dens)
    end if

    if (GL_use_dielectric) then
       interpolation_order_field = 1
       interpolation_order_to_density = 1

       call af_add_cc_variable(tree, "eps", ix=i_eps)
       call af_set_cc_methods(tree, i_eps, af_bc_neumann_zero, &
            af_gc_prolong_copy, af_prolong_zeroth)
    end if

    if (GL_use_electrode) then
      call af_add_cc_variable(tree, "lsf", .true., ix=i_lsf)
      bc_species => bc_species_dirichlet_zero ! TODO why this instead of neumann zero?
    end if

    call GL_rng%set_random_seed()

  end subroutine GL_initialize

  !> Convert particle coordinates to desired format
  pure function get_coordinates(my_part) result(coord)
    type(PC_part_t), intent(in) :: my_part
    real(dp)                    :: coord(NDIM)
    coord = get_coordinates_x(my_part%x)
  end function get_coordinates

  !> Convert particle coordinates to desired format
  pure function get_coordinates_x(x) result(coord)
    real(dp), intent(in) :: x(3)
    real(dp)             :: coord(NDIM)

    if (GL_cylindrical) then
       coord(1) = norm2(x([1, 3])) ! radius
       coord(2) = x(2)
    else
       coord = x(1:NDIM)
    end if
  end function get_coordinates_x

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
