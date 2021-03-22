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
  integer, protected :: ifc_E = -1 ! Face-centered electric field
  integer, protected :: i_E = -1     ! norm(E) (face-centered)
  integer, protected :: i_rhs      = -1 ! Source term Poisson
  integer, protected :: i_residual = -1 ! Multigrid residual
  integer, protected :: i_ppc      = -1 ! Particles per cell
  integer, protected :: i_energy   = -1 ! Energy density
  integer, protected :: i_eps      = -1 ! Dielectric permittivity
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

  ! Interval for writing to .dat file
  real(dp), protected :: GL_write_to_dat_interval(2) = 1e100_dp

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
    call af_add_cc_variable(tree, "energy", .false., ix=i_energy)

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
    
  subroutine output_log(tree, filename, out_cnt, wc_time)
    !use m_field
    use m_user_methods
    use m_pc_all
    
    type(af_t), intent(in)       :: tree
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: out_cnt !< Output number
    real(dp), intent(in)         :: wc_time !< Wallclock time
    character(len=50), save      :: fmt
    integer                      :: my_unit, n
    real(dp)                     :: velocity, dt
    real(dp), save               :: prev_pos(NDIM) = 0
    real(dp)                     :: sum_elec, sum_pos_ion
    real(dp)                     :: max_elec, max_field, max_Er, min_Er
    real(dp)                     :: sum_elem_charge, tmp, ne_zminmax(2)
    real(dp)                     :: elecdens_threshold, max_field_tip, field_voltage
    real(dp)                     :: r0(NDIM), r1(NDIM)
    type(af_loc_t)               :: loc_elec, loc_field, loc_Er, loc_tip
    integer                      :: i, n_reals, n_user_vars
    character(len=name_len)      :: var_names(user_max_log_vars)
    real(dp)                     :: var_values(user_max_log_vars)
    logical, save                :: first_time     = .true.

    if (associated(user_log_variables)) then
       ! Set default names for the user variables
       do i = 1, user_max_log_vars
          write(var_names, "(A,I0)") "uservar_", i
       end do
       call user_log_variables(tree, n_user_vars, var_names, var_values)
    else
       n_user_vars = 0
    end if

    call af_tree_sum_cc(tree, i_electron, sum_elec)
    call af_tree_sum_cc(tree, i_pos_ion, sum_pos_ion)
    call af_tree_max_cc(tree, i_electron, max_elec, loc_elec)
    call af_tree_max_cc(tree, i_E, max_field, loc_field)
    call af_tree_max_fc(tree, 1, ifc_E, max_Er, loc_Er)
    call af_tree_min_fc(tree, 1, ifc_E, min_Er)

    ! Scale threshold with gas number density
    elecdens_threshold = density_threshold * &
         (GAS_number_dens/2.414e25_dp)**2
    call analysis_zmin_zmax_threshold(tree, i_electron, elecdens_threshold, &
         [GL_domain_len(NDIM), 0.0_dp], ne_zminmax)

    ! Assume streamer tip is located at the farthest plasma density (above a
    ! threshold) away from the electrode boundaries
    r0 = 0.0
    r1 = r0 + GL_domain_len(NDIM)

    if (ne_zminmax(1)  < &
         GL_domain_len(NDIM) - ne_zminmax(2)) then
       r0(NDIM) = ne_zminmax(2) - 0.02_dp * GL_domain_len(NDIM)
       r1(NDIM) = ne_zminmax(2) + 0.02_dp * GL_domain_len(NDIM)
    else
       r0(NDIM) = ne_zminmax(1) - 0.02_dp * GL_domain_len(NDIM)
       r1(NDIM) = ne_zminmax(1) + 0.02_dp * GL_domain_len(NDIM)
    end if

    call analysis_max_var_region(tree, i_E, r0, r1, &
         max_field_tip, loc_tip)

    dt = GL_dt_output

    if (first_time) then
       first_time = .false.

       open(newunit=my_unit, file=trim(filename), action="write")
#if NDIM == 1
       write(my_unit, "(A)", advance="no") "it time dt v sum(n_e) sum(n_i) &
            &max(E) x max(n_e) x voltage ne_zmin ne_zmax &
            &max(Etip) x wc_time n_cells min(dx) &
            &highest(lvl)"
#elif NDIM == 2
       write(my_unit, "(A)", advance="no") "it time dt v sum(n_e) sum(n_i) &
            & max(E) x y max(n_e) x y max(E_r) x y min(E_r) voltage &
            &ne_zmin ne_zmax max(Etip) x y wc_time n_cells min(dx) highest(lvl)"
#elif NDIM == 3
       write(my_unit, "(A)", advance="no") "it time dt v sum(n_e) sum(n_i) &
            & max(E) x y z max(n_e) x y z voltage &
            &ne_zmin ne_zmax max(Etip) x y z wc_time n_cells min(dx) highest(lvl)"
#endif
       if (associated(user_log_variables)) then
          do i = 1, n_user_vars
             write(my_unit, "(A)", advance="no") " "//trim(var_names(i))
          end do
       end if
       write(my_unit, *) ""
       close(my_unit)

       ! Start with velocity zero
       prev_pos = af_r_loc(tree, loc_field)
    end if

#if NDIM == 1
    n_reals = 15
#elif NDIM == 2
    n_reals = 22
#elif NDIM == 3
    n_reals = 21
#endif

    if (associated(user_log_variables)) then
       write(fmt, "(A,I0,A,I0,A)") "(I6,", n_reals, "E16.8,I12,1E16.8,I3,", &
            n_user_vars, "E16.8)"
    else
       write(fmt, "(A,I0,A)") "(I6,", n_reals, "E16.8,I12,1E16.8,I3)"
    end if

    velocity = norm2(af_r_loc(tree, loc_field) - prev_pos) / GL_dt_output
    prev_pos = af_r_loc(tree, loc_field)

    open(newunit=my_unit, file=trim(filename), action="write", &
         position="append")
#if NDIM == 1
    write(my_unit, fmt) out_cnt, GL_time, dt, velocity, sum_elec, &
         sum_pos_ion, &
         max_field, af_r_loc(tree, loc_field), max_elec, &
         af_r_loc(tree, loc_elec), field_voltage, ne_zminmax, &
         max_field_tip, af_r_loc(tree, loc_tip), &
         wc_time, af_num_cells_used(tree), &
         af_min_dr(tree),tree%highest_lvl, &
         var_values(1:n_user_vars)
#elif NDIM == 2
    write(my_unit, fmt) out_cnt, GL_time, dt, velocity, sum_elec, &
         sum_pos_ion, &
         max_field, af_r_loc(tree, loc_field), max_elec, &
         af_r_loc(tree, loc_elec), max_Er, af_r_loc(tree, loc_Er), min_Er, &
         field_voltage, ne_zminmax, max_field_tip, af_r_loc(tree, loc_tip), &
         wc_time, af_num_cells_used(tree), af_min_dr(tree),tree%highest_lvl, &
         var_values(1:n_user_vars)
#elif NDIM == 3
    write(my_unit, fmt) out_cnt, GL_time, dt, velocity, sum_elec, &
         sum_pos_ion,&
         max_field, af_r_loc(tree, loc_field), max_elec, &
         af_r_loc(tree, loc_elec), field_voltage, ne_zminmax, &
         max_field_tip, af_r_loc(tree, loc_tip), &
         wc_time, af_num_cells_used(tree), &
         af_min_dr(tree),tree%highest_lvl, &
         var_values(1:n_user_vars)
#endif
    close(my_unit)

  end subroutine output_log
  
  !> Find minimum and maximum z coordinate where a variable exceeds a threshold
  subroutine analysis_zmin_zmax_threshold(tree, iv, threshold, limits, z_minmax)
    type(af_t), intent(in) :: tree
    integer, intent(in)    :: iv !< Index of variable
    !> Threshold for variable
    real(dp), intent(in)   :: threshold
    !> Limits for min/max
    real(dp), intent(in)   :: limits(2)
    !> Minimum/maximum z coordinate above threshold
    real(dp), intent(out)  :: z_minmax(2)

    call af_reduction_vec(tree, box_minmax_z, reduce_minmax, &
         limits, z_minmax, 2)

  contains

    ! Find cell with min/max z coordinate that has a density exceeding a
    ! threshold
    function box_minmax_z(box, n_vals) result(vec)
      type(box_t), intent(in) :: box
      integer, intent(in)     :: n_vals
      real(dp)                :: vec(n_vals)

      integer  :: i, j, n, nc, ix(NDIM)
      logical  :: above
      real(dp) :: r(NDIM)

      nc = box%n_cell
      i = -1
      j = -1

      do n = 1, nc
#if NDIM == 1
         above = box%cc(n, iv) > threshold
#elif NDIM == 2
         above = maxval(box%cc(1:nc, n, iv)) > threshold
#elif NDIM == 3
         above = maxval(box%cc(1:nc, 1:nc, n, iv)) > threshold
#endif
         if (above) then
            if (i == -1) i = n
            j = n
         end if
      end do

      vec = [1e100_dp, -1e100_dp]
      ix = 1
      if (i /= -1) then
         ix(NDIM) = i
         r = af_r_cc(box, ix)
         vec(1) = r(NDIM)
      end if

      if (j /= -1) then
         ix(NDIM) = i
         r = af_r_cc(box, ix)
         vec(2) = r(NDIM)
      end if
    end function box_minmax_z

    !> Reduction method (e.g., min, max, sum)
    function reduce_minmax(vec_1, vec_2, n_vals) result(vec)
      integer, intent(in)  :: n_vals
      real(dp), intent(in) :: vec_1(n_vals), vec_2(n_vals)
      real(dp)             :: vec(n_vals)
      vec(1) = min(vec_1(1), vec_2(1))
      vec(2) = max(vec_1(2), vec_2(2))
    end function reduce_minmax
  end subroutine analysis_zmin_zmax_threshold
  
  !> Find maximal value for boxes that are (at least partially) in the rectangle
  !> from r0 to r1
  subroutine analysis_max_var_region(tree, iv, r0, r1, max_value, loc)
    type(af_t), intent(in)      :: tree
    integer, intent(in)         :: iv        !< Index of variable
    real(dp), intent(in)        :: r0(NDIM)  !< Minimum coordinates
    real(dp), intent(in)        :: r1(NDIM)  !< Maximum coordinates
    real(dp), intent(out)       :: max_value !< Found maximum
    type(af_loc_t), intent(out) :: loc

    call af_reduction_loc(tree, iv, box_max_region, reduce_max, &
         -1e100_dp, max_value, loc)

  contains

    subroutine box_max_region(box, iv, val, ix)
      type(box_t), intent(in) :: box
      integer, intent(in)     :: iv
      real(dp), intent(out)   :: val
      integer, intent(out)    :: ix(NDIM)
      integer                 :: nc
      real(dp)                :: r_max(NDIM)

      nc = box%n_cell
      r_max = box%r_min + box%n_cell * box%dr

      if (any(box%r_min > r1) .or. any(r_max < r0)) then
         ix = -1
         val = -1e100_dp
      else
         ix = maxloc(box%cc(DTIMES(1:nc), iv))
         val = box%cc(DINDEX(ix), iv)
      end if
    end subroutine box_max_region

  end subroutine analysis_max_var_region  
    
  real(dp) function reduce_max(a, b)
    real(dp), intent(in) :: a, b
    reduce_max = max(a, b)
  end function reduce_max  

end module m_globals
