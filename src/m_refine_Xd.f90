#include "afivo/src/cpp_macros_$Dd.h"
module m_refine_$Dd
  use m_globals_$Dd
  use m_a$D_all
  use m_lookup_table

  implicit none
  public

  ! The refinement buffer width in cells (around flagged cells)
  integer, protected :: refine_buffer_width = 2

  ! The number of steps after which the mesh is updated
  integer, protected :: refine_per_steps = 10

  ! The grid spacing will always be larger than this value
  real(dp), protected :: refine_min_dx = 1.0e-6_dp

  ! The grid spacing will always be smaller than this value
  real(dp), protected :: refine_max_dx = 1.0e-3_dp

  ! Refine if alpha*dx is larger than this value
  real(dp), protected :: refine_adx = 1.0_dp

  ! Refine if dx^2 * rhs is larger than this value (Volts)
  real(dp), protected :: refine_cphi = 1.0e50_dp

  ! Only refine if electron density is above this value
  real(dp), protected :: refine_elec_dens = 1.0e13_dp

  ! Only derefine if grid spacing if smaller than this value
  real(dp), protected :: derefine_dx = 5e-5_dp

  ! Refine around initial conditions up to this time
  real(dp), protected :: refine_init_time = 10e-9_dp

  ! Refine until dx is smaller than this factor times the seed width
  real(dp), protected :: refine_init_fac = 0.25_dp

  ! Refine a region up to this grid spacing
  real(dp), protected, allocatable :: refine_regions_dr(:)

  ! Refine regions up to this simulation time
  real(dp), protected, allocatable :: refine_regions_tstop(:)

  ! Minimum coordinate of the refinement regions
  real(dp), protected, allocatable :: refine_regions_rmin(:,:)

  ! Maximum coordinate of the refinement regions
  real(dp), protected, allocatable :: refine_regions_rmax(:,:)

  ! Table with transport data vs electric field
  type(lookup_table_t), protected :: td_tbl

  ! Index of alpha transport data
  integer, parameter :: i_td_alpha = 1

contains

  subroutine refine_init(cfg, ndim)
    use m_config
    type(CFG_t), intent(inout) :: cfg
    integer, intent(in)        :: ndim !< Number of dimensions
    integer                    :: n
    real(dp)                   :: vec(ndim)
    real(dp), allocatable      :: dbuffer(:)

    call CFG_add_get(cfg, "refine%buffer_width", refine_buffer_width, &
         "The refinement buffer width in cells (around flagged cells)")
    call CFG_add_get(cfg, "refine%per_steps", refine_per_steps, &
         "The number of steps after which the mesh is updated")
    call CFG_add_get(cfg, "refine%min_dx", refine_min_dx, &
         "The grid spacing will always be larger than this value")
    call CFG_add_get(cfg, "refine%max_dx", refine_max_dx, &
         "The grid spacing will always be smaller than this value")

    if (refine_min_dx > refine_max_dx) &
         error stop "Cannot have refine_min_dx < refine_max_dx"

    call CFG_add_get(cfg, "refine%adx", refine_adx, &
         "Refine if alpha*dx is larger than this value")
    call CFG_add_get(cfg, "refine%cphi", refine_cphi, &
         "Refine if the curvature in phi is larger than this value")
    call CFG_add_get(cfg, "refine%elec_dens", refine_elec_dens, &
         "Only refine if electron density is above this value")
    call CFG_add_get(cfg, "refine%derefine_dx", derefine_dx, &
         "Only derefine if grid spacing if smaller than this value")
    call CFG_add_get(cfg, "refine%init_time", refine_init_time, &
         "Refine around initial conditions up to this time")
    call CFG_add_get(cfg, "refine%init_fac", refine_init_fac, &
         "Refine until dx is smaller than this factor times the seed width")

    call CFG_add(cfg, "refine%regions_dr", [1.0e99_dp], &
         "Refine regions up to this grid spacing", .true.)
    call CFG_add(cfg, "refine%regions_tstop", [-1.0e99_dp], &
         "Refine regions up to this simulation time", .true.)
    vec = 0.0_dp
    call CFG_add(cfg, "refine%regions_rmin", vec, &
         "Minimum coordinate of the refinement regions", .true.)
    call CFG_add(cfg, "refine%regions_rmax", vec, &
         "Maximum coordinate of the refinement regions", .true.)

    call CFG_get_size(cfg, "refine%regions_dr", n)
    allocate(refine_regions_dr(n))
    allocate(refine_regions_tstop(n))
    allocate(refine_regions_rmin(ndim, n))
    allocate(refine_regions_rmax(ndim, n))
    allocate(dbuffer(ndim * n))

    call CFG_get(cfg, "refine%regions_dr", refine_regions_dr)
    call CFG_get(cfg, "refine%regions_tstop", refine_regions_tstop)
    call CFG_get(cfg, "refine%regions_rmin", dbuffer)
    refine_regions_rmin = reshape(dbuffer, [ndim, n])
    call CFG_get(cfg, "refine%regions_rmax", dbuffer)
    refine_regions_rmax = reshape(dbuffer, [ndim, n])

    call load_transport_data(cfg)

  end subroutine refine_init

  ! This routine sets the cell refinement flags for box
  subroutine refine_routine(box, cell_flags)
    use m_geometry
    use m_init_cond_$Dd
    type(box$D_t), intent(in) :: box
    ! Refinement flags for the cells of the box
    integer, intent(out)     :: &
         cell_flags(DTIMES(box%n_cell))
    integer                  :: IJK, n, nc
    real(dp)                 :: dx, dx2, fld, alpha, adx, cphi, elec_dens
    real(dp)                 :: rmin($D), rmax($D)

    nc = box%n_cell
    dx = box%dr
    dx2 = dx**2

    do KJI_DO(1,nc)
       fld   = box%cc(IJK, i_E)
       alpha = LT_get_col(td_tbl, i_td_alpha, fld)
       adx   = box%dr * alpha

       ! The refinement is also based on the intensity of the source term.
       ! Here we estimate the curvature of phi (given by dx**2 *
       ! Laplacian(phi))
       cphi = dx2 * abs(box%cc(IJK, i_rhs))

       elec_dens = box%cc(IJK, i_electron)

       if ((adx > refine_adx .and. elec_dens > refine_elec_dens) &
            .or. cphi > refine_cphi) then
          cell_flags(IJK) = af_do_ref
       else if (adx < 0.125_dp * refine_adx .or. &
            elec_dens < 0.125_dp * refine_elec_dens) then
          cell_flags(IJK) = af_rm_ref
       else
          cell_flags(IJK) = af_keep_ref
       end if
    end do; CLOSE_DO

    ! Check fixed refinements
    rmin = box%r_min
    rmax = box%r_min + box%dr * box%n_cell

    do n = 1, size(refine_regions_dr)
       if (ST_time <= refine_regions_tstop(n) .and. &
            dx > refine_regions_dr(n) .and. all(&
            rmax >= refine_regions_rmin(:, n) .and. &
            rmin <= refine_regions_rmax(:, n))) then
          ! Mark just the center cell to prevent refining neighbors
          cell_flags(DTIMES(nc/2)) = af_do_ref
       end if
    end do

    ! Make sure we don't have or get a too fine or too coarse grid
    if (dx > refine_max_dx) then
       cell_flags = af_do_ref
    else if (dx < 2 * refine_min_dx) then
       where(cell_flags == af_do_ref) cell_flags = af_keep_ref
    end if

  end subroutine refine_routine

  !> Initialize the transport coefficients
  subroutine load_transport_data(cfg)
    use m_transport_data
    use m_config

    type(CFG_t), intent(inout) :: cfg

    character(len=ST_slen)     :: td_file = "input/transport_data.txt"
    integer                    :: table_size       = 500
    real(dp)                   :: max_electric_fld = 3e7_dp
    real(dp), allocatable      :: x_data(:), y_data(:)
    character(len=ST_slen)     :: data_name

    call CFG_add_get(cfg, "gas%transport_data_file", td_file, &
         "Input file with transport data")
    call CFG_add_get(cfg, "lookup_table_size", table_size, &
         "The transport data table size in the fluid model")
    call CFG_add_get(cfg, "lookup_table_max_efield", max_electric_fld, &
         "The maximum electric field in the fluid model coefficients")

    ! Create a lookup table for the model coefficients
    td_tbl = LT_create(0.0_dp, max_electric_fld, table_size, 1)

    ! Fill table with data
    data_name = "efield[V/m]_vs_alpha[1/m]"
    call CFG_add_get(cfg, "td_alpha_name", data_name, &
         "The name of the eff. ionization coeff.")
    call TD_get_from_file(td_file, ST_gas_name, &
         trim(data_name), x_data, y_data)
    call LT_set_col(td_tbl, i_td_alpha, x_data, y_data)

  end subroutine load_transport_data

end module m_refine_$Dd
