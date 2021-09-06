#include "../afivo/src/cpp_macros.h"
module m_refine
  use m_globals
  use m_af_all
  use m_lookup_table

  implicit none
  public

  ! The refinement buffer width in cells (around flagged cells)
  integer, protected :: refine_buffer_width = 2

  ! The number of steps after which the mesh is updated
  integer, protected :: refine_per_steps = 10

  ! The grid spacing will always be larger than this value
  real(dp), protected :: refine_min_dx = 1.0e-6_dp

  ! The grid spacing on the surface will always be smaller than this value
  real(dp), protected :: refine_surface_max_dx = 1.0e-3_dp

  ! The grid spacing will always be smaller than this value
  real(dp), protected :: refine_max_dx = 1.0e-3_dp

  ! Refine if alpha*dx is larger than this value
  real(dp), protected :: refine_adx = 1.0_dp

  ! Only refine if electron density is above this value
  real(dp), protected :: refine_elec_dens = 1.0e13_dp

  ! Only derefine if grid spacing if smaller than this value
  real(dp), protected :: derefine_dx = 5e-5_dp

  ! Ensure grid spacing around electrode is less than this value
  real(dp), protected :: refine_electrode_dx = 1e99_dp

  ! Refine a region up to this grid spacing
  real(dp), protected, allocatable :: refine_regions_dr(:)

  ! Refine regions up to this simulation time
  real(dp), protected, allocatable :: refine_regions_tstop(:)

  ! Minimum coordinate of the refinement regions
  real(dp), protected, allocatable :: refine_regions_rmin(:,:)

  ! Maximum coordinate of the refinement regions
  real(dp), protected, allocatable :: refine_regions_rmax(:,:)

  ! Table with transport data vs electric field
  type(LT_t), protected :: td_tbl

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
    call CFG_add_get(cfg, "refine%surface_max_dx", refine_surface_max_dx, &
         "The grid spacing for boxes on a dielectric surface will always be smaller than this value")
    call CFG_add_get(cfg, "refine%electrode_dx", refine_electrode_dx, &
         "Ensure grid spacing around electrode is less than this value")

    if (refine_min_dx > refine_surface_max_dx) &
     error stop "Cannot have refine_min_dx < refine_surface_max_dx"

    if (refine_min_dx > refine_max_dx) &
         error stop "Cannot have refine_min_dx < refine_max_dx"

    call CFG_add_get(cfg, "refine%adx", refine_adx, &
         "Refine if alpha*dx is larger than this value")
    call CFG_add_get(cfg, "refine%elec_dens", refine_elec_dens, &
         "Only refine if electron density is above this value")
    call CFG_add_get(cfg, "refine%derefine_dx", derefine_dx, &
         "Only derefine if grid spacing is smaller than this value")

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
    type(box_t), intent(in) :: box
    ! Refinement flags for the cells of the box
    integer, intent(out)     :: &
         cell_flags(DTIMES(box%n_cell))
    integer                  :: IJK, n, nc, id
    real(dp)                 :: max_dx, fld, alpha, adx, elec_dens
    real(dp)                 :: rmin(NDIM), rmax(NDIM)

    nc = box%n_cell
    max_dx = maxval(box%dr)

    do KJI_DO(1,nc)
       fld   = box%cc(IJK, i_E)
       alpha = LT_get_col(td_tbl, i_td_alpha, fld)
       adx   = max_dx * alpha

       elec_dens = box%cc(IJK, i_electron)

       if (adx > refine_adx .and. elec_dens > refine_elec_dens) then
          cell_flags(IJK) = af_do_ref
       else if ((adx < 0.33_dp * refine_adx .or. &
            elec_dens < 0.125_dp * refine_elec_dens) .and. &
            max_dx < derefine_dx) then
          cell_flags(IJK) = af_rm_ref
       else
          cell_flags(IJK) = af_keep_ref
       end if

       ! Refine around electrode
       if (box%tag == mg_lsf_box .and. max_dx > refine_electrode_dx) then
          if (box%cc(IJK, i_lsf) < 0) then
             cell_flags(IJK) = af_do_ref
          end if
       end if
    end do; CLOSE_DO

    ! Check fixed refinements
    rmin = box%r_min
    rmax = box%r_min + box%dr * box%n_cell

    do n = 1, size(refine_regions_dr)
       if (GL_time <= refine_regions_tstop(n) .and. &
            max_dx > refine_regions_dr(n) .and. all(&
            rmax >= refine_regions_rmin(:, n) .and. &
            rmin <= refine_regions_rmax(:, n))) then
          ! Mark just the center cell to prevent refining neighbors
          cell_flags(DTIMES(nc/2)) = af_do_ref
       end if
    end do

    ! Ensure that the minimum refinement for boxes on the surface are satisfied
    if (GL_use_dielectric .and. max_dx > refine_surface_max_dx) then
       id = box%neighbor_mat(DTIMES(0))
       if (diel%box_id_in_to_surface_ix(id) /= surface_none .or. &
            diel%box_id_out_to_surface_ix(id) /= surface_none) then
          ! Mark just the center cell to prevent refining neighbors
          cell_flags(DTIMES(nc/2)) = af_do_ref
       end if
    end if

    ! Make sure we don't have or get a too fine or too coarse grid
    if (max_dx > refine_max_dx) then
       cell_flags = af_do_ref
    else if (max_dx < 2 * refine_min_dx) then
       where(cell_flags == af_do_ref) cell_flags = af_keep_ref
    end if

  end subroutine refine_routine

  !> Initialize the transport coefficients
  subroutine load_transport_data(cfg)
    use m_config
    use m_gas

    type(CFG_t), intent(inout) :: cfg

    character(len=GL_slen)     :: td_file = "input/transport_data.txt"
    logical                    :: td_old_style = .true.
    integer                    :: table_size       = 500
    real(dp)                   :: max_electric_fld = 3e7_dp
    real(dp), allocatable      :: x_data(:), y_data(:)
    character(len=GL_slen)     :: data_name

    call CFG_add_get(cfg, "gas%transport_data_file", td_file, &
         "Input file with transport data")
    call CFG_add_get(cfg, "lookup_table_size", table_size, &
         "The transport data table size in the fluid model")
    call CFG_add_get(cfg, "lookup_table_max_efield", max_electric_fld, &
         "The maximum electric field in the fluid model coefficients")

    ! Create a lookup table for the model coefficients
    td_tbl = LT_create(0.0_dp, max_electric_fld, table_size, 1)

    call CFG_add_get(cfg, "td_old_style", td_old_style, &
            "Whether to use old or new style transport data")

    ! Fill table with data
    if (td_old_style) then
       data_name = "efield[V/m]_vs_alpha[1/m]"
       call CFG_add_get(cfg, "td_alpha_name", data_name, &
            "The name of the eff. ionization coeff. (V/m vs alpha)")
       call table_from_file(td_file, trim(data_name), x_data, y_data)
    else
       data_name = "Townsend ioniz. coef. alpha/N (m2)"
       call CFG_add_get(cfg, "td_alpha_name", data_name, &
            "The name of the eff. ionization coeff. (Td vs alpha/N)")
       call table_from_file(td_file, trim(data_name), x_data, y_data)

       ! Convert Td to V/m and alpha to 1/m
       x_data = x_data * gas_number_dens / 1e21_dp
       y_data = y_data * gas_number_dens
    end if

    call LT_set_col(td_tbl, i_td_alpha, x_data, y_data)

  end subroutine load_transport_data

  !> Routine to read in tabulated data from a file
  subroutine table_from_file(file_name, data_name, x_data, y_data)
    character(len=*), intent(in)       :: file_name, data_name
    real(dp), allocatable, intent(out) :: x_data(:), y_data(:)

    ! The maximum number of rows per entry
    integer, parameter :: table_max_rows   = 1500
    integer, parameter :: string_len = 100

    ! Temporary variables
    integer                   :: ioState, nL
    integer                   :: n_rows
    integer                   :: my_unit
    character(LEN=40)         :: line_fmt
    character(LEN=string_len) :: line
    real(dp)                  :: temp_table(2, table_max_rows)
    real(dp)                  :: factor

    nL = 0 ! Set the number of lines to 0

    ! Set the line format to read, only depends on string_len currently
    write(line_fmt, FMT = "(I6)") string_len
    line_fmt = "(A" // trim(adjustl(line_fmt)) // ")"

    ! Open 'file_name' (with error checking)
    open(newunit=my_unit, file = trim(file_name), action = "read", &
         err = 999, iostat = ioState, status="old")

    ! Table format

    !     table_name
    !     FACTOR: 1.0                   [optional: multiply with this factor]
    !     [other lines]
    !     ------------------            [at least 5 dashes]
    !     xxx       xxx                 [data in two column format]
    !     ...       ...
    !     xxx       xxx
    !     ------------------

    ! The outer DO loop, running until the end of the file is reached
    do
       ! Search for 'data_name' in the file
       do
          read(my_unit, FMT = line_fmt, ERR = 999, end = 888) line; nL = nL+1
          if (line == data_name) exit
       end do

       factor = 1.0_dp

       ! Now we can check whether there is a comment, while scanning lines until
       ! dashes are found, which indicate the start of the data
       do
          read(my_unit, FMT = line_fmt, ERR = 999, end = 777) line; nL = nL+1
          line = adjustl(line)
          if ( line(1:5) == "-----" ) then
             exit
          else if (line(1:7) == "FACTOR:") then
             read(line(8:), *) factor
          else if (line(1:8) == "COMMENT:") then
             continue
          else
             print *, "In file ", trim(file_name), " at line", nL
             print *, trim(line)
             error stop "Unknown statement in input file"
          end if
       end do

       ! Read the data into a temporary array
       n_rows = 0
       do
          read(my_unit, FMT = line_fmt, ERR = 999, end = 777) line; nL = nL+1
          line = adjustl(line)
          if ( line(1:5) == "-----" ) then
             exit  ! Dashes mark the end of the data
          else if (trim(line) == "" .or. line(1:1) == "#") then
             cycle ! Ignore whitespace or comments
          else if (n_rows < table_max_rows) then
             n_rows = n_rows + 1
             read(line, FMT = *, ERR = 999, end = 777) temp_table(:, n_rows)
          else
             print *, "CS_read_file error: too many rows in ", &
                  file_name, " at line ", nL
          end if
       end do

       ! Store the data in the actual table
       if (allocated(x_data)) deallocate(x_data)
       if (allocated(y_data)) deallocate(y_data)
       allocate(x_data(n_rows))
       allocate(y_data(n_rows))

       x_data = temp_table(1, 1:n_rows)
       y_data = factor * temp_table(2, 1:n_rows)

       exit                   ! Done
    end do

    close(my_unit)
    return

777 continue ! If the end of the file is reached after finding data
    print *, "table_from_file unexpectedly reached end of " // trim(file_name)
    print *, "searching '" // trim(data_name) // "'"
    error stop

888 continue ! If the end of the file is reached without finding data
    print *, "table_from_file: no data in " // trim(file_name)
    print *, "searching '" // trim(data_name) // "'"
    error stop

999 continue ! If there was an input error, the routine will end here
    print *, "table_from_file error at line", nL
    print *, "ioState = ", ioState, " in ", trim(file_name)
    print *, "searching '" // trim(data_name) // "'"
    error stop

  end subroutine table_from_file

end module m_refine
