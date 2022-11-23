#include "../afivo/src/cpp_macros.h"
!> Module for writing output
module m_output
  use m_globals

  private

  public :: output_log

contains

  subroutine output_log(tree, filename, out_cnt, wc_time)
    !use m_field
    use m_user_methods
    use m_pc_all
    use m_field
    type(af_t), intent(in)       :: tree
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: out_cnt !< Output number
    real(dp), intent(in)         :: wc_time !< Wallclock time
    character(len=50), save      :: fmt
    integer                      :: my_unit
    real(dp)                     :: velocity, dt
    real(dp), save               :: prev_pos(NDIM) = 0
    real(dp)                     :: sum_elec, sum_pos_ion, sum_neg_ion
    real(dp)                     :: max_elec, max_field, max_Er, min_Er
    real(dp)                     :: ne_zminmax(2)
    real(dp)                     :: elecdens_threshold, max_field_tip
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
    call af_tree_sum_cc(tree, i_neg_ion, sum_neg_ion)
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
    n_reals = 16
#elif NDIM == 2
    n_reals = 23
#elif NDIM == 3
    n_reals = 22
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
         sum_pos_ion, sum_neg_ion, &
         max_field, af_r_loc(tree, loc_field), max_elec, &
         af_r_loc(tree, loc_elec), field_voltage, ne_zminmax, &
         max_field_tip, af_r_loc(tree, loc_tip), &
         wc_time, af_num_cells_used(tree), &
         af_min_dr(tree),tree%highest_lvl, &
         var_values(1:n_user_vars)
#elif NDIM == 2
    write(my_unit, fmt) out_cnt, GL_time, dt, velocity, sum_elec, &
         sum_pos_ion, sum_neg_ion, &
         max_field, af_r_loc(tree, loc_field), max_elec, &
         af_r_loc(tree, loc_elec), max_Er, af_r_loc(tree, loc_Er), min_Er, &
         field_voltage, ne_zminmax, max_field_tip, af_r_loc(tree, loc_tip), &
         wc_time, af_num_cells_used(tree), af_min_dr(tree),tree%highest_lvl, &
         var_values(1:n_user_vars)
#elif NDIM == 3
    write(my_unit, fmt) out_cnt, GL_time, dt, velocity, sum_elec, &
         sum_pos_ion, sum_neg_ion, &
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

end module m_output
