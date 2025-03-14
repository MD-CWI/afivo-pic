#include "../afivo/src/cpp_macros.h"
!> Module for writing output
module m_output
  use m_globals

  implicit none
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
    real(dp)                     :: sum_elec, sum_pos_ion
    real(dp)                     :: max_elec, max_field, max_Er, min_Er
    type(af_loc_t)               :: loc_elec, loc_field, loc_Er
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

    dt = GL_dt_output

    if (first_time) then
       first_time = .false.

       open(newunit=my_unit, file=trim(filename), action="write")
#if NDIM == 1
       write(my_unit, "(A)", advance="no") "it time dt v sum(n_e) sum(n_i) &
            &max(E) x max(n_e) x voltage &
            &wc_time n_cells min(dx) &
            &highest(lvl)"
#elif NDIM == 2
       write(my_unit, "(A)", advance="no") "it time dt v sum(n_e) sum(n_i) &
            &max(E) x y max(n_e) x y max(E_r) x y min(E_r) voltage &
            &wc_time n_cells min(dx) highest(lvl)"
#elif NDIM == 3
       write(my_unit, "(A)", advance="no") "it time dt v sum(n_e) sum(n_i) &
            &max(E) x y z max(n_e) x y z voltage &
            &wc_time n_cells min(dx) highest(lvl)"
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
    n_reals = 11
#elif NDIM == 2
    n_reals = 17
#elif NDIM == 3
    n_reals = 15
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
         af_r_loc(tree, loc_elec), field_voltage, &
         wc_time, af_num_cells_used(tree), &
         af_min_dr(tree),tree%highest_lvl, &
         var_values(1:n_user_vars)
#elif NDIM == 2
    write(my_unit, fmt) out_cnt, GL_time, dt, velocity, sum_elec, &
         sum_pos_ion, max_field, af_r_loc(tree, loc_field), max_elec, &
         af_r_loc(tree, loc_elec), max_Er, af_r_loc(tree, loc_Er), min_Er, &
         field_voltage, &
         wc_time, af_num_cells_used(tree), af_min_dr(tree),tree%highest_lvl, &
         var_values(1:n_user_vars)
#elif NDIM == 3
    write(my_unit, fmt) out_cnt, GL_time, dt, velocity, sum_elec, &
         sum_pos_ion,&
         max_field, af_r_loc(tree, loc_field), max_elec, &
         af_r_loc(tree, loc_elec), field_voltage, &
         wc_time, af_num_cells_used(tree), &
         af_min_dr(tree),tree%highest_lvl, &
         var_values(1:n_user_vars)
#endif
    close(my_unit)

  end subroutine output_log

end module m_output
