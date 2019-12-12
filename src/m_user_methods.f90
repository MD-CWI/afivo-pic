#include "../afivo/src/cpp_macros.h"
!> This module contains methods that users can customize
module m_user_methods
  use m_af_all
  use m_particle_core

  implicit none
  public

  !> User-defined refinement routine
  procedure(af_subr_ref), pointer :: user_refine => null()

  !> If defined, call this routine after setting initial conditions
  procedure(init_part), pointer :: user_initial_particles => null()

  !> To set custom boundary conditions for the electric potential
  procedure(af_subr_bc), pointer :: user_potential_bc => null()

  procedure(log_subr), pointer :: user_write_log => null()

  interface
     subroutine init_part(pc)
       import
       type(PC_t), intent(inout) :: pc
     end subroutine init_part

     subroutine log_subr(tree, filename, out_cnt)
       import
       type(af_t), intent(in)       :: tree
       character(len=*), intent(in) :: filename
       integer, intent(in)          :: out_cnt
     end subroutine log_subr
  end interface

end module m_user_methods
