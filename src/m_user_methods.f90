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

  !> If defined, call this routine to set initial ion densities
  procedure(init_ions), pointer :: user_initial_ion_density => null()

  !> Call this routine every time step to generate particles
  procedure(gen_part), pointer :: user_generate_particles => null()

  !> If defined, call this routine to set the dielectric permittivity
  procedure(af_subr), pointer :: user_set_dielectric_eps => null()

  !> If defined, set the surface charge on a dielectric surface
  procedure(sigma_func), pointer :: user_set_surface_charge => null()

  !> To set custom boundary conditions for the electric potential
  procedure(af_subr_bc), pointer :: user_potential_bc => null()

  procedure(log_subr), pointer :: user_write_log => null()

  interface
     subroutine init_part(pc)
       import
       type(PC_t), intent(inout) :: pc
     end subroutine init_part

     subroutine init_ions(box)
       import
       type(box_t), intent(inout) :: box
     end subroutine init_ions

     subroutine gen_part(pc, time, time_elapsed)
       import
       type(PC_t), intent(inout) :: pc
       real(dp), intent(in)      :: time         !< Current time
       real(dp), intent(in)      :: time_elapsed !< Time since last call
     end subroutine gen_part

     subroutine log_subr(tree, filename, out_cnt)
       import
       type(af_t), intent(in)       :: tree
       character(len=*), intent(in) :: filename
       integer, intent(in)          :: out_cnt
     end subroutine log_subr

     real(dp) function sigma_func(r)
       import
       real(dp), intent(in) :: r(NDIM)
     end function sigma_func
  end interface

end module m_user_methods
