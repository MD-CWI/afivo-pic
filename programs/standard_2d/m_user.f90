!> Template for user code, this one simply uses the default routines
module m_user
  use m_config
  use m_user_methods
  use m_globals
  use m_domain

  implicit none
  private

  ! Public methods
  public :: user_initialize

contains

  subroutine user_initialize(cfg)
    type(CFG_t), intent(inout) :: cfg

    user_initial_particles => init_particles
    user_set_dielectric_eps => set_epsilon
  end subroutine user_initialize

  subroutine init_particles(pctest)
    use m_particle_core
    type(PC_t), intent(inout) :: pctest
    integer                   :: n
    real(dp)                  :: pos(3)
    type(PC_part_t)           :: part

    part%v      = 0.0_dp
    part%a      = 0.0_dp
    part%t_left = 0.0_dp

    do n = 1, 100
       pos(1:2) = [0.5_dp, 0.75_dp] * domain_len
       pos(3)   = 0.0_dp
       part%w   = 1.0_dp
       part%x(1:2) = pos(1:2) + GL_rng%two_normals() * 1e-5_dp

       if (outside_check(part) <= 0) then
          call pctest%add_part(part)
       end if
    end do
  end subroutine init_particles

  subroutine set_epsilon(box)
    type(box_t), intent(inout) :: box
    real(dp)                   :: r(2)
    integer                    :: i, j

    do j = 0, box%n_cell+1
       do i = 0, box%n_cell+1
          r = af_r_cc(box, [i, j])

          if (r(2)/domain_len(2) < 0.5_dp) then
             box%cc(i, j, i_eps) = 100.0_dp
          else
             box%cc(i, j, i_eps) = 1.0_dp
          end if
       end do
    end do
  end subroutine set_epsilon

end module m_user
