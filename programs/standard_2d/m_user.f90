!> Template for user code, this one simply uses the default routines
module m_user
  use m_config
  use m_user_methods

  implicit none
  private

  ! Public methods
  public :: user_initialize

contains

  subroutine user_initialize(cfg)
    type(CFG_t), intent(inout) :: cfg

    user_initial_particles => init_particles
  end subroutine user_initialize

  subroutine init_particles(pctest)
    use m_particle_core
    use m_domain
    use m_globals
    type(PC_t), intent(inout) :: pctest
    integer                   :: n
    real(dp)                  :: pos(3)
    type(PC_part_t)           :: part

    part%v      = 0.0_dp
    part%a      = 0.0_dp
    part%t_left = 0.0_dp

    do n = 1, 100
       pos(1:2) = 0.5_dp * domain_len
       pos(3)   = 0.0_dp
       part%w   = 1.0_dp
       part%x(1:2) = pos(1:2) + GL_rng%two_normals() * 1e-3_dp

       if (outside_check(part) <= 0) then
          call pctest%add_part(part)
       end if
    end do
  end subroutine init_particles

end module m_user
