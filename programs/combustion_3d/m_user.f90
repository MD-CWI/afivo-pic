!> Template for user code, this one simply uses the default routines
module m_user
  use m_config
  use m_user_methods
  use m_globals
  use m_domain

  implicit none
  private

  real(dp)  :: seed_location(3) = [0.5, 0.5, 0.825]
  real(dp)  :: seed_variance = 2.5e-4
  integer   :: seed_particles = 10000

  ! Public methods
  public :: user_initialize

contains

  subroutine user_initialize(cfg)
    type(CFG_t), intent(inout) :: cfg

    call CFG_add_get(cfg, "seed%location", seed_location, &
        "Location of the initial seed")
    call CFG_add_get(cfg, "seed%variance", seed_variance, &
        "Variance of the (Gaussian) initial seed")
    call CFG_add_get(cfg, "seed%num_particles", seed_particles, &
        "Number of particles in the initial seed")

    user_initial_particles => init_particles
  end subroutine user_initialize

  subroutine init_particles(pc)
    use m_particle_core
    type(PC_t), intent(inout) :: pc
    integer                   :: n
    real(dp)                  :: pos(3)
    type(PC_part_t)           :: part

    part%v      = 0.0_dp
    part%a      = 0.0_dp
    part%t_left = 0.0_dp

    do n = 1, seed_particles
       pos(1:3) = seed_location * domain_len
       part%w   = 1.0_dp
       part%x(1:3) = pos(1:3) + GL_rng%three_normals() * seed_variance

       if (outside_check(part) <= 0) then
          call pc%add_part(part)
       end if
    end do
  end subroutine init_particles

end module m_user
