!> Template for user code, this one simply uses the default routines
module m_user
  use m_config
  use m_user_methods

  implicit none
  private

  real(dp)  :: seed_pos(3) = [0.5_dp, 0.5_dp, 0.5_dp]
  integer   :: seed_num_particles = 10000
  real(dp)  :: seed_particle_weight = 1e4
  real(dp)  :: seed_sigma = 1e-4_dp

  ! Public methods
  public :: user_initialize

contains

  subroutine user_initialize(cfg)
    type(CFG_t), intent(inout) :: cfg

    call CFG_add_get(cfg, "seed_pos", seed_pos, &
         "relative position of initial seed")
    call CFG_add_get(cfg, "seed_sigma", seed_sigma, &
         "characteristic size of the initial seed")
    call CFG_add_get(cfg, "seed_num_particles", seed_num_particles, &
         "number of particles in the seed")
    call CFG_add_get(cfg, "seed_particle_weight", seed_particle_weight, &
         "weight of the particles in the seed")

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

    do n = 1, seed_num_particles
       pos(1:3) = seed_pos * domain_len
       part%w   = seed_particle_weight
       part%x(1:2) = pos(1:2) + GL_rng%two_normals() * seed_sigma
       part%x(2:3) = pos(2:3) + GL_rng%two_normals() * seed_sigma

       if (outside_check(part) <= 0) then
          call pctest%add_part(part)
       end if
    end do
  end subroutine init_particles

end module m_user
