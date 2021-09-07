!> Template for user code, this one simply uses the default routines
module m_user
  use m_config
  use m_user_methods
  use m_globals
  use m_domain

  implicit none
  private

  real(dp) :: seed_pos(2) = [0.5_dp, 0.5_dp]
  real(dp) :: seed_sigma = 1e-4_dp
  integer :: seed_num_particles = 10000
  real(dp) :: seed_particle_weight = 1e4

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
    ! user_potential_bc => my_potential
  end subroutine user_initialize

  subroutine init_particles(pc)
    use m_particle_core
    type(PC_t), intent(inout) :: pc
    integer                   :: n
    real(dp)                  :: tmp_vec(2)
    type(PC_part_t)           :: part

    part%v      = 0.0_dp
    part%a      = 0.0_dp
    part%t_left = 0.0_dp

    do n = 1, seed_num_particles
       if (GL_cylindrical) then
          tmp_vec = GL_rng%two_normals() * seed_sigma
          part%x(1) = norm2(tmp_vec)
          tmp_vec = GL_rng%two_normals() * seed_sigma
          part%x(2) = tmp_vec(1)
          part%x(1:2) = part%x(1:2) + seed_pos * domain_len
          part%x(3) = 0.0_dp
       else
          part%x(1:2) = GL_rng%two_normals() * seed_sigma
          part%x(1:2) = part%x(1:2) + seed_pos * domain_len
          part%x(3) = 0.0_dp
       end if

       part%w = seed_particle_weight

       if (outside_check(part) <= 0) then
          call pc%add_part(part)
       end if
    end do
  end subroutine init_particles

  subroutine my_potential(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type
    integer                 :: ii

    select case (nb)
    case (af_neighb_highy)
        bc_type = af_bc_dirichlet
        bc_val = -2.5e4_dp
      case (af_neighb_lowy)
        bc_type = af_bc_dirichlet
        do ii = 1, box%n_cell**(NDIM-1)
          bc_val(ii) = 1.75e3_dp * exp( - (coords(1, ii)/domain_len(1) - 0.5_dp)**2 / 5.0e-4_dp)
        end do
      case default
        bc_type = af_bc_neumann
        bc_val = 0.0_dp
    end select

  end subroutine my_potential

end module m_user
