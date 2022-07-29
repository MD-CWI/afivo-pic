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
  real(dp) :: dielectric_eps = 3.0_dp
  character(len=20) :: dielectric_type = "left"
  logical  :: user_init_pc = .true.

  ! Public methods
  public :: user_initialize

contains

  subroutine user_initialize(cfg)
    type(CFG_t), intent(inout) :: cfg

    call CFG_add_get(cfg, "pc%seed_pos", seed_pos, &
         "relative position of initial seed")
    call CFG_add_get(cfg, "pc%seed_sigma", seed_sigma, &
         "characteristic size of the initial seed")
    call CFG_add_get(cfg, "pc%seed_num_particles", seed_num_particles, &
         "number of particles in the seed")
    call CFG_add_get(cfg, "pc%seed_particle_weight", seed_particle_weight, &
         "weight of the particles in the seed")
    call CFG_add_get(cfg, "pc%user_init", user_init_pc, &
         "use user_initial_particles or not")
    if (user_init_pc) then
      user_initial_particles => init_particles
    end if

    call CFG_add_get(cfg, "dielectric_type", dielectric_type, &
         "What kind of dielectric to use")
    call CFG_add_get(cfg, "dielectric_eps", dielectric_eps, &
         "relative permittivity of the dielectric")
    ! user_set_surface_charge => init_surface_charge
    user_set_dielectric_eps => set_epsilon
    !user_potential_bc => my_potential
  end subroutine user_initialize

  subroutine init_particles(pctest)
    use m_particle_core
    type(PC_t), intent(inout) :: pctest
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
          call pctest%add_part(part)
       end if
    end do
  end subroutine init_particles

  real(dp) function init_surface_charge(r)
    real(dp), intent(in) :: r(2)
    real(dp)             :: z_dist
    real(dp), parameter :: rel_z_pos = 0.8_dp
    real(dp), parameter :: surface_density = 1.5e15_dp
    real(dp), parameter :: decay_length = 0.4e-3_dp

    z_dist = r(2) - rel_z_pos * domain_len(2)
    init_surface_charge = surface_density * exp(-z_dist**2/decay_length**2)
  end function init_surface_charge

  subroutine set_epsilon(box)
    type(box_t), intent(inout) :: box
    real(dp)                   :: r(2)
    integer                    :: i, j

    select case (dielectric_type)
    case ("left")
      do j = 0, box%n_cell+1
       do i = 0, box%n_cell+1
          r = af_r_cc(box, [i, j])

          if (r(1)/domain_len(1) < 0.25_dp) then
             box%cc(i, j, i_eps) = dielectric_eps
          else
             box%cc(i, j, i_eps) = 1.0_dp
          end if
       end do
      end do
    case ("left_right")
      do j = 0, box%n_cell+1
       do i = 0, box%n_cell+1
          r = af_r_cc(box, [i, j])

          if (r(1)/domain_len(1) < 0.25_dp .or. r(1)/domain_len(1) > 0.75_dp) then
             box%cc(i, j, i_eps) = dielectric_eps
          else
             box%cc(i, j, i_eps) = 1.0_dp
          end if
       end do
      end do
    case default
      error stop "Unknown dielectric_type"
    end select
  end subroutine set_epsilon

  subroutine my_potential(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type

    select case (nb)
    case (af_neighb_lowy)
        bc_type = af_bc_dirichlet
        bc_val = 0.0_dp
      case (af_neighb_highy)
        bc_type = af_bc_dirichlet
        bc_val = 12e3_dp
      case default
        bc_type = af_bc_neumann
        bc_val = 0.0_dp
    end select

  end subroutine my_potential


end module m_user
