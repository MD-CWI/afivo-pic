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
    ! user_initial_ion_density => init_ion_density
    user_set_surface_charge => init_surface_charge
    user_set_dielectric_eps => set_epsilon
    user_potential_bc => my_potential
  end subroutine user_initialize

  subroutine init_particles(pc)
    use m_particle_core
    type(PC_t), intent(inout) :: pc
    integer                   :: n
    real(dp)                  :: pos(3)
    real(dp)                  :: rand_seed(2)
    type(PC_part_t)           :: part

    part%v      = 0.0_dp
    part%a      = 0.0_dp
    part%t_left = 0.0_dp

    do n = 1, 1000
       pos(1:2) = [0.25_dp, 0.5_dp] * domain_len
       pos(3)   = 0.0_dp
       part%w   = 1.0_dp

       rand_seed = GL_rng%two_normals()
       part%x(1:2) = pos(1:2) + [rand_seed(1) * 1e-3_dp, rand_seed(2) * 1e-3_dp]
       ! part%x(1:2) = pos(1:2) + GL_rng%two_normals() * 1e-5_dp

       if (outside_check(part) <= 0) then
          call pc%add_part(part)
       end if
    end do
  end subroutine init_particles

  real(dp) function init_surface_charge(r)
    real(dp), intent(in) :: r(2)
    real(dp)             :: z_dist
    real(dp), parameter :: rel_z_pos = 0.5_dp
    real(dp), parameter :: surface_density = 1.5e15_dp
    real(dp), parameter :: decay_length = 4e-4_dp

    z_dist = r(2) - rel_z_pos * domain_len(2)
    init_surface_charge = surface_density * exp(-z_dist**2/decay_length**2)
  end function init_surface_charge

  subroutine init_ion_density(box)
    use m_geometry
    type(box_t), intent(inout) :: box
    integer                    :: i, j, nc
    real(dp)                   :: density, rr(2)

    nc = box%n_cell

    do j = 0, nc+1
       do i = 0, nc+1
          rr = af_r_cc(box, [i, j])
          density = 2e18_dp * GM_density_line(rr, &
               [0.26_dp, 0.95_dp] * domain_len, &
               [0.26_dp, 0.9_dp] * domain_len, &
               NDIM, 5e-5_dp, "smoothstep")
          box%cc(i, j, i_pos_ion) = box%cc(i, j, i_pos_ion) + &
               density
       end do
    end do
  end subroutine init_ion_density

  subroutine set_epsilon(box)
    type(box_t), intent(inout) :: box
    real(dp)                   :: r(2)
    integer                    :: i, j

    do j = 0, box%n_cell+1
       do i = 0, box%n_cell+1
          r = af_r_cc(box, [i, j])

          if (r(1)/domain_len(1) < 0.25_dp) then
             box%cc(i, j, i_eps) = 2.0_dp
          else
             box%cc(i, j, i_eps) = 1.0_dp
          end if
       end do
    end do
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
        bc_val = 3.0e4_dp
      case default
        bc_type = af_bc_neumann
        bc_val = 0.0_dp
    end select

  end subroutine my_potential


end module m_user
