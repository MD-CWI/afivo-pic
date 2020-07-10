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
    user_potential_bc => my_potential
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

    do n = 1, 100
       pos(1:3) = [0.5_dp, 0.5_dp, 0.7_dp] * domain_len
       part%w   = 1.0_dp
       part%x(1:2) = pos(1:2) + GL_rng%two_normals() * 1e-5_dp
       part%x(2:3) = pos(2:3) + GL_rng%two_normals() * 1e-5_dp

       if (outside_check(part) <= 0) then
          call pc%add_part(part)
       end if
    end do
  end subroutine init_particles

  subroutine set_epsilon(box)
    type(box_t), intent(inout) :: box
    real(dp)                   :: r(3)
    integer                    :: i, j, k

    do j = 0, box%n_cell+1
       do i = 0, box%n_cell+1
         do k = 0, box%n_cell+1
            r = af_r_cc(box, [i, j, k])

            if (r(3)/domain_len(3) < 0.25_dp) then
               box%cc(i, j, k, i_eps) = 4.5_dp
            else
               box%cc(i, j, k, i_eps) = 1.0_dp
            end if
         end do
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
    case (af_neighb_lowz)
        bc_type = af_bc_dirichlet
        bc_val = 0.0_dp
      case (af_neighb_highz)
        bc_type = af_bc_dirichlet
        bc_val = - 1.65e4_dp
      case default
        bc_type = af_bc_neumann
        bc_val = 0.0_dp
    end select

  end subroutine my_potential


end module m_user
