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
    ! user_generate_particles => init_electrons_only
    user_set_dielectric_eps => set_epsilon
    ! user_set_surface_charge => set_surface_charge
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

    do n = 1, 50000
       pos(1:2) = [0.5_dp, 0.15_dp] * domain_len
       pos(3)   = 0.0_dp
       part%w   = 1.0_dp
       part%x(1:2) = pos(1:2) + GL_rng%two_normals() * 1e-4_dp

       if (outside_check(part) <= 0) then
          call pc%add_part(part)
       end if
    end do
  end subroutine init_particles

  subroutine init_electrons_only(pc, time, time_elapsed)
    ! This model can be used to generate electrons (without ions)
    use m_particle_core
    type(PC_t), intent(inout)  :: pc
    real(dp), intent(in)      :: time         !< Current time
    real(dp), intent(in)      :: time_elapsed !< Time since last call

    integer                   :: n
    real(dp)                  :: pos(3)
    type(PC_part_t)           :: part

    if (time <= 0.0_dp) then
          part%v      = 0.0_dp
          part%a      = 0.0_dp
          part%t_left = 0.0_dp

          do n = 1, 1000
             pos(1:2) = [0.5_dp, 0.15_dp] * domain_len
             pos(3)   = 0.0_dp
             part%w   = 1.0_dp
             part%x(1:2) = pos(1:2) + GL_rng%two_normals() * 1.0e-5_dp

             if (outside_check(part) <= 0) then
                call pc%add_part(part)
             end if
          end do
    end if
  end subroutine

  subroutine set_epsilon(box)
    type(box_t), intent(inout) :: box
    real(dp)                   :: r(2)
    integer                    :: i, j

    do j = 0, box%n_cell+1
       do i = 0, box%n_cell+1
          r = af_r_cc(box, [i, j])

          if (r(2)/domain_len(2) < 0.25_dp) then
             box%cc(i, j, i_eps) = 4.5_dp
          else
             box%cc(i, j, i_eps) = 1.0_dp
          end if
       end do
    end do
  end subroutine set_epsilon

  real(dp) function set_surface_charge(r)
    ! Set the surface charge on all dielectric surfaces as a function of r
    real(dp), intent(in) :: r(NDIM)
    real(dp)  :: coord, var, amplitude


    coord     = r(1)/domain_len(1) - 0.5_dp
    var       = 5.0e-5_dp
    amplitude = 3.5e16_dp

    set_surface_charge = amplitude * exp(- coord**2 / var)
  end function set_surface_charge

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
        bc_val = 0.0_dp
      case (af_neighb_lowy)
        bc_type = af_bc_dirichlet
        do ii = 1, box%n_cell**(NDIM-1)
          bc_val(ii) = 1.5e4_dp + 6.0e3_dp * exp( - (coords(1, ii)/domain_len(1) - 0.5_dp)**2 / 2.0e-2_dp)
        end do
      case default
        bc_type = af_bc_neumann
        bc_val = 0.0_dp
    end select

  end subroutine my_potential


end module m_user
