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
  ! elaspe time to generate new source particles
  real(dp)  :: elaps_set_time = 5e-10

contains

  subroutine user_initialize(cfg)
    type(CFG_t), intent(inout) :: cfg

    user_initial_particles => init_particles
    ! user_generate_particles => background_charge
    user_set_dielectric_eps => set_epsilon
    user_set_dielectric_charge => set_initial_charge
    user_potential_bc => my_potential

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


    do n = 1,100
      pos(1:2) = [0.4_dp, 0.3_dp] * domain_len
       pos(3)   = 0.0_dp
       part%w   = 1.0_dp
       part%x(1:2) = pos(1:2) + GL_rng%two_normals() * 1e-5_dp

       if (outside_check(part) <= 0) then
          call pctest%add_part(part)
       end if
    end do
  end subroutine init_particles

  subroutine background_charge(particle, current_time, elaps_time)
    use m_particle_core
    type(PC_t), intent(inout) :: particle
    real(dp), intent(in)      :: current_time        !< Current time
    real(dp), intent(in)      :: elaps_time !< Time since last call
    integer                   :: n
    real(dp)                  :: pos(3) , temp
    type(PC_part_t)           :: part

    part%v      = 0.0_dp
    part%a      = 0.0_dp
    part%t_left = 0.0_dp

    temp = abs(elaps_time - elaps_set_time)
    if ( temp <= 1e-12) then

      do n = 1, 5000
        pos(1:2) = [0.6_dp, 0.3_dp] * domain_len
        pos(3)   = 0.0_dp
        part%w   = 1.0_dp
        part%x(1:2) = pos(1:2) + GL_rng%two_normals() * 1e-5_dp

          if (outside_check(part) <= 0) then
            call particle%add_part(part)
          end if
      end do

       do n = 1, 100
          pos(1:2) = [0.8_dp, 0.3_dp] * domain_len
          pos(3)   = 0.0_dp
          part%w   = 1.0_dp
          part%x(1:2) = pos(1:2) + GL_rng%two_normals() * 1e-5_dp

          if (outside_check(part) <= 0) then
             call particle%add_part(part)
          end if
        end do

        do n = 1, 10
           pos(1:2) = [1.2_dp, 0.3_dp] * domain_len
           pos(3)   = 0.0_dp
           part%w   = 1.0_dp
           part%x(1:2) = pos(1:2) + GL_rng%two_normals() * 1e-5_dp

           if (outside_check(part) <= 0) then
              call particle%add_part(part)
           end if
         end do

        ! do n = 1, 10000
        !    pos(1:2) = [0.8_dp, 0.5_dp] * domain_len
        !    pos(3)   = 0.0_dp
        !    part%w   = 1.0_dp
        !    part%x(1:2) = pos(1:2) + GL_rng%two_normals() * 1e-5_dp
        !
        !    if (outside_check(part) <= 0) then
        !       call particle%add_part(part)
        !    end if
        !  end do
    end if

  end subroutine background_charge

  subroutine set_epsilon(box)
    type(box_t), intent(inout) :: box
    real(dp)                   :: r(2)
    integer                    :: i, j

    do j = 0, box%n_cell+1
       do i = 0, box%n_cell+1
          r = af_r_cc(box, [i, j])
          if (r(2)/domain_len(2) <= 0.25_dp) then
             box%cc(i, j, i_eps) = 10.0_dp
          else if (r(2)/domain_len(2) >= 0.75_dp) then
             box%cc(i, j, i_eps) = 10.0_dp
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

    if (af_neighb_dim(nb) == NDIM) then
       if (af_neighb_low(nb)) then
          bc_type = af_bc_dirichlet
          bc_val = 0.0_dp
       else
          bc_type = af_bc_dirichlet
          bc_val  = 6e3_dp!coords(1, :) / domain_len(1) * (-3e3_dp)
       end if
    else
       bc_type = af_bc_neumann
       bc_val = 0.0_dp
    end if
  end subroutine my_potential

  real(dp) function set_initial_charge(r)
    real(dp), intent(in)    :: r(NDIM)  !< the corrdinates of the cell in dielectric surface

    if (r(1)/domain_len(1) <= 0.3_dp) then
        if (r(2)/domain_len(2) <= 0.5_dp)  then
          set_initial_charge = (1-r(1)/domain_len(1) )*(-4.0e15_dp)
        else
          set_initial_charge = (1-r(1)/domain_len(1) )*(4.0e15_dp)
        end if
    else
        set_initial_charge = 0.0_dp
    end if



    ! Screen electric field outside dielectric
    ! sigma_function = -epsilon_high/interface_location

    ! Screen electric field inside dielectric
    ! sigma_function = 1/(1 - interface_location)

    ! Equal electric field on left and right
    ! sigma_function = 1 - epsilon_high

    ! No surface charge
    ! set_initial_charge = 1.0e12_dp
  end function set_initial_charge


end module m_user
