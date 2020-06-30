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
    user_initial_particles_and_ions => init_particles_and_ions
    ! user_generate_particles => background_charge
    user_set_dielectric_eps => set_epsilon!set_epsilon_boxwise!
    ! user_set_dielectric_charge => set_surface_charge
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

    do n = 1,100
       pos(1:2) = [0.2_dp, 0.35_dp] * domain_len
       pos(3)   = 0.0_dp
       part%w   = 1.0_dp
       part%x(1:2) = pos(1:2) + GL_rng%two_normals() * 1e-5_dp

       if (outside_check(part) <= 0) then
          call pc%add_part(part)
       end if
    end do

  end subroutine init_particles

  subroutine init_particles_and_ions(pc_elec, pc_ions)
    use m_particle_core
    type(PC_t), intent(inout) :: pc_elec, pc_ions
    integer                   :: n
    real(dp)                  :: pos(3)
    type(PC_part_t)           :: part

    part%v      = 0.0_dp
    part%a      = 0.0_dp
    part%t_left = 0.0_dp

    do n = 1,100
       pos(1:2) = [0.2_dp, 0.35_dp] * domain_len
       pos(3)   = 0.0_dp
       part%w   = 1.0_dp
       part%x(1:2) = pos(1:2) + GL_rng%two_normals() * 1e-5_dp

       if (outside_check(part) <= 0) then
          call pc_elec%add_part(part)
          call pc_ions%add_part(part)
       end if
    end do

  end subroutine init_particles_and_ions

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
          if (r(2)/domain_len(2) < 0.3_dp) then
             box%cc(i, j, i_eps) = 10.0_dp
          else if (r(2)/domain_len(2) > 0.7_dp) then
             box%cc(i, j, i_eps) = 10.0_dp
          else
             box%cc(i, j, i_eps) = 1.0_dp
          end if
       end do
    end do
  end subroutine set_epsilon

  subroutine set_epsilon_boxwise(box)
    type(box_t), intent(inout) :: box
    real(dp)                   :: r(2)
    integer                    :: i, j
    real(dp)  :: die_pos_l, die_pos_h, temp

    temp = box%dr(2) * box%n_cell

    die_pos_l = domain_len(2) * 0.2_dp
    die_pos_h = domain_len(2) * 0.8_dp

    if (die_pos_l > box%r_min(2) .and. die_pos_l < box%r_min(2) + temp) &
      die_pos_l = box%r_min(2)
    if (die_pos_h > box%r_min(2) .and. die_pos_h < box%r_min(2) + temp) &
      die_pos_h = box%r_min(2)

      ! print *, "the low dielectric position", die_pos_l
      ! print *, "the high dielectric position", die_pos_h
    do j = 0, box%n_cell+1
       do i = 0, box%n_cell+1
          r = af_r_cc(box, [i, j])
          if (r(2) < die_pos_l) then
             box%cc(i, j, i_eps) = 10.0_dp
          else if (r(2) > die_pos_h) then
             box%cc(i, j, i_eps) = 10.0_dp
          else
             box%cc(i, j, i_eps) = 1.0_dp
          end if
       end do
    end do
  end subroutine set_epsilon_boxwise


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
        bc_val = 5.0e2_dp !1.5e3_dp
        bc_val = 2e3_dp !1.5e3_dp
      case default
        bc_type = af_bc_neumann
        bc_val = 0.0_dp
    end select
  end subroutine my_potential

  real(dp) function set_surface_charge(r)
    real(dp), intent(in)    :: r(NDIM)  !< the corrdinates of the cell in dielectric surface

    if (r(1)/domain_len(1) <= 0.1_dp) then
        if (r(2)/domain_len(2) <= 0.5_dp)  then
          set_surface_charge = (1-r(1)/domain_len(1) )*(3.1e15_dp) ! !2.50e15_dp
        else
          set_surface_charge = (1-r(1)/domain_len(1) )*(-3.1e15_dp)
        end if
    else
        set_surface_charge = 0.0_dp
    end if

    ! Screen electric field outside dielectric
    ! sigma_function = -epsilon_high/interface_location

    ! Screen electric field inside dielectric
    ! sigma_function = 1/(1 - interface_location)

    ! Equal electric field on left and right
    ! sigma_function = 1 - epsilon_high

    ! No surface charge
    ! set_initial_charge = 1.0e12_dp
  end function set_surface_charge


end module m_user
