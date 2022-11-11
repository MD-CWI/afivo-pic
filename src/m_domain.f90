!> Module for the computational domain
module m_domain
  use m_globals
  use m_particle_core
  use m_af_all

  implicit none
  private

  ! The length of the (square) domain
  real(dp), protected, public :: domain_len(NDIM) = 4e-3_dp

  ! The volume of the domain
  real(dp), protected, public :: domain_volume = -1.0_dp

  ! The coarse grid size (in number of cells)
  integer, protected, public :: coarse_grid_size(NDIM) = -1

  ! The size of the boxes that we use to construct our mesh
  integer, protected, public :: box_size = 8

  public :: domain_init
  public :: outside_check

contains

  subroutine domain_init(cfg)
    use m_config
    type(CFG_t), intent(inout) :: cfg
    real(dp), parameter        :: pi = acos(-1.0_dp)

    call CFG_add_get(cfg, "domain_len", domain_len, &
         "The length of the domain (m)")
    call CFG_add_get(cfg, "box_size", box_size, &
         "The number of grid cells per coordinate in a box")
    call CFG_add_get(cfg, "coarse_grid_size", coarse_grid_size, &
         "The size of the coarse grid (relevant with electrode)")

    if (all(coarse_grid_size == -1)) then
       ! Automatically set coarse grid size
       coarse_grid_size = nint(domain_len/minval(domain_len)) * box_size
    end if

    if (GL_cylindrical) then
       domain_volume = pi * domain_len(1)**2 * domain_len(2)
    else
       domain_volume = product(domain_len)
    end if
  end subroutine domain_init

  integer function outside_check(my_part)
    type(PC_part_t), intent(inout) :: my_part
    real(dp)                       :: x(NDIM)

    x = get_coordinates(my_part)
    outside_check = 0

    if (any(x < 0.0_dp .or. x > domain_len)) then
       outside_check = outside_domain
       return
    end if

    if (GL_use_dielectric) then
       if (is_in_dielectric(my_part)) then
          ! The particle is IN the domain and IN the dielectric
          outside_check = inside_dielectric
          return
       end if
    end if

    if (GL_use_electrode) then
      if (is_in_electrode(my_part)) then
        outside_check = outside_domain
        return
      end if
    end if
  end function outside_check

  logical function is_in_dielectric(my_part)
    type(PC_part_t), intent(inout) :: my_part
    real(dp)                       :: eps(1)
    logical                        :: success

    eps = af_interp0(tree, get_coordinates(my_part), &
         [i_eps], success, my_part%id)
    if (.not. success) then
       print *, "coordinates: ", get_coordinates(my_part)
       error stop "unexpected particle outside domain"
    end if
    is_in_dielectric = (eps(1) > 1)
  end function is_in_dielectric

  logical function is_in_electrode(my_part)
    type(PC_part_t), intent(inout) :: my_part
    real(dp)                       :: lsf(1)
    logical                        :: success

    lsf = af_interp1(tree, get_coordinates(my_part), &
         [i_lsf], success, my_part%id)
    if (.not. success) then
       print *, "coordinates: ", get_coordinates(my_part)
       error stop "unexpected particle outside domain"
    end if
    is_in_electrode = (lsf(1) < 0)
  end function is_in_electrode
end module m_domain
