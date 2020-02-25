!> Module for the computational domain
module m_domain
  use m_globals
  use m_particle_core
  use m_af_all

  implicit none
  private

  ! The length of the (square) domain
  real(dp), protected, public :: domain_len(NDIM) = 4e-3_dp

  ! The coarse grid size (in number of cells)
  integer, protected, public :: coarse_grid_size(NDIM)

  ! The size of the boxes that we use to construct our mesh
  integer, protected, public :: box_size = 8

  public :: domain_init
  public :: outside_check

contains

  subroutine domain_init(cfg)
    use m_config
    type(CFG_t), intent(inout) :: cfg

    call CFG_add_get(cfg, "domain_len", domain_len, &
         "The length of the domain (m)")
    call CFG_add_get(cfg, "box_size", box_size, &
         "The number of grid cells per coordinate in a box")

    coarse_grid_size = nint(domain_len/minval(domain_len)) * box_size

  end subroutine domain_init

  integer function outside_check(my_part)
    type(PC_part_t), intent(inout) :: my_part
    real(dp)                       :: x(NDIM)

    x = my_part%x(1:NDIM)
    outside_check = 0

    if (any(x < 0.0_dp .or. x > domain_len)) then
       outside_check = outside_domain
    else if (GL_use_dielectric) then
       if (is_in_dielectric(my_part)) then
          ! The particle is IN the domain and IN the dielectric
          outside_check = inside_dielectric
       end if
    end if
  end function outside_check

  logical function is_in_dielectric(my_part)
    type(PC_part_t), intent(inout) :: my_part
    real(dp)                       :: eps(1)
    logical                        :: success

    eps = af_interp0(tree, my_part%x(1:NDIM), [i_eps], success, my_part%id)
    if (.not. success) error stop "unexpected particle outside domain"
    is_in_dielectric = (eps(1) > 1)
  end function is_in_dielectric

!! =============================== !!
! These alternatives are created for more functionality, (eventually should be merged)
  integer function outside_check_x(x)
    real(dp), intent(in)  :: x(NDIM)

    outside_check = 0

    if (any(x < 0.0_dp .or. x > domain_len)) then
       outside_check = outside_domain
    else if (GL_use_dielectric) then
       if (is_in_dielectric_x(x)) then
          ! The particle is IN the domain and IN the dielectric
          outside_check = inside_dielectric
       end if
    end if
  end function outside_check_x

  logical function is_in_dielectric_x(x)
    real(dp), intent(in)           :: x(NDIM)
    real(dp)                       :: eps(1)
    logical                        :: success

    eps = af_interp0(tree, x(1:NDIM), [i_eps], success)
    if (.not. success) error stop "unexpected coordinate outside domain"
    is_in_dielectric = (eps(1) > 1)
  end function is_in_dielectric_x

!! ============================== !!

end module m_domain
