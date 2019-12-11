module m_domain
  implicit none
  public

  integer, parameter, private :: dp = kind(0.0d0)

  ! The length of the (square) domain
  real(dp), protected :: domain_len = 4e-3_dp

  ! The size of the boxes that we use to construct our mesh
  integer, protected :: box_size = 8

contains

  subroutine domain_init(cfg)
    use m_config
    type(CFG_t), intent(inout) :: cfg

    call CFG_add_get(cfg, "domain_len", domain_len, &
         "The length of the domain (m)")
    call CFG_add_get(cfg, "box_size", box_size, &
         "The number of grid cells per coordinate in a box")

  end subroutine domain_init

  pure integer function outside_check_pos(x)
    use m_particle_core
    real(dp), intent(in) :: x(NDIM)

    if (any(x(1:NDIM) < 0.0_dp .or. x(1:NDIM) > domain_len)) then
       outside_check_pos = 1
    else
       outside_check_pos = 0
    end if
  end function outside_check_pos

  integer function outside_check(my_part)
    use m_particle_core
    type(PC_part_t), intent(inout) :: my_part

    outside_check = outside_check_pos(my_part%x(1:NDIM))
  end function outside_check

end module m_domain
