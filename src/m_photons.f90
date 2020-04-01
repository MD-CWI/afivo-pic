module m_photons
! This module contains the functionality for
! photon interaction (ionization and electron emission from dielectric surface)
use m_globals
use m_af_all
use m_particle_core
use m_cross_sec
use m_units_constants
use m_particles, only: get_accel

implicit none

private

real(dp)              :: pi_quench_fac
real(dp)              :: pi_min_inv_abs_len, pi_max_inv_abs_len
! real(dp), allocatable :: pi_photo_eff_table1(:), pi_photo_eff_table2(:)

logical, public, protected :: photoi_enabled = .true. !TODO fix this in config file
logical, public, protected :: photoe_enabled = .false. !TODO fix this in config file
real(dp) :: photoe_coefficient = 1.0e-2 !TODO fix this in config file

! procedure(subr_mover), pointer :: photoionization => null()
! procedure(subr_mover), pointer :: photoemission => null()

! Public methods
public  :: photons_initialize
public  :: photoi_Zheleznyak
public  :: photoe_Zheleznyak

contains

  subroutine photons_initialize(cfg)
    use m_gas
    use m_config
    type(CFG_t), intent(inout) :: cfg
    ! integer                    :: t_size, t_size_2
    real(dp)                   :: frac_O2, temp_vec(2)

    ! Photoemission parameters
      ! TODO Include this here as well. Now it is always off!

    ! Photoionization parameters for AIR
    call CFG_add_get(cfg, "photoi%enabled", photoi_enabled, & !TODO This doesnt work
         "Whether photoionization is used")
    !call CFG_add(cfg, "photoi%efield_table", &
    !     [0.0D0, 0.25D7, 0.4D7, 0.75D7, 1.5D7, 3.75D7], &
    !     "Tabulated values of the electric field (for the photo-efficiency)")
    !call CFG_add(cfg, "photoi%efficiency_table", &
    !     [0.0D0, 0.05D0, 0.12D0, 0.08D0, 0.06D0, 0.04D0], &
    !     "The tabulated values of the the photo-efficiency")
    !call CFG_add(cfg, "photoi%frequencies", [2.925D12, 3.059D12], &
    !     "The lower/upper bound for the frequency of the photons, not currently used")
    call CFG_add(cfg, "photoi%absorp_inv_lengths", &
         [3.5D0 / UC_torr_to_bar, 200D0 / UC_torr_to_bar], &
         "The inverse min/max absorption length, will be scaled by pO2")

    if (photoi_enabled .or. photoe_enabled) then
       frac_O2 = GAS_get_fraction("O2")
       if (frac_O2 <= epsilon(1.0_dp)) then
          error stop "There is no oxygen, you should disable photoionzation"
       end if

       call CFG_get(cfg, "photoi%absorp_inv_lengths", temp_vec)
       pi_min_inv_abs_len = temp_vec(1) * frac_O2 * GAS_pressure
       pi_max_inv_abs_len = temp_vec(2) * frac_O2 * GAS_pressure

       pi_quench_fac = (40.0D0 * UC_torr_to_bar) / &
            (GAS_pressure + (40.0D0 * UC_torr_to_bar))

       ! call CFG_get_size(cfg, "photoi%efficiency_table", t_size)
       ! call CFG_get_size(cfg, "photoi%efield_table", t_size_2)
       ! if (t_size_2 /= t_size) then
       !    print *, "size(photoi_efield_table) /= size(photoi_efficiency_table)"
       !    stop
       ! end if

       !allocate(pi_photo_eff_table1(t_size))
       !allocate(pi_photo_eff_table2(t_size))
       !call CFG_get(cfg, "photoi%efield_table", pi_photo_eff_table1)
       !call CFG_get(cfg, "photoi%efficiency_table", pi_photo_eff_table2)
    end if
  end subroutine photons_initialize

  subroutine photoi_Zheleznyak(tree, pc, photo_pos, photo_w, n_photons)
    ! Perform photon generation and ionization according to the Zheleznyak model
    ! use m_domain
    ! use m_particles, only: get_accel

    type(af_t), intent(inout)     :: tree
    type(PC_t), intent(inout)     :: pc
    real(dp), intent(inout)       :: photo_pos(:, :)
    real(dp), intent(inout)       :: photo_w(:)
    integer, intent(out)          :: n_photons

    ! type(PC_part_t) :: new_part
    integer         :: i, n, m, n_uv
    real(dp)        :: x_start(3), x_stop(3)

    i = 0
    do n = 1, pc%n_events
       if (pc%event_list(n)%ctype == CS_ionize_t) then
          n_uv = GL_rng%poisson(get_mean_n_photons(pc%event_list(n)%part))

          do m = 1, n_uv
             x_start = pc%event_list(n)%part%x
             x_stop  = get_x_stop(x_start)
             if (is_in_gas(tree, x_stop)) then
               call single_photoionization_event(tree, pc, i, photo_pos, photo_w, x_stop)
             end if
          end do
       end if
    end do

    n_photons = i
  end subroutine photoi_Zheleznyak

  subroutine photoe_Zheleznyak(tree, pc)
    ! Perform photoemission based on the Zheleznyak model for air
    type(af_t), intent(inout)     :: tree
    type(PC_t), intent(inout)     :: pc
    logical       :: on_surface
    integer       :: n, m, n_uv, n_low_en, n_high_en
    real(dp)      :: x_start(3), x_stop(3)

    do n = 1, pc%n_events
       if (pc%event_list(n)%ctype == CS_ionize_t) then
         n_low_en = GL_rng%poisson(get_mean_n_photons(pc%event_list(n)%part))
         n_high_en = GL_rng%poisson(get_mean_n_photons(pc%event_list(n)%part))
         n_uv = n_low_en + n_high_en ! low-energy photons come from O2 and are generated at the same rate as high energy photons from N2

          do m = 1, n_uv
            x_start = pc%event_list(n)%part%x
            x_stop  = get_x_stop(x_start)

            call photon_diel_absorbtion(tree, x_start, x_stop, on_surface)
            if (on_surface) then
              call single_photoemission_event(tree, pc, particle_min_weight, x_start, x_stop)
            end if
          end do
      end if
    end do
  end subroutine photoe_Zheleznyak

  real(dp) function get_mean_n_photons(part)
    ! Calculate the mean number of generate photons according to Zheleznyak model
    type(PC_part_t), intent(in)  :: part
    real(dp)                     :: fld
    fld         = norm2(part%a / UC_elec_q_over_m)
    get_mean_n_photons = get_photoi_eff(fld) * pi_quench_fac * &
         part%w / particle_min_weight
  end function

  function get_x_stop(x_start) result(x_stop)
    ! Calculate the coordinates of photon absorption according to Zheleznyak model
    real(dp)  :: en_frac, fly_len
    real(dp), intent(in)  :: x_start(3)
    real(dp)              :: x_stop(3)

    en_frac = GL_rng%unif_01()
    fly_len = -log(1.0_dp - GL_rng%unif_01()) / get_photoi_lambda(en_frac)
    x_stop  = x_start + GL_rng%sphere(fly_len)
  end function

  real(dp) function get_photoi_eff(fld)
    ! Returns the photo-efficiency (for ionization) coefficient corresponding to an electric
    ! field of strength fld, according to Zheleznyak model
    use m_lookup_table
    real(dp), intent(in) :: fld
    get_photoi_eff = 0.075_dp
    ! call LT_lin_interp_list(pi_photo_eff_table1, &
         ! pi_photo_eff_table2, fld, get_photoi_eff)
  end function get_photoi_eff

  real(dp) function get_photoi_lambda(en_frac)
    ! Returns the inverse mean free path for a photon according to Zheleznyak model
    real(dp), intent(in) :: en_frac
    get_photoi_lambda = pi_min_inv_abs_len * &
         (pi_max_inv_abs_len/pi_min_inv_abs_len)**en_frac
  end function get_photoi_lambda

! =====

  subroutine single_photoemission_event(tree, pc, photon_w, x_gas, x_diel)
    use m_particles, only: surface_charge_to_particle
    ! Generate photo-emitted electrons on the surface of a dielectric
    type(af_t), intent(inout)        :: tree
    type(PC_t), intent(inout)        :: pc
    type(PC_part_t)                  :: new_part
    real(dp), intent(in)             :: photon_w, x_gas(3), x_diel(3)

    if (GL_rng%unif_01() < photoe_coefficient) then
      ! Create photo-emitted electron
      new_part%x(:) = x_gas
      new_part%v(:) = 0.0_dp
      new_part%a(:) = get_accel(new_part)
      new_part%w    = photon_w
      new_part%id   = af_get_id_at(tree, x_gas(1:NDIM))

      call pc%add_part(new_part)
      ! Update charge density on the dielectric
      call surface_charge_to_particle(tree, new_part, i_surf_elec)
    end if
  end subroutine single_photoemission_event

  subroutine single_photoionization_event(tree, pc, i, photo_pos, photo_w, x_stop)
    use m_particles, only: get_accel

    type(af_t), intent(inout)     :: tree
    type(PC_t), intent(inout)     :: pc
    real(dp), intent(inout)       :: photo_pos(:, :)
    real(dp), intent(inout)       :: photo_w(:)

    type(PC_part_t)         :: new_part
    integer, intent(inout)  :: i
    integer                 :: i_cpy
    real(dp)                :: x_stop(3)

    !$omp critical
    i = i + 1
    i_cpy = i
    !$omp end critical
    if (i_cpy > size(photo_w)) error stop "Too many photons were generated"
    ! Return coordinates and weights for ion production
    photo_pos(:, i_cpy) = x_stop(1:NDIM)
    photo_w(i_cpy)      = particle_min_weight

    ! Create new particle
    new_part%x = 0 !TODO is this neccesary?
    new_part%x(1:NDIM) = x_stop(1:NDIM)
    new_part%v = 0
    new_part%a = get_accel(new_part)
    new_part%w = particle_min_weight
    new_part%id = af_get_id_at(tree, x_stop(1:NDIM))

    call pc%add_part(new_part)
  end subroutine single_photoionization_event
! ==================

  logical function is_in_gas(tree, x_stop)
    ! Check if the photon is absorbed in the gas
    type(af_t), intent(in)         :: tree
    real(dp), intent(in)           :: x_stop(NDIM)
    real(dp)                       :: m_eps(1)
    logical                        :: success

    m_eps = af_interp0(tree, x_stop, [i_eps], success)
    if (m_eps(1) > 1.0_dp .or. .not. success) then
      is_in_gas     = .false.
    else
      is_in_gas     = .true.
    end if
  end function is_in_gas

  subroutine photon_diel_absorbtion(tree, x_start, x_stop, on_surface)
    ! Determine if an emitted photon is absorbed by the gas or dielectric surface
    ! Return coordinates of gas/diel points nearest to intersection
    type(af_t), intent(in)         :: tree
    real(dp), intent(inout)        :: x_start(NDIM), x_stop(NDIM) !< Coordinates of photon event
    real(dp)                       :: m_eps(1)
    logical, intent(inout)         :: on_surface
    logical                        :: success

    m_eps = af_interp0(tree, x_stop, [i_eps], success)
    ! Perform bisection and determine on_surface
    if (m_eps(1) > 1.0_dp .or. .not. success) then
      call bisect_line(tree, x_start, x_stop, i_eps)
      m_eps = af_interp0(tree, x_stop, [i_eps], success)
      if (.not. success) then ! Photon travelled outside the domain
         on_surface = .false.
      else ! Photon is absorbed by the dielectric
         on_surface = .true.
      end if
    else ! The photon is absorbed by the gas
      on_surface = .false.
    end if
  end subroutine photon_diel_absorbtion

  subroutine bisect_line(tree, x_start, x_stop, i_eps)
    ! given start (in gas) and stop (not in gas), this method finds the possible transition.
    ! x_start and x_stop will be moved to corresponding points
    type(af_t), intent(in)  :: tree
    real(dp), intent(inout) :: x_start(NDIM), x_stop(NDIM) !< Coordinates of interval
    real(dp)                :: distance
    real(dp)                :: x_mid(NDIM), m_eps(1) ! middle point variables
    integer                 :: n, n_steps
    integer, intent(in)     :: i_eps
    logical                 :: success

    distance = norm2(x_start - x_stop)
    n_steps = -ceiling(log(af_min_dr(tree)/distance) / log(2.0_dp))

    do n = 1, n_steps
      x_mid = 0.5_dp * (x_start + x_stop)
      m_eps = af_interp0(tree, x_mid, [i_eps], success)
      if (m_eps(1) > 1.0_dp .or. .not. success) then
        x_stop = x_mid ! Move the end to the middle
      else
        x_start = x_mid ! Move the start to the middle
      end if
    end do

  end subroutine bisect_line

end module m_photons
