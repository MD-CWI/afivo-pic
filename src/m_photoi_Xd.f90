module m_photoi_$Dd
  use m_globals_$Dd

  implicit none
  private

  real(dp)              :: pi_quench_fac
  real(dp)              :: pi_min_inv_abs_len, pi_max_inv_abs_len
  real(dp), allocatable :: pi_photo_eff_table1(:), pi_photo_eff_table2(:)
  type(RNG_t)           :: pi_rng

  logical, public :: photoi_enabled = .false.

  ! Public methods
  public :: pi_initialize
  public :: get_photoionization

contains

  subroutine pi_initialize(cfg)
    use m_gas
    use m_units_constants
    use m_config
    type(CFG_t), intent(inout) :: cfg
    integer                    :: t_size, t_size_2
    real(dp)                   :: frac_O2, temp_vec(2)

    ! Photoionization parameters
    call CFG_add_get(cfg, "photoi%enabled", photoi_enabled, &
         "Whether photoionization is used")
    call CFG_add(cfg, "photoi%efield_table", &
         [0.0D0, 0.25D7, 0.4D7, 0.75D7, 1.5D7, 3.75D7], &
         "Tabulated values of the electric field (for the photo-efficiency)")
    call CFG_add(cfg, "photoi%efficiency_table", &
         [0.0D0, 0.05D0, 0.12D0, 0.08D0, 0.06D0, 0.04D0], &
         "The tabulated values of the the photo-efficiency")
    CALL CFG_add(cfg, "photoi%frequencies", [2.925D12, 3.059D12], &
         "The lower/upper bound for the frequency of the photons, not currently used")
    call CFG_add(cfg, "photoi%absorp_inv_lengths", &
         [3.5D0 / UC_torr_to_bar, 200D0 / UC_torr_to_bar], &
         "The inverse min/max absorption length, will be scaled by pO2")

    frac_O2 = GAS_get_fraction("O2")
    if (frac_O2 <= epsilon(1.0_dp)) then
       error stop "There is no oxygen, you should disable photoionzation"
    end if

    call CFG_get(cfg, "photoi%absorp_inv_lengths", temp_vec)
    pi_min_inv_abs_len = temp_vec(1) * frac_O2 * GAS_pressure
    pi_max_inv_abs_len = temp_vec(2) * frac_O2 * GAS_pressure

    ! print *, "Max abs. length photoi.", 1.0d3 / pi_min_inv_abs_len, "mm"
    ! print *, "Min abs. length photoi.", 1.0d3 / pi_max_inv_abs_len, "mm"

    pi_quench_fac = (30.0D0 * UC_torr_to_bar) / &
         (GAS_pressure + (30.0D0 * UC_torr_to_bar))

    call CFG_get_size(cfg, "photoi%efficiency_table", t_size)
    call CFG_get_size(cfg, "photoi%efield_table", t_size_2)
    if (t_size_2 /= t_size) then
       print *, "size(photoi_efield_table) /= size(photoi_efficiency_table)"
       stop
    end if

    allocate(pi_photo_eff_table1(t_size))
    allocate(pi_photo_eff_table2(t_size))
    call CFG_get(cfg, "photoi%efield_table", pi_photo_eff_table1)
    call CFG_get(cfg, "photoi%efficiency_table", pi_photo_eff_table2)
  end subroutine pi_initialize

  subroutine get_photoionization(events, photo_pos, photo_w, n_photons)
    use m_cross_sec
    use m_particle_core
    use m_units_constants
    use m_domain_$Dd

    type(PC_events_t), intent(in) :: events
    real(dp), intent(inout)       :: photo_pos(:, :)
    real(dp), intent(inout)       :: photo_w(:)
    integer, intent(out)          :: n_photons

    integer  :: i, i_cpy, n, m, n_events, n_uv
    real(dp) :: en_frac, fld, fly_len, mean_gammas, x(3)

    n_events = events%n_stored

    i = 0
    do n = 1, n_events
       if (events%list(n)%ctype == CS_ionize_t) then

          fld         = norm2(events%list(n)%part%a / UC_elec_q_over_m)
          mean_gammas = get_photoi_eff(fld) * pi_quench_fac * &
               events%list(n)%part%w / particle_min_weight
          n_uv        = ST_rng%poisson(mean_gammas)

          do m = 1, n_uv
             en_frac = ST_rng%unif_01()
             fly_len = -log(1.0_dp - ST_rng%unif_01()) / get_photoi_lambda(en_frac)
             x       = events%list(n)%part%x + ST_rng%sphere(fly_len)

             if (outside_check_pos(x) == 0) then
                !$omp critical
                i = i + 1
                i_cpy = i
                !$omp end critical
                if (i_cpy > size(photo_w)) error stop "Too many photons were generated"
                photo_pos(:, i_cpy) = x(1:$D)
                photo_w(i_cpy)      = particle_min_weight
             end if
          end do
       end if
    end do

    n_photons = i

  end subroutine get_photoionization

  ! Returns the photo-efficiency coefficient corresponding to an electric
  ! field of strength fld
  real(dp) function get_photoi_eff(fld)
    use m_lookup_table
    real(dp), intent(in) :: fld
    call LT_lin_interp_list(pi_photo_eff_table1, &
         pi_photo_eff_table2, fld, get_photoi_eff)
  end function get_photoi_eff

  ! Returns the inverse mean free path for a photon.
  real(dp) function get_photoi_lambda(en_frac)
    real(dp), intent(in) :: en_frac
    get_photoi_lambda = pi_min_inv_abs_len * &
         (pi_max_inv_abs_len/pi_min_inv_abs_len)**en_frac
  end function get_photoi_lambda

end module m_photoi_$Dd
