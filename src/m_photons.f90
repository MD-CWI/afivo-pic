#include "../afivo/src/cpp_macros.h"
!> This module contains the functionality for photon interaction (ionization and
!> electron emission from dielectric surface)
!>
!> Authors: Xiaoran Li, Jannis Teunissen
!>
!> Note: a prefix pfi_ is used for variables and functions that compute
!> photoionization proportional to ionization.
module m_photons
  use m_globals
  use m_af_all
  use m_particle_core
  use m_cross_sec
  use m_units_constants
  use m_domain

  implicit none
  private

  !> Whether photoionization is enabled
  logical, public, protected :: photoionization_enabled = .false.
  !> Whether photon-induced secondary emission from surfaces is enabled
  logical, public, protected :: photoemission_enabled   = .false.

  !> Whether to include photoionization from collisions
  logical :: photoemission_from_collisions   = .false.
  !> Whether to include photoemission from collisions
  logical :: photoionization_from_collisions = .false.

  !> Whether photoionization from (or proportional to) ionization is enabled, as
  !> is for example used in in Zheleznyak's model
  logical  :: pfi_enabled                            = .false.
  logical  :: quenching_enabled                            = .true.
  !> Relative photoemission probability (compared to photoionization), if both
  !> are assumed proportional to ionization
  real(dp) :: pfi_relative_photoemission_probability = 1.0e-2_dp
  !> Quenching factor for photoionization from ionization
  real(dp) :: pfi_quench_factor                      = 1.0_dp
  !> Minimal inverse absorption length for photoionization from ionization
  real(dp) :: pfi_min_inv_abs_len
  !> Maximal inverse absorption length for photoionization from ionization
  real(dp) :: pfi_max_inv_abs_len

  !> Tabulated field strengths for photoionization from ionization
  real(dp), allocatable :: pfi_photo_eff_table1(:)
  !> Tabulated photoionization efficiencies for photoionization from ionization
  real(dp), allocatable :: pfi_photo_eff_table2(:)

  !> Ignore photoionization below this coordinate
  real(dp) :: photoionization_rmin(NDIM) = -1e10_dp
  !> Ignore photoionization above this coordinate
  real(dp) :: photoionization_rmax(NDIM) = 1e10_dp

  !> Pressure of the gas that is absorbing photons
  real(dp) :: photons_absorbing_gas_pressure = 0.0_dp

  !> Quenching pressure for photoionization (bar)
  real(dp) :: photon_quenching_pressure = 40e-3_dp

  ! Public methods
  public  :: photons_initialize
  public  :: photons_photoionization
  public  :: photons_photoemission

contains

  subroutine photons_initialize(cfg)
    use m_gas
    use m_config
    type(CFG_t), intent(inout) :: cfg
    real(dp)                   :: frac_gas, absorp_inv_lengths(2)
    real(dp)                   :: temp_vec(2), dummy_vec(0)
    integer                    :: t_size, t_size_2
    character(CFG_name_len)    :: model = "none"

    call CFG_add_get(cfg, "photon%model", model, &
         "The model that is used for photoionization")
    call CFG_add_get(cfg, "photon%photoemission_enabled", &
         photoemission_enabled, &
         "Whether photoemission is used")
    call CFG_add_get(cfg, "photon%quenching_enabled", &
         quenching_enabled, &
         "Whether quenching is used")
    call CFG_add_get(cfg, "photon%relative_photoemission_probability", &
         pfi_relative_photoemission_probability, &
         "The relative probability of photoemission compared to photoionization "&
         "(when both are proportional to ionization)")
    call CFG_add_get(cfg, "photon%rmin", photoionization_rmin, &
         "Ignore photoionization below this coordinate")
    call CFG_add_get(cfg, "photon%rmax", photoionization_rmax, &
         "Ignore photoionization above this coordinate")

    call CFG_add_get(cfg, "photon%photoionization_from_collisions", &
         photoionization_from_collisions, &
         "Whether to include photoionization from collisions")
    call CFG_add_get(cfg, "photon%photoemission_from_collisions", &
         photoemission_from_collisions, &
         "Whether to include photoemission from collisions")

    call CFG_add_get(cfg, "photon%quenching_pressure", &
         photon_quenching_pressure, &
         "Quenching pressure for photoionization")

    pfi_photo_eff_table1 = dummy_vec
    pfi_photo_eff_table2 = dummy_vec
    absorp_inv_lengths = [0.0_dp, 0.0_dp]

    ! Initialize parameters and pointers according to the selected model
    select case (model)
    case("Zheleznyak")
       ! Zheleznyak photoionization model for air. See "Photoionization of
       ! nitrogen and oxygen mixtures by radiation from a gas discharge" by
       ! Zheleznyak et al., 1982
       photoionization_enabled = .true.
       pfi_enabled = .true.
       !TODO change to TD
       pfi_photo_eff_table1 = [0.0D0, 100D0, 160D0, 300D0, 600D0, 1500D0]
       pfi_photo_eff_table2 = [0.0D0, 0.05D0, 0.12D0, 0.08D0, 0.06D0, 0.04D0]
       absorp_inv_lengths = [3.5D0 / UC_torr_to_bar, 200D0 / UC_torr_to_bar]

       frac_gas = GAS_get_fraction("O2")
       if (frac_gas <= epsilon(1.0_dp)) &
            error stop "Zheleznyak model requires oxygen"
       photons_absorbing_gas_pressure = frac_gas * GAS_pressure
       pfi_quench_factor = photon_quenching_pressure / &
            (GAS_pressure + photon_quenching_pressure)

       pfi_enabled = .true.
    case("CO2-experimental")
       ! CO2 photoionization model, see "Photoionization produced by
       ! low-current discharges in O2, air, N2 and > CO2", Pancheshnyi, 2015
       photoionization_enabled = .true.
       pfi_enabled = .true.
       pfi_photo_eff_table1 = [200D0, 259D0, 334D0, 534D0, 1000D0]
       pfi_photo_eff_table2 = [0.0D0, 0.6D-4, 1.2D-4, 2.8D-4, 4.8D-4]
       absorp_inv_lengths = [34D0 / UC_torr_to_bar, 220D0 / UC_torr_to_bar]

       frac_gas = GAS_get_fraction("CO2")
       if (frac_gas <= epsilon(1.0_dp)) &
            error stop "CO2-experimental model requires CO2"
       photons_absorbing_gas_pressure = frac_gas * GAS_pressure
       if (quenching_enabled) then
         pfi_quench_factor = photon_quenching_pressure / &
            (GAS_pressure + photon_quenching_pressure)
       else
         pfi_quench_factor = 1
       end if
    case ("none")
       photoionization_enabled = .false.
       photoemission_enabled = .false.
    case default
       error stop "Unrecognized photon%model"
    end select

    call CFG_add(cfg, "photon%efield_table", pfi_photo_eff_table1, &
         "Tabulated electric fields (for the photo-efficiency)", &
         .true.)
    call CFG_add(cfg, "photon%efficiency_table", pfi_photo_eff_table2, &
         "Tabulated photoionization efficiencies (one value means constant)", &
         .true.)
    call CFG_add(cfg, "photon%absorp_inv_lengths", absorp_inv_lengths, &
         "The inverse min/max absorption length, will be scaled by gas density")

    ! Check if arrays have the same size
    call CFG_get_size(cfg, "photon%efficiency_table", t_size)
    call CFG_get_size(cfg, "photon%efield_table", t_size_2)
    if (t_size_2 /= t_size) &
         error stop "size(photon%efield_table) /= size(photon%efficiency_table)"

    deallocate(pfi_photo_eff_table1, pfi_photo_eff_table2)
    allocate(pfi_photo_eff_table1(t_size))
    allocate(pfi_photo_eff_table2(t_size))
    call CFG_get(cfg, "photon%efield_table", pfi_photo_eff_table1)
    call CFG_get(cfg, "photon%efficiency_table", pfi_photo_eff_table2)

    call CFG_get(cfg, "photon%absorp_inv_lengths", temp_vec)
    pfi_min_inv_abs_len = temp_vec(1) * photons_absorbing_gas_pressure
    pfi_max_inv_abs_len = temp_vec(2) * photons_absorbing_gas_pressure

  end subroutine photons_initialize

  !> Compute photoionization due to a list of particle events
  subroutine photons_photoionization(tree, pc, n_photons)
    use omp_lib, only: omp_get_max_threads, omp_get_thread_num
    type(af_t), intent(inout) :: tree
    type(PC_t), intent(inout) :: pc
    integer, intent(out)      :: n_photons
    type(prng_t)              :: prng
    integer                   :: n, m, n_uv, tid, cix
    integer                   :: n_part_before
    real(dp)                  :: x_start(3), x_stop(3)
    real(dp)                  :: en_frac, lambda, mean_production

    call prng%init_parallel(omp_get_max_threads(), GL_rng)

    n_part_before = pc%n_part

    !$omp parallel private(n, n_uv, x_start, x_stop, m, tid, mean_production, &
    !$omp& en_frac, cix, lambda)
    tid = omp_get_thread_num() + 1
    !$omp do
    do n = 1, pc%n_events
       mean_production = mean_photoionization_production(pc%event_list(n))

       if (mean_production > 0) then
          n_uv = prng%rngs(tid)%poisson(mean_production)

          do m = 1, n_uv
             x_start = pc%event_list(n)%part%x

             if (pc%event_list(n)%ctype == CS_ionize_t) then
                en_frac = prng%rngs(tid)%unif_01()
                lambda = pfi_photoionization_lambda(en_frac)
             else
                cix = pc%event_list(n)%cix
                lambda = pc%colls(cix)%c3 * photons_absorbing_gas_pressure
             end if

             x_stop = get_x_stop(x_start, lambda, prng%rngs(tid))
             x_stop(1:NDIM) = get_coordinates_x(x_stop)

             if (is_in_gas(tree, x_stop) .and. .not. &
                  ignore_photon(x_stop)) then
                call single_photoionization_event(tree, pc, x_stop)
             end if
          end do
       end if
    end do
    !$omp end do
    !$omp end parallel

    n_photons = pc%n_part - n_part_before
    call prng%update_seed(GL_rng)
  end subroutine photons_photoionization

  !> Compute photoemission due to a list of particle events
  subroutine photons_photoemission(tree, pc)
    use omp_lib, only: omp_get_max_threads, omp_get_thread_num
    type(af_t), intent(inout) :: tree
    type(PC_t), intent(inout) :: pc
    logical                   :: on_surface
    integer                   :: n, m, n_uv, tid, cix
    real(dp)                  :: x_start(3), x_stop(3)
    real(dp)                  :: en_frac, lambda, mean_production
    type(prng_t)              :: prng

    call prng%init_parallel(omp_get_max_threads(), GL_rng)

    !$omp parallel private(n, n_uv, x_start, x_stop, on_surface, m, tid, &
    !$omp& mean_production, en_frac, cix)
    tid = omp_get_thread_num() + 1
    !$omp do
    do n = 1, pc%n_events
       mean_production = mean_photoemission_production(pc%event_list(n))

       if (mean_production > 0) then
          n_uv = prng%rngs(tid)%poisson(mean_production)

          do m = 1, n_uv
             x_start = pc%event_list(n)%part%x

             if (pc%event_list(n)%ctype == CS_ionize_t) then
                en_frac = prng%rngs(tid)%unif_01()
                lambda = pfi_photoionization_lambda(en_frac)
             else
                cix = pc%event_list(n)%cix
                lambda = pc%colls(cix)%c3 * photons_absorbing_gas_pressure
             end if

             x_stop = get_x_stop(x_start, lambda, prng%rngs(tid))
             x_stop(1:NDIM) = get_coordinates_x(x_stop)

             call photon_diel_absorbtion(tree, x_start, x_stop, on_surface)

             if (on_surface) then
                !$omp critical
                call single_photoemission_event(tree, pc, particle_min_weight, &
                     x_start, x_stop)
                !$omp end critical
             end if
          end do
       end if
    end do
    !$omp end do
    !$omp end parallel
    call prng%update_seed(GL_rng)
  end subroutine photons_photoemission

  !> Calculate the mean number of photo-ionizations for an event
  real(dp) function mean_photoionization_production(event)
    use m_lookup_table
    type(PC_event_t), intent(in) :: event

    mean_photoionization_production = 0.0_dp

    if (pfi_enabled .and. event%ctype == CS_ionize_t) then
       mean_photoionization_production = pfi_photoionization_efficiency(event) * &
            pfi_quench_factor * event%part%w / particle_min_weight
    else if (photoionization_from_collisions .and. &
         event%ctype == CS_photoemission_t) then
       ! The coefficient c2 should hold the photoionization probability
       mean_photoionization_production = pc%colls(event%cix)%c2 * &
            event%part%w / particle_min_weight
    end if

  end function mean_photoionization_production

  !> Calculate the mean number of generated photo-emissions for an event
  real(dp) function mean_photoemission_production(event)
    use m_lookup_table
    type(PC_event_t), intent(in) :: event

    mean_photoemission_production = 0.0_dp

    if (pfi_enabled .and. event%ctype == CS_ionize_t) then
       ! Same number of photons as for photoionization, but we multiply with
       ! pfi_relative_photoemission_probability
       mean_photoemission_production = pfi_photoionization_efficiency(event) * &
            pfi_quench_factor * event%part%w / particle_min_weight * &
            pfi_relative_photoemission_probability
    else if (photoemission_from_collisions .and. &
         event%ctype == CS_photoemission_t) then
       ! The coefficient c1 should hold the photoemission probability
       mean_photoemission_production = pc%colls(event%cix)%c1 * &
            event%part%w / particle_min_weight
    end if

  end function mean_photoemission_production

  !> Determine photoionization efficiency coefficient (per ionization) for the
  !> particles' electric field strength
  real(dp) function pfi_photoionization_efficiency(event)
    use m_lookup_table
    use m_gas
    type(PC_event_t), intent(in) :: event
    real(dp)                     :: reduced_fld

    if (size(pfi_photo_eff_table2) == 1) then
       ! If the efficiency table contains a single element, use it as a constant
       pfi_photoionization_efficiency = pfi_photo_eff_table2(1)
    else
      reduced_fld = norm2(event%part%a / UC_elec_q_over_m) / GAS_number_dens / 1e-21
       call LT_lin_interp_list(pfi_photo_eff_table1, &
            pfi_photo_eff_table2, reduced_fld, pfi_photoionization_efficiency)
    end if
  end function pfi_photoionization_efficiency

  ! Calculate the coordinates of photon absorption according to Zheleznyak model
  function get_x_stop(x_start, lambda, rng) result(x_stop)
    real(dp), intent(in)       :: x_start(3) ! Start coordinate
    real(dp), intent(in)       :: lambda ! Inverse absorption length
    type(rng_t), intent(inout) :: rng
    real(dp)                   :: x_stop(3)
    real(dp)                   :: fly_len

    if (lambda > 0) then
       fly_len = -log(1.0_dp - rng%unif_01()) / lambda
    else
       ! No absorption, let the photon fly out of the domain
       fly_len = norm2(domain_len)
    end if

    x_stop  = x_start + rng%sphere(fly_len)
  end function get_x_stop

  ! Returns the inverse mean free path for a photon according to Zheleznyak model
  real(dp) function pfi_photoionization_lambda(en_frac)
    real(dp), intent(in) :: en_frac
    pfi_photoionization_lambda = pfi_min_inv_abs_len * &
         (pfi_max_inv_abs_len/pfi_min_inv_abs_len)**en_frac
  end function pfi_photoionization_lambda

  ! Generate photo-emitted electrons on the surface of a dielectric
  subroutine single_photoemission_event(tree, pc, photon_w, x_gas, x_diel)
    type(af_t), intent(inout)        :: tree
    type(PC_t), intent(inout)        :: pc
    type(PC_part_t)                  :: new_part
    real(dp), intent(in)             :: photon_w, x_gas(3), x_diel(3)

    new_part%x(:) = x_gas
    new_part%v(:) = 0.0_dp
    new_part%a(:) = pc%accel_function(new_part)

    new_part%w    = photon_w
    new_part%id   = af_get_id_at(tree, x_gas(1:NDIM))

    call pc%add_part(new_part)
    ! Update charge density on the dielectric
    call surface_charge_to_particle(tree, new_part, i_surf_elec)
  end subroutine single_photoemission_event

  ! Input: a particle that is ejected from the dielectric.
  ! The surface charge is altered by the charge leaving
  subroutine surface_charge_to_particle(tree, my_part, i_surf)
    use m_units_constants

    type(af_t), intent(in)      :: tree
    type(PC_part_t), intent(in) :: my_part
    integer, intent(in)         :: i_surf !< Surface variable
    integer                     :: ix_surf, ix_cell(NDIM-1)

    call surface_get_surface_cell(tree, diel, my_part%x(1:NDIM), &
         ix_surf, ix_cell)

    ! Update the charge in the surface cell
#if NDIM == 2
    diel%surfaces(ix_surf)%sd(ix_cell(1), i_surf) = &
         diel%surfaces(ix_surf)%sd(ix_cell(1), i_surf) - &
         my_part%w / diel%surfaces(ix_surf)%dr(1)
#elif NDIM == 3
    diel%surfaces(ix_surf)%sd(ix_cell(1), ix_cell(2), i_surf) = &
         diel%surfaces(ix_surf)%sd(ix_cell(1), ix_cell(2), i_surf) - &
         my_part%w / product(diel%surfaces(ix_surf)%dr)
#endif
  end subroutine surface_charge_to_particle

  subroutine single_photoionization_event(tree, pc, x_stop)
    type(af_t), intent(inout) :: tree !
    type(PC_t), intent(inout) :: pc
    real(dp), intent(in)      :: x_stop(3)
    type(PC_part_t)           :: new_part

    ! Create new particle
    new_part%x = 0
    new_part%x(1:NDIM) = x_stop(1:NDIM)
    new_part%v = 0
    new_part%a = pc%accel_function(new_part)
    new_part%w = particle_min_weight
    new_part%id = af_get_id_at(tree, x_stop(1:NDIM))

    ! It is key that the photons are added at the end of the list
    ! Note: OMP critical is present in add_part
    call pc%add_part(new_part)
  end subroutine single_photoionization_event

  ! Check if the photon is absorbed in the gas
  logical function is_in_gas(tree, x_stop)
    type(af_t), intent(in)         :: tree
    real(dp), intent(in)           :: x_stop(NDIM)
    real(dp)                       :: m_eps(1)
    logical                        :: success

    if (GL_use_dielectric) then
       m_eps = af_interp0(tree, x_stop, [i_eps], success)
       if (m_eps(1) > 1.0_dp .or. .not. success) then
          is_in_gas     = .false.
       else
          is_in_gas     = .true.
       end if
    else
       if (any(x_stop < 0.0_dp .or. x_stop > domain_len)) then
          is_in_gas     = .false.
       else
          is_in_gas     = .true.
       end if
    end if
  end function is_in_gas

  ! Check if the photon is in the ignore area
  logical function ignore_photon(x_stop)
    real(dp), intent(in) :: x_stop(NDIM)
    ignore_photon = any(x_stop < photoionization_rmin) .or. &
         any(x_stop > photoionization_rmax)
  end function ignore_photon

  ! Determine if an emitted photon is absorbed by the gas or dielectric surface
  ! Return coordinates of gas/diel points nearest to intersection
  subroutine photon_diel_absorbtion(tree, x_start, x_stop, on_surface)
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
       on_surface = success
    else ! The photon is absorbed by the gas
       on_surface = .false.
    end if
  end subroutine photon_diel_absorbtion

  ! given start (in gas) and stop (not in gas), this method finds the possible transition.
  ! x_start and x_stop will be moved to corresponding points
  subroutine bisect_line(tree, x_start, x_stop, i_eps)
    type(af_t), intent(in)  :: tree
    real(dp), intent(inout) :: x_start(NDIM), x_stop(NDIM) !< Coordinates of interval
    real(dp)                :: distance
    real(dp)                :: x_mid(NDIM), m_eps(1) ! middle point variables
    integer                 :: n, n_steps
    integer, intent(in)     :: i_eps
    logical                 :: success

    distance = norm2(x_start - x_stop)
    ! Determine the number of steps to ensure a minimum error smaller than the cell size
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
