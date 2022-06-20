#include "../afivo/src/cpp_macros.h"
module m_photons
! This module contains the functionality for
! photon interaction (ionization and electron emission from dielectric surface)
use m_globals
use m_af_all
use m_particle_core
use m_cross_sec
use m_units_constants

use m_time_step
use m_domain

implicit none

private

logical, public         :: photoi_enabled
logical, public         :: photoe_enabled
real(dp)                :: photoe_probability = 1.0e-2_dp
character(CFG_name_len) :: model

real(dp)              :: pi_quench_fac
real(dp)              :: pi_min_inv_abs_len, pi_max_inv_abs_len
real(dp), allocatable :: pi_photo_eff_table1(:), pi_photo_eff_table2(:)

!> Ignore photoionization below this coordinate
real(dp)              :: pi_photon_rmin(NDIM) = -1e10_dp
!> Ignore photoionization above this coordinate
real(dp)              :: pi_photon_rmax(NDIM) = 1e10_dp

! the proportional factor for photoionization
real(dp)              :: pi_eff = 0.075_dp

procedure(photoi), pointer :: photoionization => null()
procedure(photoe), pointer :: photoemission => null()

! Public methods
public  :: photons_initialize
public  :: photoionization
public  :: photoemission

abstract interface
  subroutine photoi(tree, pc, n_photons)
    import
    type(af_t), intent(inout)     :: tree
    type(PC_t), intent(inout)     :: pc
    integer, intent(out)          :: n_photons
  end subroutine photoi

  subroutine photoe(tree, pc)
    import
    type(af_t), intent(inout)     :: tree
    type(PC_t), intent(inout)     :: pc
  end subroutine photoe
end interface

contains

  subroutine photons_initialize(cfg)
    use m_gas
    use m_config
    type(CFG_t), intent(inout) :: cfg

    call CFG_add_get(cfg, "photon%model", model, &
         "The model that is used for photon interactions")

    call CFG_add_get(cfg, "photon%em_enabled", photoe_enabled, &
         "Whether photoionization is used")
    call CFG_add_get(cfg, "photon%em_probability", photoe_probability, &
        "The probability for a single photoemission event")
    call CFG_add_get(cfg, "photon%ion_enabled", photoi_enabled, &
        "Whether photoionization is used")
    call CFG_add_get(cfg, "photon%rmin", pi_photon_rmin, &
         "Ignore photoionization below this coordinate")
    call CFG_add_get(cfg, "photon%rmax", pi_photon_rmax, &
        "Ignore photoionization above this coordinate")
    call CFG_add_get(cfg, "photon%ion_eff", pi_eff, &
        "the proportional factor for photoionization")

    ! Initialize parameters and pointers according to the selected model
    select case (model)
      case("Zheleznyak")
        call Zheleznyak_initialize(cfg)
        photoionization => photoi_Zheleznyak
        photoemission => photoe_Zheleznyak
      case default
        error stop "Unrecognized photon model"
    end select

  end subroutine photons_initialize

  ! ==== Now the modules for Zheleznyak air model
  subroutine Zheleznyak_initialize(cfg)
    use m_gas
    use m_config
    type(CFG_t), intent(inout) :: cfg
    integer                    :: t_size, t_size_2
    real(dp)                   :: frac_O2, temp_vec(2)

    ! Photoionization parameters for AIR
    call CFG_add(cfg, "photon%efield_table", &
        [0.0D0, 0.25D7, 0.4D7, 0.75D7, 1.5D7, 3.75D7], &
        "Tabulated values of the electric field (for the photo-efficiency)")
    call CFG_add(cfg, "photon%efficiency_table", &
        [0.0D0, 0.05D0, 0.12D0, 0.08D0, 0.06D0, 0.04D0], &
        "The tabulated values of the the photo-efficiency")
    call CFG_add(cfg, "photon%frequencies", [2.925D12, 3.059D12], &
        "The lower/upper bound for the frequency of the photons, not currently used")
    call CFG_add(cfg, "photon%absorp_inv_lengths", &
         [3.5D0 / UC_torr_to_bar, 200D0 / UC_torr_to_bar], &
         "The inverse min/max absorption length, will be scaled by pO2")

    if (photoi_enabled) then
       frac_O2 = GAS_get_fraction("O2")
       if (frac_O2 <= epsilon(1.0_dp)) then
          error stop "There is no oxygen, you should disable photoionzation"
       end if

       call CFG_get(cfg, "photon%absorp_inv_lengths", temp_vec)
       pi_min_inv_abs_len = temp_vec(1) * frac_O2 * GAS_pressure
       pi_max_inv_abs_len = temp_vec(2) * frac_O2 * GAS_pressure

       pi_quench_fac = (40.0D0 * UC_torr_to_bar) / &
            (GAS_pressure + (40.0D0 * UC_torr_to_bar))

       call CFG_get_size(cfg, "photon%efficiency_table", t_size)
       call CFG_get_size(cfg, "photon%efield_table", t_size_2)
       if (t_size_2 /= t_size) then
          print *, "size(photoi_efield_table) /= size(photoi_efficiency_table)"
          stop
       end if

       allocate(pi_photo_eff_table1(t_size))
       allocate(pi_photo_eff_table2(t_size))
       call CFG_get(cfg, "photon%efield_table", pi_photo_eff_table1)
       call CFG_get(cfg, "photon%efficiency_table", pi_photo_eff_table2)
    end if
  end subroutine Zheleznyak_initialize

  ! Perform photon generation and ionization according to the Zheleznyak model
  subroutine photoi_Zheleznyak(tree, pc, n_photons)
    use omp_lib, only: omp_get_max_threads, omp_get_thread_num
    type(af_t), intent(inout)     :: tree
    type(PC_t), intent(inout)     :: pc
    integer, intent(out)          :: n_photons
    type(prng_t)                  :: prng

    integer  :: n, m, n_uv, tid
    integer  :: n_part_before
    real(dp) :: x_start(3), x_stop(3)

    call prng%init_parallel(omp_get_max_threads(), GL_rng)

    n_part_before = pc%n_part

    !$omp parallel private(n, n_uv, x_start, x_stop, m, tid)
    tid = omp_get_thread_num() + 1
    !$omp do
    do n = 1, pc%n_events
       if (pc%event_list(n)%ctype == CS_ionize_t) then
          n_uv = prng%rngs(tid)%poisson(get_mean_n_photons(pc%event_list(n)%part))

          do m = 1, n_uv
             x_start = pc%event_list(n)%part%x
             x_stop  = get_x_stop(x_start, prng%rngs(tid))
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
  end subroutine photoi_Zheleznyak

  ! Perform photoemission based on the Zheleznyak model for air
  subroutine photoe_Zheleznyak(tree, pc)
    use omp_lib, only: omp_get_max_threads, omp_get_thread_num
    type(af_t), intent(inout)     :: tree
    type(PC_t), intent(inout)     :: pc
    logical                       :: on_surface
    integer                       :: n, m, n_uv, tid
    real(dp)                      :: x_start(3), x_stop(3)
    type(prng_t)                  :: prng

    call prng%init_parallel(omp_get_max_threads(), GL_rng)

    !$omp parallel private(n, n_uv, x_start, x_stop, on_surface, m, tid)
    tid = omp_get_thread_num() + 1
    !$omp do
    do n = 1, pc%n_events
       if (pc%event_list(n)%ctype == CS_ionize_t) then
         n_uv = prng%rngs(tid)%poisson(get_mean_n_photons(pc%event_list(n)%part))

          do m = 1, n_uv
            x_start = pc%event_list(n)%part%x
            x_stop  = get_x_stop(x_start, prng%rngs(tid))
            call photon_diel_absorbtion(tree, x_start, x_stop, on_surface)
            if (on_surface) then
              if ((prng%rngs(tid)%unif_01() < photoe_probability)) then
                !$omp critical
                call single_photoemission_event(tree, pc, particle_min_weight, x_start, x_stop)
                !$omp end critical
              end if
            end if
          end do
      end if
    end do
    !$omp end do
    !$omp end parallel
    call prng%update_seed(GL_rng)
  end subroutine photoe_Zheleznyak

  ! Calculate the mean number of generate photons according to Zheleznyak model
  real(dp) function get_mean_n_photons(part)
    type(PC_part_t), intent(in)  :: part
    real(dp)                     :: fld

    fld         = norm2(part%a / UC_elec_q_over_m)
    get_mean_n_photons = get_photoi_eff(fld) * pi_quench_fac * &
         part%w / particle_min_weight
  end function

  ! Calculate the coordinates of photon absorption according to Zheleznyak model
  function get_x_stop(x_start, rng) result(x_stop)
    type(rng_t), intent(inout)  :: rng
    real(dp), intent(in)        :: x_start(3)
    real(dp)                    :: x_stop(3)
    real(dp)                    :: en_frac, fly_len

    en_frac = rng%unif_01()
    fly_len = -log(1.0_dp - rng%unif_01()) / get_photoi_lambda(en_frac)
    x_stop  = x_start + rng%sphere(fly_len)
  end function

  ! Returns the photo-efficiency (for ionization) coefficient corresponding to an electric
  ! field of strength fld, according to Zheleznyak model
  real(dp) function get_photoi_eff(fld)
    use m_lookup_table
    real(dp), intent(in) :: fld
    get_photoi_eff = pi_eff
    ! call LT_lin_interp_list(pi_photo_eff_table1, &
         ! pi_photo_eff_table2, fld, get_photoi_eff)
  end function get_photoi_eff

  ! Returns the inverse mean free path for a photon according to Zheleznyak model
  real(dp) function get_photoi_lambda(en_frac)
    real(dp), intent(in) :: en_frac
    get_photoi_lambda = pi_min_inv_abs_len * &
         (pi_max_inv_abs_len/pi_min_inv_abs_len)**en_frac
  end function get_photoi_lambda

  ! ==== General modules

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
    ignore_photon = any(x_stop < pi_photon_rmin) .or. any(x_stop > pi_photon_rmax)
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
