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
real(dp)              :: cs_O2(2), cs_CH4, frac_O2, frac_CH4, frac_N2
real(dp)              :: pi_min_inv_abs_len, pi_max_inv_abs_len
real(dp), allocatable :: pi_photo_eff_table1(:), pi_photo_eff_table2(:)

!Default reaction rates for excited Ar and Ar2 pools
real(dp)  :: k_Ar_decay_rad  = 0.0_dp
real(dp)  :: k_Ar_quench     = 0.0_dp
real(dp)  :: k_Ar2_prod_rate = 3.0e-46_dp
real(dp)  :: k_Ar2_decay_rad = 6.0e7_dp

integer, protected, public  :: i_Ar_pool = -1 !Index of cell-centered pool of excited argon species
integer, protected, public  :: i_Ar2_pool = -1 !Index of cell-centered Argon excimer

procedure(photoi), pointer :: photoionization => null()
procedure(photoe), pointer :: photoemission => null()

! Public methods
public  :: photons_initialize
public  :: photoionization
public  :: photoemission

abstract interface
  subroutine photoi(tree, pc, photo_pos, photo_w, n_photons)
    import
    type(af_t), intent(inout)     :: tree
    type(PC_t), intent(inout)     :: pc
    real(dp), intent(inout)       :: photo_pos(:, :)
    real(dp), intent(inout)       :: photo_w(:)
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
        "Whether photoionization is used")
        call CFG_add_get(cfg, "photon%ion_enabled", photoi_enabled, &
        "Whether photoionization is used")

    ! Initialize parameters and pointers according to the selected model
    select case (model)
    case("Zheleznyak")
        call Zheleznyak_initialize(cfg)
        photoionization => photoi_Zheleznyak
        photoemission => photoe_Zheleznyak
      case("Argon")
        call Argon_initialize(cfg)
        photoionization => null() ! No photoionization is considered
        photoemission => photoe_Argon
      case default
        if (photoi_enabled .or. photoe_enabled) then
          error stop "Unrecognized photon model"
        else
          photoionization => null()
          photoemission => null()
        end if
    end select

  end subroutine photons_initialize

! ===== Model photon interactions in argon

  subroutine Argon_initialize(cfg)
    use m_config
    type(CFG_t), intent(inout) :: cfg

    call CFG_add_get(cfg, "photon%k_Ar_decay_rad", k_Ar_decay_rad, &
         "Rate of radiative decay of excited Ar")
    call CFG_add_get(cfg, "photon%k_Ar_quench", k_Ar_quench, &
         "Rate of quenching of excited Ar")
    call CFG_add_get(cfg, "photon%k_Ar2_prod_rate", k_Ar2_prod_rate, &
         "Rate of production of excited Argon2")
    call CFG_add_get(cfg, "photon%k_Ar2_decay_rad", k_Ar2_decay_rad, &
         "Rate of radiative decay of excited Ar2")

    if (GL_use_dielectric .and. photoe_enabled) then
      where(pc%colls(:)%type == CS_excite_t)
        pc%coll_is_event(:) = .true.
      end where

      !Create empty argon pool and excimer
      call af_add_cc_variable(tree, "Ar_pool", .true., ix=i_Ar_pool)
      call af_set_cc_methods(tree, i_Ar_pool, af_bc_neumann_zero, &
           prolong=af_prolong_limit)

      call af_add_cc_variable(tree, "Ar2_pool", .true., ix=i_Ar2_pool)
      call af_set_cc_methods(tree, i_Ar2_pool, af_bc_neumann_zero, &
                prolong=af_prolong_limit)
    end if

  end subroutine Argon_initialize

  subroutine photoe_Argon(tree, pc)
    ! Perform photoemission due to pool of excited argon species
    type(af_t), intent(inout)     :: tree
    type(PC_t), intent(inout)     :: pc

    real(dp), allocatable, save :: coords(:, :)
    real(dp), allocatable, save :: weights(:)
    integer, allocatable, save  :: id_guess(:)

    integer       :: n, jj

    if (.not. allocated(weights)) then
       allocate(coords(NDIM, pc%n_events))
       allocate(weights(pc%n_events))
       allocate(id_guess(pc%n_events))
    else if (size(weights) < pc%n_events) then
       deallocate(coords)
       deallocate(weights)
       deallocate(id_guess)
       allocate(coords(NDIM, pc%n_events))
       allocate(weights(pc%n_events))
       allocate(id_guess(pc%n_events))
    end if

    call Ar2_radiative_decay(tree, pc) ! Do photoemission events (updates the Ar2 pool)
    call af_loop_box(tree, perform_argon_reactions, .true.)

    ! Calculate production to the pool of excited Argon species
    jj = 0
    do n = 1, pc%n_events
       if (pc%event_list(n)%ctype == CS_excite_t) then
         jj = jj + 1
         coords(:, jj) = pc%event_list(n)%part%x(1:NDIM)
         weights(jj) = pc%event_list(n)%part%w
         id_guess(jj) = pc%event_list(n)%part%id
         ! end if
       end if
     end do
     call af_particles_to_grid(tree, i_Ar_pool, coords(:, 1:jj), &
          weights(1:jj), jj, 0, id_guess(1:jj)) ! Use zeroth order interpolation for simplicity
end subroutine photoe_Argon

subroutine Ar2_radiative_decay(tree, pc)
  ! Routine that performs Ar2_exc radiative decay for every cell in a box
  ! This routine also performs photoemission

  ! TODO Code from af_loop_box is copied (because of particle creation and non-local effects). Thats why this code is nasty
  type(af_t), intent(inout)  :: tree
  type(PC_t), intent(inout)  :: pc

  integer    :: ii, jj, nn
  integer    :: lvl, i, id
  integer    :: n_uv
  real(dp)   :: mean_ph
  real(dp)   :: cell_size
  real(dp)   :: x_start(3), x_cc(3), x_stop(3)
  logical    :: on_surface

  if (.not. tree%ready) stop "Ar2_radiative_decay: set_base has not been called"

#if NDIM == 2
  !!$omp parallel private(lvl, i, id, ii, jj, nn)
  do lvl = 1, tree%highest_lvl
    !!$omp do
    do i = 1, size(tree%lvls(lvl)%leaves)
      id = tree%lvls(lvl)%leaves(i)
      cell_size = product(tree%boxes(id)%dr)
      do ii = 1, tree%boxes(id)%n_cell
        do jj = 1, tree%boxes(id)%n_cell
          mean_ph = cell_size * GL_dt * k_Ar2_decay_rad &
            * tree%boxes(id)%cc(ii, jj, i_Ar2_pool) !/ particle_min_weight
          n_uv = GL_rng%poisson(mean_ph)

          x_cc(1:NDIM) = af_r_cc(tree%boxes(id), [ii, jj])
          do nn = 1, n_uv
            x_start = x_cc
            x_stop(1:NDIM) = x_start(1:NDIM) + GL_rng%circle(norm2(domain_len)) ! Always scatter out of the domain
            call photon_diel_absorbtion(tree, x_start, x_stop, on_surface)
            if (on_surface) then
              ! call single_photoemission_event(tree, pc, particle_min_weight, x_start, x_stop)
              call single_photoemission_event(tree, pc, 1.0_dp, x_start, x_stop)
              ! print *, "photoemission event has occurred"
            end if
          end do
        end do
      end do
      !!$omp end do
    end do
    !!$omp end do
  end do
  !!$omp end parallel
#elif NDIM == 3
    error stop
#endif
  end subroutine Ar2_radiative_decay

  subroutine perform_argon_reactions(box)
    use m_gas
    ! Per cell, do explicit Euler to update the Ar-Ar2 pools due to reaction mechanism (excluding Ar* production and Ar2 radiative decay)
    ! type(af_t), intent(inout) :: tree
    type(box_t), intent(inout)  :: box
    integer    :: ii, jj

#if NDIM == 2
    do ii = 1, box%n_cell
      do jj = 1, box%n_cell
        ! Production of Ar2
        box%cc(ii, jj, i_Ar2_pool) = box%cc(ii, jj, i_Ar2_pool) + &
          GL_dt * k_Ar2_prod_rate * GAS_number_dens**2 * box%cc(ii, jj, i_Ar_pool)
        ! decay, quenching and Ar2-prod for Ar
        box%cc(ii, jj, i_Ar_pool) = box%cc(ii, jj, i_Ar_pool) - GL_dt * &
          (k_Ar_decay_rad + k_Ar_quench * GAS_number_dens + k_Ar2_prod_rate * GAS_number_dens**2) * box%cc(ii, jj, i_Ar_pool)
        if (box%cc(ii, jj, i_Ar2_pool) < 0.0_dp .or. box%cc(ii, jj, i_Ar_pool) < 0.0_dp) then
          error stop "Negative density for excited states of Argon or Argon2 found after performing chemical reactions."
        end if
      end do
    end do
#elif NDIM == 3
      error stop
#endif
  end subroutine perform_argon_reactions

  ! ==== Now the modules for Zheleznyak model

  subroutine Zheleznyak_initialize(cfg)
    use m_gas
    use m_config
    type(CFG_t), intent(inout) :: cfg
    integer                    :: t_size, t_size_2
    real(dp)                   :: temp_vec(2)
    real(dp)                   :: pi_quench_effective

    ! Photoionization parameters for AIR (possible mixed with CH4)
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

    cs_O2 = [1.05e-22_dp, 5.9e-21_dp] ! Min/max radiative cross sections for oxygen [m^2]
    cs_CH4 = 3.0e-21_dp ! units of m^2

    if (photoi_enabled .or. photoe_enabled) then
       frac_O2 = GAS_get_fraction("O2")
       frac_N2 = GAS_get_fraction("N2")
       if (frac_O2 <= epsilon(1.0_dp)) then
          error stop "There is no oxygen, you should disable photoionzation"
       end if

       frac_CH4 = GAS_get_fraction("CH4")
       if (frac_CH4 <= epsilon(1.0_dp)) then
         call CFG_get(cfg, "photon%absorp_inv_lengths", temp_vec)
         pi_min_inv_abs_len = temp_vec(1) * frac_O2 * GAS_pressure
         pi_max_inv_abs_len = temp_vec(2) * frac_O2 * GAS_pressure

         pi_quench_fac = (40.0D0 * UC_torr_to_bar) / &
             (GAS_pressure + (40.0D0 * UC_torr_to_bar))
       else
         ! CH4 is present and can absorbs photons
         pi_min_inv_abs_len = GAS_number_dens * (frac_O2 * cs_O2(1) + frac_CH4 * cs_CH4)
         pi_max_inv_abs_len = GAS_number_dens * (frac_O2 * cs_O2(2) + frac_CH4 * cs_CH4)

         pi_quench_effective = (frac_N2 / 40.0_dp + frac_O2 / 3.5_dp + frac_CH4 / 1.8_dp) ** (-1) * UC_torr_to_bar
         pi_quench_fac = pi_quench_effective / (GAS_pressure + pi_quench_effective)
       end if

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
  subroutine photoi_Zheleznyak(tree, pc, photo_pos, photo_w, n_photons)
    use omp_lib, only: omp_get_max_threads, omp_get_thread_num
    type(af_t), intent(inout)     :: tree
    type(PC_t), intent(inout)     :: pc
    real(dp), intent(inout)       :: photo_pos(:, :)
    real(dp), intent(inout)       :: photo_w(:)
    integer, intent(out)          :: n_photons
    type(prng_t)                  :: prng

    integer         :: i, n, m, n_uv, tid
    real(dp)        :: x_start(3), x_stop(3)
    real(dp)        :: en_frac


    call prng%init_parallel(omp_get_max_threads(), GL_rng)

    i = 0
    !$omp parallel private(n, n_uv, x_start, x_stop, m, tid, en_frac)
    tid = omp_get_thread_num() + 1
    !omp do
    do n = 1, pc%n_events
       if (pc%event_list(n)%ctype == CS_ionize_t) then
          n_uv = prng%rngs(tid)%poisson(get_mean_n_photons(pc%event_list(n)%part))

          do m = 1, n_uv
             en_frac = prng%rngs(tid)%unif_01()
             if (frac_CH4 > epsilon(1.0_dp) .and. is_absorbed_by_CH4(en_frac)) cycle
             x_start = pc%event_list(n)%part%x
             x_stop  = get_x_stop(x_start, en_frac, prng%rngs(tid))

             if (is_in_gas(tree, x_stop)) then
               !$omp critical
               call single_photoionization_event(tree, pc, i, photo_pos, photo_w, x_stop)
               !$omp end critical
             end if
          end do
       end if
    end do
    !omp end do
    !$omp end parallel
    n_photons = i
    call prng%update_seed(GL_rng)
  end subroutine photoi_Zheleznyak

  ! Perform photoemission based on the Zheleznyak model
  subroutine photoe_Zheleznyak(tree, pc)
    use omp_lib, only: omp_get_max_threads, omp_get_thread_num
    type(af_t), intent(inout)     :: tree
    type(PC_t), intent(inout)     :: pc
    logical                       :: on_surface
    integer                       :: n, m, n_uv, tid
    real(dp)                      :: x_start(3), x_stop(3)
    real(dp)                      :: en_frac
    type(prng_t)                  :: prng

    call prng%init_parallel(omp_get_max_threads(), GL_rng)

    !$omp parallel private(n, n_uv, x_start, x_stop, on_surface, m, tid)
    tid = omp_get_thread_num() + 1
    !omp do
    do n = 1, pc%n_events
       if (pc%event_list(n)%ctype == CS_ionize_t) then
         n_uv = prng%rngs(tid)%poisson(get_mean_n_photons(pc%event_list(n)%part))

          do m = 1, n_uv ! High&low en photons are handled equally
            en_frac = GL_rng%unif_01()
            x_start = pc%event_list(n)%part%x
            x_stop  = get_x_stop(x_start, en_frac, prng%rngs(tid))

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
    !omp end do
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
  function get_x_stop(x_start, en_frac, rng) result(x_stop)
    type(rng_t), intent(inout)  :: rng
    real(dp), intent(in)        :: en_frac
    real(dp), intent(in)        :: x_start(3)
    real(dp)                    :: x_stop(3)
    real(dp)                    :: fly_len

    fly_len = -log(1.0_dp - rng%unif_01()) / get_photoi_lambda(en_frac)
    x_stop  = x_start + rng%sphere(fly_len)
  end function

  ! Returns the photo-efficiency (for ionization) coefficient corresponding to an electric
  ! field of strength fld, according to Zheleznyak model
  real(dp) function get_photoi_eff(fld)
    use m_lookup_table
    real(dp), intent(in) :: fld
    get_photoi_eff = 0.075_dp
    ! call LT_lin_interp_list(pi_photo_eff_table1, &
         ! pi_photo_eff_table2, fld, get_photoi_eff)
  end function get_photoi_eff

  ! Returns the inverse mean free path for a photon according to Zheleznyak model
  real(dp) function get_photoi_lambda(en_frac)
    real(dp), intent(in) :: en_frac
    get_photoi_lambda = pi_min_inv_abs_len * &
         (pi_max_inv_abs_len/pi_min_inv_abs_len)**en_frac
  end function get_photoi_lambda

  real(dp) function get_cs_O2(en_frac)
    ! Returns the radiative cross section of O2 based on a power-law approximation
    real(dp), intent(in)  ::  en_frac
    get_cs_O2 = (cs_O2(1)**(1-en_frac)) * (cs_O2(2) **en_frac)
  end function

  logical function is_absorbed_by_CH4(en_frac)
    ! Determine if a photon is absorbed by CH4
    real(dp), intent(in)  ::  en_frac
    is_absorbed_by_CH4 = (GL_rng%unif_01() < &
      get_cs_O2(en_frac) * frac_O2 / (get_cs_O2(en_frac) * frac_O2 + cs_CH4 * frac_CH4))
  end function

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

    call dielectric_get_surface_cell(tree, diel, my_part%x(1:NDIM), &
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

  subroutine single_photoionization_event(tree, pc, i, photo_pos, photo_w, x_stop)
    type(af_t), intent(inout)     :: tree!
    type(PC_t), intent(inout)     :: pc
    real(dp), intent(inout)       :: photo_pos(:, :)
    real(dp), intent(inout)       :: photo_w(:)

    type(PC_part_t)         :: new_part
    integer, intent(inout)  :: i
    integer                 :: i_cpy
    real(dp)                :: x_stop(3)

    !!$omp critical
    i = i + 1
    i_cpy = i
    !!$omp end critical
    if (i_cpy > size(photo_w)) error stop "Too many photons were generated"
    ! Return coordinates and weights for ion production
    photo_pos(:, i_cpy) = x_stop(1:NDIM)
    photo_w(i_cpy)      = particle_min_weight

    ! Create new particle
    new_part%x = 0
    new_part%x(1:NDIM) = x_stop(1:NDIM)
    new_part%v = 0
    new_part%a = pc%accel_function(new_part)
    new_part%w = particle_min_weight
    new_part%id = af_get_id_at(tree, x_stop(1:NDIM))

    call pc%add_part(new_part)
  end subroutine single_photoionization_event

  ! Check if the photon is absorbed in the gas
  logical function is_in_gas(tree, x_stop)
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
