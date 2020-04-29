#include "../afivo/src/cpp_macros.h"
module m_photons
! This module contains the functionality for
! photon interaction (ionization and electron emission from dielectric surface)
use m_globals
use m_af_all
use m_particle_core
use m_cross_sec
use m_units_constants
use m_af_all

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

!Default reaction rates for excited Ar and Ar2 pools
real(dp)  :: k_Ar_decay_rad  = 0.0_dp
real(dp)  :: k_Ar_quench     = 0.0_dp
real(dp)  :: k_Ar2_prod_rate = 3.0e-46_dp
real(dp)  :: k_Ar2_decay_rad = 6.0e7_dp


real(dp)   :: k_Ar4p_Ar4s = 1.0e7_dp
real(dp)   :: k_Ar4d_Ar4s = 1.6e5_dp
real(dp)   :: k_Ar_e_que  = 4.5e-16_dp
real(dp)   :: k_Ar_e_ion  = 1.0e-12_dp
! real(dp) :: k_Ar_Ar_ion = 6.0e-16_dp
real(dp)   :: k_Ar1s5_Ar2 = 1.0e-44_dp
real(dp)   :: k_Ar1s4_Ar2 = 3.0e-46_dp
real(dp)   :: k_Ar2_Ar2_1 = 5.0e-14_dp
real(dp)   :: k_Ar2_Ar    = 1.0e-14_dp
real(dp)   :: k_Ar2_Ar2_p = 1.0e-12_dp
real(dp)   :: k_Ar2_3_hv  = 3.5e5_dp
real(dp)   :: k_Ar2_1_hv  = 2.14e8_dp
real(dp)   :: k_Ar1s4_hv  = 1.2e8_dp



integer :: ground_gas = -3 !< make sure can not be reached by other variable
integer, protected, public  :: i_Ar_pool = -1 !Index of cell-centered pool of excited argon species
integer, protected, public  :: i_Ar_s5_pool = -1 !Index of cell-centered Argon excimer
integer, protected, public  :: i_Ar_s4_pool = -1 !Index of cell-centered Argon excimer
integer, protected, public  :: i_Ar_4p_pool = -1 !Index of cell-centered Argon excimer
integer, protected, public  :: i_Ar_4d_pool = -1 !Index of cell-centered Argon excimer

integer, protected, public  :: i_Ar2_pool = -1 !Index of cell-centered Argon excimer
integer, protected, public  :: i_Ar2_1_pool = -1 !Index of cell-centered Argon excimer



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
        error stop "Unrecognized photon model"
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

      call af_add_cc_variable(tree, "Ar_s5_pool", .true., ix=i_Ar_s5_pool)
      call af_set_cc_methods(tree, i_Ar_s5_pool, af_bc_neumann_zero, &
      prolong=af_prolong_limit)
      call af_add_cc_variable(tree, "Ar_s4_pool", .true., ix=i_Ar_s4_pool)
      call af_set_cc_methods(tree, i_Ar_s4_pool, af_bc_neumann_zero, &
      prolong=af_prolong_limit)
      call af_add_cc_variable(tree, "Ar_4p_pool", .true., ix=i_Ar_4p_pool)
      call af_set_cc_methods(tree, i_Ar_4p_pool, af_bc_neumann_zero, &
      prolong=af_prolong_limit)
      call af_add_cc_variable(tree, "Ar_4d_pool", .true., ix=i_Ar_4d_pool)
      call af_set_cc_methods(tree, i_Ar_4d_pool, af_bc_neumann_zero, &
      prolong=af_prolong_limit)

      call af_add_cc_variable(tree, "Ar2_pool", .true., ix=i_Ar2_pool)
      call af_set_cc_methods(tree, i_Ar2_pool, af_bc_neumann_zero, &
      prolong=af_prolong_limit)
      call af_add_cc_variable(tree, "Ar2_1_pool", .true., ix=i_Ar2_1_pool)
      call af_set_cc_methods(tree, i_Ar2_1_pool, af_bc_neumann_zero, &
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
    real, allocatable,save :: filter(:,:)
    integer :: species = 5
    integer       :: n, jj

    if (.not. allocated(weights)) then
       allocate(coords(NDIM, pc%n_events))
       allocate(weights(pc%n_events))
       allocate(id_guess(pc%n_events))
       allocate(filter(species,pc%n_events))
    else if (size(weights) < pc%n_events) then
       deallocate(coords)
       deallocate(weights)
       deallocate(id_guess)
       deallocate(filter)
       allocate(coords(NDIM, pc%n_events))
       allocate(weights(pc%n_events))
       allocate(id_guess(pc%n_events))
       allocate(filter(species,pc%n_events))

    end if

    ! Calculate production to the pool of excited Argon species
    jj = 0
    filter = 0

    do n = 1, pc%n_events
       if (pc%event_list(n)%ctype == CS_excite_t) then
         jj = jj + 1
         coords(:, jj) = pc%event_list(n)%part%x(1:NDIM)
         weights(jj) = pc%event_list(n)%part%w
         id_guess(jj) = pc%event_list(n)%part%id
         !> filter the 1s5
         if (pc%event_list(n)%cIx == 2) then ! TODO do this more generally
           filter(1,jj) = 1.0_dp
         !> filter the 1s4
         else if (pc%event_list(n)%cIx == 3) then
           filter(2,jj) = 1.0_dp
         !> filter the 4p
         else if (pc%event_list(n)%cIx >= 6 .and. &
                  pc%event_list(n)%cIx <= 14) then
           filter(3,jj) = 1.0_dp
         !> filter the 4d (3d, 5s, 5p, 4d)
        else if (pc%event_list(n)%cIx >= 15 .and. &
                  pc%event_list(n)%cIx <= 26) then
           filter(4,jj) = 1.0_dp
         !> filter the 4s(1s3, 1s2)
         else
           filter(5,jj) = 1.0_dp
         end if
       end if
     end do
     call af_particles_to_grid(tree, i_Ar_s5_pool, coords(:, 1:jj), &
              filter(1,1:jj) * weights(1:jj), jj, 0, id_guess(1:jj)) ! Use zeroth order interpolation for simplicity
     call af_particles_to_grid(tree, i_Ar_s4_pool, coords(:, 1:jj), &
              filter(2,1:jj) * weights(1:jj), jj, 0, id_guess(1:jj)) ! Use zeroth order interpolation for simplicity
     call af_particles_to_grid(tree, i_Ar_4p_pool, coords(:, 1:jj), &
              filter(3,1:jj) * weights(1:jj), jj, 0, id_guess(1:jj)) ! Use zeroth order interpolation for simplicity
     call af_particles_to_grid(tree, i_Ar_4d_pool, coords(:, 1:jj), &
              filter(4,1:jj) * weights(1:jj), jj, 0, id_guess(1:jj)) ! Use zeroth order interpolation for simplicity
     call af_particles_to_grid(tree, i_Ar_pool, coords(:, 1:jj), &
              filter(5,1:jj) * weights(1:jj), jj, 0, id_guess(1:jj)) ! Use zeroth order interpolation for simplicity



     call af_loop_box(tree, perform_argon_reaction3, .true.)
     call af_loop_box(tree, perform_argon_reaction2, .true.)
     call af_loop_box(tree, perform_argon_reaction1, .true.)



     call Ar2_radiative_decay(tree, pc) ! Do photoemission events (updates the Ar2 pool)

end subroutine photoe_Argon

!TODO reorganize with more generic methods
! subroutine react_3body(box, s_in1, s_in2, s_in3, rate_constant, s_out)
!       use m_gas
!       type(box_t), intent(inout)  :: box
!       integer, intent(in) ::  s_in1, s_in2, s_in3, s_out !< the index of the variables (species)
!       real(dp) :: rate_constant !< the reaction rate constant with unit /s
!       integer    :: ii, jj, k, n, nc
!       integer   :: combs(3)
!
!       ! ! TODO auto judge if s1, s2, s3 is ground gas
!       ! equivalence (s_in1, combs(1))
!       ! equivalence (s_in2, combs(2))
!       ! equivalence (s_in3, combs(3))
!       !
!       ! n = 0
!       ! do k = 1, 3
!       !   if (combs(k) == ground_gas) &
!       !     n = n + 1
!       ! end do
!       nc  = box%n_cell
!
!   ! #if NDIM == 2
!         ! if (n == 2) then
!           box%cc(1:nc, 1:nc, s_out) = box%cc(1:nc, 1:nc, s_out) + &
!             GL_dt * rate_constant * GAS_number_dens**n * box%cc(1:nc, 1:nc, s_in1)
!           ! decay, quenching and Ar2-prod for Ar
!           box%cc(1:nc, 1:nc, s_in1) = box%cc(1:nc, 1:nc, s_in1) * (1- &
!             GL_dt * rate_constant * GAS_number_dens**n)
!         ! else
!         !   error stop
!         ! end if
!     ! #elif NDIM == 3
!     !     error stop
! end subroutine react_3body

subroutine Ar2_radiative_decay(tree, pc)
  ! Routine that performs Ar2_exc radiative decay for every cell in a box
  ! This routine also performs photoemission

  ! TODO Code from af_loop_box is copied (because of particle creation and non-local effects). Thats why this code is nasty
  type(af_t), intent(inout)  :: tree
  type(PC_t), intent(inout)  :: pc

  integer    :: ii, jj, nn
  integer    :: lvl, i, id
  integer    :: n_uv, n_uv1, n_uv2
  real(dp)   :: mean_ph1, mean_ph2
  real(dp)   :: cell_size
  real(dp)   :: x_start(3), x_cc(3), x_stop(3)
  logical    :: on_surface

  if (.not. tree%ready) stop "Ar2_radiative_decay: set_base has not been called"

  !TODO Find out how to work with this omp stuff
  !!$omp parallel private(lvl, i, id, ii, jj, nn)
  do lvl = 1, tree%highest_lvl
    !!$omp do
    do i = 1, size(tree%lvls(lvl)%leaves)
      id = tree%lvls(lvl)%leaves(i)
      cell_size = product(tree%boxes(id)%dr)
      do ii = 1, tree%boxes(id)%n_cell
        do jj = 1, tree%boxes(id)%n_cell

          mean_ph1 = cell_size * GL_dt * k_Ar2_1_hv * &
            tree%boxes(id)%cc(ii, jj, i_Ar2_1_pool)
          mean_ph2 = cell_size * GL_dt * k_Ar1s4_hv * &
            tree%boxes(id)%cc(ii, jj, i_Ar_s4_pool)  !/ particle_min_weight

          if (mean_ph1 > 1_dp) then
              n_uv1 = GL_rng%poisson(mean_ph1)   !*(GL_dt/(4*1e-9_dp))
            else
              n_uv1 = 0
            end if
            if (mean_ph2 > 1_dp) then
              n_uv2 = GL_rng%poisson(mean_ph2)   !*(GL_dt/(4*1e-9_dp))
            else
              n_uv2 = 0
            end if

            n_uv = n_uv1

            ! if (n_uv > 0) &
            !   print *, "the number photon is",n_uv

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

          tree%boxes(id)%cc(ii, jj, i_Ar2_1_pool) = tree%boxes(id)%cc(ii, jj, i_Ar2_1_pool) - &
          ! n_uv * particle_min_weight / cell_size !Update Ar2 density due to radiative decay
          n_uv1 * 1.0_dp / cell_size !Update Ar2 density due to radiative decay
          tree%boxes(id)%cc(ii, jj, i_Ar_s4_pool) = tree%boxes(id)%cc(ii, jj, i_Ar_s4_pool) - &
          ! n_uv * particle_min_weight / cell_size !Update Ar2 density due to radiative decay
          n_uv2 * 1.0_dp / cell_size !Update Ar2 density due to radiative decay

          if (tree%boxes(id)%cc(ii, jj, i_Ar2_1_pool) < 0.0_dp .or. &
              tree%boxes(id)%cc(ii, jj, i_Ar_s4_pool) < 0.0_dp) then
            write(*, "(A20,E12.4)") "ar excimer density : ", &
            tree%boxes(id)%cc(ii, jj, i_Ar2_1_pool)
            write(*, "(A20,E12.4)") "ar atomic density : ", &
            tree%boxes(id)%cc(ii, jj, i_Ar_s4_pool)
            ! print *, "photons from excimer : ", n_uv1
            ! print *, "photons from atom : ", n_uv2
            ! print *, "mean photons from excimer: ",mean_ph1
            ! print *, "mean photons from atom: ",mean_ph2

            tree%boxes(id)%cc(ii, jj, i_Ar2_1_pool) = max(0.0_dp, tree%boxes(id)%cc(ii, jj, i_Ar2_1_pool))
            tree%boxes(id)%cc(ii, jj, i_Ar2_1_pool) = max(0.0_dp, tree%boxes(id)%cc(ii, jj, i_Ar_s4_pool))
            error stop "Negative density for excited states of Argon2 found after radiative decay."
          end if

        end do
      end do
    end do
    !!$omp end do
  end do
  !!$omp end parallel
  end subroutine Ar2_radiative_decay
!
  subroutine perform_argon_reaction1(box)
    use m_gas
    ! Per cell, do explicit Euler to update the Ar-Ar2 pools due to reaction mechanism (excluding Ar* production and Ar2 radiative decay)
    ! type(af_t), intent(inout) :: tree
    type(box_t), intent(inout)  :: box
    integer   :: nc
    !> test debug
    real(dp) :: min1_den , min2_den
    !< test debug

    nc = box%n_cell

    ! the density of Ar* (1s4)
    box%cc(1:nc, 1:nc, i_Ar_s4_pool) = box%cc(1:nc, 1:nc, i_Ar_s4_pool)  + &
    !> Ar4p, Ar4d -> Ar4s
      GL_dt * k_Ar4p_Ar4s * box%cc(1:nc, 1:nc, i_Ar_4p_pool) + &
      GL_dt * k_Ar4d_Ar4s * box%cc(1:nc, 1:nc, i_Ar_4d_pool) - &
      ! >Ar1s4 -> Ar2_siglet (3 body reaction)
      GL_dt * k_Ar1s4_Ar2 * GAS_number_dens**2 * box%cc(1:nc, 1:nc, i_Ar_s4_pool) - &
      !> Ar* -> Ar and Ar+
      GL_dt * (k_Ar_e_ion + k_Ar_e_que) * box%cc(1:nc, 1:nc, i_electron) * box%cc(1:nc, 1:nc, i_Ar_s4_pool)

    ! Production of Ar2* singlet
    box%cc(1:nc, 1:nc, i_Ar2_1_pool) = box%cc(1:nc, 1:nc, i_Ar2_1_pool)  + &
      !> Ar1s4 -> Ar2_siglet (3 body reaction)
      GL_dt * k_Ar1s4_Ar2 * GAS_number_dens**2 * box%cc(1:nc, 1:nc, i_Ar_s4_pool) + &
     !> Ar2_triplet -> Ar2_siglet (electron collision)
      GL_dt * k_Ar2_Ar2_1 * box%cc(1:nc, 1:nc, i_electron) * box%cc(1:nc, 1:nc, i_Ar2_pool) - &
      !> Ar2_siglet -> Ar* and Ar2+
      GL_dt * (k_Ar2_Ar + k_Ar2_Ar2_p) * box%cc(1:nc, 1:nc, i_electron) * box%cc(1:nc, 1:nc, i_Ar2_1_pool)

      !> test debug
      call af_tree_min_cc(tree, i_Ar_s4_pool, min1_den)
      call af_tree_min_cc(tree, i_Ar2_1_pool, min2_den)
      if (min1_den < 0 .or. min2_den < 0) then
      ! print *, "WARNING: the density of 1s4 is : ", min1_den
      ! print *, "WARNING: the density of Ar2 singlet is : ", min2_den
      box%cc(1:nc, 1:nc, i_Ar_s4_pool) = max(0.0_dp, box%cc(1:nc, 1:nc, i_Ar_s4_pool))
      box%cc(1:nc, 1:nc, i_Ar2_1_pool) = max(00.0_dp, box%cc(1:nc, 1:nc, i_Ar2_1_pool))
      ! error stop "negative density"
      end if
      !< test debug
  end subroutine perform_argon_reaction1

  subroutine perform_argon_reaction2(box)
    use m_gas
    ! Per cell, do explicit Euler to update the Ar-Ar2 pools due to reaction mechanism (excluding Ar* production and Ar2 radiative decay)
    ! type(af_t), intent(inout) :: tree
    type(box_t), intent(inout)  :: box
    integer   :: nc

    !> test debug
    real(dp) :: min1_den , min2_den
    !< test debug

    nc = box%n_cell

    ! the density of Ar* (1s5)
    box%cc(1:nc, 1:nc, i_Ar_s5_pool) = box%cc(1:nc, 1:nc, i_Ar_s5_pool)  + &
    !> Ar4p, Ar4d -> Ar4s
      GL_dt * k_Ar4p_Ar4s * box%cc(1:nc, 1:nc, i_Ar_4p_pool) + &
      GL_dt * k_Ar4d_Ar4s * box%cc(1:nc, 1:nc, i_Ar_4d_pool) - &
      ! >Ar1s5 -> Ar2_triplet (3 body reaction)
      GL_dt * k_Ar1s5_Ar2 * GAS_number_dens**2 * box%cc(1:nc, 1:nc, i_Ar_s5_pool) - &
      !> Ar* -> Ar and Ar+
      GL_dt * (k_Ar_e_ion + k_Ar_e_que) * box%cc(1:nc, 1:nc, i_electron) * box%cc(1:nc, 1:nc, i_Ar_s5_pool)

     ! Production of Ar2* triplet
     box%cc(1:nc, 1:nc, i_Ar2_pool) = box%cc(1:nc, 1:nc, i_Ar2_pool)  + &
      !> Ar1s5 -> Ar2_triplet (3 body reaction)
      GL_dt * k_Ar1s5_Ar2 * GAS_number_dens**2 * box%cc(1:nc, 1:nc, i_Ar_s5_pool) - &
     !> Ar2_triplet -> Ar2_siglet (electron collision)
      GL_dt * k_Ar2_Ar2_1 * box%cc(1:nc, 1:nc, i_electron) * box%cc(1:nc, 1:nc, i_Ar2_pool) - &
      !> Ar2_triplet -> Ar* and Ar2+
      GL_dt * (k_Ar2_Ar + k_Ar2_Ar2_p) * box%cc(1:nc, 1:nc, i_electron) * box%cc(1:nc, 1:nc, i_Ar2_pool)

      !> test debug
      call af_tree_min_cc(tree, i_Ar_s5_pool, min1_den)
      call af_tree_min_cc(tree, i_Ar2_pool, min2_den)
      if (min1_den < 0 .or. min2_den < 0) then
        ! print *, "WARNING: the density of 1s5 is : ", min1_den
        ! print *, "WARNING: the density of Ar2 triplet is : ", min2_den
        box%cc(1:nc, 1:nc, i_Ar_s5_pool) = max(0.0_dp, box%cc(1:nc, 1:nc, i_Ar_s5_pool))
        box%cc(1:nc, 1:nc, i_Ar2_pool) = max(00.0_dp, box%cc(1:nc, 1:nc, i_Ar2_pool))
        ! error stop "negative density"
      end if
      !< test debug
  end subroutine perform_argon_reaction2

  subroutine perform_argon_reaction3(box)
    use m_gas
    ! Per cell, do explicit Euler to update the Ar-Ar2 pools due to reaction mechanism (excluding Ar* production and Ar2 radiative decay)
    ! type(af_t), intent(inout) :: tree
    type(box_t), intent(inout)  :: box
    integer   :: nc
    !> test debug
    real(dp) :: min1_den , min2_den
    !< test debug

    ! the density of high energy atomic state Ar* (4p , 4d)
    nc = box%n_cell
        ! the density of Ar*(4p)
        box%cc(1:nc, 1:nc, i_Ar_4p_pool) = box%cc(1:nc, 1:nc, i_Ar_4p_pool) - &
          !> Ar4p -> Ar4s
          GL_dt * k_Ar4p_Ar4s*2 * box%cc(1:nc, 1:nc, i_Ar_4p_pool) - &
          !> Ar4p -> Ar and Ar4p -> Ar+ (2body reaction with e)
          GL_dt * (k_Ar_e_ion + k_Ar_e_que) * box%cc(1:nc, 1:nc, i_electron) * box%cc(1:nc, 1:nc, i_Ar_4p_pool)
        ! the density of Ar*(4d)
        box%cc(1:nc, 1:nc, i_Ar_4d_pool) = box%cc(1:nc, 1:nc, i_Ar_4d_pool) - &
          !> Ar4d -> Ar4s
          GL_dt * k_Ar4d_Ar4s*2 * box%cc(1:nc, 1:nc, i_Ar_4d_pool) - &
          !> Ar4d -> Ar and Ar4d -> Ar+ (2body reaction with e)
          GL_dt * (k_Ar_e_ion + k_Ar_e_que) * box%cc(1:nc, 1:nc, i_electron) * box%cc(1:nc, 1:nc, i_Ar_4d_pool)

          !> test debug
          call af_tree_min_cc(tree, i_Ar_4p_pool, min1_den)
          call af_tree_min_cc(tree, i_Ar_4d_pool, min2_den)
          if (min1_den < 0 .or. min2_den < 0) then
          ! print *, "WARNING: the density of 1s5 is : ", min1_den
          ! print *, "WARNING: the density of Ar2 triplet is : ", min2_den
          box%cc(1:nc, 1:nc, i_Ar_4p_pool) = max(0.0_dp, box%cc(1:nc, 1:nc, i_Ar_4p_pool))
          box%cc(1:nc, 1:nc, i_Ar_4d_pool) = max(0.0_dp, box%cc(1:nc, 1:nc, i_Ar_4d_pool))
          ! error stop "negative density"
          end if
          !< test debug

  end subroutine perform_argon_reaction3

  ! subroutine perform_argon_reactions(box)
  !   use m_gas
  !   ! Per cell, do explicit Euler to update the Ar-Ar2 pools due to reaction mechanism (excluding Ar* production and Ar2 radiative decay)
  !   ! type(af_t), intent(inout) :: tree
  !   type(box_t), intent(inout)  :: box
  !   integer   :: nc
  !   integer    :: ii, jj
  !
  !   nc = box%n_cell
  !       ! Production of Ar2
  !       box%cc(1:nc, 1:nc, i_Ar2_pool) = box%cc(1:nc, 1:nc, i_Ar2_pool)  + &
  !         GL_dt * k_Ar2_prod_rate * GAS_number_dens**2 * box%cc(ii, jj, i_Ar_pool)
  !       ! decay, quenching and Ar2-prod for Ar
  !       box%cc(ii, jj, i_Ar_pool) = box%cc(ii, jj, i_Ar_pool) - GL_dt * &
  !         (k_Ar_decay_rad + k_Ar_quench * GAS_number_dens + k_Ar2_prod_rate * GAS_number_dens**2) * box%cc(ii, jj, i_Ar_pool)
  !       if (box%cc(ii, jj, i_Ar2_pool) < 0.0_dp .or. box%cc(ii, jj, i_Ar_pool) < 0.0_dp) then
  !         error stop "Negative density for excited states of Argon or Argon2 found after performing chemical reactions."
  !       end if
  !     end do
  !   end do
  ! end subroutine perform_argon_reactions


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

    if (photoi_enabled .or. photoe_enabled) then
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

  subroutine photoi_Zheleznyak(tree, pc, photo_pos, photo_w, n_photons)
    ! Perform photon generation and ionization according to the Zheleznyak model
    ! use m_domain

    type(af_t), intent(inout)     :: tree
    type(PC_t), intent(inout)     :: pc
    real(dp), intent(inout)       :: photo_pos(:, :)
    real(dp), intent(inout)       :: photo_w(:)
    integer, intent(out)          :: n_photons

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

  ! ==== General modules

  subroutine single_photoemission_event(tree, pc, photon_w, x_gas, x_diel)
    ! use m_particles, only: surface_charge_to_particle
    ! Generate photo-emitted electrons on the surface of a dielectric
    type(af_t), intent(inout)        :: tree
    type(PC_t), intent(inout)        :: pc
    type(PC_part_t)                  :: new_part
    real(dp), intent(in)             :: photon_w, x_gas(3), x_diel(3)

    if (GL_rng%unif_01() < photoe_probability) then
      ! Create photo-emitted electron
      new_part%x(:) = x_gas
      new_part%v(:) = 0.0_dp
      new_part%a(:) = pc%accel_function(new_part)
      new_part%w    = photon_w
      new_part%id   = af_get_id_at(tree, x_gas(1:NDIM))

      call pc%add_part(new_part)
      ! Update charge density on the dielectric
      call surface_charge_to_particle(tree, new_part, i_surf_elec)
    end if
  end subroutine single_photoemission_event

  subroutine surface_charge_to_particle(tree, my_part, i_surf)
    ! Input: a particle that is ejected from the dielectric.
    ! The surface charge is altered by the charge leaving
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
    error stop
#endif
  end subroutine surface_charge_to_particle

  subroutine single_photoionization_event(tree, pc, i, photo_pos, photo_w, x_stop)
    ! use m_particles, only: get_accel

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
    new_part%x = 0
    new_part%x(1:NDIM) = x_stop(1:NDIM)
    new_part%v = 0
    new_part%a = pc%accel_function(new_part)
    new_part%w = particle_min_weight
    new_part%id = af_get_id_at(tree, x_stop(1:NDIM))

    call pc%add_part(new_part)
  end subroutine single_photoionization_event

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
