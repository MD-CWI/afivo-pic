#include "../afivo/src/cpp_macros.h"
!> Module for the electric field
module m_field
  use m_af_all
  use m_globals
  use m_domain

  implicit none
  private

  !> Start modifying the vertical background field after this time
  real(dp) :: field_mod_t0 = 1e99_dp

  !> Amplitude of sinusoidal modification
  real(dp) :: field_sin_amplitude = 0.0_dp

  !> Frequency (Hz) of sinusoidal modification
  real(dp) :: field_sin_freq = 0.0_dp

  !> Linear derivative of background field
  real(dp) :: field_lin_deriv = 0.0_dp

  !> Decay time of background field
  real(dp) :: field_decay_time = huge(1.0_dp)

  !> The applied electric field (vertical direction)
  real(dp) :: field_amplitude = 1.0e6_dp

  !> The applied voltage (vertical direction)
  real(dp) :: field_voltage

  !> Drop-off radius
  real(dp) :: field_dropoff_radius = 1e-3_dp

  !> Relative width over which the potential drops
  real(dp) :: field_dropoff_relwidth = 0.5_dp

  ! Number of V-cycles to perform per time step
  integer, protected :: multigrid_num_vcycles = 2

  character(GL_slen) :: field_bc_type = "homogeneous"

  public :: field_initialize
  public :: field_compute
  public :: field_from_potential
  public :: field_get_amplitude
  public :: field_set_voltage
  public :: field_bc_homogeneous

contains

  !> Set boundary conditions for the electric potential/field
  subroutine field_initialize(cfg, mg)
    use m_user_methods
    type(CFG_t), intent(inout)  :: cfg !< Settings
    type(mg_t), intent(inout) :: mg  !< Multigrid option struct

    call CFG_add_get(cfg, "field%mod_t0", field_mod_t0, &
         "Modify electric field after this time (s)")
    call CFG_add_get(cfg, "field%sin_amplitude", field_sin_amplitude, &
         "Amplitude of sinusoidal modification (V/m)")
    call CFG_add_get(cfg, "field%sin_freq", field_sin_freq, &
         "Frequency of sinusoidal modification (Hz)")
    call CFG_add_get(cfg, "field%lin_deriv", field_lin_deriv, &
         "Linear derivative of field [V/(ms)]")
    call CFG_add_get(cfg, "field%decay_time", field_decay_time, &
         "Decay time of field (s)")
    call CFG_add_get(cfg, "field%amplitude", field_amplitude, &
         "The applied electric field (V/m) (vertical)")
    call CFG_add_get(cfg, "field%bc_type", field_bc_type, &
         "Type of boundary condition to use (homogeneous, ...)")

    call CFG_add_get(cfg, "field%dropoff_radius", field_dropoff_radius, &
         "Potential stays constant up to this radius")
    call CFG_add_get(cfg, "field%dropoff_relwidth", field_dropoff_relwidth, &
         "Relative width over which the potential drops")
    call CFG_add_get(cfg, "multigrid_num_vcycles", multigrid_num_vcycles, &
         "Number of V-cycles to perform per time step")

    call field_set_voltage(0.0_dp)

    if (associated(user_potential_bc)) then
       mg%sides_bc => user_potential_bc
    else
       mg%sides_bc => field_bc_homogeneous
    end if

  end subroutine field_initialize

  !> Compute electric field on the tree. First perform multigrid to get electric
  !> potential, then take numerical gradient to get field.
  subroutine field_compute(tree, mg, have_guess)
    use m_units_constants
    type(af_t), intent(inout) :: tree
    type(mg_t), intent(inout) :: mg ! Multigrid option struct
    logical, intent(in)       :: have_guess
    real(dp), parameter       :: fac = -UC_elem_charge / UC_eps0
    integer                   :: lvl, i, id, nc

    nc = tree%n_cell

    ! Set the source term (rhs)
    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          tree%boxes(id)%cc(DTIMES(:), i_rhs) = fac * (&
               tree%boxes(id)%cc(DTIMES(:), i_pos_ion) - &
               tree%boxes(id)%cc(DTIMES(:), i_electron))
       end do
       !$omp end do nowait
    end do
    !$omp end parallel

    if (GL_use_dielectric) then
       ! Map surface charge to the right-hand side
       call dielectric_surface_charge_to_rhs(tree, diel, i_surf_sum_dens, &
            i_rhs, fac)
    end if

    call field_set_voltage(GL_time)

    if (.not. have_guess) then
       ! Perform a FMG cycle when we have no guess
       call mg_fas_fmg(tree, mg, .false., have_guess)
    end if

    ! Perform cheaper V-cycles
    do i = 1, multigrid_num_vcycles
       call mg_fas_vcycle(tree, mg, .false.)
    end do

    ! Compute field from potential
    call af_loop_box(tree, field_from_potential)

    ! Set ghost cells for the field components
    call af_gc_tree(tree, [i_E_all])

    if (GL_use_dielectric) then
       ! call dielectric_correct_field_cc(tree, diel, i_surf_sum_dens, &
            ! i_E_all, i_phi, -fac)
       call dielectric_correct_field_fc(tree, diel, i_surf_sum_dens, &
            ifc_E, i_phi, -fac)
    end if

    call af_loop_box(tree, compute_field_norm)

    ! Set the field norm also in ghost cells
    call af_gc_tree(tree, [i_E])

  end subroutine field_compute

  !> Compute the electric field at a given time
  function field_get_amplitude(time) result(electric_fld)
    use m_units_constants
    real(dp), intent(in) :: time
    real(dp)             :: electric_fld, t_rel

    t_rel = time - field_mod_t0
    if (t_rel > 0) then
       electric_fld = field_amplitude * exp(-t_rel/field_decay_time) + &
            t_rel * field_lin_deriv + &
            field_sin_amplitude * &
            sin(t_rel * field_sin_freq * 2 * UC_pi)
    else
       electric_fld = field_amplitude
    end if
  end function field_get_amplitude

  !> Compute the voltage at a given time
  subroutine field_set_voltage(time)
    real(dp), intent(in) :: time

    field_voltage = -domain_len(NDIM) * field_get_amplitude(time)
  end subroutine field_set_voltage

  !> This fills ghost cells near physical boundaries for the potential
  subroutine field_bc_homogeneous(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type

    if (af_neighb_dim(nb) == NDIM) then
       if (af_neighb_low(nb)) then
          bc_type = af_bc_dirichlet
          bc_val = 0.0_dp
       else
          bc_type = af_bc_dirichlet
          bc_val  = field_voltage
       end if
    else
       bc_type = af_bc_neumann
       bc_val = 0.0_dp
    end if
  end subroutine field_bc_homogeneous

  !> Compute electric field from electrical potential
  subroutine field_from_potential(box)
    type(box_t), intent(inout) :: box
    integer                     :: nc
    real(dp)                    :: inv_dr(NDIM)

    nc     = box%n_cell
    inv_dr = 1 / box%dr

#if NDIM == 2
    box%cc(1:nc, 1:nc, i_Ex) = 0.5_dp * inv_dr(1) * &
         (box%cc(0:nc-1, 1:nc, i_phi) - box%cc(2:nc+1, 1:nc, i_phi))
    box%cc(1:nc, 1:nc, i_Ey) = 0.5_dp * inv_dr(2) * &
         (box%cc(1:nc, 0:nc-1, i_phi) - box%cc(1:nc, 2:nc+1, i_phi))
    box%fc(1:nc+1, 1:nc, 1, ifc_E) = inv_dr(1) * &
         (box%cc(0:nc, 1:nc, i_phi) - box%cc(1:nc+1, 1:nc, i_phi))
    box%fc(1:nc, 1:nc+1, 2, ifc_E) = inv_dr(2) * &
         (box%cc(1:nc, 0:nc, i_phi) - box%cc(1:nc, 1:nc+1, i_phi))
#elif NDIM == 3
    box%cc(1:nc, 1:nc, 1:nc, i_Ex) = 0.5_dp * inv_dr(1) * &
         (box%cc(0:nc-1, 1:nc, 1:nc, i_phi) - box%cc(2:nc+1, 1:nc, 1:nc, i_phi))
    box%cc(1:nc, 1:nc, 1:nc, i_Ey) = 0.5_dp * inv_dr(2) * &
         (box%cc(1:nc, 0:nc-1, 1:nc, i_phi) - box%cc(1:nc, 2:nc+1, 1:nc, i_phi))
    box%cc(1:nc, 1:nc, 1:nc, i_Ez) = 0.5_dp * inv_dr(3) * &
         (box%cc(1:nc, 1:nc, 0:nc-1, i_phi) - box%cc(1:nc, 1:nc, 2:nc+1, i_phi))
    box%fc(1:nc+1, 1:nc, 1:nc, 1, ifc_E) = inv_dr(1) * &
         (box%cc(0:nc, 1:nc, 1:nc, i_phi) - &
         box%cc(1:nc+1, 1:nc, 1:nc, i_phi))
    box%fc(1:nc, 1:nc+1, 1:nc, 2, ifc_E) = inv_dr(2) * &
         (box%cc(1:nc, 0:nc, 1:nc, i_phi) - &
         box%cc(1:nc, 1:nc+1, 1:nc, i_phi))
    box%fc(1:nc, 1:nc, 1:nc+1, 3, ifc_E) = inv_dr(3) * &
         (box%cc(1:nc, 1:nc, 0:nc, i_phi) - &
         box%cc(1:nc, 1:nc, 1:nc+1, i_phi))
#endif

  end subroutine field_from_potential

  subroutine compute_field_norm(box)
    type(box_t), intent(inout) :: box
    integer                    :: nc

    nc = box%n_cell
#if NDIM == 2
    box%cc(1:nc, 1:nc, i_E) = sqrt(box%cc(1:nc, 1:nc, i_Ex)**2 + &
         box%cc(1:nc, 1:nc, i_Ey)**2)
    box%cc(1:nc, 1:nc, i_E_v2) = 0.5_dp * sqrt(&
         (box%fc(1:nc, 1:nc, 1, ifc_E) + &
         box%fc(2:nc+1, 1:nc, 1, ifc_E))**2 + &
         (box%fc(1:nc, 1:nc, 2, ifc_E) + &
         box%fc(1:nc, 2:nc+1, 2, ifc_E))**2)
#elif NDIM == 3
    box%cc(1:nc, 1:nc, 1:nc, i_E) = sqrt(&
         box%cc(1:nc, 1:nc, 1:nc, i_Ex)**2 + &
         box%cc(1:nc, 1:nc, 1:nc, i_Ey)**2 + &
         box%cc(1:nc, 1:nc, 1:nc, i_Ez)**2)
    box%cc(1:nc, 1:nc, 1:nc, i_E_v2) = 0.5_dp * sqrt(&
             (box%fc(1:nc, 1:nc, 1:nc, 1, ifc_E) + &
              box%fc(2:nc+1, 1:nc, 1:nc, 1, ifc_E))**2 + &
             (box%fc(1:nc, 1:nc, 1:nc, 2, ifc_E) + &
              box%fc(1:nc, 2:nc+1, 1:nc, 2, ifc_E))**2 + &
             (box%fc(1:nc, 1:nc, 1:nc, 3, ifc_E) + &
              box%fc(1:nc, 1:nc, 2:nc+1, 3, ifc_E))**2)
#endif
  end subroutine compute_field_norm

end module m_field
