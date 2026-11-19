#include "../afivo/src/cpp_macros.h"
!> Module for the electric field
module m_field
  use m_af_all
  use m_globals
  use m_domain
  use m_geometry

  implicit none
  private

  !> Magnetic field
  real(dp), public, protected :: magnetic_field(3) = 0.0_dp

  !> Whether there is a magnetic field
  logical, public, protected :: magnetic_field_used = .false.

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
  real(dp), public, protected :: field_voltage

  !> Drop-off radius
  real(dp) :: field_dropoff_radius = 1e-3_dp

  !> Relative width over which the potential drops
  real(dp) :: field_dropoff_relwidth = 0.5_dp

  ! Number of V-cycles to perform per time step
  integer, protected :: multigrid_num_vcycles = 2

  character(GL_slen) :: field_bc_type = "homogeneous"

  !> Whether the electrode is grounded or at the applied voltage
  logical, public, protected :: field_electrode_grounded = .false.

  !> Electrode relative start position (for standard rod electrode)
  real(dp), public, protected :: field_rod_r0(NDIM) = -1.0_dp

  !> Electrode relative end position (for standard rod electrode)
  real(dp), public, protected :: field_rod_r1(NDIM) = -1.0_dp

  !> Electrode radius (in m, for standard rod electrode)
  real(dp), public, protected :: field_rod_radius = -1.0_dp

  !> ===========================================================================
  !> Electrode 1: fraction of conical part (if conical)
  real(dp) :: cone_length_frac = -1.0e100_dp
  !> Electrode 1: tip radius (if conical)
  real(dp) :: cone_tip_radius = -1.0e100_dp
  !> Electrode 1: 'origin' of spherical tip (if conical)
  real(dp) :: cone_tip_center(NDIM)
  !> Electrode 1: radius of curvature of spherical tip (if conical)
  real(dp) :: cone_tip_r_curvature
  character(GL_slen) :: electrode_type

  public :: field_initialize
  public :: field_compute
  public :: field_from_potential
  public :: field_get_amplitude
  public :: field_set_voltage
  public :: field_bc_homogeneous
  public :: field_set_lsf_box
  public :: field_get_lsf

contains

  !> Set boundary conditions for the electric potential/field
  subroutine field_initialize(cfg, mg)
    use m_user_methods
    type(CFG_t), intent(inout)  :: cfg !< Settings
    type(mg_t), intent(inout) :: mg  !< Multigrid option struct

    call CFG_add_get(cfg, "field%magnetic", magnetic_field, &
         "Background magnetic field vector (T)")
    magnetic_field_used = any(abs(magnetic_field) > 0)

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

    call CFG_add_get(cfg, "field%electrode_grounded", field_electrode_grounded, &
         "Whether the electrode is grounded or at the applied voltage")
    call CFG_add_get(cfg, "field%rod_r0", field_rod_r0, &
         "Electrode relative start position (for standard rod electrode)")
    call CFG_add_get(cfg, "field%rod_r1", field_rod_r1, &
         "Electrode relative end position (for standard rod electrode)")
    call CFG_add_get(cfg, "field%rod_radius", field_rod_radius, &
         "Electrode radius (in m, for standard rod electrode)")
    
    !> ===========================================================================
    call CFG_add_get(cfg, "field%cone_tip_radius", cone_tip_radius, &
         "Electrode 1: tip radius (if conical)")
    call CFG_add_get(cfg, "field%cone_length_frac", cone_length_frac, &
         "Electrode 1: fraction of conical part (if conical)")
    electrode_type = "rod"
    call CFG_add_get(cfg, "field%field_electrode_type", electrode_type, &
         "Type of electrode (sphere, rod, rod_cone_top, rod_rod, user)")
    !< [electrode_settings]
    select case (electrode_type)
    case ("rod")
       ! Standard cylindrical rod.
       ! No additional initialization required.
    case ("rod_cone_top")
    if (any(field_rod_r0 < 0.0_dp) .or. &
        any(field_rod_r1 < 0.0_dp)) then
        error stop "field rod coordinates not correctly specified"
    end if

    if (field_rod_radius <= 0.0_dp) then
       error stop "field_rod_radius must be positive"
    end if

    if (cone_tip_radius <= 0.0_dp .or. &
        cone_tip_radius > field_rod_radius) then
        error stop "Invalid cone_tip_radius"
    end if

    if (cone_length_frac <= 0.0_dp .or. &
        cone_length_frac > 1.0_dp) then
        error stop "Invalid cone_length_frac"
    end if
    call get_conical_rod_properties( &
          field_rod_r0 * domain_len, &
          field_rod_r1 * domain_len, &
          field_rod_radius, &
          cone_tip_radius, &
          cone_tip_center, &
          cone_tip_r_curvature)
   case default
           error stop "Unknown electrode type"
   end select

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
    real(dp)                  :: max_rhs, residual_threshold, conv_fac
    real(dp)                  :: residual_ratio
    integer, parameter        :: max_initial_iterations = 30
    real(dp), parameter       :: max_residual = 1e8_dp
    real(dp), parameter       :: min_residual = 1e-6_dp
    real(dp)                  :: residuals(max_initial_iterations)
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
          if (GL_use_electrode) then
             where (tree%boxes(id)%cc(DTIMES(:), i_lsf) <= 0)
                tree%boxes(id)%cc(DTIMES(:), i_rhs) = 0.0_dp
             end where
          end if
       end do
       !$omp end do nowait
    end do
    !$omp end parallel

    if (GL_use_dielectric) then
       ! Map surface charge to the right-hand side
       call surface_surface_charge_to_rhs(tree, diel, i_surf_sum_dens, &
            i_rhs, fac)
    end if

    call field_set_voltage(GL_time)

    call af_tree_maxabs_cc(tree, mg%i_rhs, max_rhs)

    ! With an electrode, the initial convergence testing should be less strict
    if (GL_use_electrode) then
       conv_fac = 1e-8_dp
       if (field_electrode_grounded) then
          mg%lsf_boundary_value = 0.0_dp
       else
          mg%lsf_boundary_value = field_voltage
       end if
    else
       conv_fac = 1e-10_dp
    end if

    ! Set threshold based on rhs and on estimate of round-off error, given by
    ! delta phi / dx^2 = (phi/L * dx)/dx^2
    ! Note that we use min_residual in case max_rhs and field_voltage are zero
    residual_threshold = max(min_residual, &
         max_rhs * GL_multigrid_max_rel_residual, &
         conv_fac * abs(field_voltage)/(domain_len(NDIM) * af_min_dr(tree)))

    ! Perform a FMG cycle when we have no guess
    if (.not. have_guess) then
      do i = 1, max_initial_iterations
        call mg_fas_fmg(tree, mg, .true., .true.)
        call af_tree_maxabs_cc(tree, mg%i_tmp, residuals(i))

        if (residuals(i) < residual_threshold) then
           exit
        else if (i > 2) then
           ! Check if the residual is not changing much anymore, and if it is
           ! small enough, in which case we assume convergence
           residual_ratio = minval(residuals(i-2:i)) / &
                maxval(residuals(i-2:i))
           if (residual_ratio < 2.0_dp .and. residual_ratio > 0.5_dp &
                .and. residuals(i) < max_residual) exit
        end if
      end do

       ! Check for convergence
       if (i == max_initial_iterations + 1) then
          print *, "Iteration    residual"
          do i = 1, max_initial_iterations
             write(*, "(I4,E18.2)") i, residuals(i)
          end do
          print *, "Maybe increase multigrid_max_rel_residual?"
          error stop "No convergence in initial field computation"
       end if
    end if

    ! Perform cheaper V-cycles
    do i = 1, multigrid_num_vcycles
       call mg_fas_vcycle(tree, mg, .true.)
       call af_tree_maxabs_cc(tree, mg%i_tmp, residuals(i))
       if (residuals(i) < residual_threshold) exit
    end do

    ! Compute field from potential
    if (GL_use_electrode) then
      call mg_compute_phi_gradient(tree, mg, ifc_E, -1.0_dp, i_E)
    else
      ! Keep this routine to ensure compatibility with trilinear interpolation
      ! If tril-interp will be used, incorporate this in `mg_compute_phi_gradient`
      call af_loop_box(tree, field_from_potential)
    end if

    ! Set ghost cells for the field components
    call af_gc_tree(tree, [i_E])

    if (GL_use_dielectric) then
       ! call surface_correct_field_cc(tree, diel, i_surf_sum_dens, &
            ! i_E_all, i_phi, -fac)
       call surface_correct_field_fc(tree, diel, i_surf_sum_dens, &
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

  ! This routine sets the level set function for a simple rod
  subroutine field_set_lsf_box(box, iv)
    use m_geometry
    type(box_t), intent(inout) :: box
    integer, intent(in)        :: iv
    integer                    :: IJK, nc
    real(dp)                   :: rr(NDIM)

    nc = box%n_cell

    do KJI_DO(0,nc+1)
       rr = af_r_cc(box, [IJK])
       box%cc(IJK, iv) = field_get_lsf(rr)
    end do; CLOSE_DO
  end subroutine field_set_lsf_box

  real(dp) function field_get_lsf(r)
    real(dp), intent(in) :: r(NDIM)
    select case (electrode_type)
    case ("rod")
       field_get_lsf = GM_dist_line( &
            r, &
            field_rod_r0 * domain_len, &
            field_rod_r1 * domain_len, &
            NDIM) - field_rod_radius

    case ("rod_cone_top")
       field_get_lsf = conical_rod_lsf(r)

    case default
       error stop "Unkown electrode type in field_get_lsf"
    end select
  end function field_get_lsf

  !> ===========================================================================
  !> Add routine to include conical electrode as is done in Afivo-streamer
  !> Compute several parameters for a conical rod
  subroutine get_conical_rod_properties(r0, r1, &
       rod_radius, tip_radius, cone_tip_center, cone_tip_r_curvature)
    real(dp), intent(in)  :: r0(NDIM) !< Beginning of rod
    real(dp), intent(in)  :: r1(NDIM) !< End of conical part
    real(dp), intent(in)  :: rod_radius   !< Radius of rod
    real(dp), intent(in)  :: tip_radius   !< Radius of curvature of tip
    real(dp), intent(out) :: cone_tip_center(NDIM)
    real(dp), intent(out) :: cone_tip_r_curvature
    real(dp)              :: cone_angle, cone_length

    ! Determine (half) the opening angle of the top cone, which goes from
    ! rod_radius to tip_radius over cone_length
    cone_length = cone_length_frac * norm2(r1 - r0)
    cone_angle  = atan((rod_radius - tip_radius) / cone_length)

    ! We have a point on a sphere with coordinates of the form (R*cos(a),
    ! R*sin(a)) = (tip_radius, y), so we can get R and subtract R sin(a) to
    ! obtain the center of the sphere
    cone_tip_r_curvature = tip_radius/cos(cone_angle)
    cone_tip_center = r1 - sin(cone_angle) * &
         cone_tip_r_curvature * (r1 - r0)/norm2(r1 - r0)
  end subroutine get_conical_rod_properties

  !> Helper function to compute lsf for a rod with a conical tip
  pure function conical_rod_lsf_arg(r, r0, r1, cone_tip_center, &
       rod_radius, tip_radius, cone_length_frac, r_curvature) result (lsf)
    use m_geometry
    real(dp), intent(in) :: r(NDIM)
    real(dp), intent(in) :: r0(NDIM), r1(NDIM), cone_tip_center(NDIM)
    real(dp), intent(in) :: rod_radius, tip_radius, cone_length_frac
    real(dp), intent(in) :: r_curvature
    real(dp)             :: lsf
    real(dp)             :: dist_vec(NDIM), frac, radius_at_height, tmp

    ! Project onto line from r0 to r1
    call GM_dist_vec_line1(r, r0, r1, NDIM, dist_vec, frac)

    if (frac <= 1 - cone_length_frac) then
       ! Cylindrical part
       lsf = norm2(dist_vec) - rod_radius
    else if (frac < 1.0_dp) then
       ! Conical part
       tmp = (1 - frac) / cone_length_frac ! between 0 and 1
       radius_at_height = tip_radius + tmp * (rod_radius - tip_radius)
       lsf = norm2(dist_vec) - radius_at_height
    else
       ! Spherical tip
       lsf = norm2(r - cone_tip_center) - r_curvature
    end if
  end function conical_rod_lsf_arg

  real(dp) function conical_rod_lsf(r)
    real(dp), intent(in) :: r(NDIM)

    conical_rod_lsf = conical_rod_lsf_arg( &
         r, &
         field_rod_r0 * domain_len, &
         field_rod_r1 * domain_len, &
         cone_tip_center, &
         field_rod_radius, &
         cone_tip_radius, &
         cone_length_frac, &
         cone_tip_r_curvature)

  end function conical_rod_lsf
  !> ===========================================================================

  !> Compute electric field from electrical potential
  subroutine field_from_potential(box)
    type(box_t), intent(inout) :: box
    integer                     :: nc
    real(dp)                    :: inv_dr(NDIM)

    nc     = box%n_cell
    inv_dr = 1 / box%dr

#if NDIM == 2
    box%fc(1:nc+1, 1:nc, 1, ifc_E) = inv_dr(1) * &
         (box%cc(0:nc, 1:nc, i_phi) - box%cc(1:nc+1, 1:nc, i_phi))
    box%fc(1:nc, 1:nc+1, 2, ifc_E) = inv_dr(2) * &
         (box%cc(1:nc, 0:nc, i_phi) - box%cc(1:nc, 1:nc+1, i_phi))
#elif NDIM == 3
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
    box%cc(1:nc, 1:nc, i_E) = 0.5_dp * sqrt(&
         (box%fc(1:nc, 1:nc, 1, ifc_E) + &
         box%fc(2:nc+1, 1:nc, 1, ifc_E))**2 + &
         (box%fc(1:nc, 1:nc, 2, ifc_E) + &
         box%fc(1:nc, 2:nc+1, 2, ifc_E))**2)
#elif NDIM == 3
    box%cc(1:nc, 1:nc, 1:nc, i_E) = 0.5_dp * sqrt(&
             (box%fc(1:nc, 1:nc, 1:nc, 1, ifc_E) + &
              box%fc(2:nc+1, 1:nc, 1:nc, 1, ifc_E))**2 + &
             (box%fc(1:nc, 1:nc, 1:nc, 2, ifc_E) + &
              box%fc(1:nc, 2:nc+1, 1:nc, 2, ifc_E))**2 + &
             (box%fc(1:nc, 1:nc, 1:nc, 3, ifc_E) + &
              box%fc(1:nc, 1:nc, 2:nc+1, 3, ifc_E))**2)
#endif
  end subroutine compute_field_norm

end module m_field
