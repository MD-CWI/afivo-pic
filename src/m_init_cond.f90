#include "afivo/src/cpp_macros.h"
module m_init_cond
  use m_globals
  use m_af_all
  use m_domain

  implicit none
  private

  ! Type to store initial conditions in
  type initcnd_t
     real(dp)                        :: background_density
     real(dp)                        :: stochastic_density
     integer                         :: n_cond
     real(dp), allocatable           :: seed_r0(:, :)
     real(dp), allocatable           :: seed_r1(:, :)
     real(dp), allocatable           :: seed_density(:)
     integer, allocatable            :: seed_charge_type(:)
     real(dp), allocatable           :: seed_width(:)
     character(ST_slen), allocatable :: seed_falloff(:)
  end type initcnd_t

  ! This will contain the initial conditions
  type(initcnd_t), protected, public :: init_conds

  integer               :: num_particle_seeds
  integer, allocatable  :: seed_counts(:)
  real(dp), allocatable :: seed_weights(:)
  real(dp), allocatable :: seed_sigmas(:)
  real(dp), allocatable :: seed_pos(:, :)
  real(dp)              :: background_density = 0.0_dp

  public :: init_cond_initialize
  public :: init_cond_particles
  public :: init_cond_set_box

contains

  ! Set the initial conditions from the configuration
  subroutine init_cond_initialize(cfg, n_dim)
    type(CFG_t), intent(inout) :: cfg !< Settings
    integer, intent(in)        :: n_dim

    integer                    :: n_cond, varsize
    real(dp), allocatable      :: tmp_vec(:)
    type(initcnd_t)            :: ic

    call CFG_add(cfg, "background_density", 0.0_dp, &
         "The background ion and electron density (1/m3)")
    call CFG_add(cfg, "seed%density", [0.0e15_dp], &
         "Initial density of the seed (1/m3)", .true.)
    call CFG_add(cfg, "seed%rel_r0", [DTIMES(0.5_dp)], &
         "The relative start position of the initial seed", .true.)
    call CFG_add(cfg, "seed%rel_r1", [DTIMES(0.5_dp)], &
         "The relative end position of the initial seed", .true.)
    call CFG_add(cfg, "seed%charge_type", [0], &
         "Type of seed: neutral (0), ions (1) or electrons (-1)", .true.)
    call CFG_add(cfg, "seed%width", [0.25d-3], &
         "Seed width (m)", .true.)
    call CFG_add(cfg, "seed%falloff", ['gaussian'], &
         "Fall-off type for seed (gaussian, step), default=gaussian", .true.)

    call CFG_get_size(cfg, "seed%density", n_cond)
    ic%n_cond = n_cond

    call CFG_get_size(cfg, "seed%rel_r0", varsize)
    if (varsize /= n_dim * n_cond) &
         stop "seed%rel_r0 variable has incompatible size"

    call CFG_get_size(cfg, "seed%rel_r1", varsize)
    if (varsize /= n_dim * n_cond) &
         stop "seed%rel_r1 variable has incompatible size"

    call CFG_get_size(cfg, "seed%charge_type", varsize)
    if (varsize /= n_cond) &
         stop "seed%charge_type variable has incompatible size"

    call CFG_get_size(cfg, "seed%width", varsize)
    if (varsize /= n_cond) &
         stop "seed%width variable has incompatible size"

    allocate(ic%seed_density(n_cond))
    allocate(ic%seed_charge_type(n_cond))
    allocate(ic%seed_r0(n_dim, n_cond))
    allocate(ic%seed_r1(n_dim, n_cond))
    allocate(ic%seed_width(n_cond))
    allocate(ic%seed_falloff(n_cond))

    allocate(tmp_vec(n_dim * n_cond))
    call CFG_get(cfg, "seed%rel_r0", tmp_vec)
    ic%seed_r0 = domain_len * reshape(tmp_vec, [n_dim, n_cond])
    call CFG_get(cfg, "seed%rel_r1", tmp_vec)
    ic%seed_r1 = domain_len * reshape(tmp_vec, [n_dim, n_cond])

    call CFG_get(cfg, "background_density", ic%background_density)
    call CFG_get(cfg, "seed%density", ic%seed_density)
    call CFG_get(cfg, "seed%charge_type", ic%seed_charge_type)
    call CFG_get(cfg, "seed%width", ic%seed_width)
    call CFG_get(cfg, "seed%falloff", ic%seed_falloff)

    init_conds = ic

    call CFG_add_get(cfg, "pseed%background_density", background_density, &
         "Background density of electrons")
    call CFG_add(cfg, "pseed%counts", [100], &
         "Number of particles per seed", .true.)
    call CFG_get_size(cfg, "pseed%counts", n_cond)
    num_particle_seeds = n_cond

    allocate(seed_weights(n_cond))
    seed_weights = 1.0_dp
    allocate(seed_counts(n_cond))
    seed_counts = 100
    allocate(seed_sigmas(n_cond))
    seed_sigmas = 1e-4_dp
    allocate(seed_pos(n_dim, n_cond))
    deallocate(tmp_vec)
    allocate(tmp_vec(n_dim * n_cond))
    tmp_vec = 0.5_dp

    call CFG_get(cfg, "pseed%counts", seed_counts)
    call CFG_add_get(cfg, "pseed%weights", seed_weights, &
         "Weight of electrons for the Gaussian seeds", .true.)
    call CFG_add_get(cfg, "pseed%sigmas", seed_sigmas, &
         "Width of the Gaussian seeds", .true.)
    call CFG_add_get(cfg, "pseed%pos", tmp_vec, &
         "Relative location of the Gaussian seeds", .true.)
    seed_pos = reshape(tmp_vec * domain_len, [n_dim, n_cond])

  end subroutine init_cond_initialize

  subroutine init_cond_particles(tree, pc)
    use m_particle_core
    type(af_t), intent(inout) :: tree
    type(PC_t), intent(inout) :: pc
    integer                   :: n, i
    real(dp)                  :: pos(3)
    type(PC_part_t)           :: part

    part%v      = 0.0_dp
    part%a      = 0.0_dp
    part%t_left = 0.0_dp

    ! n_background = background_density * domain_len**2
    ! do i = 1, 

    ! end do

    do n = 1, num_particle_seeds
#if NDIM == 2
       pos(1:2)  = tree%r_base + seed_pos(1:2, n)
       pos(3)    = 0.0_dp
       part%x(3) = pos(3)
#elif NDIM == 3
       pos(1:3)  = tree%r_base + seed_pos(1:3, n)
#endif
       part%w    = seed_weights(n)

       do i = 1, seed_counts(n)
#if NDIM == 2
          part%x(1:2) = pos(1:2) + ST_rng%two_normals() * seed_sigmas(n)
#elif NDIM == 3
          ! TODO: avoid setting x(2) twice
          part%x(1:2) = pos(1:2) + ST_rng%two_normals() * seed_sigmas(n)
          part%x(2:3) = pos(2:3) + ST_rng%two_normals() * seed_sigmas(n)
#endif

          if (outside_check(part) <= 0) then
             call pc%add_part(part)
          end if
       end do
    end do
  end subroutine init_cond_particles

  !> Sets the initial condition
  subroutine init_cond_set_box(box)
    use m_geometry
    type(box_t), intent(inout) :: box
    integer                     :: IJK, n, nc
    real(dp)                    :: rr(NDIM)
    real(dp)                    :: density

    nc = box%n_cell
    box%cc(DTIMES(:), i_electron) = init_conds%background_density
    box%cc(DTIMES(:), i_pos_ion)  = init_conds%background_density
    box%cc(DTIMES(:), i_phi)      = 0 ! Inital potential set to zero

    do KJI_DO(0,nc+1)
       rr   = af_r_cc(box, [IJK])

       do n = 1, init_conds%n_cond
          density = init_conds%seed_density(n) * &
               GM_density_line(rr, init_conds%seed_r0(:, n), &
               init_conds%seed_r1(:, n), NDIM, &
               init_conds%seed_width(n), &
               init_conds%seed_falloff(n))

          ! Add electrons and/or ions depending on the seed charge type
          ! (positive, negative or neutral)
          if (init_conds%seed_charge_type(n) <= 0) then
             box%cc(IJK, i_electron) = box%cc(IJK, i_electron) + density
          end if

          if (init_conds%seed_charge_type(n) >= 0) then
             box%cc(IJK, i_pos_ion) = box%cc(IJK, i_pos_ion) + density
          end if
       end do
    end do; CLOSE_DO

  end subroutine init_cond_set_box

end module m_init_cond
