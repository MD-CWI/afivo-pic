#include "../afivo/src/cpp_macros.h"
!> Module to help setting up initial conditions
module m_init_cond
  use m_globals
  use m_domain
  use m_af_all

  implicit none
  private

  ! Type to store initial conditions in
  type initcnd_t
     real(dp)                        :: background_density
     integer                         :: n_cond
     real(dp), allocatable           :: seed_r0(:, :)
     real(dp), allocatable           :: seed_r1(:, :)
     real(dp), allocatable           :: seed_density(:)
     real(dp), allocatable           :: seed_width(:)
     character(64), allocatable      :: seed_falloff(:)
  end type initcnd_t

  ! This will contain the initial conditions
  type(initcnd_t), protected, public :: init_conds

  public :: init_cond_initialize
  public :: init_cond_set_box

contains

  ! Set the initial conditions from the configuration
  subroutine init_cond_initialize(cfg)
    use m_config
    type(CFG_t), intent(inout) :: cfg !< Settings

    integer                    :: n, n_cond, varsize
    real(dp)                   :: empty_real(0)
    real(dp), allocatable      :: tmp_vec(:)
    character(len=name_len)    :: empty_string(0)
    type(initcnd_t)            :: ic

    call CFG_add(cfg, "background_density", 0.0_dp, &
         "The background ion and electron density (1/m3)")
    call CFG_add(cfg, "seed_density", empty_real, &
         "Initial density of the seed (1/m3)", .true.)
    call CFG_add(cfg, "seed_rel_r0", empty_real, &
         "The relative start position of the initial seed", .true.)
    call CFG_add(cfg, "seed_rel_r1", empty_real, &
         "The relative end position of the initial seed", .true.)
    call CFG_add(cfg, "seed_width", empty_real, &
         "Seed width (m)", .true.)
    call CFG_add(cfg, "seed_falloff", empty_string, &
         "Fall-off type for seed (sigmoid, gaussian, smoothstep, step, laser)", .true.)

    call CFG_get_size(cfg, "seed_density", n_cond)
    ic%n_cond = n_cond

    call CFG_get_size(cfg, "seed_rel_r0", varsize)
    if (varsize /= NDIM * n_cond) &
         stop "seed_rel_r0 variable has incompatible size"

    call CFG_get_size(cfg, "seed_rel_r1", varsize)
    if (varsize /= NDIM * n_cond) &
         stop "seed_rel_r1 variable has incompatible size"

    call CFG_get_size(cfg, "seed_width", varsize)
    if (varsize /= n_cond) &
         stop "seed_width variable has incompatible size"

    allocate(ic%seed_density(n_cond))
    allocate(ic%seed_r0(NDIM, n_cond))
    allocate(ic%seed_r1(NDIM, n_cond))
    allocate(ic%seed_width(n_cond))
    allocate(ic%seed_falloff(n_cond))

    allocate(tmp_vec(NDIM * n_cond))
    call CFG_get(cfg, "seed_rel_r0", tmp_vec)
    ic%seed_r0 = reshape(tmp_vec, [NDIM, n_cond])
    call CFG_get(cfg, "seed_rel_r1", tmp_vec)
    ic%seed_r1 = reshape(tmp_vec, [NDIM, n_cond])

    do n = 1, n_cond
       ic%seed_r0(:, n) = ic%seed_r0(:, n) * domain_len
       ic%seed_r1(:, n) = ic%seed_r1(:, n) * domain_len
    end do

    call CFG_get(cfg, "background_density", ic%background_density)
    call CFG_get(cfg, "seed_density", ic%seed_density)
    call CFG_get(cfg, "seed_width", ic%seed_width)
    call CFG_get(cfg, "seed_falloff", ic%seed_falloff)

    init_conds = ic

  end subroutine init_cond_initialize

  !> Sets the initial condition
  subroutine init_cond_set_box(box)
    use m_geometry
    use m_gas
    use m_user_methods
    type(box_t), intent(inout) :: box
    integer                    :: IJK, n, nc, ix, n_particles
    real(dp)                   :: rr(NDIM), w
    real(dp)                   :: density, pos(NDIM, ceiling(particle_per_cell))
    real(dp)                   :: volume
    type(PC_part_t)            :: new_part

    nc = box%n_cell
    volume = product(box%dr)

    do KJI_DO(1,nc)
       rr = af_r_cc(box, [IJK])
#if NDIM == 2
       if (GL_cylindrical) volume = af_cyl_volume_cc(box, i)
#endif

       density = init_conds%background_density
       do n = 1, init_conds%n_cond
          density = density + init_conds%seed_density(n) * &
               GM_density_line(rr, init_conds%seed_r0(:, n), &
               init_conds%seed_r1(:, n), NDIM, &
               init_conds%seed_width(n), &
               init_conds%seed_falloff(n))
       end do

       ! Determine number of particles to generate
       n_particles = nint(min(density * volume / particle_min_weight, &
            particle_per_cell))
       ! Determine weight
       w = density * volume / n_particles

       ! Sample particle positions within the cell
       if (GL_cylindrical) then
          ! @todo fix this with some complex formula in the future
          do ix = 1, n_particles
             pos(1, ix) = rr(1) + box%dr(1) * &
                  (sqrt(GL_rng%unif_01()) - 0.5) ! Not fully correct!
             pos(2, ix) = rr(2) + box%dr(2) * &
                  (sqrt(GL_rng%unif_01()) - 0.5_dp)
          end do
       else
          do ix = 1, n_particles
             pos(:, ix) = rr + box%dr * &
                  [DTIMES(GL_rng%unif_01() - 0.5_dp)]
          end do
       end if

       new_part%a = 0
       new_part%w = w

       do ix = 1, n_particles
          new_part%x(1:NDIM) = pos(:, ix)
          if (outside_check(new_part) == 0) then
             call pc%add_part(new_part)
          end if
       end do

    end do; CLOSE_DO

  end subroutine init_cond_set_box

end module m_init_cond
