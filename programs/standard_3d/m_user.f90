!> Template for user code, this one simply uses the default routines
module m_user
  use m_config
  use m_user_methods

  implicit none
  private

  real(dp)  :: seed_pos(3) = [0.5_dp, 0.5_dp, 0.5_dp]
  integer   :: seed_num_particles = 10000
  real(dp)  :: seed_particle_weight = 1e4
  real(dp)  :: seed_sigma = 1e-4_dp

  logical   :: beam_enabled = .false.
  real(dp)  :: beam_current_A = 0.0_dp
  real(dp)  :: beam_energy_ev = 1000.0_dp
  real(dp)  :: beam_beta_z = 0.0_dp
  real(dp)  :: beam_particle_weight = 1e4_dp
  real(dp)  :: beam_radius = 1e-4_dp
  real(dp)  :: beam_center_rel(2) = [0.5_dp, 0.5_dp]
  real(dp)  :: beam_z0_rel = 0.01_dp
  real(dp)  :: beam_divergence_rad = 0.0_dp

  ! Public methods
  public :: user_initialize

contains

  subroutine user_initialize(cfg)
    type(CFG_t), intent(inout) :: cfg

    call CFG_add_get(cfg, "seed_pos", seed_pos, &
         "relative position of initial seed")
    call CFG_add_get(cfg, "seed_sigma", seed_sigma, &
         "characteristic size of the initial seed")
    call CFG_add_get(cfg, "seed_num_particles", seed_num_particles, &
         "number of particles in the seed")
    call CFG_add_get(cfg, "seed_particle_weight", seed_particle_weight, &
         "weight of the particles in the seed")

    user_initial_particles => init_particles

    call CFG_add_get(cfg, "beam%enabled", beam_enabled, &
         "Whether to inject an electron beam")
    call CFG_add_get(cfg, "beam%current_A", beam_current_A, &
         "Beam current in Ampere (electron current magnitude)")
    call CFG_add_get(cfg, "beam%energy_ev", beam_energy_ev, &
         "Beam kinetic energy in eV (directed +z)")
    call CFG_add_get(cfg, "beam%beta_z", beam_beta_z, &
         "Directed z-velocity as fraction of c (overrides energy_ev if nonzero)")
    call CFG_add_get(cfg, "beam%particle_weight", beam_particle_weight, &
         "Macro-particle weight for injected beam electrons")
    call CFG_add_get(cfg, "beam%radius", beam_radius, &
         "Beam radius in meters (uniform disk)")
    call CFG_add_get(cfg, "beam%center_rel", beam_center_rel, &
         "Beam center in relative x,y coordinates")
    call CFG_add_get(cfg, "beam%z0_rel", beam_z0_rel, &
         "Injection plane z in relative coordinate")
    call CFG_add_get(cfg, "beam%divergence_rad", beam_divergence_rad, &
         "Transverse RMS divergence in radians")

    if (beam_enabled) user_generate_particles => generate_beam
  end subroutine user_initialize

  subroutine init_particles(pctest)
    use m_particle_core
    use m_domain
    use m_globals
    use m_af_all
    type(PC_t), intent(inout) :: pctest
    integer                   :: n
    real(dp)                  :: pos(3)
    type(PC_part_t)           :: part

    part%v      = 0.0_dp
    part%a      = 0.0_dp
    part%t_left = 0.0_dp

    do n = 1, seed_num_particles
       pos(1:3) = seed_pos * domain_len
       part%w   = seed_particle_weight
       part%x(1:2) = pos(1:2) + GL_rng%two_normals() * seed_sigma
       part%x(3) = pos(3) + GL_rng%normal() * seed_sigma
       part%id = af_get_id_at(tree, part%x(1:3))
       if (associated(pctest%accel_function)) part%a = pctest%accel_function(part)

       if (outside_check(part) <= 0) then
          call pctest%add_part(part)
       end if
    end do
  end subroutine init_particles

  subroutine generate_beam(pc, time, time_elapsed)
    use m_units_constants
    use m_particle_core, only: PC_t, PC_part_t
    use m_globals,       only: GL_rng, tree
    use m_time_step,     only: GL_dt
    use m_domain,        only: domain_len, outside_check
    use m_af_all,        only: af_get_id_at
    type(PC_t), intent(inout) :: pc
    real(dp), intent(in)      :: time
    real(dp), intent(in)      :: time_elapsed

    integer         :: n, n_new
    real(dp)        :: dt_inj, mean_macro
    real(dp)        :: center(3), v0, sigma_v
    real(dp)        :: u1, u2, r, theta, c, s
    real(dp)        :: xy(2)
    type(PC_part_t) :: part

    if (.not. beam_enabled) return
    if (beam_current_A <= 0.0_dp) return
    if (beam_particle_weight <= 0.0_dp) return
    if (abs(beam_beta_z) <= 0.0_dp) then
       if (beam_energy_ev <= 0.0_dp) return
    else
       if (abs(beam_beta_z) >= 1.0_dp) error stop "beam%beta_z must satisfy |beta_z| < 1"
    end if

    dt_inj = merge(time_elapsed, GL_dt, time_elapsed > 0.0_dp)
    if (dt_inj <= 0.0_dp) return
    mean_macro = (beam_current_A * dt_inj) / (UC_elem_charge * beam_particle_weight)
    if (mean_macro <= 0.0_dp) return

    if (mean_macro < 1.0e6_dp) then
       n_new = GL_rng%poisson(mean_macro)
    else
       n_new = nint(mean_macro)
    end if
    if (n_new <= 0) return

    center(1:2) = beam_center_rel * domain_len(1:2)
    center(3) = beam_z0_rel * domain_len(3)

    if (abs(beam_beta_z) > 0.0_dp) then
       v0 = beam_beta_z * UC_lightspeed
    else
       v0 = sqrt(2.0_dp * beam_energy_ev * UC_elec_volt / UC_elec_mass)
    end if
    sigma_v = abs(v0) * max(0.0_dp, beam_divergence_rad)

    part%a      = 0.0_dp
    part%t_left = 0.0_dp
    part%w      = beam_particle_weight

    do n = 1, n_new
       u1 = GL_rng%unif_01()
       u2 = GL_rng%unif_01()
       r = beam_radius * sqrt(u1)
       theta = 2.0_dp * UC_pi * u2
       c = cos(theta)
       s = sin(theta)
       part%x(1) = center(1) + r * c
       part%x(2) = center(2) + r * s
       part%x(3) = center(3)

       if (sigma_v > 0.0_dp) then
          xy = GL_rng%two_normals()
          part%v(1) = sigma_v * xy(1)
          part%v(2) = sigma_v * xy(2)
       else
          part%v(1:2) = 0.0_dp
       end if
       part%v(3) = v0

       part%id = af_get_id_at(tree, part%x(1:3))
       if (associated(pc%accel_function)) part%a = pc%accel_function(part)

       if (outside_check(part) <= 0) then
          call pc%add_part(part)
       end if
    end do
  end subroutine generate_beam

end module m_user
