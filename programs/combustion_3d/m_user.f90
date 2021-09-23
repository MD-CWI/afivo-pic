!> Template for user code, this one simply uses the default routines
module m_user
  use m_config
  use m_user_methods
  use m_globals
  use m_domain

  implicit none
  private

  real(dp)  :: seed_pos(3) = [0.5, 0.5, 0.825]
  real(dp) :: seed_sigma = 1e-4_dp
  integer :: seed_num_particles = 10000
  real(dp) :: seed_particle_weight = 1e4

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
    user_init_CAS_array => init_CAS_array_IBK
  end subroutine user_initialize

  subroutine init_particles(pc)
    use m_particle_core
    type(PC_t), intent(inout) :: pc
    integer                   :: n
    ! real(dp)                  :: tmp_vec(2)
    type(PC_part_t)           :: part

    part%v      = 0.0_dp
    part%a      = 0.0_dp
    part%t_left = 0.0_dp

    n = seed_num_particles

    do while(n > 0)
       part%x(1:3) = GL_rng%three_normals() * seed_sigma + seed_pos * domain_len
       part%w = seed_particle_weight

       if (outside_check(part) <= 0) then
          call pc%add_part(part)
          n = n-1
       end if
    end do
  end subroutine init_particles

  subroutine init_CAS_array_IBK(cIx_groups, cIx_labels, n_species)
    real(dp), allocatable, intent(inout) :: cIx_groups(:, :)
    character(len=GL_slen), allocatable, intent(inout) :: cIx_labels(:)
    integer, intent(in)  :: n_species ! maximum number of species
    integer :: n_groups

    n_groups = 6
    print *, n_species
    allocate(cIx_groups(n_groups, n_species), cIx_labels(n_groups))
    cIx_groups = 0.0_dp ! Initialize zero weight

    !> INSERT GROUP NAMES HERE
    cIx_labels(1) = "N2_triplets"
    cIx_labels(2) = "O2_singlets"
    cIx_labels(3) = "H"
    cIx_labels(4) = "H2"
    cIx_labels(5) = "N"
    cIx_labels(6) = "O"

    !> INSERT WEIGHTS HERE (First column => group, second column => cIx)
    cIx_groups(1, 1:8)   = 1.0_dp ! Nitrogen triplets

    cIx_groups(2, 42:43) = 1.0_dp ! Oxygen singlets

    cIx_groups(3, 52) = 1.0_dp ! Hydrogen atoms
    cIx_groups(3, 55) = 1.0_dp ! Hydrogen atoms
    cIx_groups(3, 59) = 1.0_dp ! Hydrogen atoms
    cIx_groups(3, 62) = 1.0_dp ! Hydrogen atoms
    cIx_groups(3, 66) = 1.0_dp ! Hydrogen atoms
    cIx_groups(3, 67) = 1.0_dp ! Hydrogen atoms

    cIx_groups(4, 53) = 1.0_dp ! Hydrogen molecules
    cIx_groups(4, 54) = 1.0_dp ! Hydrogen molecules
    cIx_groups(4, 59) = 1.0_dp ! Hydrogen molecules
    cIx_groups(4, 60) = 2.0_dp ! Hydrogen molecules
    cIx_groups(4, 61) = 1.0_dp ! Hydrogen molecules
    cIx_groups(4, 64) = 1.0_dp ! Hydrogen molecules
    cIx_groups(4, 65) = 2.0_dp ! Hydrogen molecules
    cIx_groups(4, 67) = 1.0_dp ! Hydrogen molecules

    cIx_groups(5, 30) = 2.0_dp ! Dissociated nitrogen
    cIx_groups(5, 34) = 2.0_dp ! Dissociated nitrogen
    cIx_groups(5, 35) = 2.0_dp ! Dissociated nitrogen

    cIx_groups(6, 48) = 2.0_dp ! Dissociated oxygen
    cIx_groups(6, 50) = 2.0_dp ! Dissociated oxygen
    cIx_groups(6, 51) = 2.0_dp ! Dissociated oxygen


  end subroutine init_CAS_array_IBK

end module m_user
