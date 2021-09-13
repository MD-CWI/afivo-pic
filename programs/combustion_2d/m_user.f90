!> Template for user code, this one simply uses the default routines
module m_user
  use m_config
  use m_user_methods
  use m_globals
  use m_domain

  implicit none
  private

  real(dp) :: seed_pos(2) = [0.5_dp, 0.5_dp]
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
    user_init_CAS_array => init_CAS_array
  end subroutine user_initialize

  ! subroutine init_CAS_array_old(cIx_to_track, cIx_weights, cIx_labels, num_cIx)
  !   integer, allocatable, intent(inout) :: cIx_to_track(:, :)
  !   real(dp), allocatable, intent(inout) :: cIx_weights(:, :)
  !   character(len=GL_slen), allocatable, intent(inout) :: cIx_labels(:)
  !   integer, intent(inout)  :: num_cIx
  !   integer :: max_proc
  !
  !   num_cIx = 2 ! The number of grouped CAS densities we include
  !   max_proc = 6 ! The maximum number of processes in one group (needed for allocation)
  !   allocate(cIx_to_track(num_cIx, max_proc), cIx_weights(num_cIx, max_proc), cIx_labels(num_cIx))
  !   cIx_to_track(:,:) = 0
  !   cIx_weights(:,:)  = 1 ! Default weight equals one
  !
  !   cIx_labels(1) = "N2_triplets"
  !   cIx_to_track(1, :) = [1, 2, 3, 4, 5, 6]
  !
  !   cIx_labels(2) = "O2_singlets"
  !   cIx_to_track(2, 1:2) = [23, 24]
  !
  !   print *, "CAS_weights not implemented."
  !   ! TODO Consider switching to a dictionary-type structure using pointers.
  ! end subroutine init_CAS_array_old

  subroutine init_CAS_array(cIx_groups, cIx_labels, n_species)
    real(dp), allocatable, intent(inout) :: cIx_groups(:, :)
    character(len=GL_slen), allocatable, intent(inout) :: cIx_labels(:)
    integer, intent(in)  :: n_species ! maximum number of species
    integer :: n_groups

    n_groups = 4
    print *, n_species
    allocate(cIx_groups(n_groups, n_species), cIx_labels(n_groups))
    cIx_groups = 0.0_dp ! Initialize zero weight

    !> INSERT GROUP NAMES HERE
    cIx_labels(1) = "N2_triplets"
    cIx_labels(2) = "O2_singlets"
    cIx_labels(3) = "H2"
    cIx_labels(4) = "H"

    !> INSERT WEIGHTS HERE (First column => group, second column => cIx)
    cIx_groups(1, 1:6)   = 1.0_dp ! Nitrogen triplets
    cIx_groups(2, 23:24) = 1.0_dp ! Oxygen singlets

    cIx_groups(3, 34) = 1.0_dp ! Hydrogen molecules
    cIx_groups(3, 37) = 1.0_dp ! Hydrogen molecules
    cIx_groups(3, 41) = 1.0_dp ! Hydrogen molecules
    cIx_groups(3, 44) = 1.0_dp ! Hydrogen molecules
    cIx_groups(3, 48) = 1.0_dp ! Hydrogen molecules
    cIx_groups(3, 49) = 1.0_dp ! Hydrogen molecules

    cIx_groups(4, 35) = 1.0_dp ! Hydrogen atoms
    cIx_groups(4, 36) = 1.0_dp ! Hydrogen atoms
    cIx_groups(4, 41) = 1.0_dp ! Hydrogen atoms
    cIx_groups(4, 42) = 2.0_dp ! Hydrogen atoms
    cIx_groups(4, 43) = 1.0_dp ! Hydrogen atoms
    cIx_groups(4, 46) = 1.0_dp ! Hydrogen atoms
    cIx_groups(4, 47) = 2.0_dp ! Hydrogen atoms
    cIx_groups(4, 49) = 1.0_dp ! Hydrogen atoms


  end subroutine init_CAS_array

  subroutine init_particles(pc)
    use m_particle_core
    type(PC_t), intent(inout) :: pc
    integer                   :: n
    real(dp)                  :: tmp_vec(2)
    type(PC_part_t)           :: part

    part%v      = 0.0_dp
    part%a      = 0.0_dp
    part%t_left = 0.0_dp

    do n = 1, seed_num_particles
       if (GL_cylindrical) then
          tmp_vec = GL_rng%two_normals() * seed_sigma
          part%x(1) = norm2(tmp_vec)
          tmp_vec = GL_rng%two_normals() * seed_sigma
          part%x(2) = tmp_vec(1)
          part%x(1:2) = part%x(1:2) + seed_pos * domain_len
          part%x(3) = 0.0_dp
       else
          part%x(1:2) = GL_rng%two_normals() * seed_sigma
          part%x(1:2) = part%x(1:2) + seed_pos * domain_len
          part%x(3) = 0.0_dp
       end if

       part%w = seed_particle_weight

       if (outside_check(part) <= 0) then
          call pc%add_part(part)
       end if
    end do
  end subroutine init_particles

end module m_user
