#include "../afivo/src/cpp_macros.h"
!> Program to perform discharge simulations with particle-in-cell method

program test_photoemission


  use omp_lib
  use m_af_all
  use m_globals
  use m_domain
  use m_refine
  use m_user
  use m_user_methods

  implicit none

integer, parameter     :: ndim = NDIM

  ! Read command line arguments and configuration files
  call CFG_update_from_arguments(cfg)

  ! Initialize the user's code
  call user_initialize(cfg)

  ! Initialize other modules
  call domain_init(cfg)
  call refine_init(cfg, ndim)
  call GL_initialize(cfg, ndim)

  ! Initialize the tree (which contains all the mesh information)
  call init_tree(tree)

  ! Set the dielectric permittivity before initializing multigrid
  if (GL_use_dielectric) then
     if (.not. associated(user_set_dielectric_eps)) &
          error stop "user_set_dielectric_eps not defined"
     call af_loop_box(tree, user_set_dielectric_eps)
  end if

  if (GL_use_dielectric) then

     ! Initialize dielectric surfaces at the fifth refinement level
     call af_refine_up_to_lvl(tree, 5)
     call dielectric_initialize(tree, i_eps, diel, 1)
  end if

call random_photoemission_event()

contains

    !> Initialize the AMR tree
    subroutine init_tree(tree)
      type(af_t), intent(inout) :: tree

      ! Initialize tree
      if (GL_cylindrical) then
         call af_init(tree, box_size, domain_len, &
              coarse_grid_size, coord=af_cyl)
      else
         call af_init(tree, box_size, domain_len, &
              coarse_grid_size)
      end if

    end subroutine init_tree


  subroutine random_photoemission_event()

    real(dp) :: start(2), final(2)
    real(dp) :: x_start(2), x_stop(2)
    logical  :: on_surface

    ! Generate two random numbers between [0.25, 0.75]

    start(1) = rand(0) * 0.5 + 0.25
    start(2) = rand(0) * 0.5 + 0.25
    x_start = start * domain_len

    ! Generate two random numbers between [-2, 2]
    final(1) = rand(0) * 4 - 2
    final(2) = rand(0) * 4 - 2

    x_stop = final * domain_len

    print *, "=========================="
    print *, " "
    print *, "start coordinates: ", start
    print *, "stop coordinates: ", final
    print *, " "


    call bisect_line(tree, x_start, x_stop, on_surface, i_eps)

    print *, "Calculated x_start: ", x_start / domain_len
    print *, "Calculated x_stop: ", x_stop / domain_len
    print *, "on_surface = ", on_surface
    print *, " "
    print *, "=========================="

  end subroutine


end program
