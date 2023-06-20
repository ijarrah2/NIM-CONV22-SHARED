include 'include.f90'
! =========================================================================================
! ================================ MAIN PROGRAM ===========================================
! =========================================================================================
program main
  use class_general
  use class_mesh
  use class_curvilinear
! use class_read_mesh    Deleted for simplicity
  use class_create_geometry
  use class_create_mesh
  use class_prepare_mesh
  use class_parameters
  use class_geometry
  use class_nim
  use class_bc
  use class_nse
  use omp_lib
  implicit none
  character(len = 40)  :: mesh_file
  integer              :: mesh_id
  real :: start, finish
  real :: initialization_time
!  mesh_file="mesh//file.msh"

  timing_start=omp_get_wtime()


! Read msh file format created using ICEM CFD
! Deleted from this version of code
!  call read_mesh_file(mesh_file)
!  call create_cube


! Create a hollow cylindrical mesh
  call create_hollow(R_i,R_o,L,0.,0.,L/2.,NMESH_i,NMESH_j,NMESH_k)

! Create links between elements based on NIM requirements
  call prepare_mesh


! Fill global arrays
  call fill_global_arrays


  timing_initialize_1=omp_get_wtime()

! Solve Navier-stokes equations
  call solve_nse



  timing_all=omp_get_wtime()

! Report time at end of calculations

  initialization_1_time=timing_initialize_1-timing_start
  initialization_2_time=timing_initialize_2-timing_initialize_1
  integration_time=timing_integration-timing_initialize_2
  simulation_time=timing_all-timing_start
  print '("============================================================================")'
  print '("========================== TIMING SUMMARY ==================================")'
  print '("============================================================================")'
  print '("Preparing mesh time = ",f10.5," seconds.")',initialization_1_time
  print '("Initialization time = ",f10.5," seconds.")',initialization_2_time
  print '("Integration    time = ",f10.5," seconds.")',integration_time
  print '("Simulation     time = ",f10.5," seconds.")',simulation_time
  print '("============================================================================")'
  print '("============================================================================")'


! Deallocate arrays allocated during simulation
  call deallocate_arrays

end program main
