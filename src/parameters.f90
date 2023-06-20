module class_parameters
  use class_mesh
  real, dimension(1:8,3)  :: nodes2_coordinates=reshape( & 
                             (/ 1, 1, 1, &
                               -1, 1, 1, &
                               -1,-1, 1, &
                                1,-1, 1, &
                                1, 1,-1, &
                               -1, 1,-1, &
                               -1,-1,-1, &
                                1,-1,-1 /), &
                            shape(nodes2_coordinates), order=(/2,1/) )
! ++++++++++++++++ GAUSS INTEGRATION VARIABLES ++++++++++++++++++++++++++++++++++++++
!  integer,parameter              :: Gauss_int_p=3
!  real, dimension(1:Gauss_int_p) :: Gauss_x_i=(/0.7745966692414834,0.,-0.7745966692414834/)
!  real, dimension(1:Gauss_int_p) :: Gauss_w_i=(/0.5555555555555556,0.8888888888888888, 0.5555555555555556/)

  integer,parameter              :: Gauss_int_p=5
  real, dimension(1:Gauss_int_p) :: Gauss_x_i=(/-1./3.*sqrt(5.-2.*sqrt(10./7.)), 1./3.*sqrt(5.-2.*sqrt(10./7.)), 0. , &
                                                -1./3.*sqrt(5.+2.*sqrt(10./7.)), 1./3.*sqrt(5.+2.*sqrt(10./7.))/)
  real, dimension(1:Gauss_int_p) :: Gauss_w_i=(/(322.+13.*sqrt(70.))/900.,(322.+13.*sqrt(70.))/900., 128./225. , &
                                                (322.-13.*sqrt(70.))/900.,(322.-13.*sqrt(70.))/900. /)

!  integer,parameter            :: Gauss_int_p=15
!  real, dimension(1:Gauss_int_p) :: Gauss_w_i= &
!                            (/0.03075324,0.07036605,0.10715922,0.13957068,0.16626921,0.18616100,0.19843149 &
!                              ,0.20257824,0.19843149,0.18616100,0.16626921,0.13957068,0.10715922,0.07036605,0.03075324/)
!  real, dimension(1:Gauss_int_p) :: Gauss_x_i= &
!                            (/-0.98799252,-0.93727339,-0.84820658,-0.72441773,-0.57097217,-0.39415135,-0.20119409 &
!                              ,0.00000000,0.20119409,0.39415135,0.57097217,0.72441773,0.84820658,0.93727339,0.98799252/)
  real,parameter                       :: MATH_PI=3.14159265359
! =======================================================================================
  integer,parameter:: left=1
  integer,parameter:: right=2
  integer,parameter:: back=3
  integer,parameter:: front=4
  integer,parameter:: bottom=5
  integer,parameter:: top=6

! =======================================================================================

  integer,parameter:: velx1=1
  integer,parameter:: velx2=2
  integer,parameter:: velx3=3
  integer,parameter:: PPE=4
  integer,parameter:: CDE=5


! =======================================================================================
  integer,parameter:: no_slip=0
  integer,parameter:: constant_bc=1
  integer,parameter:: dirichlet=2
  integer,parameter:: zero_gradient=3
  integer,parameter:: neumann=4


! =======================================================================================

  integer,parameter:: analytical_velx1=1
  integer,parameter:: analytical_velx2=2
  integer,parameter:: analytical_velx3=3
  integer,parameter:: analytical_p=4
  integer,parameter:: analytical_CDE=101
  integer,parameter:: IC_velx1=5
  integer,parameter:: IC_velx2=6
  integer,parameter:: IC_velx3=7
  integer,parameter:: IC_PPE=8
  integer,parameter:: IC_CDE=102
  integer,parameter:: bx_fun=10
  integer,parameter:: by_fun=11
  integer,parameter:: bz_fun=12
  integer,parameter:: q_fun=103

! =======================================================================================
	integer,parameter:: bc_dvdn   = 10
	integer,parameter:: bc_dvdnt  = 11
	integer,parameter:: bc_dvdns  = 12
	integer,parameter:: bc_dvdntt = 13
	integer,parameter:: bc_dvdnss = 14
	integer,parameter:: bc_dvdnts = 15

! =======================================================================================
! Geometry variables
!  integer,parameter                    :: NMESH_i=64
!  integer,parameter                    :: NMESH_j=NMESH_i
!  integer,parameter                    :: NMESH_k=NMESH_i
  integer,parameter                    :: NMESH_j=6
  integer,parameter                    :: NMESH_i=6*NMESH_j
  integer,parameter                    :: NMESH_k=2*NMESH_j
  real,parameter                       :: R_i=1.0*(MATH_PI/3.0)!*sqrt(.5)
  real,parameter                       :: R_o=2.0
  real,parameter                       :: L=2.
  real,parameter                       :: Lx=1.
  real,parameter                       :: Ly=Lx
  real,parameter                       :: Lz=Lx
  real,parameter                       :: startx=0.
  real,parameter                       :: starty=startx
  real,parameter                       :: startz=startx

! For hollow: NMESH_i=nTheta, NMESH_j=nR, NMESH_k=Nz

  
!  integer,parameter                        :: NMESH=NMESH_i*NMESH_j*NMESH_k
!  type(mesh),dimension(0:NMESH)            :: MESHES_ARRAY
!  integer,parameter                        :: NFACE=(NMESH_i+1)*NMESH_j*NMESH_k+NMESH_i*(NMESH_j+1)*NMESH_k+NMESH_i*NMESH_j*(NMESH_k+1)
!  type(face),dimension(0:NFACE)            :: FACES_ARRAY
!  integer,parameter                        :: NNODE=(NMESH_i+1)*(NMESH_j+1)*(NMESH_k+1)
!  type(node),dimension(0:NNODE)            :: NODES_ARRAY
!  integer,parameter                        :: NBC=7
!  type(boundary_condition),dimension(1:NBC):: BCS_ARRAY
!  integer,dimension(1:NBC)                 :: BC_id


! Hollow ---< uncomment
    integer,parameter                        :: NMESH=NMESH_i*NMESH_j*NMESH_k
    type(mesh),dimension(0:NMESH)            :: MESHES_ARRAY
    integer,parameter                        :: NFACE=3*NMESH_i*NMESH_j*NMESH_k+(NMESH_i*NMESH_k)+(NMESH_i*NMESH_j)
    type(face),dimension(0:NFACE)            :: FACES_ARRAY
    integer,parameter                        :: NNODE=(NMESH_i)*(NMESH_j+1)*(NMESH_k+1)
    type(node),dimension(0:NNODE)            :: NODES_ARRAY
    integer,parameter                        :: NBC=5
    type(boundary_condition),dimension(1:NBC):: BCS_ARRAY
    integer,dimension(1:NBC)                 :: BC_id


! =======================================================================================

  real,parameter                           :: v_limit=1E-2   ! Limit where v is considered zero
  real,parameter                           :: f_limit=0.0 ! limit where F is considered zero

! =======================================================================================

! Time-dependent parameters
  real,parameter                 :: DT=0.01         ! This global variable is the global time
  real                           :: TI=0.
  integer,parameter              :: NT=20000
  integer                        :: STEP=0
  integer                        :: STEP_OUT=100 ! dump every STEP_OUT steps

! RMS variables
  real                           :: RMS_time=0.
  integer                        :: RMS_time_ctr=0

! Solver paramters

  real,parameter                 :: Tolerance_u=1E-6
  real,parameter                 :: Tolerance_v=Tolerance_u
  real,parameter                 :: Tolerance_w=Tolerance_u
  real,parameter                 :: Tolerance_p=5E-3
  real,parameter                 :: Tolerance_T=1E-6
  real,parameter                 :: dudt_epsilon=1E-5
  real,parameter                 :: w_u=.4
  real,parameter                 :: w_v=w_u
  real,parameter                 :: w_w=w_u
  real,parameter                 :: w_p=w_u
  real,parameter                 :: w_p2=w_p
  real,parameter                 :: w_T=w_u
  integer,parameter              :: sweep_u=1
  integer,parameter              :: sweep_p=1
  integer,parameter              :: sweep_T=1

! =============================== Problem specification ===============================
  real,parameter                           :: rho=1.0
!  real,parameter                           :: Re=100.0
!  real,parameter                           :: nu=1./Re
  real,parameter                           :: Pr=0.7				! Pr number
  real,parameter                           :: Ra=1.0E5
  
  real,parameter                           :: th_tc=1.
  real,parameter                           :: g_beta=1.
  real,parameter                           :: nu=1./sqrt(Ra/(g_beta*th_tc*Pr))
  real,parameter                           :: Re=1./nu
  real,parameter                           :: alpha=nu/Pr
! =======================================================================================

! ================================ Variables may save computational time ================
  logical,parameter :: steady_state=.false.
  logical,parameter :: is_energy=.true.
  
  logical,parameter :: is_bx=.false.
  logical,parameter :: is_by=.false.
  logical,parameter :: is_bz=.false.
  logical,parameter :: is_q =.false.

  logical,parameter :: is_bx_TD=.false.     ! if bx time-dependent
  logical,parameter :: is_by_TD=.false.     ! if by time-dependent
  logical,parameter :: is_bz_TD=.false.     ! if bz time-dependent
  
  logical,parameter :: is_analytical=.false. ! if analytical solution is provided
  logical,parameter :: is_analytical_TD=.false.
  
  logical,parameter :: is_dpdx_tolerance=.false.
  logical,parameter :: is_check_steady_state=.true.
  integer,parameter :: nim_formulation=2
  
! ======================================================================================


! =============================== TIMING PARAMETERS ====================================
  real*8              :: timing_start
  real*8              :: timing_initialize_1
  real*8              :: timing_initialize_2
  real*8              :: timing_integration
  real*8              :: timing_integration_step, timing_integration_step_start
  real*8              :: timing_all
  real*8              :: initialization_1_time           ! Reading and preparing meshes
  real*8              :: initialization_2_time           ! Initialize the problem: bc, ic, .....
  real*8              :: integration_time                ! NSE integration step timing
  real*8              :: simulation_time                 ! Overall Simulation time
  
   ! integer,parameter :: NMESH2=4096
   ! integer,parameter :: NFACE2=13056
   ! integer,parameter :: NNODE2=4913
  ! ! integer,parameter :: NMESH2=32768
  ! ! integer,parameter :: NFACE2=101376
  ! ! integer,parameter :: NNODE2=35937
  ! integer,parameter :: NBC2=6
  
! =============================== Global Arrays Definitions =============================  
! Arrays related to meshes
! Unknowns variables
  real,dimension(1:NMESH)              :: meshes_uxyz_old, meshes_uxyz, meshes_uxyzt, & 	! Volume averaged variables
                                           meshes_vxyz_old, meshes_vxyz, meshes_vxyzt, &
										   meshes_wxyz_old, meshes_wxyz, meshes_wxyzt, &
										   meshes_u, meshes_v, meshes_w, meshes_p,	&		! Previous time step variables -> converted to curvilinear
										   meshes_Txyz_old, meshes_Txyz, meshes_Txyzt
  real,dimension(1:NMESH)              :: meshes_uxyz_old_con,meshes_vxyz_old_con,meshes_wxyz_old_con, meshes_Txyz_old_con
! Sources in NSE -> RC, RD, RB
  real,dimension(1:NMESH)              :: meshes_bx, meshes_by, meshes_bz,meshes_q, & 
                                           meshes_R_C_u, meshes_R_B_u, meshes_R_D_u, &
                                           meshes_R_C_v, meshes_R_B_v, meshes_R_D_v, & 
										   meshes_R_C_w, meshes_R_B_w, meshes_R_D_w, &
										   meshes_R_B_p, meshes_R_D_p,meshes_R_B_T, meshes_R_D_T
										   
  real,dimension(1:NMESH)              :: meshes_uT, meshes_vT, meshes_wT, meshes_pT, meshes_TT 	 		! Analytical values of vel and p

  real,dimension(1:NMESH)              :: meshes_dila, meshes_dila_old						! residual from the continuity equation


! NIM formulation 2
  real,dimension(1:3,1:6,1:NMESH)      :: M_A, MT_A
  real,dimension(1:5,1:NMESH)          :: M_F1, M_F2, M_F3, MT_F1, MT_F2, MT_F3
  real,dimension(1:8,1:NMESH)          :: M_F4, MT_F4


  
! Discrete equations related variables
  real,dimension(1:15,1:NMESH)         :: AA_vx, AA_vy, AA_vz	 	! Arrays of the discrete coefficients - v
  real,dimension(1:15,1:NMESH)         :: AA_px, AA_py, AA_pz		! Arrays of the discrete coefficients - p
  real,dimension(1:9,1:NMESH)          :: a_t, aT_t					! Array of the time-dependent velocity

! Curvilinear related variables:
  real,dimension(1:3,1:3,1:NMESH)      :: meshes_xi_x, meshes_gij  
  real,dimension(1:3,1:3,1:6,1:NMESH)  :: meshes_gij_faces
  real,dimension(1:NMESH)              :: meshes_G1,meshes_G2,meshes_G3
  real,dimension(1:NMESH)              :: meshes_volume
  real,dimension(1:3,1:6,1:NMESH)      :: meshes_normal_faces, meshes_tangential_faces, meshes_binormal_faces, meshes_mid_faces

! Geometry related variables
  integer,dimension(1:6,1:NMESH)       :: meshes_neighbors
  integer,dimension(1:6,1:NMESH)       :: meshes_faces
  real,dimension(1:NMESH)              :: dpdx,dpdx_old,dpdy,dpdy_old,dpdz,dpdz_old
!  real,dimension(1:8,1:NMESH)          :: meshes_nodes_x,meshes_nodes_y,meshes_nodes_z
! ==========================================================================================
! Arrays related to faces

! Unknowns
  real,dimension(1:NFACE)             :: faces_u, faces_v, faces_w, faces_p, faces_T				! Area averaged variables
  real,dimension(1:NFACE)             :: faces_uT, faces_vT, faces_wT, faces_pT, faces_TT			! Analytical solution over faces
  real,dimension(1:NFACE)             :: faces_u_old, faces_v_old, faces_w_old, faces_p_old, faces_T_old
  
  real,dimension(1:3,1:NFACE)         :: faces_normal				! Normal vectors
  real,dimension(1:3,1:NFACE)         :: faces_tangential, faces_bitangential
  real,dimension(1:3,1:NFACE)         :: faces_mid_coordinates	! mid-point coordinates
  real,dimension(1:NFACE)             :: faces_bx, faces_by, faces_bz, faces_q
  real,dimension(1:NFACE)             :: faces_area  
  
  integer,dimension(1:2,1:NFACE)      :: faces_shared_meshes
  integer,dimension(1:NFACE)          :: faces_BC					! Boundary condition ID
  
  
  integer,dimension(:),allocatable     :: BCs_faces, BCs_faces_id   	! All faces that have boundary conditions, id:BC id
  integer,dimension(:),allocatable     :: BCs_faces_position, BCs_faces_mesh_num			! Position of the face in the mesh and number of mesh
  integer,dimension(:,:),allocatable   :: BCs_faces_points									! Store faces numbers for all points to be used in FDM
  real,dimension(:,:) ,allocatable     :: BCs_faces_ap, BCs_faces_a,BCs_faces_aT		! ap=pressure discrete constants, a=finite difference constants
  
  integer                              :: tot_faces_bcs
  real,dimension(:),allocatable        :: BCs_faces_v_n, BCs_faces_v_nt, BCs_faces_v_ns,BCs_faces_v_ntt,BCs_faces_v_nss,BCs_faces_v_nts	! Derivatives with respect to the tangential on the boundaries
! ==========================================================================================
  real                                 :: Dilatation
  real                                 :: dudt_max,dTdt_max				! to check the steady-state solution
  real                                 :: RMS_meshes_u, RMS_meshes_v, RMS_meshes_w, RMS_faces_u, RMS_faces_v, RMS_faces_w
! ================================ Restart =================================================
	! Restart file specifications
	logical,parameter               :: is_write_restart=.true.
	logical,parameter               :: is_restart=.false.
	integer,parameter               :: write_restart_step=1000
	integer,parameter               :: read_restart_step=20000
end module class_parameters
