! ====================================== CLASSES DEFINITIONS ==========================================
module class_mesh
  implicit none

! ==========                            1- Curvilinear CLASS                                    =======
  type curvi3d
  ! g_i   --> covariant vectors
  ! gi    --> contavariant vectors
  ! g_ij  --> Covariant metric matrix
  ! gij   --> contravariant metric matrix
  ! gi_j  --> gi,j --> dgi/dxj
    real, dimension(1:3)             :: g_1,g_2,g_3,g1,g2,g3
    real                             :: detG,detJ,x,y,z,x_i,y_i,z_i
    real, dimension(1:3,1:3)         :: jacmat,jacinv,g_ij,gij
    real, dimension(1:3,1:3,1:3)     :: gamma_tensor,gi_j
  end type curvi3d


  type node
! --------------------------------------------
    integer                 :: n              ! Node number
    real                    :: x,y,z          ! Coordinates
    integer, dimension(1:8) :: shared_meshes  ! eight meshes are the max number
    integer, dimension(1:12):: shared_faces=& ! eight meshes are the max number
                               (/0,0,0,0,0,0,0,0,0,0,0,0/)
    real                    :: u=0, v=0, w=0  ! velocity components
    real                    :: p=0            ! Pressure
    real                    :: T=0            ! Temperature
    integer                 :: added_faces=0
    real                    :: Tt             ! True value of T
! --------------------------------------------
  end type node

    type face
! --------------------------------------------
    integer                            :: n             ! edge number
    integer, dimension (1:4)           :: nodes         ! Nodes forming the face (We can store numbers only)
    integer, dimension (1:4)           :: arranged_nodes! Nodes forming the face (We can store numbers only)
    integer, dimension (1:4)           :: edges         ! Edges forming the face (only the edges numbers)
    real                               :: area          ! Area of the face
    integer, dimension (1:2)           :: shared_meshes ! two meshes are the max number

    integer                            :: BC            ! The BC id which correspond to class_BC
    real                               :: x_mid,y_mid,z_mid ! the mid-point coordinate (used to arrange faces)
    real,dimension(1:3)                :: normal        ! vector normal
    logical                            :: is_bc_TD=.false.
	logical                            :: added_to_BC=.false.
! --------------------------------------------
  end type face



    type mesh
! --------------------------------------------
    integer                         :: n=0           ! mesh number
    integer, dimension (1:8)        :: nodes         ! Nodes forming the mesh (We can store numbers only)
    integer, dimension (1:12)       :: edges         ! Edges forming the mesh (only the edges numbers)
    integer, dimension (1:6)        :: faces &       ! faces forming the mesh (only the faces numbers)
                                      =(/0,0,0,0,0,0/)
    integer, dimension (1:6)        :: unarranged_faces & ! Faces arranged based on arrange_faces function
                                      =(/0,0,0,0,0,0/)
    integer, dimension (1:8)        :: unarranged_nodes & ! Faces arranged based on arrange_faces function
                                      =(/0,0,0,0,0,0,0,0/)
    real                            :: volume        ! volume of the mesh
    integer, dimension (1:6)        :: neighbors     ! The six neighbor meshes (0: if boundary)
                                                     ! ==================================================
                                                     ! ======= "left,right,back,front,bottom,top" =======
                                                     ! ======= "x-  ,x+   ,y-  ,y+   ,z-    ,z+ " =======
                                                     ! ==================================================
    integer                         :: added_faces=0, added_neighbors=0, added_nodes=0
    real, dimension(0:7):: C,D,E         ! Coefficient arrays (used in mapping subroutine)
    logical                         :: added_to_queue=.false., arranged=.false.
    integer                         :: face_counter=0
    integer                         :: BC=1              !! 0 if it has face with BC, 1 if it does not contain a face with BC
    real,dimension(1:7,1:3,1:3)     :: gij
    real,dimension(1:3,1:3)         :: Jacinv, xi_x      !! xi_x = transpose(Jacinv)
	real                            :: G1,G2,G3
! --------------------------------------------
  end type mesh


  type boundary_condition
! --------------------------------------------
    integer                         :: n,id          ! BC number
    integer                         :: n_faces       ! number of faces that have the bc
    integer                         :: first,last    ! number of first face and last face of the boundary
    real                :: area          ! Area of the surface boundary
    character(len=3)                :: BC_t,BC_m     ! Thermal bc type and the momentum type
    integer                         :: kind=1        ! 1: Drichlet, 2: Neumann, 3: Robin
! --------------------------------------------
  end type boundary_condition

  contains

end module class_mesh
