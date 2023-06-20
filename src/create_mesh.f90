module class_create_mesh
use class_general
use class_parameters
use class_mesh
use class_curvilinear
implicit none

contains
! =============================================================================================================
! =============================================================================================================
! =========================================      NODES OPERATIONS    ==========================================
! =============================================================================================================
! =============================================================================================================

! *****************  SET COORDINATE ******************
! Input: Node, x,y,z
  subroutine node_create (Vor,n,x,y,z)
! -----------------------------------
    type(node), intent (inout) :: Vor
    real, intent (in)            :: x,y
    real, intent (in), optional  :: z       ! if omitted ==> 2D
    integer, intent (in)         :: n

! -----------------------------------
    Vor%x = x
    Vor%y = y
    Vor%z = z
    Vor%n = n
  end subroutine node_create
! *****************  SET COORDINATE ******************
! Input: Node, x,y,z
  subroutine node_set_coordinate (Vor,x,y,z)
! -----------------------------------
    type(node), intent (inout) :: Vor
    real, intent (in)            :: x,y
    real, intent (in), optional  :: z       ! if omitted ==> 2D
! -----------------------------------
    Vor%x = x
    Vor%y = y
    Vor%z = z

  end subroutine node_set_coordinate

! *****************  SET COORDINATE ******************
! Input: Node, x,y,z
  subroutine node_set_number (Vor,n)
! -----------------------------------
    type(node), intent (inout) :: Vor
    integer, intent (in)         :: n
! -----------------------------------
    Vor%n = n
  end subroutine node_set_number

! *****************  INITIALIZE VARIABLES ******************
! Input: node
  subroutine node_initialize_variables (Vor)
! -----------------------------------
    type(node), intent (inout) :: Vor
! -----------------------------------
    Vor%u = 0
    Vor%v = 0
    Vor%w = 0
    Vor%p = 0
    Vor%T = 0
  end subroutine node_initialize_variables

! ************************** DISTANCE BETWEEN VERTICES *****************
! distance function is to calculate the distance between two coordinates
! Input Two vertices, output: distance between two coordinates
  function node_distance (V1,V2) result (x)
! ---------------------------------------
    type (node), intent(in) :: V1, V2
    real                      :: x
! ---------------------------------------
    x=sqrt((V1%x-V2%x)**2+(V1%y-V2%y)**2+(V1%z-V2%z)**2)
  end function node_distance


! *****************  FIND THE SHARE ELEMENTS ******************
! find all elements that have the same node
! Input: node, share meshes array
  subroutine node_find_shared (Vor,shared)
! -----------------------------------
    type(node), intent (inout)         :: Vor
    integer, dimension(1:8), intent (in) :: shared
! -----------------------------------
    Print*, 'NEED TO BE CREATED LATER!!!!!'
  end subroutine node_find_shared
! =============================================================================================================
! =============================================================================================================
! =========================================       BCS OPERATIONS     ==========================================
! =============================================================================================================
! =============================================================================================================


! =======================================================================================
! *****************  CREATE BC from faces ******************
! Input: BC, n_faces
  subroutine BC_create(BC,n,id,n_faces,first,last,BC_m,BC_t)
! -----------------------------------
    type(boundary_condition), intent (inout)     :: BC
    integer, intent (in)                         :: n,n_faces,first,last,id
    character(len=*), intent (in), optional      :: BC_t, BC_m
! -----------------------------------
    BC%n=n
    BC%n_faces=n_faces
    BC%first=first
    BC%last=last
    BC%BC_m='def'
    BC%BC_t='def'
    BC%id=id
    if (present (BC_m))     BC%BC_m=BC_m
    if (present (BC_t))     BC%BC_t=BC_t

!    call BC_area(BC)
  end subroutine BC_create

  subroutine BC_area(BC)
! ---------------------------------------
    type (boundary_condition), intent(inout) :: BC
    integer                                  :: i
! ---------------------------------------
    BC%area =0;
    do i=BC%first,BC%last
       BC%area=BC%area+FACES_ARRAY(i)%area
    enddo
! ********************* The area need to be calculated later ****************
  end subroutine BC_area


! ****************** Print BC information ********************************
  subroutine BC_print (F)
! ------------------------------------
    type (boundary_condition), intent (in) :: F
! ------------------------------------
    write(6,1) F%n,F%n_faces,F%first ,F%last, F%id
1 format('BC # (',i0,') has ',i0,' faces from (',i0,' to 'i0,') and has id #',i10)
    write(6,4) F%n,F%area
4 format('BC # (',i0,') has area of (',f0.5,')')
  end subroutine BC_print


! ==========================================================================================================================
! ==========================================================================================================================
! ================================ M   E   S   H   E   S           O   P  E  R  A  T  I  O  N  =============================
! ==========================================================================================================================
! ==========================================================================================================================



! *****************  CREATE mesh FROM VERICES ******************
! Input: mesh, V1,V2,V3,V4
  subroutine mesh_create (M,F1,F2,F3,F4,F5,F6,n)
! -----------------------------------
    type(mesh), intent (inout)   :: M
    integer, intent (in)         :: F1,F2,F3,F4,F5,F6
    integer, intent (in)         :: n
! -----------------------------------
    M%faces(1)=F1
    M%faces(2)=F2
    M%faces(3)=F3
    M%faces(4)=F4
    M%faces(5)=F5
    M%faces(6)=F6
    M%n=n
  end subroutine mesh_create


! add faces to the mesh arrays
  subroutine mesh_add_face (M,F1,n)
! -----------------------------------
    type(mesh), intent (inout)   :: M
    integer, intent (in)         :: F1
    integer, intent (in)         :: n
! -----------------------------------
    if (n /= 0) Then
    if (F1 /=0) Then
      M%added_faces=M%added_faces+1
      M%n=n
      if (M%added_faces > 6) Then
        write (6,1) M%n
        1 format('ERR110: Mesh number (',i0,') has more than six faces!! (check the mesh file)')
      else
        M%faces(M%added_faces)=F1
      endif
    endif
    endif
  end subroutine mesh_add_face

  subroutine mesh_check_added_faces (M)
! ====================================
    type(mesh), intent (inout)   :: M
! ====================================
    if (M%added_faces /= 6) Then
      write (6,1) M%n,M%added_faces
      1 format('ERR332: Number of faces in mesh number (',i0,') is not six (nfaces = ',i0, ') Please check mesh file')
      return
    endif
  end subroutine mesh_check_added_faces

  subroutine mesh_add_neighbor (M,F1)
! -----------------------------------
    type(mesh), intent (inout)   :: M
    integer, intent (in)         :: F1
! -----------------------------------
    if (M%n /= 0) Then
      M%added_neighbors=M%added_neighbors+1
      if (M%added_neighbors > 6) Then
      write (6,1) M%n
      1 format('ERR111: Mesh number (',i0,') has more than six neighbors!! (check the mesh file)')
      else
      M%neighbors(M%added_neighbors)=F1
      endif
    endif
!    call Face_area(F)
  end subroutine mesh_add_neighbor

! This subroutine extracts the node numbers from the faces. It helps finding the nodes of each mesh
  subroutine mesh_add_nodes (M)
! -----------------------------------
    type(mesh), intent (inout)   :: M
    integer                      :: i,j,k,node
    logical                      :: exist
! -----------------------------------
    do i=1,6
       do j=1,4
          exist=.FALSE.
          node=FACES_ARRAY(M%faces(i))%nodes(j)
!          print*, node
          do k=1,8
             if (node == M%nodes(k)) exist=.TRUE.
          enddo
          if (.NOT. exist) call mesh_add_one_node(M,node)
       enddo
    enddo
  end subroutine mesh_add_nodes
 subroutine mesh_add_one_node (M,node)
! -----------------------------------
    type(mesh), intent (inout)   :: M
    integer                      :: node
! -----------------------------------
    M%added_nodes=M%added_nodes+1
!    print*, 'called'
    if (M%added_nodes > 8) Then
    write (6,1) M%n
    1 format('ERR112: Mesh number (',i0,') has more than eight nodes!!! (check the mesh file)')
    else
    M%nodes(M%added_nodes)=node
    endif

  end subroutine mesh_add_one_node


 subroutine mesh_calculate_mapping_coefficients(M)
! Should be called after arranging the nodes and faces
! -----------------------------------
    type(mesh), intent (inout)   :: M
    real, dimension(1:8)      :: x,y,z    ! Location of the nodes
    integer                      :: i
! -----------------------------------

    do i=1,8
       x(i)=NODES_ARRAY(M%nodes(i))%x
       y(i)=NODES_ARRAY(M%nodes(i))%y
       z(i)=NODES_ARRAY(M%nodes(i))%z
    enddo
! -----------------------------------------------------
!      8------7        z  y
!     /|     /|        | /
!    5------6 |        |/
!    | |    | |        0-------x
!    | 4----|-3
!    |/     |/
!    1------2
!
! -----------------------------------------------------
    M%C(0)=(x(7)+x(8)+x(5)+x(6)+x(3)+x(4)+x(1)+x(2))/8.
    M%C(1)=(x(7)-x(8)-x(5)+x(6)+x(3)-x(4)-x(1)+x(2))/8.
    M%C(2)=(x(7)+x(8)-x(5)-x(6)+x(3)+x(4)-x(1)-x(2))/8.
    M%C(3)=(x(7)+x(8)+x(5)+x(6)-x(3)-x(4)-x(1)-x(2))/8.
    M%C(4)=(x(7)-x(8)+x(5)-x(6)+x(3)-x(4)+x(1)-x(2))/8.
    M%C(5)=(x(7)-x(8)-x(5)+x(6)-x(3)+x(4)+x(1)-x(2))/8.
    M%C(6)=(x(7)+x(8)-x(5)-x(6)-x(3)-x(4)+x(1)+x(2))/8.
    M%C(7)=(x(7)-x(8)+x(5)-x(6)-x(3)+x(4)-x(1)+x(2))/8.
! -----------------------------------------------------
    M%D(0)=(y(7)+y(8)+y(5)+y(6)+y(3)+y(4)+y(1)+y(2))/8.
    M%D(1)=(y(7)-y(8)-y(5)+y(6)+y(3)-y(4)-y(1)+y(2))/8.
    M%D(2)=(y(7)+y(8)-y(5)-y(6)+y(3)+y(4)-y(1)-y(2))/8.
    M%D(3)=(y(7)+y(8)+y(5)+y(6)-y(3)-y(4)-y(1)-y(2))/8.
    M%D(4)=(y(7)-y(8)+y(5)-y(6)+y(3)-y(4)+y(1)-y(2))/8.
    M%D(5)=(y(7)-y(8)-y(5)+y(6)-y(3)+y(4)+y(1)-y(2))/8.
    M%D(6)=(y(7)+y(8)-y(5)-y(6)-y(3)-y(4)+y(1)+y(2))/8.
    M%D(7)=(y(7)-y(8)+y(5)-y(6)-y(3)+y(4)-y(1)+y(2))/8.
! -----------------------------------------------------
    M%E(0)=(z(7)+z(8)+z(5)+z(6)+z(3)+z(4)+z(1)+z(2))/8.
    M%E(1)=(z(7)-z(8)-z(5)+z(6)+z(3)-z(4)-z(1)+z(2))/8.
    M%E(2)=(z(7)+z(8)-z(5)-z(6)+z(3)+z(4)-z(1)-z(2))/8.
    M%E(3)=(z(7)+z(8)+z(5)+z(6)-z(3)-z(4)-z(1)-z(2))/8.
    M%E(4)=(z(7)-z(8)+z(5)-z(6)+z(3)-z(4)+z(1)-z(2))/8.
    M%E(5)=(z(7)-z(8)-z(5)+z(6)-z(3)+z(4)+z(1)-z(2))/8.
    M%E(6)=(z(7)+z(8)-z(5)-z(6)-z(3)-z(4)+z(1)+z(2))/8.
    M%E(7)=(z(7)-z(8)+z(5)-z(6)-z(3)+z(4)-z(1)+z(2))/8.
  end subroutine mesh_calculate_mapping_coefficients

! *****************  Calculating average u,v,w ******************
  subroutine volume_calculate (M)
! -----------------------------------
    type(mesh), intent (in)         :: M
    real, dimension (1:Gauss_int_p,1:Gauss_int_p,1:Gauss_int_p) :: v_i
    type(curvi3d)                   :: C
    integer                         :: i,j,k
! -----------------------------------
    do i=1,Gauss_int_p
       do j=1,Gauss_int_p
          do k=1,Gauss_int_p
             call evaluate_curvilinear(M,C,Gauss_x_i(i),Gauss_x_i(j),Gauss_x_i(k))
             v_i(i,j,k)=1.*C%detJ
          enddo
       enddo
    enddo
    meshes_volume(M%n)=Gauss_integration(v_i)
  end subroutine volume_calculate


! ---------------------------------------------------------------------------------------------------------
! The Following function will arrange the nodes as the following figure!!
!          4
!      8------7        z  y
!     /| 6   /|        | /
!    5------6 |        |/
!    |1| 3  |2|        0-------x
!    | 4----|-3
!    |/  5  |/
!    1------2
  subroutine mesh_arrange_nodes (M)
! -----------------------------------
    type(mesh), intent (inout)      :: M
    integer, dimension(1:8)         :: nodes,nodes_old
    real, dimension(1:8,1:4)        :: nodes_all
    real, dimension(1:4,1:4)        :: nodes_z1,nodes_z2
    integer, dimension(1:2)         :: nodes_z1_y1,nodes_z1_y2,nodes_z2_y1,node_z2_y2 ! We may not use these
    integer                         :: i
    real, dimension(1:4)            :: temp
! NOTES: The arrangment is not very general here!!! it needs to be updated when finding a better way
!        The first thing is to seprate the nodes based on z dirction (1,2,3,4) are together
!        This can be found by finding the four points that have the lowest z_values
!        Therefore, we need to arrange the array based on z-increase
!        Then we seperate the array into two arrays (1,2,3,4) and (5,6,7,8)
!        Then we arrange each of the two arrays based on y coordinate and seperate each array into
!        two arrays ----> (1,2),(3,4),(5,6),(7,8)
!        Finally, we arrange each of the four arrays based on x coordinate

    M%unarranged_nodes=M%nodes
goto 102
    do i=1,8
       nodes_all(i,1)=REAL(M%nodes(i))
       nodes_all(i,2)=NODES_ARRAY(M%nodes(i))%x
       nodes_all(i,3)=NODES_ARRAY(M%nodes(i))%y
       nodes_all(i,4)=NODES_ARRAY(M%nodes(i))%z
    enddo
    call array_ascending2(nodes_all,8,4,4)

    call array_ascending2(nodes_all(1:4,1:4),4,4,3)
    call array_ascending2(nodes_all(5:8,1:4),4,4,3)


    call array_ascending2(nodes_all(1:2,1:4),2,4,2)
    call array_ascending2(nodes_all(3:4,1:4),2,4,2)
    call array_ascending2(nodes_all(5:6,1:4),2,4,2)
    call array_ascending2(nodes_all(7:8,1:4),2,4,2)


    M%nodes(1)=REAL(nodes_all(1,1))
    M%nodes(2)=REAL(nodes_all(2,1))
    M%nodes(3)=REAL(nodes_all(4,1))
    M%nodes(4)=REAL(nodes_all(3,1))
    M%nodes(5)=REAL(nodes_all(5,1))
    M%nodes(6)=REAL(nodes_all(6,1))
    M%nodes(7)=REAL(nodes_all(8,1))
    M%nodes(8)=REAL(nodes_all(7,1))
    print*,'NODES:',M%nodes(1),M%nodes(2),M%nodes(3),M%nodes(4),M%nodes(5),M%nodes(6)

!    print*, '============================'
102 continue
  end subroutine mesh_arrange_nodes
! The following function will be used to arrange the faces order in faces array, the desired order will be:
! 'left right back front bottom top'
! This function needs more works to be done
! The arragment is based on the node arragments (WE need to arrange nodes before)
! The node arrangment is as the following:
!          4
!      8------7        z  y
!     /| 6   /|        | /
!    5------6 |        |/
!    |1| 3  |2|        0-------x
!    | 4----|-3
!    |/  5  |/
!    1------2
!
!
!    5------8    6------7    5------6
!    |      |    |      |    |      |
!    |  1   |    |  2   |    |  3   |
!    1------4    2------3    1------2
!
!    8------7    4------3    8------7
!    |      |    |      |    |      |
!    |  4   |    |  5   |    |  6   |
!    4------3    1------2    5------6
!
! The function arranges the faces and the nodes in the faces

  subroutine mesh_arrange_faces (M)
! -----------------------------------
    type(mesh), intent (inout)      :: M
    integer, dimension(1:6)         :: faces_old
    integer, dimension(1:4)         :: nodes,arranged_nodes
    integer                         :: i
    real, dimension(1:6,1:4)        :: faces_all
    real, dimension(1:4,1:4)        :: nodes_all

! -----------------------------------
    M%unarranged_faces=M%faces

! ===========================================================
! ========= The following arrangment is based on the normal =
! ========= We choose the direction based on the largest component of normal vector
! ==========================================================

    do i=1,6
       faces_all(i,1)=REAL(M%faces(i))
       faces_all(i,2)=abs(FACES_ARRAY(M%faces(i))%normal(1))
       faces_all(i,3)=abs(FACES_ARRAY(M%faces(i))%normal(2))
       faces_all(i,4)=abs(FACES_ARRAY(M%faces(i))%normal(3))
    enddo
    call array_ascending2(faces_all,6,4,4)
    call array_ascending2(faces_all(1:4,1:4),4,4,3)
    call array_ascending2(faces_all(1:2,1:4),2,4,2)

    if (FACES_ARRAY(INT(faces_all(6,1)))%z_mid>FACES_ARRAY(INT(faces_all(5,1)))%z_mid) THEN
       M%faces(6)=faces_all(6,1)
       M%faces(5)=faces_all(5,1)
    else
       M%faces(5)=faces_all(6,1)
       M%faces(6)=faces_all(5,1)
    endif
    if (FACES_ARRAY(INT(faces_all(4,1)))%y_mid>FACES_ARRAY(INT(faces_all(3,1)))%y_mid) THEN
       M%faces(4)=faces_all(4,1)
       M%faces(3)=faces_all(3,1)
    else
       M%faces(3)=faces_all(4,1)
       M%faces(4)=faces_all(3,1)
    endif

    if (FACES_ARRAY(INT(faces_all(2,1)))%x_mid>FACES_ARRAY(INT(faces_all(1,1)))%x_mid) THEN
       M%faces(2)=faces_all(2,1)
       M%faces(1)=faces_all(1,1)
    else
       M%faces(1)=faces_all(2,1)
       M%faces(2)=faces_all(1,1)
    endif

    do i=1,4
       nodes_all(i,1)=REAL(FACES_ARRAY(M%faces(3))%nodes(i))
       nodes_all(i,2)=NODES_ARRAY(FACES_ARRAY(M%faces(3))%nodes(i))%x
       nodes_all(i,3)=NODES_ARRAY(FACES_ARRAY(M%faces(3))%nodes(i))%y
       nodes_all(i,4)=NODES_ARRAY(FACES_ARRAY(M%faces(3))%nodes(i))%z
    enddo
    call array_ascending2(nodes_all(1:4,1:4),4,4,4)

    call array_ascending2(nodes_all(1:2,1:4),2,4,2)
    call array_ascending2(nodes_all(3:4,1:4),2,4,2)
    M%nodes(1)=REAL(nodes_all(1,1))
    M%nodes(2)=REAL(nodes_all(2,1))
    M%nodes(5)=REAL(nodes_all(3,1))
    M%nodes(6)=REAL(nodes_all(4,1))


    do i=1,4
       nodes_all(i,1)=REAL(FACES_ARRAY(M%faces(4))%nodes(i))
       nodes_all(i,2)=NODES_ARRAY(FACES_ARRAY(M%faces(4))%nodes(i))%x
       nodes_all(i,3)=NODES_ARRAY(FACES_ARRAY(M%faces(4))%nodes(i))%y
       nodes_all(i,4)=NODES_ARRAY(FACES_ARRAY(M%faces(4))%nodes(i))%z

    enddo
    call array_ascending2(nodes_all(1:4,1:4),4,4,4)

    call array_ascending2(nodes_all(1:2,1:4),2,4,2)
    call array_ascending2(nodes_all(3:4,1:4),2,4,2)
    M%nodes(4)=REAL(nodes_all(1,1))
    M%nodes(3)=REAL(nodes_all(2,1))
    M%nodes(8)=REAL(nodes_all(3,1))
    M%nodes(7)=REAL(nodes_all(4,1))

  end subroutine mesh_arrange_faces

! This function is to find the face number that has the nodes array
!    8------7
!    |      |
!    |  6   |
!    5------6
  function mesh_find_face(M,nodes) result(x)
    type(mesh), intent (inout)      :: M
!    integer, dimension(1:6)         :: faces_old
    integer, dimension(1:4)         :: nodes
    integer                         :: i,j,k,x
    x=0


    do i=1,6
         if ( is_exist (FACES_ARRAY(M%unarranged_faces(i))%nodes,4,nodes(1)) .AND. &
              is_exist (FACES_ARRAY(M%unarranged_faces(i))%nodes,4,nodes(2)) .AND. &
              is_exist (FACES_ARRAY(M%unarranged_faces(i))%nodes,4,nodes(3)) .AND. &
              is_exist (FACES_ARRAY(M%unarranged_faces(i))%nodes,4,nodes(4))) THEN
         x=M%unarranged_faces(i)
         return
         else
!                    FACES_ARRAY(M%faces(i))%nodes(3),FACES_ARRAY(M%faces(i))%nodes(4)
         endif
    enddo
    print*, 'DOES NOT EXIST',M%n,nodes(1),nodes(2),nodes(3),nodes(4), M%faces(i)
   M%face_counter=M%face_counter+1;
  end function mesh_find_face

! *****************  CREATE mesh FROM VERICES ******************
! Input: mesh, V1,V2,V3,V4
  subroutine mesh_initialize (M)
! -----------------------------------
    type(mesh), intent (inout)   :: M
! -----------------------------------
    call mesh_create(M,0,0,0,0,0,0,0)
  end subroutine mesh_initialize


! This function is used to calculate the average of T at point from two points
! Used to calculate T at the center of edges
  function TP(M,p1,p2)         result(X)
    type(mesh), intent (in)      :: M
    integer, intent (in)         :: p1,p2
    real                         :: X
    X=(NODES_ARRAY(M%nodes(p1))%T+NODES_ARRAY(M%nodes(p2))%T)/2.
    X=NODES_ARRAY(M%nodes(p1))%T
   end function TP

  function opposite_face(face)  result(i)
    implicit none
    integer, intent(in)  :: face
    integer              :: i
    if (face == 1) then
       i=2
    else if (face == 2) then
       i=1
    else if (face == 3) then
       i=4
    else if (face == 4) then
       i=3
    else if (face == 5) then
       i=6
    else if (face == 6) then
       i=5
    else
       print*, 'ERR303: FACE NUMBER ', face, 'DOES NOT EXIST: PLEASE CHECK OPPOSITE_FACES ROUTINE (INTEGRATION)'
    endif

  end function opposite_face
  subroutine other_faces(i,i2,i3)
     integer, intent (in)    :: i
     integer, intent (inout) :: i2,i3
     if (i == 1 .OR. i == 2) THEN
       i2=4
       i3=6
     else if (i == 3 .OR. i == 4) Then
       i2=2
       i3=6
     else if (i == 5 .OR. i==6) Then
       i2=2
       i3=4
     endif
    if (i == 1) then
       i2=3
       i3=5
    else if (i == 2) then
       i2=4
       i3=6
    else if (i == 3) then
       i2=1
       i3=5
    else if (i == 4) then
       i2=2
       i3=6
    else if (i == 5) then
       i2=1
       i3=3
    else if (i == 6) then
       i2=2
       i3=4
    endif
  end subroutine other_faces

! *****************  CREATE FACE FROM VERICES ******************
! Input: face, V1,V2,V3,V4
  subroutine face_create (F,V1,V2,V3,V4,n,BC)
! -----------------------------------
    type(face), intent (inout)              :: F
    integer,  intent (in)                   :: V1,V2,V3,V4
    integer, intent (in)                    :: n
    integer, intent (in)                    :: BC
! -----------------------------------
    F%nodes(1)=V1
    F%nodes(2)=V2
    F%nodes(3)=V3
    F%nodes(4)=V4
    F%n=n
    F%BC=BC
    call face_mid_point(F)
    call face_find_normal(F)
  end subroutine face_create


! *****************  SET Nodes ******************
! Input: face, V1,V2,V3,V4
  subroutine face_set_nodes (F,V1,V2,V3,V4)
! -----------------------------------
    type(face), intent (inout)   :: F
    integer, intent (in)         :: V1,V2,V3,V4
! -----------------------------------
    F%nodes(1)=V1
    F%nodes(2)=V2
    F%nodes(3)=V3
    F%nodes(4)=V4
  end subroutine face_set_nodes
! *****************  SET Nodes ******************
! Input: face, V1,V2,V3,V4
  subroutine face_set_edges (F,E1,E2,E3,E4)
! -----------------------------------
    type(face), intent (inout)   :: F
    integer, intent (in)         :: E1,E2,E3,E4
! -----------------------------------
    F%edges(1)=E1
    F%edges(2)=E2
    F%edges(3)=E3
    F%edges(4)=E4
  end subroutine face_set_edges


! *****************  SET EDGE NUMBER ******************
! Input: Node, x,y,z
  subroutine face_set_number (F,n)
! -----------------------------------
    type(face), intent (inout) :: F
    integer, intent (in)       :: n
! -----------------------------------
    F%n = n
  end subroutine face_set_number

! This subroutine extracts the node numbers from the faces. It helps finding the nodes of each mesh
  subroutine node_add_faces (F,n)
! -----------------------------------
    integer,dimension(1:4), intent (in) :: F ! the node array of the faces
    integer, intent (in)                :: n ! Number of the faces
    integer                             :: i,j,k,node
    logical                             :: exist
! -----------------------------------
     do i=1,4
        exist=.FALSE.
        node=F(i)
        do j=1,12
           if (n == NODES_ARRAY(node)%shared_faces(j)) exist=.TRUE.
        enddo
        if (.NOT. exist) call node_add_one_face(NODES_ARRAY(node),n)
       enddo
  end subroutine node_add_faces


! *****************  Add the shared faces to the node ******************
! Input: mesh, V1,V2,V3,V4
  subroutine node_add_one_face (N,F1)
! -----------------------------------
    type(node), intent (inout)   :: N
    integer, intent (in)         :: F1  ! face number
    integer                      :: i
! -----------------------------------
      N%added_faces=N%added_faces+1
      if (N%added_faces > 12) Then
      write (6,1) N%n
      1 format('ERR113: node number (',i0,') has more than 12 faces!! (check the mesh file) -> class_Face')
      do i=1,12
         print*,N%shared_faces(i), F1
      enddo
      else
      N%shared_faces(N%added_faces)=F1
      endif
  end subroutine node_add_one_face

  subroutine face_set_shared (F,shared)
! Maybe we don't need the shared matix. Also we probably need the total number of elements and the element array
! -----------------------------------
    type(face), intent (inout)           :: F
    integer, dimension(1:2), intent (in) :: shared
! -----------------------------------
    F%shared_meshes=shared
  end subroutine face_set_shared
  subroutine face_find_normal (F)
! check https://www.ma.utexas.edu/users/m408m/Display12-5-4.shtml
! -----------------------------------
    type(face), intent (inout)           :: F
! -----------------------------------
    real, dimension(1:3)                 :: r_b,s_b,n
    real                                 :: temp
    integer                              :: i,p1,p2,p3
    p1=F%nodes(1)
    p2=F%nodes(2)
    p3=F%nodes(3)
    r_b=(/NODES_ARRAY(p1)%x-NODES_ARRAY(p2)%x, &
          NODES_ARRAY(p1)%y-NODES_ARRAY(p2)%y, &
          NODES_ARRAY(p1)%z-NODES_ARRAY(p2)%z/)
    s_b=(/NODES_ARRAY(p3)%x-NODES_ARRAY(p2)%x, &
          NODES_ARRAY(p3)%y-NODES_ARRAY(p2)%y, &
          NODES_ARRAY(p3)%z-NODES_ARRAY(p2)%z/)
    n=cross(s_b,r_b)
! Normalize the vector
    temp=dot(n,n)**.5
    do i=1,3
       F%normal(i)=n(i)/temp
    enddo
  end subroutine face_find_normal

! Calculates the mid-point of a face
  subroutine face_mid_point(F)
! ---------------------------------------
    type (face), intent(inout) :: F
! ---------------------------------------
    F%x_mid=(NODES_ARRAY(F%nodes(1))%x+NODES_ARRAY(F%nodes(2))%x+NODES_ARRAY(F%nodes(3))%x+NODES_ARRAY(F%nodes(4))%x)/4
    F%y_mid=(NODES_ARRAY(F%nodes(1))%y+NODES_ARRAY(F%nodes(2))%y+NODES_ARRAY(F%nodes(3))%y+NODES_ARRAY(F%nodes(4))%y)/4
    F%z_mid=(NODES_ARRAY(F%nodes(1))%z+NODES_ARRAY(F%nodes(2))%z+NODES_ARRAY(F%nodes(3))%z+NODES_ARRAY(F%nodes(4))%z)/4
  end subroutine face_mid_point

end module class_create_mesh

