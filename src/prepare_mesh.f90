module class_prepare_mesh
use class_general
use class_parameters
use class_mesh
use class_curvilinear
use class_create_mesh
implicit none

contains

subroutine prepare_mesh
  implicit none
  integer ::i
  integer :: start, finish
  INTEGER :: &
     nb_ticks_initial, & ! initial value of the clock tick counter
     nb_ticks_final,   & ! final value of the clock tick counter
     nb_ticks_max,     & ! maximum value of the clock counter
     nb_ticks_sec,     & ! number of clock ticks per second
     nb_ticks           ! number of clock ticks of the code
  REAL :: elapsed_time  ! real time in seconds
    
!    call deform_nodes
! The following loop is to find the nodes of each mesh
! The following loops can be added to class_Mesh as a subroutine and can be called there!!
    do i=1,NMESH
       call mesh_add_nodes(MESHES_ARRAY(i))
    enddo

    call link_mesh()

    call define_neighbors()


    do i=1,NFACE
       call node_add_faces(FACES_ARRAY(i)%nodes,i)
    enddo
    do i=1,NBC
       call BC_area(BCS_ARRAY(i))
    enddo

end subroutine prepare_mesh

! The purpose of this function is to link the meshes together. This function will determine
! the relation between all meshes. The left, right, top,... faces and neighbor meshes will
! be determined at the end of this function.
subroutine link_mesh()
    implicit none
    integer                      :: i,j,k,start_point
    integer, dimension (1:NMESH) :: Queue
    integer                      :: opposite1,counter,M1,M2
    integer,dimension(1:6)       :: faces2

! For the faces array in mesh class, if -1: BC, 0: not defined yet, otherwise: already defined.
! First, we start with the first mesh, we arrange faces and nodes of the mesh using mesh_arrange_faces subroutine
!    call mesh_arrange_faces(MESHES_ARRAY(1))
! Mesh 1 is ready to be used to link the other meshes.
! we need outer loop over the top and bottom queues.
!
    counter=1
    start_point=1
    call mesh_arrange_faces(MESHES_ARRAY(start_point))



    Queue(1)=start_point
    MESHES_ARRAY(start_point)%added_to_queue=.TRUE.
    MESHES_ARRAY(start_point)%arranged=.TRUE.
! Fill the queue
    do i=1,NMESH
       M1=Queue(i)
       do j=1,6
          M2=adjacent_mesh(M1,j)
          if(M2 /=0) then
            if (.NOT. MESHES_ARRAY(M2)%arranged) then
               call match_next_mesh_faces(M1,M2,j)
               MESHES_ARRAY(M2)%arranged=.TRUE.
            endif
            if (.NOT. MESHES_ARRAY(M2)%added_to_queue) then
               counter=counter+1
               Queue(counter)=M2
               MESHES_ARRAY(M2)%added_to_queue=.TRUE.
            endif

          endif
       enddo
    enddo
end subroutine link_mesh



! This function define the face numbers based on a pre-defined mesh.
! M1 is the predefined mesh, M2 is the mesh that we need to its
! faces to be defined. F is the common faces between the meshes.
subroutine match_next_mesh_faces(M1,M2,F)
    integer, intent(in)     :: M1,M2
    integer, intent(in)     :: F
    integer                 :: F2,F2_old,opp ! index of the F in mesh2
    integer,dimension(1:6)  :: fn2
    integer,dimension(1:2)  :: edge

! First, the face number F in M1 will be the opposite face in M2.
! Then we need to find the opposite face of M2(F)

!    print*,'faces==',MESHES_ARRAY(M2)%faces

    F2=opposite_face(F)
    fn2(F2)=MESHES_ARRAY(M1)%faces(F) ! The same face is M1 is opposite in M2.
    F2_old=position_in_array(MESHES_ARRAY(M2)%faces,6,fn2(F2))
    fn2(F)=opposite_face_from_nodes(M2,F2_old)
    if (F /= 5 .AND. F/= 6) then
       edge=get_top_edge(M1,F)
       fn2(6)=find_face_from_two_nodes(M2,F2_old,edge)
       edge=get_bottom_edge(M1,F)
       fn2(5)=find_face_from_two_nodes(M2,F2_old,edge)
    endif
    if (F /= 4 .AND. F/= 3) then
       edge=get_front_edge(M1,F)
       fn2(4)=find_face_from_two_nodes(M2,F2_old,edge)
       edge=get_back_edge(M1,F)
       fn2(3)=find_face_from_two_nodes(M2,F2_old,edge)
    endif
    if (F /= 2 .AND. F/= 1) then
       edge=get_right_edge(M1,F)
       fn2(2)=find_face_from_two_nodes(M2,F2_old,edge)
       edge=get_left_edge(M1,F)
       fn2(1)=find_face_from_two_nodes(M2,F2_old,edge)
    endif

    MESHES_ARRAY(M2)%faces=fn2
    call arrange_nodes_from_faces(M2)
  end subroutine match_next_mesh_faces

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

  subroutine arrange_nodes_from_faces(M)
     integer, intent(in)     :: M
     integer, dimension(1:8) :: nodes2

     nodes2(1)=node_from_three_faces(M,1,3,5)
     nodes2(2)=node_from_three_faces(M,2,3,5)
     nodes2(3)=node_from_three_faces(M,2,4,5)
     nodes2(4)=node_from_three_faces(M,1,4,5)
     nodes2(5)=node_from_three_faces(M,1,3,6)
     nodes2(6)=node_from_three_faces(M,2,3,6)
     nodes2(7)=node_from_three_faces(M,2,4,6)
     nodes2(8)=node_from_three_faces(M,1,4,6)

     MESHES_ARRAY(M)%nodes=nodes2
!     print*, nodes2
  end subroutine arrange_nodes_from_faces

  subroutine define_neighbors()
    integer               :: i

    do i=1,NMESH
       MESHES_ARRAY(i)%neighbors(1)=adjacent_mesh(i,1)
       MESHES_ARRAY(i)%neighbors(2)=adjacent_mesh(i,2)
       MESHES_ARRAY(i)%neighbors(3)=adjacent_mesh(i,3)
       MESHES_ARRAY(i)%neighbors(4)=adjacent_mesh(i,4)
       MESHES_ARRAY(i)%neighbors(5)=adjacent_mesh(i,5)
       MESHES_ARRAY(i)%neighbors(6)=adjacent_mesh(i,6)
    enddo

  end subroutine define_neighbors


  function face_is_defined(M,F)  result(is_defined)
    integer, intent(in)    :: M
    integer, intent(in)    :: F
    logical                :: is_defined
    if (MESHES_ARRAY(M)%faces(F) == 0) then
        is_defined=.FALSE.
    else
        is_defined=.TRUE.
    endif
  end function face_is_defined

!  Inputs are the mesh number and the common face which the adjacent needs to be found.
  function adjacent_mesh(M,F) result(adjacent)
    integer, intent(in)    :: M
    integer, intent(in)    :: F
    integer                :: adjacent
    adjacent=0
    if (FACES_ARRAY(MESHES_ARRAY(M)%FACES(F))%shared_meshes(1) == M) then
        adjacent=FACES_ARRAY(MESHES_ARRAY(M)%FACES(F))%shared_meshes(2)
    else
        adjacent=FACES_ARRAY(MESHES_ARRAY(M)%FACES(F))%shared_meshes(1)
    endif
!    if (adjacent == 0) print*, 'IT is on the boundary'
  end function adjacent_mesh
  function node_from_three_faces(M,F1,F2,F3) result(node)
     integer, intent(in)     :: M,F1,F2,F3
     integer                 :: i,j,k,counter
     integer                 :: node
     counter=counter+1
     do i=1,4
        counter=0
        do j=1,4
           if(FACES_ARRAY(MESHES_ARRAY(M)%faces(F1))%nodes(i)==FACES_ARRAY(MESHES_ARRAY(M)%faces(F2))%nodes(j)) &
             counter=counter+1
        enddo
        do k=1,4
           if(FACES_ARRAY(MESHES_ARRAY(M)%faces(F1))%nodes(i)==FACES_ARRAY(MESHES_ARRAY(M)%faces(F3))%nodes(k)) &
             counter=counter+1
        enddo
        if (counter==2) then
           node=FACES_ARRAY(MESHES_ARRAY(M)%faces(F1))%nodes(i)
!     print*,"----->>>>",node,i,F1,M

           return
        endif
     enddo
!     print*,"----->>>>",node
     print*,"----->>>>",node,i,F1,M,MESHES_ARRAY(M)%faces

  end function node_from_three_faces

  function opposite_face_from_nodes(M,F) result(opposite)
    integer, intent(in)     :: M
    integer, intent(in)     :: F
    integer                 :: opposite
    integer                 :: adjacent
    integer, dimension(1:4) :: face_nodes!, opposite_face_nodes
    integer                 :: i,j
    logical                 :: found=.FALSE.
    integer                 :: counter=1
    counter=1
    do i=1,8
       found=.FALSE.
       do j=1,4
          if (MESHES_ARRAY(M)%nodes(i)==FACES_ARRAY(MESHES_ARRAY(M)%faces(F))%nodes(j)) found=.TRUE.
       enddo
       if (.NOT. found) then
          face_nodes(counter)=MESHES_ARRAY(M)%nodes(i)
          counter=counter+1
!          print*, 'counter = ',counter
       endif
    enddo
    opposite=find_face_from_nodes(M,face_nodes)
  end function opposite_face_from_nodes

  function find_face_from_nodes(M,nodes)  result(F)
    integer, intent(in)                 :: M
    integer, dimension(1:4), intent(in) :: nodes
    integer                             :: F
    integer                             :: c,i,j
    integer, dimension(1:4)             :: tmp
    call quicksort(nodes,1,4)
    do i=1,6
       tmp=FACES_ARRAY(MESHES_ARRAY(M)%faces(i))%nodes
       call quicksort(tmp,1,4)
       if (are_nodes_same(nodes,tmp)) then
         F=MESHES_ARRAY(M)%faces(i)
         return
       endif
    enddo
  end function find_face_from_nodes

  function are_nodes_same(nodes1,nodes2)  result(res)
    integer, dimension(1:4), intent(in) :: nodes1,nodes2
! Both are sorted arrays
    integer                             :: i
    logical                             :: res

    do i=1,4
       if (nodes1(i) /= nodes2(i)) then
          res=.FALSE.
          return
       endif
    enddo
    res=.TRUE.
  end function are_nodes_same


  function find_face_from_two_nodes(M,F,nodes) result (FFn)
    integer, intent(in)                 :: M,F
    integer, dimension(1:2), intent(in) :: nodes
    integer                             :: i,FFn
    integer,dimension(1:2)              :: fn!,nodes2
    FFn=0
    do i=1,6
       if (is_edge_exist(FACES_ARRAY(MESHES_ARRAY(M)%faces(i))%nodes,nodes)) then
          if (F/=i) then
             FFn=MESHES_ARRAY(M)%faces(i)
             return
          endif
       endif
    enddo
  end function find_face_from_two_nodes

  function is_edge_exist(face_nodes,edge_nodes)   result(exist)
    integer, dimension(1:2), intent(in) :: edge_nodes
    integer, dimension(1:4), intent(in) :: face_nodes
    integer                             :: i,k,j,counter=0
    logical                             :: exist
    counter=0
    do j=1,2
      do k=1,4
         if(face_nodes(k)==edge_nodes(j)) counter=counter+1
      enddo
    enddo
    if (counter==2) then
      exist=.TRUE.
      return
    else
      exist=.FALSE.
      return
    endif
  end function is_edge_exist


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
  function get_top_edge(M,F) result (edge2)
    integer, intent(in)     :: M
    integer, intent(in)     :: F
    integer, dimension(1:2) :: edge1,edge2
    if (F==1) then
       edge1=(/5,8/)
    else if (F==2) then
       edge1=(/6,7/)
    else if (F==3) then
       edge1=(/5,6/)
    else if (F==4) then
       edge1=(/7,8/)
    endif
    edge2=(/MESHES_ARRAY(M)%nodes(edge1(1)),MESHES_ARRAY(M)%nodes(edge1(2))/)
  end function get_top_edge
  function get_bottom_edge(M,F) result (edge2)
    integer, intent(in)     :: M
    integer, intent(in)     :: F
    integer, dimension(1:2) :: edge1,edge2
    if (F==1) then
       edge1=(/1,4/)
    else if (F==2) then
       edge1=(/2,3/)
    else if (F==3) then
       edge1=(/1,2/)
    else if (F==4) then
       edge1=(/3,4/)
    endif
    edge2=(/MESHES_ARRAY(M)%nodes(edge1(1)),MESHES_ARRAY(M)%nodes(edge1(2))/)
  end function get_bottom_edge
  function get_right_edge(M,F) result (edge2)
    integer, intent(in)     :: M
    integer, intent(in)     :: F
    integer, dimension(1:2) :: edge1,edge2
    if (F==3) then
       edge1=(/2,6/)
    else if (F==4) then
       edge1=(/3,7/)
    else if (F==5) then
       edge1=(/2,3/)
    else if (F==6) then
       edge1=(/6,7/)
    endif
    edge2=(/MESHES_ARRAY(M)%nodes(edge1(1)),MESHES_ARRAY(M)%nodes(edge1(2))/)
  end function get_right_edge
  function get_left_edge(M,F) result (edge2)
    integer, intent(in)     :: M
    integer, intent(in)     :: F
    integer, dimension(1:2) :: edge1,edge2
    if (F==3) then
       edge1=(/1,5/)
    else if (F==4) then
       edge1=(/4,8/)
    else if (F==5) then
       edge1=(/1,4/)
    else if (F==6) then
       edge1=(/5,8/)
    endif
    edge2=(/MESHES_ARRAY(M)%nodes(edge1(1)),MESHES_ARRAY(M)%nodes(edge1(2))/)
  end function get_left_edge
  function get_back_edge(M,F) result (edge2)
    integer, intent(in)     :: M
    integer, intent(in)     :: F
    integer, dimension(1:2) :: edge1,edge2
    if (F==1) then
       edge1=(/1,5/)
    else if (F==2) then
       edge1=(/2,6/)
    else if (F==5) then
       edge1=(/1,2/)
    else if (F==6) then
       edge1=(/5,6/)
    endif
    edge2=(/MESHES_ARRAY(M)%nodes(edge1(1)),MESHES_ARRAY(M)%nodes(edge1(2))/)
  end function get_back_edge
  function get_front_edge(M,F) result (edge2)
    integer, intent(in)     :: M
    integer, intent(in)     :: F
    integer, dimension(1:2) :: edge1,edge2
    if (F==1) then
       edge1=(/4,8/)
    else if (F==2) then
       edge1=(/3,7/)
    else if (F==5) then
       edge1=(/3,4/)
    else if (F==6) then
       edge1=(/7,8/)
    endif
    edge2=(/MESHES_ARRAY(M)%nodes(edge1(1)),MESHES_ARRAY(M)%nodes(edge1(2))/)
  end function get_front_edge




end module class_prepare_mesh

