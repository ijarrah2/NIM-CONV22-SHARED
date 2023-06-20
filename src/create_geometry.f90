module class_create_geometry
use class_general
use class_mesh
use class_create_mesh
use class_parameters
  type node1
! --------------------------------------------
    integer                 :: n              ! Node number
    real                    :: x,y,z          ! Coordinates
! --------------------------------------------
  end type node1

  type mesh1
! --------------------------------------------
    integer                 :: n              ! Node number
    integer, dimension(1:8) :: nodes
    integer, dimension(1:6) :: faces
! --------------------------------------------
  end type mesh1

  type face1
! --------------------------------------------
    integer                 :: n              ! Node number
    integer, dimension(1:4) :: nodes
    integer, dimension(1:2) :: shared_meshes
    logical                 :: added_to_BC=.FALSE.
! --------------------------------------------
  end type face1

  type BC1
! --------------------------------------------
    integer                           :: n              ! Node number
    integer                           :: nfaces
    integer,allocatable,dimension(:)  :: faces
! --------------------------------------------
  end type BC1
contains
  
subroutine create_cube
  type(node1),allocatable,dimension(:,:,:)  :: nodes
  type(mesh1),allocatable,dimension(:,:,:)  :: meshes
  type(face1),allocatable,dimension(:,:,:)  :: faces_x
  type(face1),allocatable,dimension(:,:,:)  :: faces_y
  type(face1),allocatable,dimension(:,:,:)  :: faces_z

  type(node1),allocatable,dimension(:)      :: nodes2
  type(face1),allocatable,dimension(:)      :: faces2
  type(mesh1),allocatable,dimension(:)      :: meshes2
  integer                                   :: NBCs=7
  integer                                   :: NFACES,NMESHES,NNODES
  type(BC1)  ,dimension(1:7)                :: BCs
  integer                                   :: i,j,k,l,m,counter_node,counter_faces,counter_meshes,counter_BCs
  real                                      :: dx,dy,dz

  dx=Lx/NMESH_i
  dy=Ly/NMESH_j
  dz=Lz/NMESH_k

  NMESHES=NMESH_i*NMESH_j*NMESH_k
  NNODES  =(NMESH_i+1)*(NMESH_j+1)*(NMESH_k+1)
  NFACES  =(NMESH_i+1)*NMESH_j*NMESH_k+NMESH_i*(NMESH_j+1)*NMESH_k+NMESH_i*NMESH_j*(NMESH_k+1)


  Allocate(nodes2(1:NNODE))
  Allocate(faces2(1:NFACE))
  Allocate(meshes2(1:NMESH))
  Allocate(nodes(1:NMESH_i+1,1:NMESH_j+1,1:NMESH_k+1))
  Allocate(meshes(1:NMESH_i,1:NMESH_j,1:NMESH_k))
  Allocate(faces_x(1:NMESH_i+1,1:NMESH_j,1:NMESH_k))
  Allocate(faces_y(1:NMESH_i,1:NMESH_j+1,1:NMESH_k))
  Allocate(faces_z(1:NMESH_i,1:NMESH_j,1:NMESH_k+1))

  BCs(2)%nfaces=NMESH_j*NMESH_k
  BCs(3)%nfaces=NMESH_j*NMESH_k
  BCs(4)%nfaces=NMESH_i*NMESH_k
  BCs(5)%nfaces=NMESH_i*NMESH_k
  BCs(6)%nfaces=NMESH_i*NMESH_j
  BCs(7)%nfaces=NMESH_i*NMESH_j
  BCs(1)%nfaces=NFACE-2*(NMESH_j*NMESH_k+NMESH_i*NMESH_k+NMESH_i*NMESH_j)


! Allocate faces arrays inside the bcs class
  do i=1,NBCs
    Allocate(BCs(i)%faces(1:BCs(i)%nfaces))
  enddo

! ================== Create nodes ===========================
  counter_node=0
  do k=1,NMESH_k+1
     do j=1,NMESH_j+1
        do i=1,NMESH_i+1
           counter_node=counter_node+1
           nodes(i,j,k)%x=startx+(i-1)*dx
           nodes(i,j,k)%y=starty+(j-1)*dy
           nodes(i,j,k)%z=startz+(k-1)*dz
           nodes(i,j,k)%n=counter_node
           nodes2(counter_node)=nodes(i,j,k)
        enddo
     enddo
  enddo

! ================= Create faces ============================
! ------- start with x faces ----> function of theta
  counter_faces=0     ! This is for all faces
  do k=1,NMESH_k
     do j=1,NMESH_j
        do i=1,NMESH_i+1
           counter_faces=counter_faces+1
           faces_x(i,j,k)%nodes(1)=nodes(i,j,k)%n
           faces_x(i,j,k)%nodes(2)=nodes(i,j+1,k)%n
           faces_x(i,j,k)%nodes(3)=nodes(i,j+1,k+1)%n
           faces_x(i,j,k)%nodes(4)=nodes(i,j,k+1)%n
           faces_x(i,j,k)%n=counter_faces
           faces2(counter_faces)=faces_x(i,j,k)
        enddo
     enddo
  enddo
! FACES y
  do k=1,NMESH_k
     do j=1,NMESH_j+1
        do i=1,NMESH_i
           counter_faces=counter_faces+1
           faces_y(i,j,k)%nodes(1)=nodes(i,j,k)%n
           faces_y(i,j,k)%nodes(2)=nodes(i+1,j,k)%n
           faces_y(i,j,k)%nodes(3)=nodes(i+1,j,k+1)%n
           faces_y(i,j,k)%nodes(4)=nodes(i,j,k+1)%n
           faces_y(i,j,k)%n=counter_faces
           faces2(counter_faces)=faces_y(i,j,k)
        enddo
     enddo
  enddo

! FACES z
  do k=1,NMESH_k+1
     do j=1,NMESH_j
        do i=1,NMESH_i
           counter_faces=counter_faces+1
           faces_z(i,j,k)%nodes(1)=nodes(i,j,k)%n
           faces_z(i,j,k)%nodes(2)=nodes(i+1,j,k)%n
           faces_z(i,j,k)%nodes(3)=nodes(i+1,j+1,k)%n
           faces_z(i,j,k)%nodes(4)=nodes(i,j+1,k)%n
           faces_z(i,j,k)%n=counter_faces
           faces2(counter_faces)=faces_z(i,j,k)
        enddo
     enddo
  enddo
! ============================= Create MESHES =======================
  counter_meshes=0
  do k=1,NMESH_k
     do j=1,NMESH_j
        do i=1,NMESH_i
           counter_meshes=counter_meshes+1
           meshes(i,j,k)%nodes(1)=nodes(i,j,k)%n
           meshes(i,j,k)%nodes(2)=nodes(i+1,j,k)%n
           meshes(i,j,k)%nodes(3)=nodes(i+1,j+1,k)%n
           meshes(i,j,k)%nodes(4)=nodes(i,j+1,k)%n
           meshes(i,j,k)%nodes(5)=nodes(i,j,k+1)%n
           meshes(i,j,k)%nodes(6)=nodes(i+1,j,k+1)%n
           meshes(i,j,k)%nodes(7)=nodes(i+1,j+1,k+1)%n
           meshes(i,j,k)%nodes(8)=nodes(i,j+1,k+1)%n
           meshes(i,j,k)%n=counter_meshes
           meshes2(counter_meshes)=meshes(i,j,k)
!        print*, meshes(1,j,k)%n,counter_meshes

        enddo
     enddo
  enddo
! ============================= Create BCs =====================
! BC 2 Bottom BC ----------> x-faces x=0
  counter_BCs=0
  BCs(2)%n=2
  do k=1,NMESH_k
     do j=1,NMESH_j
        counter_BCs=counter_BCs+1
        BCs(2)%faces(counter_BCs)=faces_x(1,j,k)%n
        faces2(faces_x(1,j,k)%n)%added_to_BC=.TRUE.
        faces2(faces_x(1,j,k)%n)%shared_meshes(1)= meshes(1,j,k)%n
        faces2(faces_x(1,j,k)%n)%shared_meshes(2)= 0
!        print*, meshes(1,j,k)%n
     enddo
  enddo

! ============================= Create BCs =====================
! BC 3 Bottom BC ----------> x-faces x=Lx
  counter_BCs=0
  BCs(3)%n=3
  do k=1,NMESH_k
     do j=1,NMESH_j
        counter_BCs=counter_BCs+1
        BCs(3)%faces(counter_BCs)=faces_x(NMESH_i+1,j,k)%n
        faces2(faces_x(NMESH_i+1,j,k)%n)%added_to_BC=.TRUE.
        faces2(faces_x(NMESH_i+1,j,k)%n)%shared_meshes(1)= meshes(NMESH_i,j,k)%n
        faces2(faces_x(NMESH_i+1,j,k)%n)%shared_meshes(2)= 0
     enddo
  enddo


! ============================= Create BCs =====================
! BC 4 Bottom BC ----------> y-faces y=0
  counter_BCs=0
  BCs(4)%n=4
  do k=1,NMESH_k
     do i=1,NMESH_i
        counter_BCs=counter_BCs+1
        BCs(4)%faces(counter_BCs)=faces_y(i,1,k)%n
        faces2(faces_y(i,1,k)%n)%added_to_BC=.TRUE.
        faces2(faces_y(i,1,k)%n)%shared_meshes(1)= meshes(i,1,k)%n
        faces2(faces_y(i,1,k)%n)%shared_meshes(2)= 0
     enddo
  enddo

! ============================= Create BCs =====================
! BC 5 Bottom BC ----------> y-faces y=Ly
  counter_BCs=0
  BCs(5)%n=5
  do k=1,NMESH_k
     do i=1,NMESH_i
        counter_BCs=counter_BCs+1
        BCs(5)%faces(counter_BCs)=faces_y(i,NMESH_j+1,k)%n
        faces2(faces_y(i,NMESH_j+1,k)%n)%added_to_BC=.TRUE.
        faces2(faces_y(i,NMESH_j+1,k)%n)%shared_meshes(1)= meshes(i,NMESH_j,k)%n
        faces2(faces_y(i,NMESH_j+1,k)%n)%shared_meshes(2)= 0
     enddo
  enddo

! ============================= Create BCs =====================
! BC 6 Bottom BC ----------> z-faces z=0
  counter_BCs=0
  BCs(6)%n=6
  do j=1,NMESH_j
     do i=1,NMESH_i
        counter_BCs=counter_BCs+1
        BCs(6)%faces(counter_BCs)=faces_z(i,j,1)%n
        faces2(faces_z(i,j,1)%n)%added_to_BC=.TRUE.
        faces2(faces_z(i,j,1)%n)%shared_meshes(1)= meshes(i,j,1)%n
        faces2(faces_z(i,j,1)%n)%shared_meshes(2)= 0
     enddo
  enddo

! ============================= Create BCs =====================
! BC 7 Bottom BC ----------> z-faces z=Lz
  counter_BCs=0
  BCs(7)%n=7
  do j=1,NMESH_j
     do i=1,NMESH_i
        counter_BCs=counter_BCs+1
        BCs(7)%faces(counter_BCs)=faces_z(i,j,NMESH_k+1)%n
        faces2(faces_z(i,j,NMESH_k+1)%n)%added_to_BC=.TRUE.
        faces2(faces_z(i,j,NMESH_k+1)%n)%shared_meshes(1)= meshes(i,j,NMESH_k)%n
        faces2(faces_z(i,j,NMESH_k+1)%n)%shared_meshes(2)= 0
     enddo
  enddo


! Everywhere else
! BC 1
  counter_BCs=0
  BCs(1)%n=1
  do k=1,NMESH_k
     do j=1,NMESH_j
        do i=2,NMESH_i
          counter_BCs=counter_BCs+1
          BCs(1)%faces(counter_BCs)=faces_x(i,j,k)%n
          faces2(faces_x(i,j,k)%n)%added_to_BC=.TRUE.
          faces2(faces_x(i,j,k)%n)%shared_meshes(1)= meshes(i,j,k)%n
          faces2(faces_x(i,j,k)%n)%shared_meshes(2)= meshes(i-1,j,k)%n
        enddo
     enddo
  enddo
  do k=1,NMESH_k
     do j=2,NMESH_j
        do i=1,NMESH_i
          counter_BCs=counter_BCs+1
          BCs(1)%faces(counter_BCs)=faces_y(i,j,k)%n
          faces2(faces_y(i,j,k)%n)%added_to_BC=.TRUE.
          faces2(faces_y(i,j,k)%n)%shared_meshes(1)= meshes(i,j,k)%n
          faces2(faces_y(i,j,k)%n)%shared_meshes(2)= meshes(i,j-1,k)%n
        enddo
     enddo
  enddo

  do k=2,NMESH_k
     do j=1,NMESH_j
        do i=1,NMESH_i
          counter_BCs=counter_BCs+1
          BCs(1)%faces(counter_BCs)=faces_z(i,j,k)%n
          faces2(faces_z(i,j,k)%n)%added_to_BC=.TRUE.
          faces2(faces_z(i,j,k)%n)%shared_meshes(1)= meshes(i,j,k)%n
          faces2(faces_z(i,j,k)%n)%shared_meshes(2)= meshes(i,j,k-1)%n
        enddo
     enddo
  enddo


! Check all BCs are defined
    do i=1,NFACE
       if(.NOT. faces2(i)%added_to_BC ) then
          print*, 'ERROR',i
       endif
    enddo



   call fill_global_geometry(nodes2,faces2,BCs,NNODES,NFACES,NBCs)
    do i=1,7
      DEALLOCATE(BCs(i)%faces)
    enddo

  deallocate(nodes2)
  deallocate(faces2)
  deallocate(meshes2)
  deallocate(nodes)
  deallocate(meshes)
  deallocate(faces_x)
  deallocate(faces_y)
  deallocate(faces_z)

end subroutine create_cube


  subroutine create_hollow(R_i,R_o,H,cx,cy,cz,n1,n2,n3)
    real, intent(in)                         :: R_i,R_o,H,cx,cy,cz
    integer,intent(in)                       :: n1,n2,n3

    type(node1),allocatable,dimension(:,:,:)     :: nodes
    type(mesh1),allocatable,dimension(:,:,:)     :: meshes
    type(face1),allocatable,dimension(:,:,:)     :: faces_x
    type(face1),allocatable,dimension(:,:,:)     :: faces_y
    type(face1),allocatable,dimension(:,:,:)     :: faces_z


    type(node1),allocatable,dimension(:)      :: nodes2
    type(face1),allocatable,dimension(:)      :: faces2
    type(mesh1),allocatable,dimension(:)      :: meshes2
! For hollow cylinderical geometries, we have 5 boundaries, 1: is no bc
    type(BC1)  ,dimension(1:5)                :: BCs
	integer                                  :: NBCs=5
    integer                                  :: NFACES,NMESHES,NNODES
    real                                     :: PI,dtheta,dr,dz,tmp_r,tmp_theta
    integer                                  :: i,j,k,l,m,counter_node,counter_faces,counter_meshes,counter_BCs
    print*, 'ddddddddddd'
! ------------------ Dimensions evaluation ------------------
    PI=4.D0*DATAN(1.D0)    
    dr=(R_o-R_i)/n2
    dtheta=2*PI/(n1)
    dz=H/n3
! -----------------------------------------------------------
    NMESHES=n1*n2*n3
    NNODES =n1*(n2+1)*(n3+1)
    NFACES=3*n1*n2*n3+(n1*n3)+(n1*n2)
! -----------------------------------------------------------
    Allocate(nodes2(1:NNODES))
    Allocate(faces2(1:NFACES))
    Allocate(meshes2(1:NMESHES))
    Allocate(nodes(1:n1,1:n2+1,1:n3+1))
    Allocate(meshes(1:n1,1:n2,1:n3))
    Allocate(faces_x(1:n1,1:n2,1:n3))
    Allocate(faces_y(1:n1,1:n2+1,1:n3))
    Allocate(faces_z(1:n1,1:n2,1:n3+1))



    BCs(1)%nfaces=NFACES-2*n1*n2-2*n1*n3
    BCs(2)%nfaces=n1*n2
    BCs(3)%nfaces=n1*n2
    BCs(4)%nfaces=n1*n3
    BCs(5)%nfaces=n1*n3
! Allocating number of faces in each boundary
    Allocate(BCs(1)%faces(1:BCs(1)%nfaces))
    Allocate(BCs(2)%faces(1:BCs(2)%nfaces)) ! bottom
    Allocate(BCs(3)%faces(1:BCs(3)%nfaces)) ! TOP
    Allocate(BCs(4)%faces(1:BCs(4)%nfaces)) ! Inside
    Allocate(BCs(5)%faces(1:BCs(5)%nfaces)) ! Outside
! -----------------------------------------------------------

! ================== Create nodes ===========================
    counter_node=0
    do k=1,n3+1
       do j=1,n2+1
          do i=1,n1
             counter_node=counter_node+1
             tmp_r=R_i+(j-1)*dr
             tmp_theta=(i-1)*dtheta
             nodes(i,j,k)%x=tmp_r*cos(tmp_theta)
             nodes(i,j,k)%y=tmp_r*sin(tmp_theta)
             nodes(i,j,k)%z=(k-1)*dz-H/2.
             nodes(i,j,k)%n=counter_node
             nodes2(counter_node)=nodes(i,j,k)
!             print*, nodes(i,j,k)%n,nodes(i,j,k)%x,nodes(i,j,k)%y,nodes(i,j,k)%z
          enddo
       enddo
    enddo

! ================= Create faces ============================
! ------- start with x faces ----> function of theta
    counter_faces=0 	! This is for all faces
    do k=1,n3
       do j=1,n2
          do i=1,n1
             counter_faces=counter_faces+1
             faces_x(i,j,k)%nodes(1)=nodes(i,j,k)%n
             faces_x(i,j,k)%nodes(2)=nodes(i,j+1,k)%n
             faces_x(i,j,k)%nodes(3)=nodes(i,j+1,k+1)%n
             faces_x(i,j,k)%nodes(4)=nodes(i,j,k+1)%n
             faces_x(i,j,k)%n=counter_faces
             faces2(counter_faces)=faces_x(i,j,k)
          enddo
       enddo
    enddo
! FACES y
    do k=1,n3
       do j=1,n2+1
          do i=1,n1
             counter_faces=counter_faces+1
             faces_y(i,j,k)%nodes(1)=nodes(i,j,k)%n
             if (i<n1) then
             faces_y(i,j,k)%nodes(2)=nodes(i+1,j,k)%n
             faces_y(i,j,k)%nodes(3)=nodes(i+1,j,k+1)%n
             else
             faces_y(i,j,k)%nodes(2)=nodes(1,j,k)%n
             faces_y(i,j,k)%nodes(3)=nodes(1,j,k+1)%n
             endif
             faces_y(i,j,k)%nodes(4)=nodes(i,j,k+1)%n
             faces_y(i,j,k)%n=counter_faces
             faces2(counter_faces)=faces_y(i,j,k)
          enddo
       enddo
    enddo
! FACES z
    do k=1,n3+1
       do j=1,n2
          do i=1,n1
             counter_faces=counter_faces+1
             faces_z(i,j,k)%nodes(1)=nodes(i,j,k)%n
             if (i<n1) then
             faces_z(i,j,k)%nodes(2)=nodes(i+1,j,k)%n
             faces_z(i,j,k)%nodes(3)=nodes(i+1,j+1,k)%n
             else
             faces_z(i,j,k)%nodes(2)=nodes(1,j,k)%n
             faces_z(i,j,k)%nodes(3)=nodes(1,j+1,k)%n
             endif
             faces_z(i,j,k)%nodes(4)=nodes(i,j+1,k)%n
             faces_z(i,j,k)%n=counter_faces
             faces2(counter_faces)=faces_z(i,j,k)
          enddo
       enddo
    enddo
! ============================= Create MESHES =======================
    counter_meshes=0
    do k=1,n3
       do j=1,n2
          do i=1,n1
             counter_meshes=counter_meshes+1
             meshes(i,j,k)%nodes(1)=nodes(i,j,k)%n
             meshes(i,j,k)%nodes(4)=nodes(i,j+1,k)%n
             meshes(i,j,k)%nodes(5)=nodes(i,j,k+1)%n
             meshes(i,j,k)%nodes(8)=nodes(i,j+1,k+1)%n
             if (i<n1) then
             meshes(i,j,k)%nodes(2)=nodes(i+1,j,k)%n
             meshes(i,j,k)%nodes(3)=nodes(i+1,j+1,k)%n
             meshes(i,j,k)%nodes(6)=nodes(i+1,j,k+1)%n
             meshes(i,j,k)%nodes(7)=nodes(i+1,j+1,k+1)%n
             else
             meshes(i,j,k)%nodes(2)=nodes(1,j,k)%n
             meshes(i,j,k)%nodes(3)=nodes(1,j+1,k)%n
             meshes(i,j,k)%nodes(6)=nodes(1,j,k+1)%n
             meshes(i,j,k)%nodes(7)=nodes(1,j+1,k+1)%n
             endif
             meshes(i,j,k)%n=counter_meshes
             meshes2(counter_meshes)=meshes(i,j,k)
          enddo
       enddo
    enddo
! ============================= Create BCs =====================
! BC 2 Bottom BC ----------> z-faces z=0
  counter_BCs=0
  BCs(2)%n=2
  do j=1,n2
     do i=1,n1
        counter_BCs=counter_BCs+1
        BCs(2)%faces(counter_BCs)=faces_z(i,j,1)%n
        faces2(faces_z(i,j,1)%n)%added_to_BC=.TRUE.
        faces2(faces_z(i,j,1)%n)%shared_meshes(1)= meshes(i,j,1)%n
        faces2(faces_z(i,j,1)%n)%shared_meshes(2)= 0
!        print*,faces2(BCs(2)%faces(counter_BCs))%nodes,faces2(BCs(2)%faces(counter_BCs))%shared_meshes
     enddo
  enddo

! BC 3 Bottom BC ----------> z-faces z=n3+1
  counter_BCs=0
  BCs(3)%n=3
  do j=1,n2
     do i=1,n1
        counter_BCs=counter_BCs+1
        BCs(3)%faces(counter_BCs)=faces_z(i,j,n3+1)%n
        faces2(faces_z(i,j,n3+1)%n)%added_to_BC=.TRUE.
        faces2(faces_z(i,j,n3+1)%n)%shared_meshes(1)= meshes(i,j,n3)%n
        faces2(faces_z(i,j,n3+1)%n)%shared_meshes(2)= 0
!        print*,faces2(BCs(3)%faces(counter_BCs))%nodes,faces2(BCs(3)%faces(counter_BCs))%shared_meshes
     enddo
  enddo

! BC 4 Bottom BC ----------> y-faces y=R_i
  counter_BCs=0
  BCs(4)%n=4
  do k=1,n3
     do i=1,n1
        counter_BCs=counter_BCs+1
        BCs(4)%faces(counter_BCs)=faces_y(i,1,k)%n
        faces2(faces_y(i,1,k)%n)%added_to_BC=.TRUE.
        faces2(faces_y(i,1,k)%n)%shared_meshes(1)= meshes(i,1,k)%n
        faces2(faces_y(i,1,k)%n)%shared_meshes(2)= 0
!        print*,faces2(BCs(4)%faces(counter_BCs))%nodes,faces2(BCs(4)%faces(counter_BCs))%shared_meshes
     enddo
  enddo


! BC 5 Bottom BC ----------> y-faces y=R_o
  counter_BCs=0
  BCs(5)%n=5
  do k=1,n3
     do i=1,n1
        counter_BCs=counter_BCs+1
        BCs(5)%faces(counter_BCs)=faces_y(i,n2+1,k)%n
        faces2(faces_y(i,n2+1,k)%n)%added_to_BC=.TRUE.
        faces2(faces_y(i,n2+1,k)%n)%shared_meshes(1)= meshes(i,n2,k)%n
        faces2(faces_y(i,n2+1,k)%n)%shared_meshes(2)= 0
!        print*,faces2(BCs(5)%faces(counter_BCs))%nodes,faces2(BCs(5)%faces(counter_BCs))%shared_meshes
     enddo
  enddo

! BC 1 everywhere else
! first x faces
  counter_BCs=0
  BCs(1)%n=1
   do k=1,n3
       do j=1,n2
          do i=1,n1
             if(.NOT. faces2(faces_x(i,j,k)%n)%added_to_BC ) then
               counter_BCs=counter_BCs+1
               BCs(1)%faces(counter_BCs)=faces_x(i,j,k)%n
               faces2(faces_x(i,j,k)%n)%added_to_BC=.TRUE.
               faces2(faces_x(i,j,k)%n)%shared_meshes(1)= meshes(i,j,k)%n
               if (i ==1) then 
                  faces2(faces_x(i,j,k)%n)%shared_meshes(2)= meshes(n1,j,k)%n
               else
                  faces2(faces_x(i,j,k)%n)%shared_meshes(2)= meshes(i-1,j,k)%n
               endif
!               print*,faces2(BCs(1)%faces(counter_BCs))%nodes,faces2(BCs(1)%faces(counter_BCs))%shared_meshes
             endif
          enddo
       enddo
    enddo

! second y faces
!  counter_BCs=0
  BCs(1)%n=1
   do k=1,n3
       do j=2,n2
          do i=1,n1
             if(.NOT. faces2(faces_y(i,j,k)%n)%added_to_BC ) then
               counter_BCs=counter_BCs+1
               BCs(1)%faces(counter_BCs)=faces_y(i,j,k)%n
               faces2(faces_y(i,j,k)%n)%added_to_BC=.TRUE.
               faces2(faces_y(i,j,k)%n)%shared_meshes(1)= meshes(i,j,k)%n
               faces2(faces_y(i,j,k)%n)%shared_meshes(2)= meshes(i,j-1,k)%n
!               print*,faces2(BCs(1)%faces(counter_BCs))%nodes,faces2(BCs(1)%faces(counter_BCs))%shared_meshes
             endif
          enddo
       enddo
    enddo
! second z faces
!  counter_BCs=0
  BCs(1)%n=1
   do k=2,n3
       do j=1,n2
          do i=1,n1
             if(.NOT. faces2(faces_z(i,j,k)%n)%added_to_BC ) then
               counter_BCs=counter_BCs+1
               BCs(1)%faces(counter_BCs)=faces_z(i,j,k)%n
               faces2(faces_z(i,j,k)%n)%added_to_BC=.TRUE.
               faces2(faces_z(i,j,k)%n)%shared_meshes(1)= meshes(i,j,k)%n
               faces2(faces_z(i,j,k)%n)%shared_meshes(2)= meshes(i,j,k-1)%n
!               print*,faces2(BCs(1)%faces(counter_BCs))%nodes,faces2(BCs(1)%faces(counter_BCs))%shared_meshes
             endif
          enddo
       enddo
    enddo
! Check all BCs are defined
    do i=1,NFACES
       if(.NOT. faces2(i)%added_to_BC ) then
          print*, 'ERROR',i
       endif    
    enddo    


    print*, 'HELLO MAN',n1,n2,n3,NFACES,NNODES,NMESHES

  do i=1,NNODES
    nodes2(i)%x=nodes2(i)%x+cx
    nodes2(i)%y=nodes2(i)%y+cy
    nodes2(i)%z=nodes2(i)%z+cz

  enddo


   call fill_global_geometry(nodes2,faces2,BCs,NNODES,NFACES,NBCs)
    do i=1,5
      DEALLOCATE(BCs(i)%faces)
    enddo

    DEALLOCATE(nodes2)
    DEAllocate(faces2)
    DEAllocate(meshes2)
    DEAllocate(nodes)
    DEAllocate(meshes)
    DEAllocate(faces_x)
    DEAllocate(faces_y)
    DEAllocate(faces_z)


  end subroutine create_hollow



subroutine fill_global_geometry(nodes2,faces2,BCs,NNODES,NFACES,NBCs)
! Fill global arrays
  type(node1),dimension(1:NNODES),intent(in) :: nodes2
  type(face1),dimension(1:NFACES),intent(in) :: faces2
  type(BC1)  ,dimension(1:NBCs)  ,intent(in) :: BCs
  integer, intent(in) :: NNODES,NFACES,NBCs
      integer,dimension(1:6)             :: temp2
	  integer                             :: cnt=1
    do i=1,NBCs
       BC_id(i)=i
       BCS_ARRAY(i)%id=i
       BCS_ARRAY(i)%n_faces=BCs(i)%nfaces
       if (i==1) then
           BCS_ARRAY(i)%first=1
           BCS_ARRAY(i)%last=BCS_ARRAY(i)%n_faces
       else
           BCS_ARRAY(i)%first=BCS_ARRAY(i-1)%last+1
           BCS_ARRAY(i)%last=BCS_ARRAY(i)%n_faces+BCS_ARRAY(i)%first-1
       endif
    enddo
    do i=1,NNODES
       call node_create(NODES_ARRAY(i),i,nodes2(i)%x,nodes2(i)%y,nodes2(i)%z)
    enddo
    call deform_nodes
	
    do i=1,NBCs
         do j=1,BCs(i)%nfaces
		  temp2(1:4)=faces2(BCs(i)%faces(j))%nodes
		  temp2(5:6)=faces2(BCs(i)%faces(j))%shared_meshes
          call face_create (FACES_ARRAY(cnt),temp2(1),temp2(2) ,temp2(3),temp2(4),cnt,BC_id(i))
          call face_set_shared(FACES_ARRAY(cnt),(/temp2(5),temp2(6)/))
          call mesh_add_face(MESHES_ARRAY(temp2(5)),cnt,temp2(5))
          call mesh_add_face(MESHES_ARRAY(temp2(6)),cnt,temp2(6))
          if (temp2(5) /= 0)  call mesh_add_neighbor(MESHES_ARRAY(temp2(5)),temp2(6))
          if (temp2(6) /= 0)  call mesh_add_neighbor(MESHES_ARRAY(temp2(6)),temp2(5))
		  cnt=cnt+1
       enddo
    enddo
end subroutine fill_global_geometry





end module class_create_geometry