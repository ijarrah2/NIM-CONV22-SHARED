! ===============================================================================
! ==================== READ MESH FILES OF nsh and msh fortmats ==================
! ===============================================================================

module class_read_mesh
use class_general
use class_mesh
use class_create_mesh
use class_parameters
implicit none
contains

! ============================== MAIN FUNCTION TO CALL ==========================
subroutine read_mesh_file(infile)
  implicit none
  CHARACTER(len=*), intent(in)     :: infile
  integer                          :: mesh_id          ! 1:.msh, 2:.nsh
  if (infile(len(TRIM(infile))-2:len(TRIM(infile))) == "nsh") call readnsh(TRIM(infile))
  if (infile(len(TRIM(infile))-2:len(TRIM(infile))) == "msh") call readmsh(TRIM(infile)) 
!  print*, "inside read_mesh"

end subroutine read_mesh_file


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++              READ MSH FILES                ++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


subroutine readmsh (infile)
    implicit none
! -----------------------------------------------------
    CHARACTER(len=*), intent(in)                 :: infile
! -----------------------------------------------------
    integer :: lines,i,io_stat,flag,pos1,pos2,n_entries,counter_bc,read_unit,BC_id
    integer,   allocatable  :: entries(:)
! flag is used to decide if the following section is nodes or faces,...etc
    character(len=200) :: line
    character(len=6),dimension(1:6)   :: tempstr
    real,dimension(1:6)                :: temp,temp_old
    integer,dimension(1:6)             :: temp2
    integer                            :: nd
    call read_numbers_mesh (infile)
    allocate(FACES_ARRAY(0:NFACE))
    allocate(NODES_ARRAY(0:NNODE))
    allocate(MESHES_ARRAY(0:NMESH))
    read_unit=10
    counter_bc=0

    open(unit=read_unit, file=infile, action='READ', status='OLD')
    do
        read(read_unit, '(A)', iostat=io_stat) line
        if (io_stat /= 0) exit
        if (line(1:1) == '(') THEN
           if(line(2:2) /= '') THEN
              flag=str2int(line(2:3))
              if (flag == 10) THEN
                 ! NODES READING
                 read(read_unit, '(A)', iostat=io_stat) line
                 read(read_unit, '(A)', iostat=io_stat) line
                 do i=1,NNODE
                   read(read_unit, '(A)', iostat=io_stat) line
                   read (line,*) temp(1),temp(2),temp(3)
                   call node_create(NODES_ARRAY(i),i,temp(1),temp(2),temp(3))
                 enddo
                 call deform_nodes
              else if (flag == 13) THEN
                 call locating_braces(line(2:len(line)),pos1,pos2)
                 read (line(pos1+1:pos2-1),*) tempstr(1),tempstr(2),tempstr(3),tempstr(4),tempstr(5)
                 read (tempstr(1),'(Z6)') BC_id
                 if (tempstr(1) /= '0') THEN
                     counter_bc=counter_bc+1
                     print*, 'ENTERED',BCS_ARRAY(counter_bc)%first
                     do i=BCS_ARRAY(counter_bc)%first,BCS_ARRAY(counter_bc)%last
                        read(read_unit, '(A)', iostat=io_stat) line
                        read (line,*) tempstr(1),tempstr(2),tempstr(3),tempstr(4),tempstr(5),tempstr(6)
                        read (tempstr(1),'(Z6)') temp2(1)
                        read (tempstr(2),'(Z6)') temp2(2)
                        read (tempstr(3),'(Z6)') temp2(3)
                        read (tempstr(4),'(Z6)') temp2(4)
                        read (tempstr(5),'(Z6)') temp2(5)
                        read (tempstr(6),'(Z6)') temp2(6)
                        call face_create (FACES_ARRAY(i),temp2(1),temp2(2) ,temp2(3),temp2(4),i,BC_id)
                        call face_set_shared(FACES_ARRAY(i),(/temp2(5),temp2(6)/))
                        call mesh_add_face(MESHES_ARRAY(temp2(5)),i,temp2(5))
                        call mesh_add_face(MESHES_ARRAY(temp2(6)),i,temp2(6))
                        if (temp2(5) /= 0)  call mesh_add_neighbor(MESHES_ARRAY(temp2(5)),temp2(6))
                        if (temp2(6) /= 0)  call mesh_add_neighbor(MESHES_ARRAY(temp2(6)),temp2(5))

                     enddo
                 endif
              endif
           endif
        endif
    end do

    close(read_unit)
end subroutine readmsh

! ===========================================================================================
! ==================================== READ NSH FILES =======================================
! ===========================================================================================
  subroutine readnsh (infile)
! reads file extension of nsh ----> temporary formate
    implicit none
! -----------------------------------------------------
    CHARACTER(len=*), intent(in)                 :: infile
! -----------------------------------------------------
    integer :: lines,i,io_stat,flag,pos1,pos2,n_entries,counter_bc,read_unit,j
! flag is used to decide if the following section is nodes or faces,...etc
    character(len=200) :: line
    character(len=6),dimension(1:6)   :: tempstr
    real,dimension(1:6)                :: temp
    integer,dimension(1:6)             :: temp2
    read_unit=19
    counter_bc=0

    open(unit=read_unit, file=infile, action='READ', status='OLD')
    read(read_unit, '(A)', iostat=io_stat) line
    read (line,*) temp(1),temp(2),temp(3),temp(4)
    NNODE=temp(1)
    NFACE=temp(2)
    NMESH=temp(3)
    NBC  =temp(4)
	print*, '=============================== Mesh information =========================='
    print*, 'NNODES =', NNODE
    print*, 'NFACES =', NFACE
    print*, 'NMESHES=', NMESH
    print*, 'NBCs   =', NBC
	print*, '==========================================================================='

    allocate(FACES_ARRAY(0:NFACE))
    allocate(NODES_ARRAY(0:NNODE))
    allocate(MESHES_ARRAY(0:NMESH))
    allocate(BCS_ARRAY(1:NBC))
    allocate(BC_id(1:NBC))
    do i=1,NBC
       read(read_unit, '(A)', iostat=io_stat) line
       read (line,*) temp(1),temp(2)
       BC_id(i)=temp(1)
       BCS_ARRAY(i)%id=i
       BCS_ARRAY(i)%n_faces=temp(2)
!       print*, 'BC # ', BC_id(i),BCS_ARRAY(i)%n_faces
       if (i==1) then
           BCS_ARRAY(i)%first=1
           BCS_ARRAY(i)%last=BCS_ARRAY(i)%n_faces
       else
           BCS_ARRAY(i)%first=BCS_ARRAY(i-1)%last+1
           BCS_ARRAY(i)%last=BCS_ARRAY(i)%n_faces+BCS_ARRAY(i)%first-1
       endif
!       print*, 'first', BCS_ARRAY(i)%first,'last', BCS_ARRAY(i)%last
    enddo


    read(read_unit, '(A)', iostat=io_stat) line
    do i=1,NNODE
       read(read_unit, '(A)', iostat=io_stat) line
       read (line,*) temp(1),temp(2),temp(3)
       call node_create(NODES_ARRAY(i),i,temp(1),temp(2),temp(3))
    enddo
    call deform_nodes
    do i=1,NBC
       read(read_unit, '(A)', iostat=io_stat) line
       do j=BCS_ARRAY(i)%first,BCS_ARRAY(i)%last
          read(read_unit, '(A)', iostat=io_stat) line
          read (line,*) temp2(1),temp2(2),temp2(3),temp2(4),temp2(5),temp2(6)

          call face_create (FACES_ARRAY(j),temp2(1),temp2(2) ,temp2(3),temp2(4),j,BC_id(i))
          call face_set_shared(FACES_ARRAY(j),(/temp2(5),temp2(6)/))
          call mesh_add_face(MESHES_ARRAY(temp2(5)),j,temp2(5))
          call mesh_add_face(MESHES_ARRAY(temp2(6)),j,temp2(6))
          if (temp2(5) /= 0)  call mesh_add_neighbor(MESHES_ARRAY(temp2(5)),temp2(6))
          if (temp2(6) /= 0)  call mesh_add_neighbor(MESHES_ARRAY(temp2(6)),temp2(5))
       enddo
    enddo


    close(read_unit)

  end subroutine readnsh



subroutine read_numbers_mesh (infile)
    implicit none
! -----------------------------------------------------
    CHARACTER(len=*), intent(in)                 :: infile
! -----------------------------------------------------
    integer :: lines,i,io_stat,flag,pos1,pos2,n_entries,temp,counter_bc,id,first,last,n_faces
    integer,   allocatable  :: entries(:)
! flag is used to decide if the following section is nodes or faces,...etc
    character(len=200) :: line
    character(len=6),dimension(1:5)   :: tempstr
    NBC=0
    counter_bc=0
    call bcs_number (infile)
    allocate(BCS_ARRAY(1:NBC))
    allocate(BC_id(1:NBC))
    lines = number_of_lines(infile)
    open(unit=2, file=infile, action='READ', status='OLD')
    do
        read(2, '(A)', iostat=io_stat) line
        if (io_stat /= 0) exit
        if (line(1:1) == '(') THEN
           if(line(2:2) /= '') THEN
              flag=str2int(line(2:3))
              if (flag == 10) THEN
                 call locating_braces(line(2:len(line)),pos1,pos2)
                   read (line(pos1+1:pos2-1),*) tempstr(1),tempstr(2),tempstr(3),tempstr(4),tempstr(5)
                   if (tempstr(1) == '0') THEN
                     read (tempstr(3),'(Z6)') NNODE
                     print*, 'n_nodes =',NNODE
                   endif
              else if (flag == 12) THEN
                 call locating_braces(line(2:len(line)),pos1,pos2)
                   read (line(pos1+1:pos2-1),*) tempstr(1),tempstr(2),tempstr(3),tempstr(4),tempstr(5)
                   if (tempstr(1) == '0') THEN
                     read (tempstr(3),'(Z6)') NMESH
                     print*, 'n_meshes = ',NMESH
                   endif
              else if (flag == 13) THEN
                 call locating_braces(line(2:len(line)),pos1,pos2)
                   read (line(pos1+1:pos2-1),*) tempstr(1),tempstr(2),tempstr(3),tempstr(4),tempstr(5)
                   if (tempstr(1) == '0') THEN
                     read (tempstr(3),'(Z6)') NFACE
                     print*, 'n_faces = ',NFACE
                  else
                     counter_bc=counter_bc+1
                     call locating_braces(line(2:len(line)),pos1,pos2)
                     read (line(pos1+1:pos2-1),*) tempstr(1),tempstr(2),tempstr(3),tempstr(4),tempstr(5)
                     read (tempstr(1),'(Z6)') id
                     read (tempstr(2),'(Z6)') first
                     read (tempstr(3),'(Z6)') last
                     call BC_create(BCS_ARRAY(counter_bc),counter_bc,id,last-first+1,first,last)
                     call BC_print(BCS_ARRAY(counter_bc))
                     BC_id(counter_bc)=id
                 endif
              endif
           endif
        endif
    enddo
    close(2)
end subroutine read_numbers_mesh

subroutine locating_braces(line,pos1,pos2)
  implicit none
    integer, intent (inout)         :: pos1, pos2
    CHARACTER(len=*), intent(in) :: line
    integer :: length,i
    length=len(line)
    do i=1,length-1
       if(line(i:i) == '(') THEN
          pos1=i+1
       elseif (line(i:i) == ')') THEN
          pos2=i+1
          EXIT
       endif
    enddo
end subroutine locating_braces
! THIS FUNCTION READS NUMBER OF BCS IN THE MESH FILE (FACE BC)
subroutine bcs_number (infile)
    implicit none
! -----------------------------------------------------
    CHARACTER(len=*), intent(in)                 :: infile
! -----------------------------------------------------
    integer :: i,io_stat,flag,n_entries,temp,pos1,pos2
    integer,   allocatable  :: entries(:)
    character(len=200) :: line
    character(len=6),dimension(1:5)   :: tempstr
    NBC=0
    open(unit=3, file=infile, action='READ', status='OLD')
    do
        read(3, '(A)', iostat=io_stat) line
        if (io_stat /= 0) exit
        if (line(1:1) == '(') THEN
           if(line(2:2) /= '') THEN
              flag=str2int(line(2:3))
               if (flag == 13) THEN
                  call locating_braces(line(2:len(line)),pos1,pos2)
                   read (line(pos1+1:pos2-1),*) tempstr(1)
                 if (tempstr(1) /= '0') THEN
                     NBC=NBC+1            ! To count number of surfaces that have BCs
                 endif
              endif
           endif
        endif
    end do
    close(3)
end subroutine bcs_number


end module class_read_mesh
