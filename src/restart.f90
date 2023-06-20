module class_restart
use class_mesh
use class_general
use class_parameters
implicit none
contains

subroutine write_restart
    implicit none
	character (35) :: restart_file
	character (20) :: rest_id

	write (rest_id, *) STEP
	!Create a dynamic file name
	restart_file = "restart/rest.step"//trim(adjustl(rest_id))
	print*,"***************************************************************************"
	print*,"Writing restart file at time = ", TI, "|     Time step = ", STEP
	print*,"***************************************************************************"
	OPEN(111,FILE=restart_file,FORM='UNFORMATTED')
	REWIND(111)
	write(111) TI
	write(111) STEP
	write(111) faces_u
	write(111) faces_v
	write(111) faces_w
	write(111) faces_p
	write(111) faces_T
	write(111) meshes_uxyz
	write(111) meshes_vxyz
	write(111) meshes_wxyz
	write(111) meshes_Txyz
	write(111) meshes_uxyz_old
	write(111) meshes_vxyz_old
	write(111) meshes_wxyz_old
	write(111) meshes_Txyz_old
	write(111) meshes_dila_old
!	print*, faces_u
	close(111)
end subroutine write_restart


subroutine read_restart
    implicit none
	character (35) :: restart_file
	character (20) :: rest_id

	write (rest_id, *) read_restart_step
	!Create a dynamic file name
	restart_file = "restart/rest.step"//trim(adjustl(rest_id))
	OPEN(111,FILE=restart_file,FORM='UNFORMATTED')
	REWIND(111)
	read(111) TI
	read(111) STEP
	read(111) faces_u
	read(111) faces_v
	read(111) faces_w
	read(111) faces_p
	read(111) faces_T
	read(111) meshes_uxyz
	read(111) meshes_vxyz
	read(111) meshes_wxyz
	read(111) meshes_Txyz
	read(111) meshes_uxyz_old
	read(111) meshes_vxyz_old
	read(111) meshes_wxyz_old
	read(111) meshes_Txyz_old
	read(111) meshes_dila
	close(111)
end subroutine read_restart



end module class_restart
