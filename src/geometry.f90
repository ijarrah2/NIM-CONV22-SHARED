! This class contains all the constants associated with the convection-diffusion equation
module class_geometry
use class_general
use class_parameters
use class_mesh
use class_curvilinear
use class_create_mesh
use class_prepare_mesh
implicit none

contains

! ====================== This should be called once ==============================
subroutine evaluate_all_mesh_constants
  implicit none
  integer :: i

!$OMP     PARALLEL DO PRIVATE (i)
  do i=1,NMESH
     call mesh_calculate_mapping_coefficients(MESHES_ARRAY(i))
     call volume_calculate(MESHES_ARRAY(i))
	 call find_normal(i)
	 call calculate_defusion_tensor_over_faces(i)
	 call convection_diffusion_constants_mesh_i(i)
!	 call find_normal2(i)
!	 call fill_normals(i)
  enddo

!$OMP     END PARALLEL DO


end subroutine evaluate_all_mesh_constants



subroutine fill_normals(i)
  implicit none
  integer,intent(in) :: i
  integer            :: fn
  
  do fn=1,6
	faces_normal(:,meshes_faces(fn,i))=meshes_normal_faces(:,fn,i)
	faces_tangential(:,meshes_faces(fn,i))=meshes_tangential_faces(:,fn,i)
	faces_bitangential(:,meshes_faces(fn,i))=meshes_binormal_faces(:,fn,i)
	faces_mid_coordinates(:,meshes_faces(fn,i))=meshes_mid_faces(:,fn,i)
  enddo

  
end subroutine fill_normals

subroutine evaluate_time_dependent_quantities
  implicit none
  integer :: i

  if(.not. steady_state .and. (is_bx_TD .or. is_by_TD .or. is_bz_TD .or. is_analytical_TD)) then
  !$OMP     PARALLEL DO PRIVATE (i)
		do i=1,NMESH
!		   call time_dependent_face_averages(i)
!		   call time_dependent_mesh_averages(i)
		enddo
  !$OMP     END PARALLEL DO
  endif

end subroutine evaluate_time_dependent_quantities


subroutine evaluate_time_dependent_quantities_i(M)
  implicit none
  type(mesh), intent(in)        :: M



end subroutine evaluate_time_dependent_quantities_i


subroutine find_normal(i)
! Outward normals
  implicit none
  integer, intent(in) :: i
  real                :: x_i,y_i,z_i,nx,ny,nz,t1x,t1y,t1z,nn
  type(curvi3d)       :: Curvidata
  real, dimension(1:3):: t_1, t_2, n, bn, n2

! ============================= IMPORTANT ----> NORMAL IS LOCATION DEPENDENT --> FACES ARE NOT COPLANAR ===================
! ======== HERE WE EVALUATE THE NORMAL AT THE CENTER OF EACH FACE ---> MAYBE WE NEED TO FIND THE AVERAGE OF THE SURFACE ===

  ! left x_i=-1.0
  x_i=-1.0
  y_i= 0.0
  z_i= 0.0
  call evaluate_curvilinear(MESHES_ARRAY(i),curvidata,x_i,y_i,z_i)

  n=curvidata%g1/sqrt(curvidata%gij(1,1))
  t_1=curvidata%g_2/sqrt(curvidata%g_ij(2,2))
  t_2=cross(n,t_1)
  t_2=abs(t_2/sqrt(dot(t_2,t_2)))

   t_1=curvidata%g_2/sqrt(curvidata%g_ij(2,2))
   t_2=curvidata%g_3/sqrt(curvidata%g_ij(3,3))
  faces_normal      (1:3,meshes_faces(1,i))=n
  faces_tangential  (1:3,meshes_faces(1,i))=t_1
  faces_bitangential(1:3,meshes_faces(1,i))=t_2
  faces_mid_coordinates(1:3,meshes_faces(1,i))=(/curvidata%x,curvidata%y,curvidata%z/)

!  print*, dot(faces_normal(:,meshes_faces(1,i)),faces_tangential(:,meshes_faces(1,i))), dot(faces_tangential(:,meshes_faces(1,i)),faces_bitangential(:,meshes_faces(1,i)))
  ! right x_i= 1.0
  x_i= 1.0
  y_i= 0.0
  z_i= 0.0
  call evaluate_curvilinear(MESHES_ARRAY(i),curvidata,x_i,y_i,z_i)

  n=curvidata%g1/sqrt(curvidata%gij(1,1))
  t_1=curvidata%g_2/sqrt(curvidata%g_ij(2,2))
  t_2=cross(n,t_1)
  t_2=abs(t_2/sqrt(dot(t_2,t_2)))
   t_1=curvidata%g_2/sqrt(curvidata%g_ij(2,2))
   t_2=curvidata%g_3/sqrt(curvidata%g_ij(3,3))
  faces_normal      (1:3,meshes_faces(2,i))=n
  faces_tangential  (1:3,meshes_faces(2,i))=t_1
  faces_bitangential(1:3,meshes_faces(2,i))=t_2
  faces_mid_coordinates(1:3,meshes_faces(2,i))=(/curvidata%x,curvidata%y,curvidata%z/)

  ! back y_i=-1.0
  x_i= 0.0
  y_i=-1.0
  z_i= 0.0
  call evaluate_curvilinear(MESHES_ARRAY(i),curvidata,x_i,y_i,z_i)

  n=curvidata%g2/sqrt(curvidata%gij(2,2))
  t_1=curvidata%g_1/sqrt(curvidata%g_ij(1,1))
  t_2=cross(n,t_1)
  t_2=abs(t_2/sqrt(dot(t_2,t_2)))
   t_1=curvidata%g_1/sqrt(curvidata%g_ij(1,1))
   t_2=curvidata%g_3/sqrt(curvidata%g_ij(3,3))
  faces_normal      (1:3,meshes_faces(3,i))=n
  faces_tangential  (1:3,meshes_faces(3,i))=t_1
  faces_bitangential(1:3,meshes_faces(3,i))=t_2
  faces_mid_coordinates(1:3,meshes_faces(3,i))=(/curvidata%x,curvidata%y,curvidata%z/)

  ! front y_i= 1.0
  x_i= 0.0
  y_i= 1.0
  z_i= 0.0
  call evaluate_curvilinear(MESHES_ARRAY(i),curvidata,x_i,y_i,z_i)

  n=curvidata%g2/sqrt(curvidata%gij(2,2))
  t_1=curvidata%g_1/sqrt(curvidata%g_ij(1,1))
  t_2=cross(n,t_1)
  t_2=abs(t_2/sqrt(dot(t_2,t_2)))
   t_1=curvidata%g_1/sqrt(curvidata%g_ij(1,1))
   t_2=curvidata%g_3/sqrt(curvidata%g_ij(3,3))
  faces_normal      (1:3,meshes_faces(4,i))=n
  faces_tangential  (1:3,meshes_faces(4,i))=t_1
  faces_bitangential(1:3,meshes_faces(4,i))=t_2
  faces_mid_coordinates(1:3,meshes_faces(4,i))=(/curvidata%x,curvidata%y,curvidata%z/)

  ! bottom z_i=-1.0
  x_i= 0.0
  y_i= 0.0
  z_i=-1.0
  call evaluate_curvilinear(MESHES_ARRAY(i),curvidata,x_i,y_i,z_i)

  n=curvidata%g3/sqrt(curvidata%gij(3,3))
  t_1=curvidata%g_1/sqrt(curvidata%g_ij(1,1))
  t_2=cross(n,t_1)
  t_2=abs(t_2/sqrt(dot(t_2,t_2)))
   t_1=curvidata%g_1/sqrt(curvidata%g_ij(1,1))
   t_2=curvidata%g_2/sqrt(curvidata%g_ij(2,2))
  faces_normal      (1:3,meshes_faces(5,i))=n
  faces_tangential  (1:3,meshes_faces(5,i))=t_1
  faces_bitangential(1:3,meshes_faces(5,i))=t_2
  faces_mid_coordinates(1:3,meshes_faces(5,i))=(/curvidata%x,curvidata%y,curvidata%z/)

  ! top z_i= 1.0
  x_i= 0.0
  y_i= 0.0
  z_i= 1.0
  call evaluate_curvilinear(MESHES_ARRAY(i),curvidata,x_i,y_i,z_i)

  n=curvidata%g3/sqrt(curvidata%gij(3,3))
  t_1=curvidata%g_1/sqrt(curvidata%g_ij(1,1))
  t_2=cross(n,t_1)
  t_2=abs(t_2/sqrt(dot(t_2,t_2)))
   t_1=curvidata%g_1/sqrt(curvidata%g_ij(1,1))
   t_2=curvidata%g_2/sqrt(curvidata%g_ij(2,2))
  faces_normal      (1:3,meshes_faces(6,i))=n
  faces_tangential  (1:3,meshes_faces(6,i))=t_1
  faces_bitangential(1:3,meshes_faces(6,i))=t_2
  faces_mid_coordinates(1:3,meshes_faces(6,i))=(/curvidata%x,curvidata%y,curvidata%z/)
end subroutine find_normal



  subroutine find_normal2(M)
    implicit none
! -----------------------------------
    integer, intent(in)                         :: M
    real, dimension (1:Gauss_int_p,1:Gauss_int_p) ::  dS
	real, dimension (1:3,1:Gauss_int_p,1:Gauss_int_p) :: n_i,t_i,s_i, mid_i
    integer   :: i,j,k,fn
    type(curvi3d)    :: curvidata
	real             :: x_i, y_i, z_i
! Calculate diffusion tensor, area of faces, and analytical over faces

! -----------------------------------


! left faces
	x_i=-1.0
	fn=1
	do j=1,Gauss_int_p
	  do k=1,Gauss_int_p
		call evaluate_curvilinear(MESHES_ARRAY(M),curvidata,x_i,Gauss_x_i(j),Gauss_x_i(k))
		dS(j,k)=curvidata%detJ*sqrt(curvidata%gij(1,1))
		n_i(:,j,k)=curvidata%g1/sqrt(curvidata%gij(1,1))*dS(j,k)
		t_i(:,j,k)=curvidata%g_2/sqrt(curvidata%g_ij(2,2))*dS(j,k)
		s_i(:,j,k)=cross(n_i(:,j,k),t_i(:,j,k))
		s_i(:,j,k)=s_i(:,j,k)/sqrt(dot(s_i(:,j,k),s_i(:,j,k)))*dS(j,k)
		mid_i(1,j,k)=curvidata%x*dS(j,k)
		mid_i(2,j,k)=curvidata%y*dS(j,k)
		mid_i(3,j,k)=curvidata%z*dS(j,k)
       enddo
    enddo
	
	meshes_normal_faces(1,fn,M)=Gauss_integration2d(n_i(1,:,:))/faces_area(meshes_faces(fn,M))
	meshes_normal_faces(2,fn,M)=Gauss_integration2d(n_i(2,:,:))/faces_area(meshes_faces(fn,M))
	meshes_normal_faces(3,fn,M)=Gauss_integration2d(n_i(3,:,:))/faces_area(meshes_faces(fn,M))
	meshes_tangential_faces(1,fn,M)=Gauss_integration2d(t_i(1,:,:))/faces_area(meshes_faces(fn,M))
	meshes_tangential_faces(2,fn,M)=Gauss_integration2d(t_i(2,:,:))/faces_area(meshes_faces(fn,M))
	meshes_tangential_faces(3,fn,M)=Gauss_integration2d(t_i(3,:,:))/faces_area(meshes_faces(fn,M))
	meshes_binormal_faces(1,fn,M)=Gauss_integration2d(s_i(1,:,:))/faces_area(meshes_faces(fn,M))
	meshes_binormal_faces(2,fn,M)=Gauss_integration2d(s_i(2,:,:))/faces_area(meshes_faces(fn,M))
	meshes_binormal_faces(3,fn,M)=Gauss_integration2d(s_i(3,:,:))/faces_area(meshes_faces(fn,M))
	meshes_mid_faces(1,fn,M)=Gauss_integration2d(mid_i(1,:,:))/faces_area(meshes_faces(fn,M))
	meshes_mid_faces(2,fn,M)=Gauss_integration2d(mid_i(2,:,:))/faces_area(meshes_faces(fn,M))
	meshes_mid_faces(3,fn,M)=Gauss_integration2d(mid_i(3,:,:))/faces_area(meshes_faces(fn,M))

	
	
! right faces
	x_i=1.0
	fn=2
	do j=1,Gauss_int_p
	  do k=1,Gauss_int_p
		call evaluate_curvilinear(MESHES_ARRAY(M),curvidata,x_i,Gauss_x_i(j),Gauss_x_i(k))
		dS(j,k)=curvidata%detJ*sqrt(curvidata%gij(1,1))
		n_i(:,j,k)=curvidata%g1/sqrt(curvidata%gij(1,1))*dS(j,k)
		t_i(:,j,k)=curvidata%g_2/sqrt(curvidata%g_ij(2,2))*dS(j,k)
		s_i(:,j,k)=cross(n_i(:,j,k),t_i(:,j,k))
		s_i(:,j,k)=s_i(:,j,k)/sqrt(dot(s_i(:,j,k),s_i(:,j,k)))*dS(j,k)
		mid_i(1,j,k)=curvidata%x*dS(j,k)
		mid_i(2,j,k)=curvidata%y*dS(j,k)
		mid_i(3,j,k)=curvidata%z*dS(j,k)
       enddo
    enddo
	
	meshes_normal_faces(1,fn,M)=Gauss_integration2d(n_i(1,:,:))/faces_area(meshes_faces(fn,M))
	meshes_normal_faces(2,fn,M)=Gauss_integration2d(n_i(2,:,:))/faces_area(meshes_faces(fn,M))
	meshes_normal_faces(3,fn,M)=Gauss_integration2d(n_i(3,:,:))/faces_area(meshes_faces(fn,M))
	meshes_tangential_faces(1,fn,M)=Gauss_integration2d(t_i(1,:,:))/faces_area(meshes_faces(fn,M))
	meshes_tangential_faces(2,fn,M)=Gauss_integration2d(t_i(2,:,:))/faces_area(meshes_faces(fn,M))
	meshes_tangential_faces(3,fn,M)=Gauss_integration2d(t_i(3,:,:))/faces_area(meshes_faces(fn,M))
	meshes_binormal_faces(1,fn,M)=Gauss_integration2d(s_i(1,:,:))/faces_area(meshes_faces(fn,M))
	meshes_binormal_faces(2,fn,M)=Gauss_integration2d(s_i(2,:,:))/faces_area(meshes_faces(fn,M))
	meshes_binormal_faces(3,fn,M)=Gauss_integration2d(s_i(3,:,:))/faces_area(meshes_faces(fn,M))
	meshes_mid_faces(1,fn,M)=Gauss_integration2d(mid_i(1,:,:))/faces_area(meshes_faces(fn,M))
	meshes_mid_faces(2,fn,M)=Gauss_integration2d(mid_i(2,:,:))/faces_area(meshes_faces(fn,M))
	meshes_mid_faces(3,fn,M)=Gauss_integration2d(mid_i(3,:,:))/faces_area(meshes_faces(fn,M))	

	
	y_i=-1.0
	fn=3
	do j=1,Gauss_int_p
	  do k=1,Gauss_int_p
		call evaluate_curvilinear(MESHES_ARRAY(M),curvidata,Gauss_x_i(j),y_i,Gauss_x_i(k))
		dS(j,k)=curvidata%detJ*sqrt(curvidata%gij(2,2))
		n_i(:,j,k)=curvidata%g2/sqrt(curvidata%gij(2,2))*dS(j,k)
		t_i(:,j,k)=curvidata%g_1/sqrt(curvidata%g_ij(1,1))*dS(j,k)
		s_i(:,j,k)=cross(n_i(:,j,k),t_i(:,j,k))
		s_i(:,j,k)=s_i(:,j,k)/sqrt(dot(s_i(:,j,k),s_i(:,j,k)))*dS(j,k)
		mid_i(1,j,k)=curvidata%x*dS(j,k)
		mid_i(2,j,k)=curvidata%y*dS(j,k)
		mid_i(3,j,k)=curvidata%z*dS(j,k)
       enddo
    enddo
	
	meshes_normal_faces(1,fn,M)=Gauss_integration2d(n_i(1,:,:))/faces_area(meshes_faces(fn,M))
	meshes_normal_faces(2,fn,M)=Gauss_integration2d(n_i(2,:,:))/faces_area(meshes_faces(fn,M))
	meshes_normal_faces(3,fn,M)=Gauss_integration2d(n_i(3,:,:))/faces_area(meshes_faces(fn,M))
	meshes_tangential_faces(1,fn,M)=Gauss_integration2d(t_i(1,:,:))/faces_area(meshes_faces(fn,M))
	meshes_tangential_faces(2,fn,M)=Gauss_integration2d(t_i(2,:,:))/faces_area(meshes_faces(fn,M))
	meshes_tangential_faces(3,fn,M)=Gauss_integration2d(t_i(3,:,:))/faces_area(meshes_faces(fn,M))
	meshes_binormal_faces(1,fn,M)=Gauss_integration2d(s_i(1,:,:))/faces_area(meshes_faces(fn,M))
	meshes_binormal_faces(2,fn,M)=Gauss_integration2d(s_i(2,:,:))/faces_area(meshes_faces(fn,M))
	meshes_binormal_faces(3,fn,M)=Gauss_integration2d(s_i(3,:,:))/faces_area(meshes_faces(fn,M))
	meshes_mid_faces(1,fn,M)=Gauss_integration2d(mid_i(1,:,:))/faces_area(meshes_faces(fn,M))
	meshes_mid_faces(2,fn,M)=Gauss_integration2d(mid_i(2,:,:))/faces_area(meshes_faces(fn,M))
	meshes_mid_faces(3,fn,M)=Gauss_integration2d(mid_i(3,:,:))/faces_area(meshes_faces(fn,M))


	y_i=1.0
	fn=4
	do j=1,Gauss_int_p
	  do k=1,Gauss_int_p
		call evaluate_curvilinear(MESHES_ARRAY(M),curvidata,Gauss_x_i(j),y_i,Gauss_x_i(k))
		dS(j,k)=curvidata%detJ*sqrt(curvidata%gij(2,2))
		n_i(:,j,k)=curvidata%g2/sqrt(curvidata%gij(2,2))*dS(j,k)
		t_i(:,j,k)=curvidata%g_1/sqrt(curvidata%g_ij(1,1))*dS(j,k)
		s_i(:,j,k)=cross(n_i(:,j,k),t_i(:,j,k))
		s_i(:,j,k)=s_i(:,j,k)/sqrt(dot(s_i(:,j,k),s_i(:,j,k)))*dS(j,k)
		mid_i(1,j,k)=curvidata%x*dS(j,k)
		mid_i(2,j,k)=curvidata%y*dS(j,k)
		mid_i(3,j,k)=curvidata%z*dS(j,k)
       enddo
    enddo
	
	meshes_normal_faces(1,fn,M)=Gauss_integration2d(n_i(1,:,:))/faces_area(meshes_faces(fn,M))
	meshes_normal_faces(2,fn,M)=Gauss_integration2d(n_i(2,:,:))/faces_area(meshes_faces(fn,M))
	meshes_normal_faces(3,fn,M)=Gauss_integration2d(n_i(3,:,:))/faces_area(meshes_faces(fn,M))
	meshes_tangential_faces(1,fn,M)=Gauss_integration2d(t_i(1,:,:))/faces_area(meshes_faces(fn,M))
	meshes_tangential_faces(2,fn,M)=Gauss_integration2d(t_i(2,:,:))/faces_area(meshes_faces(fn,M))
	meshes_tangential_faces(3,fn,M)=Gauss_integration2d(t_i(3,:,:))/faces_area(meshes_faces(fn,M))
	meshes_binormal_faces(1,fn,M)=Gauss_integration2d(s_i(1,:,:))/faces_area(meshes_faces(fn,M))
	meshes_binormal_faces(2,fn,M)=Gauss_integration2d(s_i(2,:,:))/faces_area(meshes_faces(fn,M))
	meshes_binormal_faces(3,fn,M)=Gauss_integration2d(s_i(3,:,:))/faces_area(meshes_faces(fn,M))
	meshes_mid_faces(1,fn,M)=Gauss_integration2d(mid_i(1,:,:))/faces_area(meshes_faces(fn,M))
	meshes_mid_faces(2,fn,M)=Gauss_integration2d(mid_i(2,:,:))/faces_area(meshes_faces(fn,M))
	meshes_mid_faces(3,fn,M)=Gauss_integration2d(mid_i(3,:,:))/faces_area(meshes_faces(fn,M))


	z_i=-1.0
	fn=5
	do j=1,Gauss_int_p
	  do k=1,Gauss_int_p
		call evaluate_curvilinear(MESHES_ARRAY(M),curvidata,Gauss_x_i(j),Gauss_x_i(k),z_i)
		dS(j,k)=curvidata%detJ*sqrt(curvidata%gij(3,3))
		n_i(:,j,k)=curvidata%g3/sqrt(curvidata%gij(3,3))*dS(j,k)
		t_i(:,j,k)=curvidata%g_1/sqrt(curvidata%g_ij(1,1))*dS(j,k)
		s_i(:,j,k)=cross(n_i(:,j,k),t_i(:,j,k))
		s_i(:,j,k)=s_i(:,j,k)/sqrt(dot(s_i(:,j,k),s_i(:,j,k)))*dS(j,k)
		mid_i(1,j,k)=curvidata%x*dS(j,k)
		mid_i(2,j,k)=curvidata%y*dS(j,k)
		mid_i(3,j,k)=curvidata%z*dS(j,k)
       enddo
    enddo
	
	meshes_normal_faces(1,fn,M)=Gauss_integration2d(n_i(1,:,:))/faces_area(meshes_faces(fn,M))
	meshes_normal_faces(2,fn,M)=Gauss_integration2d(n_i(2,:,:))/faces_area(meshes_faces(fn,M))
	meshes_normal_faces(3,fn,M)=Gauss_integration2d(n_i(3,:,:))/faces_area(meshes_faces(fn,M))
	meshes_tangential_faces(1,fn,M)=Gauss_integration2d(t_i(1,:,:))/faces_area(meshes_faces(fn,M))
	meshes_tangential_faces(2,fn,M)=Gauss_integration2d(t_i(2,:,:))/faces_area(meshes_faces(fn,M))
	meshes_tangential_faces(3,fn,M)=Gauss_integration2d(t_i(3,:,:))/faces_area(meshes_faces(fn,M))
	meshes_binormal_faces(1,fn,M)=Gauss_integration2d(s_i(1,:,:))/faces_area(meshes_faces(fn,M))
	meshes_binormal_faces(2,fn,M)=Gauss_integration2d(s_i(2,:,:))/faces_area(meshes_faces(fn,M))
	meshes_binormal_faces(3,fn,M)=Gauss_integration2d(s_i(3,:,:))/faces_area(meshes_faces(fn,M))
	meshes_mid_faces(1,fn,M)=Gauss_integration2d(mid_i(1,:,:))/faces_area(meshes_faces(fn,M))
	meshes_mid_faces(2,fn,M)=Gauss_integration2d(mid_i(2,:,:))/faces_area(meshes_faces(fn,M))
	meshes_mid_faces(3,fn,M)=Gauss_integration2d(mid_i(3,:,:))/faces_area(meshes_faces(fn,M))

! left faces
	z_i=1.0
	fn=6
	do j=1,Gauss_int_p
	  do k=1,Gauss_int_p
		call evaluate_curvilinear(MESHES_ARRAY(M),curvidata,Gauss_x_i(j),Gauss_x_i(k),z_i)
		dS(j,k)=curvidata%detJ*sqrt(curvidata%gij(3,3))
		n_i(:,j,k)=curvidata%g3/sqrt(curvidata%gij(3,3))*dS(j,k)
		t_i(:,j,k)=curvidata%g_1/sqrt(curvidata%g_ij(1,1))*dS(j,k)
		s_i(:,j,k)=cross(n_i(:,j,k),t_i(:,j,k))
		s_i(:,j,k)=s_i(:,j,k)/sqrt(dot(s_i(:,j,k),s_i(:,j,k)))*dS(j,k)
		mid_i(1,j,k)=curvidata%x*dS(j,k)
		mid_i(2,j,k)=curvidata%y*dS(j,k)
		mid_i(3,j,k)=curvidata%z*dS(j,k)
       enddo
    enddo
	
	meshes_normal_faces(1,fn,M)=Gauss_integration2d(n_i(1,:,:))/faces_area(meshes_faces(fn,M))
	meshes_normal_faces(2,fn,M)=Gauss_integration2d(n_i(2,:,:))/faces_area(meshes_faces(fn,M))
	meshes_normal_faces(3,fn,M)=Gauss_integration2d(n_i(3,:,:))/faces_area(meshes_faces(fn,M))
	meshes_tangential_faces(1,fn,M)=Gauss_integration2d(t_i(1,:,:))/faces_area(meshes_faces(fn,M))
	meshes_tangential_faces(2,fn,M)=Gauss_integration2d(t_i(2,:,:))/faces_area(meshes_faces(fn,M))
	meshes_tangential_faces(3,fn,M)=Gauss_integration2d(t_i(3,:,:))/faces_area(meshes_faces(fn,M))
	meshes_binormal_faces(1,fn,M)=Gauss_integration2d(s_i(1,:,:))/faces_area(meshes_faces(fn,M))
	meshes_binormal_faces(2,fn,M)=Gauss_integration2d(s_i(2,:,:))/faces_area(meshes_faces(fn,M))
	meshes_binormal_faces(3,fn,M)=Gauss_integration2d(s_i(3,:,:))/faces_area(meshes_faces(fn,M))
	meshes_mid_faces(1,fn,M)=Gauss_integration2d(mid_i(1,:,:))/faces_area(meshes_faces(fn,M))
	meshes_mid_faces(2,fn,M)=Gauss_integration2d(mid_i(2,:,:))/faces_area(meshes_faces(fn,M))
	meshes_mid_faces(3,fn,M)=Gauss_integration2d(mid_i(3,:,:))/faces_area(meshes_faces(fn,M))


  end subroutine find_normal2



  subroutine calculate_defusion_tensor_over_faces(i)
! -----------------------------------
    integer, intent(in)             :: i
    real, dimension(1:Gauss_int_p)  :: ones_plus,ones_minus
    integer                         :: j

! -----------------------------------
    do j=1,Gauss_int_p
       ones_plus(j)=1.
       ones_minus(j)=-1.
    enddo
    call calculate_diffusion_tensor_faces_i(i,ones_minus,Gauss_x_i,Gauss_x_i,1)
    call calculate_diffusion_tensor_faces_i(i,ones_plus ,Gauss_x_i,Gauss_x_i,2)
    call calculate_diffusion_tensor_faces_i(i,Gauss_x_i,ones_minus,Gauss_x_i,3)
    call calculate_diffusion_tensor_faces_i(i,Gauss_x_i,ones_plus ,Gauss_x_i,4)
    call calculate_diffusion_tensor_faces_i(i,Gauss_x_i,Gauss_x_i,ones_minus,5)
    call calculate_diffusion_tensor_faces_i(i,Gauss_x_i,Gauss_x_i,ones_plus ,6)
 
 end subroutine calculate_defusion_tensor_over_faces

  subroutine calculate_diffusion_tensor_faces_i(M,x_array,y_array,z_array,n)
    implicit none
! -----------------------------------
    integer, intent(in)                         :: M
    integer, intent(in)                         :: n
    real, intent(in), dimension(1:Gauss_int_p)  :: x_array,y_array,z_array
    real, dimension (1:Gauss_int_p,1:Gauss_int_p,1:Gauss_int_p) :: K11,K12,K13,K21,K22,K23,K31,K32,K33,dS,u_i, v_i, w_i, p_i,bx_i,by_i,bz_i, uic_i, vic_i, wic_i, pic_i, TT_i, Tic_i, q_i
	real, dimension (1:3,1:Gauss_int_p,1:Gauss_int_p,1:Gauss_int_p) :: n_i,t_i,s_i
    integer   :: i,j,k
    real :: Diff
    type(curvi3d)    :: C
! Calculate diffusion tensor, area of faces, and analytical over faces

! -----------------------------------
    do i=1,Gauss_int_p
       do j=1,Gauss_int_p
          do k=1,Gauss_int_p
                call evaluate_curvilinear(MESHES_ARRAY(M),C,x_array(i),y_array(j),z_array(k))
                if(n==1 .or. n==2) then
                   dS(i,j,k)=C%detJ*sqrt(C%gij(1,1))
                elseif(n==3 .or. n==4) then
                   dS(i,j,k)=C%detJ*sqrt(C%gij(2,2))
                elseif(n==5 .or. n==6) then
                   dS(i,j,k)=C%detJ*sqrt(C%gij(3,3))
                endif
                K11(i,j,k)=C%gij(1,1)/sqrt(C%gij(1,1))*dS(i,j,k)
                K12(i,j,k)=C%gij(1,2)/sqrt(C%gij(1,1))*dS(i,j,k)
                K13(i,j,k)=C%gij(1,3)/sqrt(C%gij(1,1))*dS(i,j,k)
                K21(i,j,k)=C%gij(2,1)/sqrt(C%gij(2,2))*dS(i,j,k)
                K22(i,j,k)=C%gij(2,2)/sqrt(C%gij(2,2))*dS(i,j,k)
                K23(i,j,k)=C%gij(2,3)/sqrt(C%gij(2,2))*dS(i,j,k)
                K31(i,j,k)=C%gij(3,1)/sqrt(C%gij(3,3))*dS(i,j,k)
                K32(i,j,k)=C%gij(3,2)/sqrt(C%gij(3,3))*dS(i,j,k)
                K33(i,j,k)=C%gij(3,3)/sqrt(C%gij(3,3))*dS(i,j,k)
                if(.not. is_analytical_TD .and. is_analytical) then
    				u_i(i,j,k)=functions_id (C%x,C%y,C%z,0.0,velx1)*dS(i,j,k)
					v_i(i,j,k)=functions_id (C%x,C%y,C%z,0.0,velx2)*dS(i,j,k)
					w_i(i,j,k)=functions_id (C%x,C%y,C%z,0.0,velx3)*dS(i,j,k)
					p_i(i,j,k)=functions_id (C%x,C%y,C%z,0.0,PPE)*dS(i,j,k)
					TT_i(i,j,k)=functions_id (C%x,C%y,C%z,0.0,CDE)*dS(i,j,k)
				endif
				uic_i(i,j,k)=functions_id (C%x,C%y,C%z,0.0,IC_velx1)*dS(i,j,k)
				vic_i(i,j,k)=functions_id (C%x,C%y,C%z,0.0,IC_velx2)*dS(i,j,k)
				wic_i(i,j,k)=functions_id (C%x,C%y,C%z,0.0,IC_velx3)*dS(i,j,k)
				pic_i(i,j,k)=functions_id (C%x,C%y,C%z,0.0,IC_PPE)*dS(i,j,k)
				Tic_i(i,j,k)=functions_id (C%x,C%y,C%z,0.0,IC_CDE)*dS(i,j,k)
                if(is_bx .and. .not. is_bx_TD) bx_i(i,j,k)=functions_id (C%x,C%y,C%z,0.0,bx_fun)*dS(i,j,k)
                if(is_by .and. .not. is_by_TD) by_i(i,j,k)=functions_id (C%x,C%y,C%z,0.0,by_fun)*dS(i,j,k)
                if(is_bz .and. .not. is_bz_TD) bz_i(i,j,k)=functions_id (C%x,C%y,C%z,0.0,bz_fun)*dS(i,j,k)
                if(is_q .and. is_energy) q_i(i,j,k)=functions_id (C%x,C%y,C%z,0.0,q_fun)*dS(i,j,k)
          enddo
       enddo
    enddo
	faces_area(meshes_faces(n,M))=Gauss_integration2d_face(dS,n)
	if(.not. is_analytical_TD .and. is_analytical) then
		faces_uT(meshes_faces(n,M))= Gauss_integration2d_face(u_i,n)/faces_area(meshes_faces(n,M))
		faces_vT(meshes_faces(n,M))= Gauss_integration2d_face(v_i,n)/faces_area(meshes_faces(n,M))
		faces_wT(meshes_faces(n,M))= Gauss_integration2d_face(w_i,n)/faces_area(meshes_faces(n,M))
		faces_pT(meshes_faces(n,M))= Gauss_integration2d_face(p_i,n)/faces_area(meshes_faces(n,M))
		faces_TT(meshes_faces(n,M))= Gauss_integration2d_face(TT_i,n)/faces_area(meshes_faces(n,M))
	endif
	if(is_bx .and. .not. is_bx_TD) faces_bx(meshes_faces(n,M))= Gauss_integration2d_face(bx_i,n)/faces_area(meshes_faces(n,M))
	if(is_by .and. .not. is_by_TD) faces_by(meshes_faces(n,M))= Gauss_integration2d_face(by_i,n)/faces_area(meshes_faces(n,M))
	if(is_bz .and. .not. is_bz_TD) faces_bz(meshes_faces(n,M))= Gauss_integration2d_face(bz_i,n)/faces_area(meshes_faces(n,M))
	if(is_q .and. is_energy) faces_q(meshes_faces(n,M))= Gauss_integration2d_face(q_i,n)/faces_area(meshes_faces(n,M))

	faces_u(meshes_faces(n,M))= Gauss_integration2d_face(uic_i,n)/faces_area(meshes_faces(n,M))
	faces_v(meshes_faces(n,M))= Gauss_integration2d_face(vic_i,n)/faces_area(meshes_faces(n,M))
	faces_w(meshes_faces(n,M))= Gauss_integration2d_face(wic_i,n)/faces_area(meshes_faces(n,M))
	faces_p(meshes_faces(n,M))= Gauss_integration2d_face(pic_i,n)/faces_area(meshes_faces(n,M))
	faces_T(meshes_faces(n,M))= Gauss_integration2d_face(Tic_i,n)/faces_area(meshes_faces(n,M))
	
	
    meshes_gij_faces(1,1,n,M)= Gauss_integration2d_face(K11,n)/faces_area(meshes_faces(n,M))
    meshes_gij_faces(1,2,n,M)= Gauss_integration2d_face(K12,n)/faces_area(meshes_faces(n,M))
    meshes_gij_faces(1,3,n,M)= Gauss_integration2d_face(K13,n)/faces_area(meshes_faces(n,M))
    meshes_gij_faces(2,1,n,M)= Gauss_integration2d_face(K21,n)/faces_area(meshes_faces(n,M))
    meshes_gij_faces(2,2,n,M)= Gauss_integration2d_face(K22,n)/faces_area(meshes_faces(n,M))
    meshes_gij_faces(2,3,n,M)= Gauss_integration2d_face(K23,n)/faces_area(meshes_faces(n,M))
    meshes_gij_faces(3,1,n,M)= Gauss_integration2d_face(K31,n)/faces_area(meshes_faces(n,M))
    meshes_gij_faces(3,2,n,M)= Gauss_integration2d_face(K32,n)/faces_area(meshes_faces(n,M))
    meshes_gij_faces(3,3,n,M)= Gauss_integration2d_face(K33,n)/faces_area(meshes_faces(n,M))
  end subroutine calculate_diffusion_tensor_faces_i


  subroutine convection_diffusion_constants_mesh_i (M)
! -----------------------------------
    integer, intent (in)      :: M
    real, dimension (1:Gauss_int_p,1:Gauss_int_p,1:Gauss_int_p) :: K11_i,K22_i,K33_i,K12_i,K21_i,K13_i,K31_i,K23_i,K32_i,G_11,G_22,G_33,bx_i,by_i,bz_i, &
																   J11_i,J22_i,J33_i,J12_i,J21_i,J13_i,J31_i,J23_i,J32_i, u_ic, v_ic, w_ic, p_ic,u_i,v_i,w_i, TT_i,T_ic, q_i
    real                                                        :: G11,G22,G33,bx,by,bz
    integer                         :: i,j,k
    type(curvi3d)    :: C

! -----------------------------------
    do i=1,Gauss_int_p
       do j=1,Gauss_int_p
          do k=1,Gauss_int_p
                call evaluate_curvilinear(MESHES_ARRAY(M),C,Gauss_x_i(i),Gauss_x_i(j),Gauss_x_i(k))
!                if(steady_state) S=-Equation_source(C%x,C%y,C%z,0.0)
                G11=(C%gamma_tensor(1,1,2)*C%gij(1,2)+C%gamma_tensor(1,2,1)*C%gij(2,1)&
                    +C%gamma_tensor(1,1,3)*C%gij(1,3)+C%gamma_tensor(1,3,1)*C%gij(3,1)&
                    +C%gamma_tensor(1,2,3)*C%gij(2,3)+C%gamma_tensor(1,3,2)*C%gij(3,2))

                G22=(C%gamma_tensor(2,1,2)*C%gij(1,2)+C%gamma_tensor(2,2,1)*C%gij(2,1)&
                    +C%gamma_tensor(2,1,3)*C%gij(1,3)+C%gamma_tensor(2,3,1)*C%gij(3,1)&
                    +C%gamma_tensor(2,2,3)*C%gij(2,3)+C%gamma_tensor(2,3,2)*C%gij(3,2))

                G33=(C%gamma_tensor(3,1,2)*C%gij(1,2)+C%gamma_tensor(3,2,1)*C%gij(2,1)&
                    +C%gamma_tensor(3,1,3)*C%gij(1,3)+C%gamma_tensor(3,3,1)*C%gij(3,1)&
                    +C%gamma_tensor(3,2,3)*C%gij(2,3)+C%gamma_tensor(3,3,2)*C%gij(3,2))

! ============================= Filling arrays ==================================
!                S_i(i,j,k)  =S      *C%detJ
                G_11(i,j,k) =G11    *C%detJ
                G_22(i,j,k) =G22    *C%detJ
                G_33(i,j,k) =G33    *C%detJ
                K11_i(i,j,k)=C%gij(1,1)*C%detJ
                K12_i(i,j,k)=C%gij(1,2)*C%detJ
                K13_i(i,j,k)=C%gij(1,3)*C%detJ
                K21_i(i,j,k)=C%gij(2,1)*C%detJ
                K22_i(i,j,k)=C%gij(2,2)*C%detJ
                K23_i(i,j,k)=C%gij(2,3)*C%detJ
                K31_i(i,j,k)=C%gij(3,1)*C%detJ
                K32_i(i,j,k)=C%gij(3,2)*C%detJ
                K33_i(i,j,k)=C%gij(3,3)*C%detJ

                J11_i(i,j,k)=C%jacinv(1,1)*C%detJ
                J12_i(i,j,k)=C%jacinv(1,2)*C%detJ
                J13_i(i,j,k)=C%jacinv(1,3)*C%detJ
                J21_i(i,j,k)=C%jacinv(2,1)*C%detJ
                J22_i(i,j,k)=C%jacinv(2,2)*C%detJ
                J23_i(i,j,k)=C%jacinv(2,3)*C%detJ
                J31_i(i,j,k)=C%jacinv(3,1)*C%detJ
                J32_i(i,j,k)=C%jacinv(3,2)*C%detJ
                J33_i(i,j,k)=C%jacinv(3,3)*C%detJ
				
				! Initial conditions
				u_ic(i,j,k)=functions_id(C%x,C%y,C%z,0.0,IC_velx1)*C%detJ
				v_ic(i,j,k)=functions_id(C%x,C%y,C%z,0.0,IC_velx2)*C%detJ
				w_ic(i,j,k)=functions_id(C%x,C%y,C%z,0.0,IC_velx3)*C%detJ
				p_ic(i,j,k)=functions_id(C%x,C%y,C%z,0.0,IC_PPE)  *C%detJ
				T_ic(i,j,k)=functions_id(C%x,C%y,C%z,0.0,IC_PPE)  *C%detJ
				
				
				if(is_bx .and. .not. is_bx_TD) bx_i(i,j,k)=functions_id (C%x,C%y,C%z,0.0,bx_fun)*C%detJ
				if(is_by .and. .not. is_bx_TD) by_i(i,j,k)=functions_id (C%x,C%y,C%z,0.0,by_fun)*C%detJ
				if(is_bz .and. .not. is_bx_TD) bz_i(i,j,k)=functions_id (C%x,C%y,C%z,0.0,bz_fun)*C%detJ
				if(is_q .and. is_energy) q_i(i,j,k)=functions_id (C%x,C%y,C%z,0.0,q_fun)*C%detJ
				if(is_analytical) then
					u_i(i,j,k)=functions_id (C%x,C%y,C%z,0.0,analytical_velx1)*C%detJ
					v_i(i,j,k)=functions_id (C%x,C%y,C%z,0.0,analytical_velx2)*C%detJ
					w_i(i,j,k)=functions_id (C%x,C%y,C%z,0.0,analytical_velx3)*C%detJ
					TT_i(i,j,k)=functions_id (C%x,C%y,C%z,0.0,analytical_CDE)*C%detJ
				endif 

          enddo
       enddo
    enddo
! ================================================================
! initial conditions
! u
    meshes_uxyz_old(M)=Gauss_integration(u_ic)/meshes_volume(M)
	meshes_uxyz(M)=meshes_uxyz_old(M)
	meshes_uxyzt(M)=meshes_uxyz_old(M)
! v	
    meshes_vxyz_old(M)=Gauss_integration(v_ic)/meshes_volume(M)
	meshes_vxyz(M)=meshes_vxyz_old(M)
	meshes_vxyzt(M)=meshes_vxyz_old(M)
! w	
	meshes_wxyz_old(M)=Gauss_integration(w_ic)/meshes_volume(M)
	meshes_wxyz(M)=meshes_wxyz_old(M)
	meshes_wxyzt(M)=meshes_wxyz_old(M)
	
	
	meshes_Txyz_old(M)=Gauss_integration(T_ic)/meshes_volume(M)
	meshes_Txyz(M)=meshes_Txyz_old(M)
	meshes_Txyzt(M)=meshes_Txyz_old(M)
	
! p -> does not time-dependent, but in case I need to initialize it with a function.	
	meshes_p(M)=Gauss_integration(p_ic)/meshes_volume(M)
	if(is_analytical) then
		meshes_uT(M)=Gauss_integration(u_i)/meshes_volume(M)
		meshes_vT(M)=Gauss_integration(v_i)/meshes_volume(M)
		meshes_wT(M)=Gauss_integration(w_i)/meshes_volume(M)
		meshes_TT(M)=Gauss_integration(TT_i)/meshes_volume(M)
	endif 
! ================================================================
    if(is_bx .and. .not. is_bx_TD) meshes_bx(M)=Gauss_integration(bx_i)/meshes_volume(M)
    if(is_by .and. .not. is_by_TD) meshes_by(M)=Gauss_integration(by_i)/meshes_volume(M)
    if(is_bz .and. .not. is_bz_TD) meshes_bz(M)=Gauss_integration(bz_i)/meshes_volume(M)
    if(is_q .and. is_energy) meshes_q(M)=Gauss_integration(q_i)/meshes_volume(M)
! ================================================================
    meshes_G1(M)   =Gauss_integration(G_11)/meshes_volume(M)
    meshes_G2(M)   =Gauss_integration(G_22)/meshes_volume(M)
    meshes_G3(M)   =Gauss_integration(G_33)/meshes_volume(M)
! ================================================================
    meshes_gij(1,1,M)=Gauss_integration(K11_i)/meshes_volume(M)
    meshes_gij(1,2,M)=Gauss_integration(K12_i)/meshes_volume(M)
    meshes_gij(1,3,M)=Gauss_integration(K13_i)/meshes_volume(M)
    meshes_gij(2,1,M)=Gauss_integration(K21_i)/meshes_volume(M)
    meshes_gij(2,2,M)=Gauss_integration(K22_i)/meshes_volume(M)
    meshes_gij(2,3,M)=Gauss_integration(K23_i)/meshes_volume(M)
    meshes_gij(3,1,M)=Gauss_integration(K31_i)/meshes_volume(M)
    meshes_gij(3,2,M)=Gauss_integration(K32_i)/meshes_volume(M)
    meshes_gij(3,3,M)=Gauss_integration(K33_i)/meshes_volume(M)
! ===============================================================
	meshes_xi_x(1,1,M)=Gauss_integration(J11_i)/meshes_volume(M)
    meshes_xi_x(2,1,M)=Gauss_integration(J12_i)/meshes_volume(M)
    meshes_xi_x(3,1,M)=Gauss_integration(J13_i)/meshes_volume(M)
    meshes_xi_x(1,2,M)=Gauss_integration(J21_i)/meshes_volume(M)
    meshes_xi_x(2,2,M)=Gauss_integration(J22_i)/meshes_volume(M)
    meshes_xi_x(3,2,M)=Gauss_integration(J23_i)/meshes_volume(M)
    meshes_xi_x(1,3,M)=Gauss_integration(J31_i)/meshes_volume(M)
    meshes_xi_x(2,3,M)=Gauss_integration(J32_i)/meshes_volume(M)
    meshes_xi_x(3,3,M)=Gauss_integration(J33_i)/meshes_volume(M)

  end subroutine convection_diffusion_constants_mesh_i



  subroutine time_dependent_face_averages(M)
    implicit none
! -----------------------------------
    integer, intent (in)                  :: M
!~    integer, intent(in)                         :: n
    real, dimension (1:Gauss_int_p,1:Gauss_int_p,1:Gauss_int_p) :: dS,bx_i,by_i,bz_i,u_i,v_i,w_i,p_i
    integer   :: i,j,k,l,fn
    real      :: t,x_i,y_i,z_i
    type(curvi3d)    :: C
! Calculate diffusion tensor, area of faces, and analytical over faces

! -----------------------------------
! left 
    fn=1
    do j=1,Gauss_int_p
      do k=1,Gauss_int_p
        do l=1,Gauss_int_p
           x_i=-1
		   call evaluate_curvilinear(MESHES_ARRAY(M),C,x_i,Gauss_x_i(j),Gauss_x_i(k))
		   call map_ti(Gauss_x_i(l),t)
		   dS(j,k,l)=C%detJ*sqrt(C%gij(1,1))
			if(is_analytical .and. is_analytical_TD) then
				u_i(j,k,l)=functions_id (C%x,C%y,C%z,t,velx1)*dS(j,k,l)
				v_i(j,k,l)=functions_id (C%x,C%y,C%z,t,velx2)*dS(j,k,l)
				w_i(j,k,l)=functions_id (C%x,C%y,C%z,t,velx3)*dS(j,k,l)
				p_i(j,k,l)=functions_id (C%x,C%y,C%z,t,PPE)*dS(j,k,l)
			endif 
			if(is_bx .and. is_bx_TD) bx_i(j,k,l)=functions_id (C%x,C%y,C%z,t,bx_fun)*dS(j,k,l)
			if(is_by .and. is_by_TD) by_i(j,k,l)=functions_id (C%x,C%y,C%z,t,by_fun)*dS(j,k,l)
			if(is_bz .and. is_bz_TD) bz_i(j,k,l)=functions_id (C%x,C%y,C%z,t,bz_fun)*dS(j,k,l)
        enddo
      enddo
    enddo
    if(is_analytical .and. is_analytical_TD) then
	  faces_uT(meshes_faces(fn,M))= Gauss_integration(u_i)/(Gauss_integration(dS))
	  faces_vT(meshes_faces(fn,M))= Gauss_integration(v_i)/(Gauss_integration(dS))
	  faces_wT(meshes_faces(fn,M))= Gauss_integration(w_i)/(Gauss_integration(dS))
	  faces_pT(meshes_faces(fn,M))= Gauss_integration(p_i)/(Gauss_integration(dS))
    endif 
	if(is_bx .and. is_bx_TD) faces_bx(meshes_faces(fn,M))= Gauss_integration(bx_i)/(Gauss_integration(dS))
	if(is_by .and. is_by_TD) faces_by(meshes_faces(fn,M))= Gauss_integration(by_i)/(Gauss_integration(dS))
	if(is_bz .and. is_bz_TD) faces_bz(meshes_faces(fn,M))= Gauss_integration(bz_i)/(Gauss_integration(dS))

! -----------------------------------
! right
    fn=2
    do j=1,Gauss_int_p
      do k=1,Gauss_int_p
        do l=1,Gauss_int_p
           x_i=1
           call evaluate_curvilinear(MESHES_ARRAY(M),C,x_i,Gauss_x_i(j),Gauss_x_i(k))
           call map_ti(Gauss_x_i(l),t)
           dS(j,k,l)=C%detJ*sqrt(C%gij(1,1))
			if(is_analytical .and. is_analytical_TD) then
				u_i(j,k,l)=functions_id (C%x,C%y,C%z,t,velx1)*dS(j,k,l)
				v_i(j,k,l)=functions_id (C%x,C%y,C%z,t,velx2)*dS(j,k,l)
				w_i(j,k,l)=functions_id (C%x,C%y,C%z,t,velx3)*dS(j,k,l)
				p_i(j,k,l)=functions_id (C%x,C%y,C%z,t,PPE)*dS(j,k,l)
			endif 
			if(is_bx .and. is_bx_TD) bx_i(j,k,l)=functions_id (C%x,C%y,C%z,t,bx_fun)*dS(j,k,l)
			if(is_by .and. is_by_TD) by_i(j,k,l)=functions_id (C%x,C%y,C%z,t,by_fun)*dS(j,k,l)
			if(is_bz .and. is_bz_TD) bz_i(j,k,l)=functions_id (C%x,C%y,C%z,t,bz_fun)*dS(j,k,l)
        enddo
      enddo
    enddo
    if(is_analytical .and. is_analytical_TD) then
	  faces_uT(meshes_faces(fn,M))= Gauss_integration(u_i)/(Gauss_integration(dS))
	  faces_vT(meshes_faces(fn,M))= Gauss_integration(v_i)/(Gauss_integration(dS))
	  faces_wT(meshes_faces(fn,M))= Gauss_integration(w_i)/(Gauss_integration(dS))
	  faces_pT(meshes_faces(fn,M))= Gauss_integration(p_i)/(Gauss_integration(dS))
    endif 
	if(is_bx .and. is_bx_TD) faces_bx(meshes_faces(fn,M))= Gauss_integration(bx_i)/(Gauss_integration(dS))
	if(is_by .and. is_by_TD) faces_by(meshes_faces(fn,M))= Gauss_integration(by_i)/(Gauss_integration(dS))
	if(is_bz .and. is_bz_TD) faces_bz(meshes_faces(fn,M))= Gauss_integration(bz_i)/(Gauss_integration(dS))

! -----------------------------------
! back
    fn=3
    do i=1,Gauss_int_p
      do k=1,Gauss_int_p
        do l=1,Gauss_int_p
           y_i=-1
           call evaluate_curvilinear(MESHES_ARRAY(M),C,Gauss_x_i(i),y_i,Gauss_x_i(k))
           call map_ti(Gauss_x_i(l),t)
           dS(i,k,l)=C%detJ*sqrt(C%gij(2,2))
			if(is_analytical .and. is_analytical_TD) then
				u_i(i,k,l)=functions_id (C%x,C%y,C%z,t,velx1)*dS(i,k,l)
				v_i(i,k,l)=functions_id (C%x,C%y,C%z,t,velx2)*dS(i,k,l)
				w_i(i,k,l)=functions_id (C%x,C%y,C%z,t,velx3)*dS(i,k,l)
				p_i(i,k,l)=functions_id (C%x,C%y,C%z,t,PPE)*dS(i,k,l)
			endif 
			if(is_bx .and. is_bx_TD) bx_i(i,k,l)=functions_id (C%x,C%y,C%z,t,bx_fun)*dS(i,k,l)
			if(is_by .and. is_by_TD) by_i(i,k,l)=functions_id (C%x,C%y,C%z,t,by_fun)*dS(i,k,l)
			if(is_bz .and. is_bz_TD) bz_i(i,k,l)=functions_id (C%x,C%y,C%z,t,bz_fun)*dS(i,k,l)
        enddo
      enddo
    enddo
    if(is_analytical .and. is_analytical_TD) then
	  faces_uT(meshes_faces(fn,M))= Gauss_integration(u_i)/(Gauss_integration(dS))
	  faces_vT(meshes_faces(fn,M))= Gauss_integration(v_i)/(Gauss_integration(dS))
	  faces_wT(meshes_faces(fn,M))= Gauss_integration(w_i)/(Gauss_integration(dS))
	  faces_pT(meshes_faces(fn,M))= Gauss_integration(p_i)/(Gauss_integration(dS))
    endif 
	if(is_bx .and. is_bx_TD) faces_bx(meshes_faces(fn,M))= Gauss_integration(bx_i)/(Gauss_integration(dS))
	if(is_by .and. is_by_TD) faces_by(meshes_faces(fn,M))= Gauss_integration(by_i)/(Gauss_integration(dS))
	if(is_bz .and. is_bz_TD) faces_bz(meshes_faces(fn,M))= Gauss_integration(bz_i)/(Gauss_integration(dS))

! -----------------------------------
! front
    fn=4
    do i=1,Gauss_int_p
      do k=1,Gauss_int_p
        do l=1,Gauss_int_p
           y_i=1
           call evaluate_curvilinear(MESHES_ARRAY(M),C,Gauss_x_i(i),y_i,Gauss_x_i(k))
           call map_ti(Gauss_x_i(l),t)
           dS(i,k,l)=C%detJ*sqrt(C%gij(2,2))
			if(is_analytical .and. is_analytical_TD) then
				u_i(i,k,l)=functions_id (C%x,C%y,C%z,t,velx1)*dS(i,k,l)
				v_i(i,k,l)=functions_id (C%x,C%y,C%z,t,velx2)*dS(i,k,l)
				w_i(i,k,l)=functions_id (C%x,C%y,C%z,t,velx3)*dS(i,k,l)
				p_i(i,k,l)=functions_id (C%x,C%y,C%z,t,PPE)*dS(i,k,l)
			endif 
			if(is_bx .and. is_bx_TD) bx_i(i,k,l)=functions_id (C%x,C%y,C%z,t,bx_fun)*dS(i,k,l)
			if(is_by .and. is_by_TD) by_i(i,k,l)=functions_id (C%x,C%y,C%z,t,by_fun)*dS(i,k,l)
			if(is_bz .and. is_bz_TD) bz_i(i,k,l)=functions_id (C%x,C%y,C%z,t,bz_fun)*dS(i,k,l)
        enddo
      enddo
    enddo
    if(is_analytical .and. is_analytical_TD) then
	  faces_uT(meshes_faces(fn,M))= Gauss_integration(u_i)/(Gauss_integration(dS))
	  faces_vT(meshes_faces(fn,M))= Gauss_integration(v_i)/(Gauss_integration(dS))
	  faces_wT(meshes_faces(fn,M))= Gauss_integration(w_i)/(Gauss_integration(dS))
	  faces_pT(meshes_faces(fn,M))= Gauss_integration(p_i)/(Gauss_integration(dS))
    endif 
	if(is_bx .and. is_bx_TD) faces_bx(meshes_faces(fn,M))= Gauss_integration(bx_i)/(Gauss_integration(dS))
	if(is_by .and. is_by_TD) faces_by(meshes_faces(fn,M))= Gauss_integration(by_i)/(Gauss_integration(dS))
	if(is_bz .and. is_bz_TD) faces_bz(meshes_faces(fn,M))= Gauss_integration(bz_i)/(Gauss_integration(dS))

! -----------------------------------
! bottom
    fn=5
    do i=1,Gauss_int_p
      do j=1,Gauss_int_p
        do l=1,Gauss_int_p
           z_i=-1
           call evaluate_curvilinear(MESHES_ARRAY(M),C,Gauss_x_i(i),Gauss_x_i(j),z_i)
           call map_ti(Gauss_x_i(l),t)
           dS(i,j,l)=C%detJ*sqrt(C%gij(3,3))
			if(is_analytical .and. is_analytical_TD) then
				u_i(i,j,l)=functions_id (C%x,C%y,C%z,t,velx1)*dS(i,j,l)
				v_i(i,j,l)=functions_id (C%x,C%y,C%z,t,velx2)*dS(i,j,l)
				w_i(i,j,l)=functions_id (C%x,C%y,C%z,t,velx3)*dS(i,j,l)
				p_i(i,j,l)=functions_id (C%x,C%y,C%z,t,PPE)*dS(i,j,l)
			endif 
			if(is_bx .and. is_bx_TD) bx_i(i,j,l)=functions_id (C%x,C%y,C%z,t,bx_fun)*dS(i,j,l)
			if(is_by .and. is_by_TD) by_i(i,j,l)=functions_id (C%x,C%y,C%z,t,by_fun)*dS(i,j,l)
			if(is_bz .and. is_bz_TD) bz_i(i,j,l)=functions_id (C%x,C%y,C%z,t,bz_fun)*dS(i,j,l)
        enddo
      enddo
    enddo
    if(is_analytical .and. is_analytical_TD) then
	  faces_uT(meshes_faces(fn,M))= Gauss_integration(u_i)/(Gauss_integration(dS))
	  faces_vT(meshes_faces(fn,M))= Gauss_integration(v_i)/(Gauss_integration(dS))
	  faces_wT(meshes_faces(fn,M))= Gauss_integration(w_i)/(Gauss_integration(dS))
	  faces_pT(meshes_faces(fn,M))= Gauss_integration(p_i)/(Gauss_integration(dS))
    endif 
	if(is_bx .and. is_bx_TD) faces_bx(meshes_faces(fn,M))= Gauss_integration(bx_i)/(Gauss_integration(dS))
	if(is_by .and. is_by_TD) faces_by(meshes_faces(fn,M))= Gauss_integration(by_i)/(Gauss_integration(dS))
	if(is_bz .and. is_bz_TD) faces_bz(meshes_faces(fn,M))= Gauss_integration(bz_i)/(Gauss_integration(dS))
! -----------------------------------
! top
    fn=6
    do i=1,Gauss_int_p
      do j=1,Gauss_int_p
        do l=1,Gauss_int_p
           z_i=1
           call evaluate_curvilinear(MESHES_ARRAY(M),C,Gauss_x_i(i),Gauss_x_i(j),z_i)
           call map_ti(Gauss_x_i(l),t)
           dS(i,j,l)=C%detJ*sqrt(C%gij(3,3))
			if(is_analytical .and. is_analytical_TD) then
				u_i(i,j,l)=functions_id (C%x,C%y,C%z,t,velx1)*dS(i,j,l)
				v_i(i,j,l)=functions_id (C%x,C%y,C%z,t,velx2)*dS(i,j,l)
				w_i(i,j,l)=functions_id (C%x,C%y,C%z,t,velx3)*dS(i,j,l)
				p_i(i,j,l)=functions_id (C%x,C%y,C%z,t,PPE)*dS(i,j,l)
			endif 
			if(is_bx .and. is_bx_TD) bx_i(i,j,l)=functions_id (C%x,C%y,C%z,t,bx_fun)*dS(i,j,l)
			if(is_by .and. is_by_TD) by_i(i,j,l)=functions_id (C%x,C%y,C%z,t,by_fun)*dS(i,j,l)
			if(is_bz .and. is_bz_TD) bz_i(i,j,l)=functions_id (C%x,C%y,C%z,t,bz_fun)*dS(i,j,l)
        enddo
      enddo
    enddo
    if(is_analytical .and. is_analytical_TD) then
	  faces_uT(meshes_faces(fn,M))= Gauss_integration(u_i)/(Gauss_integration(dS))
	  faces_vT(meshes_faces(fn,M))= Gauss_integration(v_i)/(Gauss_integration(dS))
	  faces_wT(meshes_faces(fn,M))= Gauss_integration(w_i)/(Gauss_integration(dS))
	  faces_pT(meshes_faces(fn,M))= Gauss_integration(p_i)/(Gauss_integration(dS))
    endif 
	if(is_bx .and. is_bx_TD) faces_bx(meshes_faces(fn,M))= Gauss_integration(bx_i)/(Gauss_integration(dS))
	if(is_by .and. is_by_TD) faces_by(meshes_faces(fn,M))= Gauss_integration(by_i)/(Gauss_integration(dS))
	if(is_bz .and. is_bz_TD) faces_bz(meshes_faces(fn,M))= Gauss_integration(bz_i)/(Gauss_integration(dS))
  end subroutine time_dependent_face_averages


  subroutine time_dependent_mesh_averages (M)
! -----------------------------------
    integer, intent (in)      :: M
    real, dimension (1:Gauss_int_p,1:Gauss_int_p,1:Gauss_int_p,1:Gauss_int_p) ::bx_i,by_i,bz_i,u_i,v_i,w_i
    real                            :: t
    integer                         :: i,j,k,l
    type(curvi3d)                   :: C

! -----------------------------------
    do i=1,Gauss_int_p
       do j=1,Gauss_int_p
          do k=1,Gauss_int_p
             do l=1,Gauss_int_p
                call evaluate_curvilinear(MESHES_ARRAY(M),C,Gauss_x_i(i),Gauss_x_i(j),Gauss_x_i(k))
                call map_ti(Gauss_x_i(l),t)
				if(is_analytical .and. is_analytical_TD) then
					u_i(i,j,k,l)=functions_id (C%x,C%y,C%z,t,velx1)*C%detJ
					v_i(i,j,k,l)=functions_id (C%x,C%y,C%z,t,velx2)*C%detJ
					w_i(i,j,k,l)=functions_id (C%x,C%y,C%z,t,velx3)*C%detJ
				endif 
				if(is_bx .and. is_bx_TD) bx_i(i,j,k,l)=functions_id (C%x,C%y,C%z,t,bx_fun)*C%detJ
				if(is_by .and. is_by_TD) by_i(i,j,k,l)=functions_id (C%x,C%y,C%z,t,by_fun)*C%detJ
				if(is_bz .and. is_bz_TD) bz_i(i,j,k,l)=functions_id (C%x,C%y,C%z,t,bz_fun)*C%detJ
             enddo
          enddo
       enddo
    enddo
    if(is_analytical .and. is_analytical_TD) then
	  meshes_uT(M)= Gauss_integration_time(u_i)/(2*meshes_volume(M))
	  meshes_vT(M)= Gauss_integration_time(v_i)/(2*meshes_volume(M))
	  meshes_wT(M)= Gauss_integration_time(w_i)/(2*meshes_volume(M))
    endif 
	if(is_bx .and. is_bx_TD) meshes_bx(M)= Gauss_integration_time(bx_i)/(2*meshes_volume(M))
	if(is_by .and. is_by_TD) meshes_by(M)= Gauss_integration_time(by_i)/(2*meshes_volume(M))
	if(is_bz .and. is_bz_TD) meshes_bz(M)= Gauss_integration_time(bz_i)/(2*meshes_volume(M))
  end subroutine time_dependent_mesh_averages



subroutine evaluate_pressure_bc_array
  implicit none
  integer             :: i,j
  integer             :: M1,M2,M3,M4,fn
  real,dimension(1:5) :: n,t,s         !!  distance of each points from normal (n), binormal (s) and tangential (t)
  real,dimension(1:3) :: n_vec, t_vec, s_vec
  real                :: D

	!$OMP     PARALLEL DO PRIVATE (i,j,M1,M2,M3,M4,fn,n,t,s,n_vec, t_vec, s_vec,D)
	do i=1,tot_faces_bcs
		n_vec=faces_normal      (:,BCs_faces(i))
		t_vec=faces_tangential  (:,BCs_faces(i))
		s_vec=faces_bitangential(:,BCs_faces(i))
		M1=BCs_faces_mesh_num(i)
		BCs_faces_points(1,i)=BCs_faces(i)
		! Faces on the left
		if(BCs_faces_position(i)==1) then
			fn=2
			M2=meshes_neighbors(fn,M1)
			M3=meshes_neighbors(4,M1)
			if(M3==0) M3=meshes_neighbors(3,M1)
			M4=meshes_neighbors(6,M1)
			if(M4==0) M4=meshes_neighbors(5,M1)
		elseif(BCs_faces_position(i)==2) then
			fn=1
			M2=meshes_neighbors(fn,M1)
			M3=meshes_neighbors(4,M1)
			if(M3==0) M3=meshes_neighbors(3,M1)
			M4=meshes_neighbors(6,M1)
			if(M4==0) M4=meshes_neighbors(5,M1)
		elseif(BCs_faces_position(i)==3) then
			fn=4
			M2=meshes_neighbors(fn,M1)
			M3=meshes_neighbors(2,M1)
			if(M3==0) M3=meshes_neighbors(1,M1)
			M4=meshes_neighbors(6,M1)
			if(M4==0) M4=meshes_neighbors(5,M1)
		elseif(BCs_faces_position(i)==4) then
			fn=3
			M2=meshes_neighbors(fn,M1)
			M3=meshes_neighbors(2,M1)
			if(M3==0) M3=meshes_neighbors(1,M1)
			M4=meshes_neighbors(6,M1)
			if(M4==0) M4=meshes_neighbors(5,M1)
		elseif(BCs_faces_position(i)==5) then
			fn=6
			M2=meshes_neighbors(fn,M1)
			M3=meshes_neighbors(2,M1)
			if(M3==0) M3=meshes_neighbors(1,M1)
			M4=meshes_neighbors(4,M1)
			if(M4==0) M4=meshes_neighbors(3,M1)
		elseif(BCs_faces_position(i)==6) then
			fn=5
			M2=meshes_neighbors(fn,M1)
			M3=meshes_neighbors(2,M1)
			if(M3==0) M3=meshes_neighbors(1,M1)
			M4=meshes_neighbors(4,M1)
			if(M4==0) M4=meshes_neighbors(3,M1)
		endif
		BCs_faces_points(2,i)=meshes_faces(fn,M1)
		BCs_faces_points(3,i)=meshes_faces(fn,M2)
		BCs_faces_points(4,i)=meshes_faces(fn,M3) ! j+1 or j-1
		BCs_faces_points(5,i)=meshes_faces(fn,M4)
		! find the distances
		do j=1,5
			n(j)=dot((faces_mid_coordinates(:,BCs_faces_points(j,i))-faces_mid_coordinates(:,BCs_faces_points(1,i))),n_vec)
			t(j)=dot((faces_mid_coordinates(:,BCs_faces_points(j,i))-faces_mid_coordinates(:,BCs_faces_points(1,i))),t_vec)
			s(j)=dot((faces_mid_coordinates(:,BCs_faces_points(j,i))-faces_mid_coordinates(:,BCs_faces_points(1,i))),s_vec)
		enddo
		D=(n(3)*s(4)*t(2) - n(3)*s(5)*t(2) - n(2)*s(4)*t(3) + n(2)*s(5)*t(3) - n(3)*s(2)*t(4) + n(2)*s(3)*t(4) - n(2)*s(5)*t(4) +  &
           n(3)*s(5)*t(4) + n(5)*s(3)*t(2) - n(5)*s(4)*t(2) - n(5)*s(2)*t(3) + n(5)*s(4)*t(3) + n(5)*s(2)*t(4) - n(5)*s(3)*t(4) +  &
	       n(3)*s(2)*t(5) - n(2)*s(3)*t(5) + n(2)*s(4)*t(5) - n(3)*s(4)*t(5) + n(4)*(s(5)*t(2) + s(2)*t(3) - s(5)*t(3) - s(2)*t(5) + s(3)*(-t(2) + t(5))))
		BCs_faces_a(2,i)=2*(-n(5)*s(4)*t(3) + n(4)*s(5)*t(3) + n(5)*s(3)*t(4) - n(3)*s(5)*t(4) - n(4)*s(3)*t(5) + n(3)*s(4)*t(5))/(n(2)**2.*D)
		BCs_faces_a(3,i)=2*( n(5)*s(4)*t(2) - n(4)*s(5)*t(2) - n(5)*s(2)*t(4) + n(2)*s(5)*t(4) + n(4)*s(2)*t(5) - n(2)*s(4)*t(5))/(n(3)**2.*D)
		BCs_faces_a(4,i)=2*( n(5)*s(3)*t(2) - n(3)*s(5)*t(2) - n(5)*s(2)*t(3) + n(2)*s(5)*t(3) + n(3)*s(2)*t(5) - n(2)*s(3)*t(5))/(n(4)**2.*D)
		BCs_faces_a(5,i)=2*( n(4)*s(3)*t(2) - n(3)*s(4)*t(2) - n(4)*s(2)*t(3) + n(2)*s(4)*t(3) + n(3)*s(2)*t(4) - n(2)*s(3)*t(4))/(n(5)**2.*D)
		BCs_faces_a(1 ,i)=-BCs_faces_a(2,i)-BCs_faces_a(3,i)-BCs_faces_a(4,i)-BCs_faces_a(5,i)
		BCs_faces_a(6 ,i)=-BCs_faces_a(2,i)*n(2)            -BCs_faces_a(3,i)*n(3)            -BCs_faces_a(4,i)*n(4)            -BCs_faces_a(5,i)*n(5)
		BCs_faces_a(7 ,i)=-BCs_faces_a(2,i)*t(2)            -BCs_faces_a(3,i)*t(3)            -BCs_faces_a(4,i)*t(4)            -BCs_faces_a(5,i)*t(5)
		BCs_faces_a(8 ,i)=-BCs_faces_a(2,i)*s(2)            -BCs_faces_a(3,i)*s(3)            -BCs_faces_a(4,i)*s(4)            -BCs_faces_a(5,i)*s(5)
		BCs_faces_a(9 ,i)=-BCs_faces_a(2,i)*n(2)*t(2)       -BCs_faces_a(3,i)*n(3)*t(3)       -BCs_faces_a(4,i)*n(4)*t(4)       -BCs_faces_a(5,i)*n(5)*t(5)
		BCs_faces_a(10,i)=-BCs_faces_a(2,i)*n(2)*s(2)       -BCs_faces_a(3,i)*n(3)*s(3)       -BCs_faces_a(4,i)*n(4)*s(4)       -BCs_faces_a(5,i)*n(5)*s(5)
		BCs_faces_a(11,i)=-BCs_faces_a(2,i)*s(2)*t(2)       -BCs_faces_a(3,i)*s(3)*t(3)       -BCs_faces_a(4,i)*s(4)*t(4)       -BCs_faces_a(5,i)*s(5)*t(5)
		BCs_faces_a(12,i)=-BCs_faces_a(2,i)*t(2)**2./2.     -BCs_faces_a(3,i)*t(3)**2./2.     -BCs_faces_a(4,i)*t(4)**2./2.     -BCs_faces_a(5,i)*t(5)**2./2.
		BCs_faces_a(13,i)=-BCs_faces_a(2,i)*s(2)**2./2.     -BCs_faces_a(3,i)*s(3)**2./2.     -BCs_faces_a(4,i)*s(4)**2./2.     -BCs_faces_a(5,i)*s(5)**2./2.
		BCs_faces_a(14,i)=-BCs_faces_a(2,i)*n(2)*t(2)**2./2.-BCs_faces_a(3,i)*n(3)*t(3)**2./2.-BCs_faces_a(4,i)*n(4)*t(4)**2./2.-BCs_faces_a(5,i)*n(5)*t(5)**2./2.
		BCs_faces_a(15,i)=-BCs_faces_a(2,i)*n(2)*s(2)**2./2.-BCs_faces_a(3,i)*n(3)*s(3)**2./2.-BCs_faces_a(4,i)*n(4)*s(4)**2./2.-BCs_faces_a(5,i)*n(5)*s(5)**2./2.
		BCs_faces_a(16,i)=-BCs_faces_a(2,i)*s(2)*t(2)**2./2.-BCs_faces_a(3,i)*s(3)*t(3)**2./2.-BCs_faces_a(4,i)*s(4)*t(4)**2./2.-BCs_faces_a(5,i)*s(5)*t(5)**2./2.
		BCs_faces_a(17,i)=-BCs_faces_a(2,i)*t(2)*s(2)**2./2.-BCs_faces_a(3,i)*t(3)*s(3)**2./2.-BCs_faces_a(4,i)*t(4)*s(4)**2./2.-BCs_faces_a(5,i)*t(5)*s(5)**2./2.
		BCs_faces_a(18,i)=-BCs_faces_a(2,i)*n(2)*t(2)*s(2)  -BCs_faces_a(3,i)*n(3)*t(3)*s(3)  -BCs_faces_a(4,i)*n(4)*t(4)*s(4)  -BCs_faces_a(5,i)*n(5)*t(5)*s(5)
		BCs_faces_a(19,i)=-BCs_faces_a(2,i)*t(2)**3./6.     -BCs_faces_a(3,i)*t(3)**3./6.     -BCs_faces_a(4,i)*t(4)**3./6.     -BCs_faces_a(5,i)*t(5)**3./6.
		BCs_faces_a(20,i)=-BCs_faces_a(2,i)*s(2)**3./6.     -BCs_faces_a(3,i)*s(3)**3./6.     -BCs_faces_a(4,i)*s(4)**3./6.     -BCs_faces_a(5,i)*s(5)**3./6.
		
	enddo
	!$OMP     END PARALLEL DO

end subroutine evaluate_pressure_bc_array


subroutine fill_global_arrays
    integer i,j,k,l
	do j=1,6
		meshes_neighbors(j,1:NMESH)=MESHES_ARRAY(1:NMESH)%neighbors(j)
		meshes_faces(j,1:NMESH)=MESHES_ARRAY(1:NMESH)%faces(j)
	enddo
	do i=1,2
        faces_shared_meshes(i,1:NFACE)=FACES_ARRAY(1:NFACE)%shared_meshes(i)	
	enddo
	faces_BC=FACES_ARRAY(1:NFACE)%BC

! ================================	
! allocate BCs_faces
	l=0
	do i=2,NBC
		l=l+BCS_ARRAY(i)%n_faces
	enddo
	tot_faces_bcs=l
	print*, "Total # of faces that have BCs = ", tot_faces_bcs
	allocate(BCs_faces(1:l))
	allocate(BCs_faces_id(1:l))
	allocate(BCs_faces_position(1:l))
	allocate(BCs_faces_mesh_num(1:l))
	allocate(BCs_faces_a (1:20,1:l))
	allocate(BCs_faces_ap(1:8,1:l))
        allocate(BCs_faces_aT(1:3,1:l))
	allocate(BCs_faces_points(1:5,1:l))
	allocate(BCs_faces_v_n(1:l))
	allocate(BCs_faces_v_nt(1:l))
	allocate(BCs_faces_v_ns(1:l))
	allocate(BCs_faces_v_ntt(1:l))
	allocate(BCs_faces_v_nss(1:l))
	allocate(BCs_faces_v_nts(1:l))
! Fill the arrays
	l=0
	do i=2,NBC
		do j=BCS_ARRAY(i)%first,BCS_ARRAY(i)%last
			l=l+1
			BCs_faces(l)=j
			BCs_faces_id(l)=BCS_ARRAY(i)%id
		enddo
	enddo
! Find the position of the face in each mesh
	do i=1,tot_faces_bcs
	   if(faces_shared_meshes(2,BCs_faces(i)) == 0) then
		  j=faces_shared_meshes(1,BCs_faces(i))
	   elseif(faces_shared_meshes(1,BCs_faces(i)) == 0) then
		  j=faces_shared_meshes(2,BCs_faces(i))
	   endif
	   BCs_faces_mesh_num(i)=j
	   BCs_faces_position(i)=position_in_array(meshes_faces(:,j),6,BCs_faces(i))
!	   print*, BCs_faces(i), BCs_faces_mesh_num(i), BCs_faces_position(i), BCs_faces_id(i)
	enddo
! ================================
end subroutine fill_global_arrays



subroutine deallocate_arrays
	deallocate(BCs_faces)
	deallocate(BCs_faces_id)
	deallocate(BCs_faces_ap)
	deallocate(BCs_faces_a)
	deallocate(BCs_faces_position)
	deallocate(BCs_faces_mesh_num)
	deallocate(BCs_faces_points)
end subroutine deallocate_arrays


! =========================================================================


end module class_geometry

