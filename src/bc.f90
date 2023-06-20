! Imposing boundary conditions
module class_bc
  use class_general
  use class_nim
  use class_parameters
  use class_mesh
  use class_geometry
  implicit none

  contains


  subroutine apply_all_dirichlet_bcs
  ! ---------------------------------
    call apply_dirichlet_velocity_bcs
    call apply_dirichlet_temperature_bcs
!	call apply_dirichlet_pressure_bcs

  end subroutine apply_all_dirichlet_bcs


  subroutine apply_dirichlet_velocity_bcs
    implicit none
    integer              :: i
	!$OMP     PARALLEL DO PRIVATE (i)
    do i=1,tot_faces_bcs
! Fill BCs for fields
	   faces_u(BCs_faces(i))=apply_dirichlet_bc_i(BCs_faces_mesh_num(i),BCs_faces_position(i),velx1,BCs_faces_id(i))
	   faces_v(BCs_faces(i))=apply_dirichlet_bc_i(BCs_faces_mesh_num(i),BCs_faces_position(i),velx2,BCs_faces_id(i))
	   faces_w(BCs_faces(i))=apply_dirichlet_bc_i(BCs_faces_mesh_num(i),BCs_faces_position(i),velx3,BCs_faces_id(i))
    enddo
	!$OMP     END PARALLEL DO
  end subroutine apply_dirichlet_velocity_bcs


  subroutine apply_dirichlet_pressure_bcs
    implicit none
    integer              :: i
	!$OMP     PARALLEL DO PRIVATE (i)
    do i=1,tot_faces_bcs
! Fill BCs for pressure
	   faces_p(BCs_faces(i))=apply_dirichlet_bc_i(BCs_faces_mesh_num(i),BCs_faces_position(i),PPE,BCs_faces_id(i))
    enddo
	!$OMP     END PARALLEL DO
  end subroutine apply_dirichlet_pressure_bcs


  subroutine apply_dirichlet_temperature_bcs
    implicit none
    integer              :: i
	!$OMP     PARALLEL DO PRIVATE (i)
    do i=1,tot_faces_bcs
! Fill BCs for pressure
	   faces_T(BCs_faces(i))=apply_dirichlet_bc_i(BCs_faces_mesh_num(i),BCs_faces_position(i),CDE,BCs_faces_id(i))
    enddo
	!$OMP     END PARALLEL DO
  end subroutine apply_dirichlet_temperature_bcs

  function apply_dirichlet_bc_i(M,n,field_id,bc_id) result(TT)
  implicit none
! -----------------------------------
    integer, intent (in)                          :: M
    integer, intent(in)                           :: n,field_id,bc_id
    real, dimension (1:Gauss_int_p,1:Gauss_int_p) :: TT_i, dS
	real      :: TT, area
    integer   :: i,j,k
    type(curvi3d)    :: C
! -----------------------------------
    do i=1,Gauss_int_p
       do j=1,Gauss_int_p
			if(n==1 .or. n==2) then
			   if (n==1) call evaluate_curvilinear(MESHES_ARRAY(M),C,-1.0,Gauss_x_i(i),Gauss_x_i(j))
			   if (n==2) call evaluate_curvilinear(MESHES_ARRAY(M),C, 1.0,Gauss_x_i(i),Gauss_x_i(j))
			   dS(i,j)=C%detJ*sqrt(C%gij(1,1))
			elseif(n==3 .or. n==4) then
			   if (n==3) call evaluate_curvilinear(MESHES_ARRAY(M),C,Gauss_x_i(i),-1.0,Gauss_x_i(j))
			   if (n==4) call evaluate_curvilinear(MESHES_ARRAY(M),C,Gauss_x_i(i), 1.0,Gauss_x_i(j))
			   dS(i,j)=C%detJ*sqrt(C%gij(2,2))
			elseif(n==5 .or. n==6) then
			   if (n==5) call evaluate_curvilinear(MESHES_ARRAY(M),C,Gauss_x_i(i),Gauss_x_i(j),-1.0)
			   if (n==6) call evaluate_curvilinear(MESHES_ARRAY(M),C,Gauss_x_i(i),Gauss_x_i(j), 1.0)
			   dS(i,j)=C%detJ*sqrt(C%gij(3,3))
			endif
			TT_i(i,j)=BCs_id(C%x,C%y,C%z,field_id,bc_id)*dS(i,j)
       enddo
    enddo
	TT=Gauss_integration2d(TT_i)/faces_area(meshes_faces(n,M))
 end function apply_dirichlet_bc_i

  function find_bcs_velocity_derivatives(M,n,bc_id,dvdt_id) result(TT)
  implicit none
! -----------------------------------
    integer, intent (in)                          :: M
    integer, intent(in)                           :: n,bc_id,dvdt_id
    real, dimension (1:Gauss_int_p,1:Gauss_int_p) :: TT_i, dS
	real      :: TT, area
    integer   :: i,j,k
    type(curvi3d)    :: C
! -----------------------------------
    do i=1,Gauss_int_p
       do j=1,Gauss_int_p
			if(n==1 .or. n==2) then
			   if (n==1) call evaluate_curvilinear(MESHES_ARRAY(M),C,-1.0,Gauss_x_i(i),Gauss_x_i(j))
			   if (n==2) call evaluate_curvilinear(MESHES_ARRAY(M),C, 1.0,Gauss_x_i(i),Gauss_x_i(j))
			   dS(i,j)=C%detJ*sqrt(C%gij(1,1))
			elseif(n==3 .or. n==4) then
			   if (n==3) call evaluate_curvilinear(MESHES_ARRAY(M),C,Gauss_x_i(i),-1.0,Gauss_x_i(j))
			   if (n==4) call evaluate_curvilinear(MESHES_ARRAY(M),C,Gauss_x_i(i), 1.0,Gauss_x_i(j))
			   dS(i,j)=C%detJ*sqrt(C%gij(2,2))
			elseif(n==5 .or. n==6) then
			   if (n==5) call evaluate_curvilinear(MESHES_ARRAY(M),C,Gauss_x_i(i),Gauss_x_i(j),-1.0)
			   if (n==6) call evaluate_curvilinear(MESHES_ARRAY(M),C,Gauss_x_i(i),Gauss_x_i(j), 1.0)
			   dS(i,j)=C%detJ*sqrt(C%gij(3,3))
			endif
			TT_i(i,j)=BCs_derivatives(C%x,C%y,C%z,bc_id,dvdt_id)*dS(i,j)
       enddo
    enddo
	TT=Gauss_integration2d(TT_i)/faces_area(meshes_faces(n,M))
 end function find_bcs_velocity_derivatives


subroutine find_derivatives_pressure_boundary
	integer                  :: i
	integer,dimension(1:6)  :: n1_faces

	!$OMP     PARALLEL DO PRIVATE (i,n1_faces)
    do i=1,tot_faces_bcs	
		BCs_faces_v_n(i)  =find_bcs_velocity_derivatives(BCs_faces_mesh_num(i),BCs_faces_position(i),BCs_faces_id(i),bc_dvdn)
		BCs_faces_v_nt(i) =find_bcs_velocity_derivatives(BCs_faces_mesh_num(i),BCs_faces_position(i),BCs_faces_id(i),bc_dvdnt)
		BCs_faces_v_ns(i) =find_bcs_velocity_derivatives(BCs_faces_mesh_num(i),BCs_faces_position(i),BCs_faces_id(i),bc_dvdns)
		BCs_faces_v_ntt(i)=find_bcs_velocity_derivatives(BCs_faces_mesh_num(i),BCs_faces_position(i),BCs_faces_id(i),bc_dvdntt)
		BCs_faces_v_nss(i)=find_bcs_velocity_derivatives(BCs_faces_mesh_num(i),BCs_faces_position(i),BCs_faces_id(i),bc_dvdnss)
		BCs_faces_v_nts(i)=find_bcs_velocity_derivatives(BCs_faces_mesh_num(i),BCs_faces_position(i),BCs_faces_id(i),bc_dvdnts)

	enddo

end subroutine find_derivatives_pressure_boundary



subroutine apply_pressure_bc
   integer                 :: i,j,M1,M2,M3,fn
   real,dimension(1:3)     :: v,b,n
   real,dimension(1:5)     :: v_n		! normal velocities for the 5 points
   real                    :: b_n,d2udn2, B1,B2,B3		! B1 is the pseudo source terms. B2,B3 are the second derivetives
   integer,dimension(1:6)  :: n1_faces					! faces of the mesh
	!$OMP     PARALLEL DO PRIVATE (i,j,n,v,b,b_n,v_n,d2udn2,M1,M2,M3,fn,B1,B2,B3,n1_faces)
    do i=1,tot_faces_bcs
		! Get the mesh where the face exists
		M1=BCs_faces_mesh_num(i)
		! Get the faces of the mesh
		n1_faces=meshes_faces(:,M1)
		! Get the normal to the face
		n=faces_normal(:,BCs_faces(i))
		! Get the b-normal
		b=(/faces_bx(BCs_faces(i)),faces_by(BCs_faces(i)),faces_bz(BCs_faces(i))/)
		b_n=dot(b,n)
		! Get the normal velocities
		do j=1,5
			v=(/faces_u(BCs_faces_points(j,i)),faces_v(BCs_faces_points(j,i)),faces_w(BCs_faces_points(j,i))/)
			v_n(j)=dot(v,n)
		enddo
!		vbc_dvdn=find_bcs_velocity_derivatives(BCs_faces_mesh_num(i),BCs_faces_position(i),BCs_faces_id(i),bc_dvdn)
		d2udn2=-BCs_faces_a(1,i)*v_n(1)-BCs_faces_a(2,i)*v_n(2)-BCs_faces_a(3,i)*v_n(3)-BCs_faces_a(4,i)*v_n(4)-BCs_faces_a(5,i)*v_n(5)&
			   -BCs_faces_a(6 ,i)*BCs_faces_v_n  (i)-BCs_faces_a(9 ,i)*BCs_faces_v_nt (i)-BCs_faces_a(10,i)*BCs_faces_v_ns (i)&
			   -BCs_faces_a(14,i)*BCs_faces_v_ntt(i)-BCs_faces_a(15,i)*BCs_faces_v_nss(i)-BCs_faces_a(18,i)*BCs_faces_v_nts(i)
!		d2udn2=0.
	! ---------> The integration part <----------
		if(BCs_faces_position(i)==1) then
			fn=1
			M2=meshes_neighbors(3,M1)
			M3=meshes_neighbors(4,M1)
			if(M2==0) then
				B2=(-3.*faces_p(meshes_faces(fn,M1))+4.*faces_p(meshes_faces(fn,M3))-faces_p(meshes_faces(fn,meshes_neighbors(4,M3))))/4.0*meshes_gij_faces(1,2,fn,M1)
			elseif(M3==0) then
				B2=( 3.*faces_p(meshes_faces(fn,M1))-4.*faces_p(meshes_faces(fn,M2))+faces_p(meshes_faces(fn,meshes_neighbors(3,M2))))/4.0*meshes_gij_faces(1,2,fn,M1)
			else
				B2=(faces_p(meshes_faces(fn,M3))-faces_p(meshes_faces(fn,M2)))/4.0*meshes_gij_faces(1,2,fn,M1)
			endif
			M2=meshes_neighbors(5,M1)
			M3=meshes_neighbors(6,M1)
			if(M2==0) then
				B3=(-3.*faces_p(meshes_faces(fn,M1))+4.*faces_p(meshes_faces(fn,M3))-faces_p(meshes_faces(fn,meshes_neighbors(6,M3))))/4.0*meshes_gij_faces(1,3,fn,M1)
			elseif(M3==0) then
				B3=( 3.*faces_p(meshes_faces(fn,M1))-4.*faces_p(meshes_faces(fn,M2))+faces_p(meshes_faces(fn,meshes_neighbors(5,M2))))/4.0*meshes_gij_faces(1,3,fn,M1)
			else
				B3=(faces_p(meshes_faces(fn,M3))-faces_p(meshes_faces(fn,M2)))/4.0*meshes_gij_faces(1,3,fn,M1)
			endif
			B1=BCs_faces_ap(8,i)*(meshes_R_D_p(M1)+meshes_R_B_p(M1))
			faces_p(n1_faces(1))=(1-w_p2)*faces_p(n1_faces(1)) + w_p2/(BCs_faces_ap(1,i)) * ( &
								 -BCs_faces_ap(2,i)*faces_p(n1_faces(2)) &
								 -BCs_faces_ap(3,i)*faces_p(n1_faces(3)) &
								 -BCs_faces_ap(4,i)*faces_p(n1_faces(4)) &
								 -BCs_faces_ap(5,i)*faces_p(n1_faces(5)) &
								 -BCs_faces_ap(6,i)*faces_p(n1_faces(6)) &
								 +B1-B2-B3+(nu*rho*d2udn2-b_n*rho))
		elseif(BCs_faces_position(i)==2) then
			fn=2
			M2=meshes_neighbors(3,M1)
			M3=meshes_neighbors(4,M1)
			if(M2==0) then
				B2=(-3.*faces_p(meshes_faces(fn,M1))+4.*faces_p(meshes_faces(fn,M3))-faces_p(meshes_faces(fn,meshes_neighbors(4,M3))))/4.0*meshes_gij_faces(1,2,fn,M1)
			elseif(M3==0) then
				B2=( 3.*faces_p(meshes_faces(fn,M1))-4.*faces_p(meshes_faces(fn,M2))+faces_p(meshes_faces(fn,meshes_neighbors(3,M2))))/4.0*meshes_gij_faces(1,2,fn,M1)
			else
				B2=(faces_p(meshes_faces(fn,M3))-faces_p(meshes_faces(fn,M2)))/4.0*meshes_gij_faces(1,2,fn,M1)
			endif
			M2=meshes_neighbors(5,M1)
			M3=meshes_neighbors(6,M1)
			if(M2==0) then
				B3=(-3.*faces_p(meshes_faces(fn,M1))+4.*faces_p(meshes_faces(fn,M3))-faces_p(meshes_faces(fn,meshes_neighbors(6,M3))))/4.0*meshes_gij_faces(1,3,fn,M1)
			elseif(M3==0) then
				B3=( 3.*faces_p(meshes_faces(fn,M1))-4.*faces_p(meshes_faces(fn,M2))+faces_p(meshes_faces(fn,meshes_neighbors(5,M2))))/4.0*meshes_gij_faces(1,3,fn,M1)
			else
				B3=(faces_p(meshes_faces(fn,M3))-faces_p(meshes_faces(fn,M2)))/4.0*meshes_gij_faces(1,3,fn,M1)
			endif
			B1=BCs_faces_ap(8,i)*(meshes_R_D_p(M1)+meshes_R_B_p(M1))
			faces_p(n1_faces(2))=(1-w_p2)*faces_p(n1_faces(2)) + w_p2/(BCs_faces_ap(2,i)) * ( &
								 -BCs_faces_ap(1,i)*faces_p(n1_faces(1)) &
								 -BCs_faces_ap(3,i)*faces_p(n1_faces(3)) &
								 -BCs_faces_ap(4,i)*faces_p(n1_faces(4)) &
								 -BCs_faces_ap(5,i)*faces_p(n1_faces(5)) &
								 -BCs_faces_ap(6,i)*faces_p(n1_faces(6)) &
								 +B1-B2-B3+(nu*rho*d2udn2-b_n*rho))
		elseif(BCs_faces_position(i)==3) then
			fn=3
			M2=meshes_neighbors(1,M1)
			M3=meshes_neighbors(2,M1)
			if(M2==0) then
				B2=(-3.*faces_p(meshes_faces(fn,M1))+4.*faces_p(meshes_faces(fn,M3))-faces_p(meshes_faces(fn,meshes_neighbors(2,M3))))/4.0*meshes_gij_faces(2,1,fn,M1)
			elseif(M3==0) then
				B2=( 3.*faces_p(meshes_faces(fn,M1))-4.*faces_p(meshes_faces(fn,M2))+faces_p(meshes_faces(fn,meshes_neighbors(1,M2))))/4.0*meshes_gij_faces(2,1,fn,M1)
			else
				B2=(faces_p(meshes_faces(fn,M3))-faces_p(meshes_faces(fn,M2)))/4.0*meshes_gij_faces(2,1,fn,M1)
			endif
			M2=meshes_neighbors(5,M1)
			M3=meshes_neighbors(6,M1)
			if(M2==0) then
				B3=(-3.*faces_p(meshes_faces(fn,M1))+4.*faces_p(meshes_faces(fn,M3))-faces_p(meshes_faces(fn,meshes_neighbors(6,M3))))/4.0*meshes_gij_faces(2,3,fn,M1)
			elseif(M3==0) then
				B3=( 3.*faces_p(meshes_faces(fn,M1))-4.*faces_p(meshes_faces(fn,M2))+faces_p(meshes_faces(fn,meshes_neighbors(5,M2))))/4.0*meshes_gij_faces(2,3,fn,M1)
			else
				B3=(faces_p(meshes_faces(fn,M3))-faces_p(meshes_faces(fn,M2)))/4.0*meshes_gij_faces(2,3,fn,M1)
			endif
			B1=BCs_faces_ap(8,i)*(meshes_R_D_p(M1)+meshes_R_B_p(M1))
			faces_p(n1_faces(3))=(1-w_p2)*faces_p(n1_faces(3)) + w_p2/(BCs_faces_ap(3,i)) * ( &
								 -BCs_faces_ap(1,i)*faces_p(n1_faces(1)) &
								 -BCs_faces_ap(2,i)*faces_p(n1_faces(2)) &
								 -BCs_faces_ap(4,i)*faces_p(n1_faces(4)) &
								 -BCs_faces_ap(5,i)*faces_p(n1_faces(5)) &
								 -BCs_faces_ap(6,i)*faces_p(n1_faces(6)) &
								 +B1-B2-B3+(nu*rho*d2udn2-b_n*rho))
		elseif(BCs_faces_position(i)==4) then
			fn=4
			M2=meshes_neighbors(1,M1)
			M3=meshes_neighbors(2,M1)
			if(M2==0) then
				B2=(-3.*faces_p(meshes_faces(fn,M1))+4.*faces_p(meshes_faces(fn,M3))-faces_p(meshes_faces(fn,meshes_neighbors(2,M3))))/4.0*meshes_gij_faces(2,1,fn,M1)
			elseif(M3==0) then
				B2=( 3.*faces_p(meshes_faces(fn,M1))-4.*faces_p(meshes_faces(fn,M2))+faces_p(meshes_faces(fn,meshes_neighbors(1,M2))))/4.0*meshes_gij_faces(2,1,fn,M1)
			else
				B2=(faces_p(meshes_faces(fn,M3))-faces_p(meshes_faces(fn,M2)))/4.0*meshes_gij_faces(2,1,fn,M1)
			endif
			M2=meshes_neighbors(5,M1)
			M3=meshes_neighbors(6,M1)
			if(M2==0) then
				B3=(-3.*faces_p(meshes_faces(fn,M1))+4.*faces_p(meshes_faces(fn,M3))-faces_p(meshes_faces(fn,meshes_neighbors(6,M3))))/4.0*meshes_gij_faces(2,3,fn,M1)
			elseif(M3==0) then
				B3=( 3.*faces_p(meshes_faces(fn,M1))-4.*faces_p(meshes_faces(fn,M2))+faces_p(meshes_faces(fn,meshes_neighbors(5,M2))))/4.0*meshes_gij_faces(2,3,fn,M1)
			else
				B3=(faces_p(meshes_faces(fn,M3))-faces_p(meshes_faces(fn,M2)))/4.0*meshes_gij_faces(2,3,fn,M1)
			endif
			B1=BCs_faces_ap(8,i)*(meshes_R_D_p(M1)+meshes_R_B_p(M1))
			faces_p(n1_faces(4))=(1-w_p2)*faces_p(n1_faces(4)) + w_p2/(BCs_faces_ap(4,i)) * ( &
								 -BCs_faces_ap(1,i)*faces_p(n1_faces(1)) &
								 -BCs_faces_ap(2,i)*faces_p(n1_faces(2)) &
								 -BCs_faces_ap(3,i)*faces_p(n1_faces(3)) &
								 -BCs_faces_ap(5,i)*faces_p(n1_faces(5)) &
								 -BCs_faces_ap(6,i)*faces_p(n1_faces(6)) &
								 +B1-B2-B3+(nu*rho*d2udn2-b_n*rho))
		elseif(BCs_faces_position(i)==5) then
			fn=5
			M2=meshes_neighbors(1,M1)
			M3=meshes_neighbors(2,M1)
			if(M2==0) then
				B2=(-3.*faces_p(meshes_faces(fn,M1))+4.*faces_p(meshes_faces(fn,M3))-faces_p(meshes_faces(fn,meshes_neighbors(2,M3))))/4.0*meshes_gij_faces(3,1,fn,M1)
			elseif(M3==0) then
				B2=( 3.*faces_p(meshes_faces(fn,M1))-4.*faces_p(meshes_faces(fn,M2))+faces_p(meshes_faces(fn,meshes_neighbors(1,M2))))/4.0*meshes_gij_faces(3,1,fn,M1)
			else
				B2=(faces_p(meshes_faces(fn,M3))-faces_p(meshes_faces(fn,M2)))/4.0*meshes_gij_faces(3,1,fn,M1)
			endif
			M2=meshes_neighbors(3,M1)
			M3=meshes_neighbors(4,M1)
			if(M2==0) then
				B3=(-3.*faces_p(meshes_faces(fn,M1))+4.*faces_p(meshes_faces(fn,M3))-faces_p(meshes_faces(fn,meshes_neighbors(4,M3))))/4.0*meshes_gij_faces(3,2,fn,M1)
			elseif(M3==0) then
				B3=( 3.*faces_p(meshes_faces(fn,M1))-4.*faces_p(meshes_faces(fn,M2))+faces_p(meshes_faces(fn,meshes_neighbors(3,M2))))/4.0*meshes_gij_faces(3,2,fn,M1)
			else
				B3=(faces_p(meshes_faces(fn,M3))-faces_p(meshes_faces(fn,M2)))/4.0*meshes_gij_faces(3,2,fn,M1)
			endif
			B1=BCs_faces_ap(8,i)*(meshes_R_D_p(M1)+meshes_R_B_p(M1))
			faces_p(n1_faces(5))=(1-w_p2)*faces_p(n1_faces(5)) + w_p2/(BCs_faces_ap(5,i)) * ( &
								 -BCs_faces_ap(1,i)*faces_p(n1_faces(1)) &
								 -BCs_faces_ap(2,i)*faces_p(n1_faces(2)) &
								 -BCs_faces_ap(3,i)*faces_p(n1_faces(3)) &
								 -BCs_faces_ap(4,i)*faces_p(n1_faces(4)) &
								 -BCs_faces_ap(6,i)*faces_p(n1_faces(6)) &
								 +B1-B2-B3+(nu*rho*d2udn2-b_n*rho))
		elseif(BCs_faces_position(i)==6) then
			fn=6
			M2=meshes_neighbors(1,M1)
			M3=meshes_neighbors(2,M1)
			if(M2==0) then
				B2=(-3.*faces_p(meshes_faces(fn,M1))+4.*faces_p(meshes_faces(fn,M3))-faces_p(meshes_faces(fn,meshes_neighbors(2,M3))))/4.0*meshes_gij_faces(3,1,fn,M1)
			elseif(M3==0) then
				B2=( 3.*faces_p(meshes_faces(fn,M1))-4.*faces_p(meshes_faces(fn,M2))+faces_p(meshes_faces(fn,meshes_neighbors(1,M2))))/4.0*meshes_gij_faces(3,1,fn,M1)
			else
				B2=(faces_p(meshes_faces(fn,M3))-faces_p(meshes_faces(fn,M2)))/4.0*meshes_gij_faces(3,2,fn,M1)
			endif
			M2=meshes_neighbors(3,M1)
			M3=meshes_neighbors(4,M1)
			if(M2==0) then
				B3=(-3.*faces_p(meshes_faces(fn,M1))+4.*faces_p(meshes_faces(fn,M3))-faces_p(meshes_faces(fn,meshes_neighbors(4,M3))))/4.0*meshes_gij_faces(3,2,fn,M1)
			elseif(M3==0) then
				B3=( 3.*faces_p(meshes_faces(fn,M1))-4.*faces_p(meshes_faces(fn,M2))+faces_p(meshes_faces(fn,meshes_neighbors(3,M2))))/4.0*meshes_gij_faces(3,2,fn,M1)
			else
				B3=(faces_p(meshes_faces(fn,M3))-faces_p(meshes_faces(fn,M2)))/4.0*meshes_gij_faces(3,2,fn,M1)
			endif
			B1=BCs_faces_ap(8,i)*(meshes_R_D_p(M1)+meshes_R_B_p(M1))
			faces_p(n1_faces(6))=(1-w_p2)*faces_p(n1_faces(6)) + w_p2/(BCs_faces_ap(6,i)) * ( &
								 -BCs_faces_ap(1,i)*faces_p(n1_faces(1)) &
								 -BCs_faces_ap(2,i)*faces_p(n1_faces(2)) &
								 -BCs_faces_ap(3,i)*faces_p(n1_faces(3)) &
								 -BCs_faces_ap(4,i)*faces_p(n1_faces(4)) &
								 -BCs_faces_ap(5,i)*faces_p(n1_faces(5)) &
								 +B1-B2-B3+(nu*rho*d2udn2-b_n*rho))
		endif
!		faces_p(i)=faces_pT(i)
		
    enddo
	!$OMP     END PARALLEL DO
end subroutine apply_pressure_bc


subroutine apply_zero_flux
   integer                 :: i,M,f1,f2,f3,f4
   integer,dimension(1:6)  :: n1_faces					! faces of the mesh
	real,dimension(1:4)   :: h
	!$OMP     PARALLEL DO PRIVATE (i,M,n1_faces,h,f1,f2,f3,f4)
    do i=1,tot_faces_bcs
		! Get the mesh where the face exists
		M=BCs_faces_mesh_num(i)
		! Get the faces of the mesh
		n1_faces=meshes_faces(:,M)
		if(BCs_faces_position(i)==1) then
		elseif(BCs_faces_position(i)==2) then
		elseif(BCs_faces_position(i)==3) then
		elseif(BCs_faces_position(i)==4) then
		elseif(BCs_faces_position(i)==5) then
			  f1=BCs_faces(i)
			  f2=meshes_faces(6,M)
			  f3=meshes_faces(6,meshes_neighbors(6,M))
			  f4=meshes_faces(6,meshes_neighbors(6,meshes_neighbors(6,M)))
			  h(1)=sqrt((faces_mid_coordinates(1,f1)-faces_mid_coordinates(1,f2))**2+(faces_mid_coordinates(2,f1)-faces_mid_coordinates(2,f2))**2+(faces_mid_coordinates(3,f1)-faces_mid_coordinates(3,f2))**2)
			  h(2)=sqrt((faces_mid_coordinates(1,f2)-faces_mid_coordinates(1,f3))**2+(faces_mid_coordinates(2,f2)-faces_mid_coordinates(2,f3))**2+(faces_mid_coordinates(3,f2)-faces_mid_coordinates(3,f3))**2)
			  h(3)=sqrt((faces_mid_coordinates(1,f3)-faces_mid_coordinates(1,f4))**2+(faces_mid_coordinates(2,f3)-faces_mid_coordinates(2,f4))**2+(faces_mid_coordinates(3,f3)-faces_mid_coordinates(3,f4))**2)
			  faces_T(f1)=(+((h(1)+h(2))*(h(1)+h(2)+h(3)))/(h(1)*h(2)*(h(2)+h(3)))*faces_T(f2) &
					-(h(1)*(h(1)+h(2)+h(3)))/(h(2)*h(3)*(h(1)+h(2)))*faces_T(f3) &
					+(h(1)*(h(1)+h(2)))/(h(3)*(h(2)+h(3))*(h(1)+h(2)+h(3)))*faces_T(f4))/(1./h(1)+1./(h(1)+h(2))+1./(h(1)+h(2)+h(3)))


!			faces_T(n1_faces(5))=faces_T(n1_faces(6))
		elseif(BCs_faces_position(i)==6) then
			  f1=BCs_faces(i)
			  f2=meshes_faces(5,M)
			  f3=meshes_faces(5,meshes_neighbors(5,M))
			  f4=meshes_faces(5,meshes_neighbors(5,meshes_neighbors(5,M)))
			  h(1)=sqrt((faces_mid_coordinates(1,f1)-faces_mid_coordinates(1,f2))**2+(faces_mid_coordinates(2,f1)-faces_mid_coordinates(2,f2))**2+(faces_mid_coordinates(3,f1)-faces_mid_coordinates(3,f2))**2)
			  h(2)=sqrt((faces_mid_coordinates(1,f2)-faces_mid_coordinates(1,f3))**2+(faces_mid_coordinates(2,f2)-faces_mid_coordinates(2,f3))**2+(faces_mid_coordinates(3,f2)-faces_mid_coordinates(3,f3))**2)
			  h(3)=sqrt((faces_mid_coordinates(1,f3)-faces_mid_coordinates(1,f4))**2+(faces_mid_coordinates(2,f3)-faces_mid_coordinates(2,f4))**2+(faces_mid_coordinates(3,f3)-faces_mid_coordinates(3,f4))**2)
			  faces_T(f1)=( +((h(1)+h(2))*(h(1)+h(2)+h(3)))/(h(1)*h(2)*(h(2)+h(3)))*faces_T(f2) &
					-(h(1)*(h(1)+h(2)+h(3)))/(h(2)*h(3)*(h(1)+h(2)))*faces_T(f3) &
					+(h(1)*(h(1)+h(2)))/(h(3)*(h(2)+h(3))*(h(1)+h(2)+h(3)))*faces_T(f4))/(1./h(1)+1./(h(1)+h(2))+1./(h(1)+h(2)+h(3)))

!			faces_T(n1_faces(6))=faces_T(n1_faces(5))

		endif
    enddo
	!$OMP     END PARALLEL DO
end subroutine apply_zero_flux


end module class_bc

