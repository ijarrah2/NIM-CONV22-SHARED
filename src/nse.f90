! NSE class -> solve NSE using NIM
module class_nse
  use class_general
  use class_nim
  use class_parameters
  use class_mesh
  use class_geometry
  use class_bc
  use omp_lib
  use class_restart
!  use class_prepare_mesh
  implicit none

  contains

! =========================================================
! =========================================================
! ==================== Main Algorithm =====================
! =============== Time Integration Steps ==================
! =========================================================
! =========================================================

  subroutine solve_nse
    implicit none
    integer l,i,j,k,m
! *********************************************************
! ************    INITIALIZE THE DOMAIN    ****************
! *********************************************************

! Setting the current NIM time
	TI=DT
	STEP=1
	
    call initialize_nse


! Apply all Dirichlet velocity boundary conditions - if not TD
    call apply_all_dirichlet_bcs

!    call cpu_time(timing_initialize_2)
	timing_initialize_2=omp_get_wtime()

! *********************************************************
! ************    START THE TIME INTEGRAL    **************
! *********************************************************

    do l=1,NT
! Update the time-dependent variables : bx,by,bz + analytical solution.
	   call evaluate_time_dependent_quantities
! Calculate the velocity NIM constants based on the previous time step.
	if(nim_formulation==2) then
		call calculate_velocity_nim_constants_formulation2
	else
		call calculate_velocity_nim_constants
	endif 
! Main integtation function
       call integrate_nse
		if(is_write_restart .and. mod(STEP,write_restart_step) == 0) call write_restart
		if(mod(STEP,STEP_OUT) == 0)   call write_data
		if(is_check_steady_state) then
!			call check_steady_state
			if(dudt_max<dudt_epsilon .and. dTdt_max<dudt_epsilon) then
				print*, "+++++++++++++++++++++++++++++++++ Steady state is reached ++++++++++++++++++++++++++++++"
				print*, "Simulation is stopped with dudt = ", dudt_max
				if(is_write_restart) call write_restart
				call write_data
				exit
			endif
		endif

! Fix time step and copy the arrays
	   call update_time_step
	   call find_Nu2
  if(is_analytical) call find_rms
    enddo

! *********************************************************
! ************    END THE TIME INTEGRAL    ****************
! *********************************************************

  if(is_analytical) call find_rms
!  call cpu_time(timing_integration)
  timing_integration=omp_get_wtime()

  end subroutine solve_nse


! =========================================================
  subroutine initialize_nse
  
! *********************************************************
! ************    INITIALIZE THE DOMAIN    ****************
! *********************************************************
! evaluate all constants related to mesh: normals, areas, volume, metrics ... etc
! apply initial condition, physical sources (if steady-state)
    call evaluate_all_mesh_constants

	call evaluate_pressure_bc_array

! Evaluate all PPE constans
	call calculate_pressure_nim_constants
	call find_derivatives_pressure_boundary

	if (is_restart) then
		call read_restart
		call write_data
		call update_time_step
!		call evaluate_D
	endif
! Write the initial condition
    if (.not. is_restart) call write_data
  end subroutine initialize_nse
! =========================================================



subroutine integrate_nse
  implicit none
  integer                              :: i,j,ff,cnt
  real                                 :: err_u_max,err_v_max,err_w_max,err_p_max,err_T_max
! cnt is iteration counter
  cnt=0
  err_u_max=1.
  err_v_max=1.
  err_w_max=1.
  err_p_max=1.
  err_T_max=1.

 do while (err_u_max>TOLERANCE_u .or. err_v_max>TOLERANCE_v .or. err_w_max>TOLERANCE_w .or. err_p_max>TOLERANCE_p .or. err_T_max>TOLERANCE_T)
!   do i=1,6000
   cnt=cnt+1
  


! =================================== VELOCITY INTEGRATION ================================

   do ff=1,sweep_u
   ! call calculate_velocity_nim_constants_formulation2
		call update_velocity_sources
		call integrate_velocity_field


   enddo

! =================================== PRESSURE INTEGRATION =================================

   do ff=1,sweep_p

		call update_pressure_sources
		call apply_pressure_bc
		call integrate_pressure_field

   enddo
   do ff=1,sweep_T
		call calculate_temperature_nim_constants_formulation2
		call apply_zero_flux
		call update_temperature_sources
		call integrate_temperature_field
   enddo
   call find_error(err_u_max,err_v_max,err_w_max,err_p_max,err_T_max)


!  print*, err_u_max,err_v_max,err_w_max,err_p_max,err_T_max

  enddo

  if(is_check_steady_state) call check_steady_state
  call evaluate_D
  write (*,'(A, I5, A, I4, A, 1p,E11.4, A, 1p,E11.4, A, 1p,E11.4,A, 1p,E11.4,A, 1p,E11.4,A, 1p,E11.4)') 'step #', STEP, '  | # iters. ', cnt, '  | err V=', &
		real(MAX(err_u_max,err_v_max,err_w_max),4), '  | err p=', err_p_max, '  | err T=', err_T_max & 
         , '  | D=', Dilatation, '  | dudt=', dudt_max , '  | dTdt=', dTdt_max

!  print*, "step = ", STEP, "  |  ", "iters. = ", cnt
!  print*, "err V = ", real(MAX(err_u_max,err_v_max,err_w_max),4), "Max error p = ", err_p_max
!  print*, cnt
end subroutine integrate_nse


subroutine integrate_velocity_field
    implicit none
	integer   :: i
	if(nim_formulation==2) then
		!$OMP     PARALLEL DO  PRIVATE (i)
		do i=1,NMESH
			call integrate_velocity_field_formulation2_i(i)
		enddo
		!$OMP     END PARALLEL DO
	else
		!$OMP     PARALLEL DO  PRIVATE (i)
		do i=1,NMESH
			call integrate_velocity_field_i(i)
		enddo
		!$OMP     END PARALLEL DO
	endif
end subroutine integrate_velocity_field



subroutine update_velocity_sources
	integer i
	call calculate_gforces
     !$OMP     PARALLEL DO  PRIVATE (i)
     do i=1,NMESH

		call calculate_R_b(i)
		call calculate_R_c_new(i)
		call calculate_R_D2(i)
     enddo
     !$OMP     END PARALLEL DO
	 
end subroutine update_velocity_sources




subroutine integrate_pressure_field
	integer i
     !$OMP     PARALLEL DO  PRIVATE (i)
     do i=1,NMESH
        call integrate_pressure_field_i(i)
     enddo
     !$OMP     END PARALLEL DO
	 
end subroutine integrate_pressure_field



subroutine update_pressure_sources
	integer i
     !$OMP     PARALLEL DO  PRIVATE (i)
     do i=1,NMESH
		call calculate_R_b_p(i)
		call calculate_R_D_p2(i)
     enddo
     !$OMP     END PARALLEL DO
	 
end subroutine update_pressure_sources

subroutine integrate_temperature_field
	integer i
     !$OMP     PARALLEL DO  PRIVATE (i)
     do i=1,NMESH
        call integrate_temperature_field_formulation2_i(i)
     enddo
     !$OMP     END PARALLEL DO
	 
end subroutine integrate_temperature_field

subroutine update_temperature_sources
	integer i
     !$OMP     PARALLEL DO  PRIVATE (i)
     do i=1,NMESH
		call calculate_R_b_T(i)
		call calculate_R_D_T(i)
     enddo
     !$OMP     END PARALLEL DO
	 
end subroutine update_temperature_sources

subroutine find_error(err_u_max,err_v_max,err_w_max,err_p_max,err_T_max)
	real, intent(inout)         :: err_u_max,err_v_max,err_w_max,err_p_max,err_T_max
	integer                     :: i
	real                        :: error_max_dpdx, error_max_dpdy,error_max_dpdz,error_max_vxyz,error_max_uxyz,error_max_wxyz,error_max_Txyz
	err_u_max=0.
	err_v_max=0.
	err_w_max=0.
	err_p_max=0.
	err_T_max=0.
	error_max_dpdx=0.
	error_max_dpdy=0.
	error_max_dpdz=0.
	error_max_uxyz=0.
	error_max_vxyz=0.
	error_max_wxyz=0.
	error_max_Txyz=0.
	!$OMP     PARALLEL DO  PRIVATE (i) REDUCTION(max: err_u_max,err_v_max,err_w_max,err_p_max,err_T_max)
	do i=1,NFACE
		err_u_max=MAX(err_u_max,abs(faces_u(i)-faces_u_old(i)))
		err_v_max=MAX(err_v_max,abs(faces_v(i)-faces_v_old(i)))
		err_w_max=MAX(err_w_max,abs(faces_w(i)-faces_w_old(i)))
		err_p_max=MAX(err_p_max,abs(faces_p(i)-faces_p_old(i)))
		err_T_max=MAX(err_T_max,abs(faces_T(i)-faces_T_old(i)))
		
		faces_u_old(i)=faces_u(i)
		faces_v_old(i)=faces_v(i)
		faces_w_old(i)=faces_w(i)
		faces_p_old(i)=faces_p(i)
		faces_T_old(i)=faces_T(i)
	enddo
	!$OMP     END PARALLEL DO

	if(is_dpdx_tolerance) then
		!$OMP     PARALLEL DO  PRIVATE (i) REDUCTION(max: error_max_dpdx,error_max_dpdy,error_max_dpdz)
		do i=1,NMESH
			dpdx(i)=(faces_p(meshes_faces(2,i))-faces_p(meshes_faces(1,i)))/2.*meshes_xi_x(1,1,i)
			dpdy(i)=(faces_p(meshes_faces(4,i))-faces_p(meshes_faces(3,i)))/2.*meshes_xi_x(2,2,i)
			dpdz(i)=(faces_p(meshes_faces(6,i))-faces_p(meshes_faces(5,i)))/2.*meshes_xi_x(3,3,i)
			error_max_dpdx=MAX(error_max_dpdx,abs(dpdx(i)-dpdx_old(i)))
			error_max_dpdy=MAX(error_max_dpdy,abs(dpdy(i)-dpdy_old(i)))
			error_max_dpdz=MAX(error_max_dpdz,abs(dpdz(i)-dpdz_old(i)))
			dpdx_old(i)=dpdx(i)
			dpdy_old(i)=dpdy(i)
			dpdz_old(i)=dpdz(i)
		enddo
		!$OMP     END PARALLEL DO
		err_p_max=MAX(error_max_dpdx,error_max_dpdy,error_max_dpdz)
	endif
	
		!$OMP     PARALLEL DO  PRIVATE (i) REDUCTION(max: error_max_uxyz,error_max_vxyz,error_max_wxyz,error_max_Txyz)
		do i=1,NMESH
			error_max_uxyz=MAX(error_max_uxyz,abs(meshes_uxyz(i)-meshes_uxyz_old_con(i)))
			error_max_vxyz=MAX(error_max_vxyz,abs(meshes_vxyz(i)-meshes_vxyz_old_con(i)))
			error_max_wxyz=MAX(error_max_wxyz,abs(meshes_wxyz(i)-meshes_wxyz_old_con(i)))
			error_max_Txyz=MAX(error_max_Txyz,abs(meshes_Txyz(i)-meshes_Txyz_old_con(i)))
			meshes_uxyz_old_con(i)=meshes_uxyz(i)
			meshes_vxyz_old_con(i)=meshes_vxyz(i)
			meshes_wxyz_old_con(i)=meshes_wxyz(i)
			meshes_Txyz_old_con(i)=meshes_Txyz(i)
		enddo
		!$OMP     END PARALLEL DO
		err_u_max=MAX(err_u_max,error_max_uxyz)
		err_v_max=MAX(err_v_max,error_max_vxyz)
		err_w_max=MAX(err_w_max,error_max_wxyz)
		err_T_max=MAX(err_T_max,error_max_Txyz)
end subroutine find_error




subroutine integrate_velocity_field_i(n)
	integer, intent(in)                  :: n
	real                                 :: B1,B2,B3,B4,B_ij
	integer                              :: fnn,n1,n2
	integer, dimension(1:6)              :: n1_faces, n2_faces
    integer                              :: ni1,ni2,ni3,nj1,nj2,nj3,nk1,nk2,nk3
    real                                 :: u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z

	n1=n

! =========================== INTEGRATE U ============================================
	n1_faces(:)=meshes_faces(:,n1)
	B_ij=a_t(7,n1)*(meshes_R_C_u(n1)+meshes_R_D_u(n1)+meshes_R_B_u(n1))
	meshes_uxyz(n1)=(1-w_u)* meshes_uxyz(n1)+w_u*( &
	             -a_t(1,n1)* faces_u(n1_faces(1)) &
	             -a_t(2,n1)* faces_u(n1_faces(2)) &
				 -a_t(3,n1)* faces_u(n1_faces(3)) &
				 -a_t(4,n1)* faces_u(n1_faces(4)) &
				 -a_t(5,n1)* faces_u(n1_faces(5)) &
				 -a_t(6,n1)* faces_u(n1_faces(6)) &
				 -a_t(8,n1)* meshes_uxyz_old(n1) &
				 -B_ij)/(a_t(9,n1))

! =========================== INTEGRATE V ============================================

	B_ij=a_t(7,n1)*(meshes_R_C_v(n1)+meshes_R_D_v(n1)+meshes_R_B_v(n1))
	meshes_vxyz(n1)=(1-w_v)* meshes_vxyz(n1)+w_v*( &
	             -a_t(1,n1)* faces_v(n1_faces(1)) &
	             -a_t(2,n1)* faces_v(n1_faces(2)) &
				 -a_t(3,n1)* faces_v(n1_faces(3)) &
				 -a_t(4,n1)* faces_v(n1_faces(4)) &
				 -a_t(5,n1)* faces_v(n1_faces(5)) &
				 -a_t(6,n1)* faces_v(n1_faces(6)) &
				 -a_t(8,n1)* meshes_vxyz_old(n1) &
				 -B_ij)/(a_t(9,n1))

! =========================== INTEGRATE W ============================================

	B_ij=a_t(7,n1)*(meshes_R_C_w(n1)+meshes_R_D_w(n1)+meshes_R_B_w(n1))
	meshes_wxyz(n1)=(1-w_w)* meshes_wxyz(n1)+w_w*( &
	             -a_t(1,n1)* faces_w(n1_faces(1)) &
	             -a_t(2,n1)* faces_w(n1_faces(2)) &
				 -a_t(3,n1)* faces_w(n1_faces(3)) &
				 -a_t(4,n1)* faces_w(n1_faces(4)) &
				 -a_t(5,n1)* faces_w(n1_faces(5)) &
				 -a_t(6,n1)* faces_w(n1_faces(6)) &
				 -a_t(8,n1)* meshes_wxyz_old(n1) &
				 -B_ij)/(a_t(9,n1))

	!                               ======================
	! ============================== ||F||A||C||E|| ||2|| =============================
	!                               ======================

	if (faces_BC(n1_faces(2)) == BC_id(1)) THEN
	   n2=meshes_neighbors(2,n1)
	   n2_faces(:)=meshes_faces(:,n2)
	   
       nj1=meshes_neighbors(3,n1)
       nj2=meshes_neighbors(4,n1)
       nk1=meshes_neighbors(5,n1)
       nk2=meshes_neighbors(6,n1)


       if(nj1==0) nj3=meshes_neighbors(4,nj2)
       if(nj2==0) nj3=meshes_neighbors(3,nj1)
       if(nk1==0) nk3=meshes_neighbors(6,nk2)
       if(nk2==0) nk3=meshes_neighbors(5,nk1)
	   
! =========================== INTEGRATE U ============================================
	   B1=(meshes_R_C_u(n1)+meshes_R_D_u(n1)+meshes_R_B_u(n1))*AA_vx(13,n1)*2*DT
	   B2=(meshes_R_C_u(n2)+meshes_R_D_u(n2)+meshes_R_B_u(n2))*AA_vx(15,n1)*2*DT
	   B3= (meshes_gij_faces(1,2,1,n2)-meshes_gij_faces(1,2,2,n1))*(faces_u(n1_faces(4))-faces_u(n1_faces(3)) + &
										                            faces_u(n2_faces(4))-faces_u(n2_faces(3)))/4.0
	   B4= (meshes_gij_faces(1,3,1,n2)-meshes_gij_faces(1,3,2,n1))*(faces_u(n1_faces(6))-faces_u(n1_faces(5)) + & 
									                        	    faces_u(n2_faces(6))-faces_u(n2_faces(5)))/4.0

       if(nj1==0) then
         u_y=(-3*faces_u(meshes_faces(2,n1))+4*faces_u(meshes_faces(2,nj2))-faces_u(meshes_faces(2,nj3)))/4.
       elseif(nj2==0) then
         u_y=(3*faces_u(meshes_faces(2,n1))-4*faces_u(meshes_faces(2,nj1))+faces_u(meshes_faces(2,nj3)))/4.
       else
         u_y=(faces_u(meshes_faces(2,nj2))-faces_u(meshes_faces(2,nj1)))/4.
       endif

       if(nk1==0) then
         u_z=(-3*faces_u(meshes_faces(2,n1))+4*faces_u(meshes_faces(2,nk2))-faces_u(meshes_faces(2,nk3)))/4.
       elseif(nk2==0) then
         u_z=(3*faces_u(meshes_faces(2,n1))-4*faces_u(meshes_faces(2,nk1))+faces_u(meshes_faces(2,nk3)))/4.
       else
         u_z=(faces_u(meshes_faces(2,nk2))-faces_u(meshes_faces(2,nk1)))/4.
       endif


       B3= (meshes_gij_faces(1,2,1,n2)-meshes_gij_faces(1,2,2,n1))*u_y
       B4= (meshes_gij_faces(1,3,1,n2)-meshes_gij_faces(1,3,2,n1))*u_z


	   faces_u(n1_faces(2))=(1-w_u)*faces_u(n1_faces(2)) + w_u/(AA_vx(2,n1)) * ( &
											-AA_vx(1 ,n1)*faces_u(n1_faces(1)) &
											-AA_vx(3 ,n1)*faces_u(n2_faces(2)) &
											-AA_vx(4 ,n1)*faces_u(n1_faces(3)) &
											-AA_vx(5 ,n1)*faces_u(n1_faces(4)) &
											-AA_vx(6 ,n1)*faces_u(n2_faces(3)) &
											-AA_vx(7 ,n1)*faces_u(n2_faces(4)) &
											-AA_vx(8 ,n1)*faces_u(n1_faces(5)) &
											-AA_vx(9 ,n1)*faces_u(n1_faces(6)) &
											-AA_vx(10,n1)*faces_u(n2_faces(5)) &
											-AA_vx(11,n1)*faces_u(n2_faces(6)) &
											-AA_vx(12,n1)*(meshes_uxyz(n1)-meshes_uxyz_old(n1)) &
											-AA_vx(14,n1)*(meshes_uxyz(n2)-meshes_uxyz_old(n2)) &
											+B1+B2+B3+B4)
											
! =========================== INTEGRATE v ============================================
	   B1=(meshes_R_C_v(n1)+meshes_R_D_v(n1)+meshes_R_B_v(n1))*AA_vx(13,n1)*2*DT
	   B2=(meshes_R_C_v(n2)+meshes_R_D_v(n2)+meshes_R_B_v(n2))*AA_vx(15,n1)*2*DT
	   B3= (meshes_gij_faces(1,2,1,n2)-meshes_gij_faces(1,2,2,n1))*(faces_v(n1_faces(4))-faces_v(n1_faces(3)) + &
										                            faces_v(n2_faces(4))-faces_v(n2_faces(3)))/4.0
	   B4= (meshes_gij_faces(1,3,1,n2)-meshes_gij_faces(1,3,2,n1))*(faces_v(n1_faces(6))-faces_v(n1_faces(5)) + & 
										                            faces_v(n2_faces(6))-faces_v(n2_faces(5)))/4.0
       if(nj1==0) then
         v_y=(-3*faces_v(meshes_faces(2,n1))+4*faces_v(meshes_faces(2,nj2))-faces_v(meshes_faces(2,nj3)))/4.
       elseif(nj2==0) then
         v_y=(3*faces_v(meshes_faces(2,n1))-4*faces_v(meshes_faces(2,nj1))+faces_v(meshes_faces(2,nj3)))/4.
       else
         v_y=(faces_v(meshes_faces(2,nj2))-faces_v(meshes_faces(2,nj1)))/4.
       endif

       if(nk1==0) then
         v_z=(-3*faces_v(meshes_faces(2,n1))+4*faces_v(meshes_faces(2,nk2))-faces_v(meshes_faces(2,nk3)))/4.
       elseif(nk2==0) then
         v_z=(3*faces_v(meshes_faces(2,n1))-4*faces_v(meshes_faces(2,nk1))+faces_v(meshes_faces(2,nk3)))/4.
       else
         v_z=(faces_v(meshes_faces(2,nk2))-faces_v(meshes_faces(2,nk1)))/4.
       endif


       B3= (meshes_gij_faces(1,2,1,n2)-meshes_gij_faces(1,2,2,n1))*v_y
       B4= (meshes_gij_faces(1,3,1,n2)-meshes_gij_faces(1,3,2,n1))*v_z
	   faces_v(n1_faces(2))=(1-w_v)*faces_v(n1_faces(2)) + w_v/(AA_vx(2,n1)) * ( &
											-AA_vx(1 ,n1)*faces_v(n1_faces(1)) &
											-AA_vx(3 ,n1)*faces_v(n2_faces(2)) &
											-AA_vx(4 ,n1)*faces_v(n1_faces(3)) &
											-AA_vx(5 ,n1)*faces_v(n1_faces(4)) &
											-AA_vx(6 ,n1)*faces_v(n2_faces(3)) &
											-AA_vx(7 ,n1)*faces_v(n2_faces(4)) &
											-AA_vx(8 ,n1)*faces_v(n1_faces(5)) &
											-AA_vx(9 ,n1)*faces_v(n1_faces(6)) &
											-AA_vx(10,n1)*faces_v(n2_faces(5)) &
											-AA_vx(11,n1)*faces_v(n2_faces(6)) &
											-AA_vx(12,n1)*(meshes_vxyz(n1)-meshes_vxyz_old(n1)) &
											-AA_vx(14,n1)*(meshes_vxyz(n2)-meshes_vxyz_old(n2)) &
											+B1+B2+B3+B4)
											
! =========================== INTEGRATE W ============================================
	   B1=(meshes_R_C_w(n1)+meshes_R_D_w(n1)+meshes_R_B_w(n1))*AA_vx(13,n1)*2*DT
	   B2=(meshes_R_C_w(n2)+meshes_R_D_w(n2)+meshes_R_B_w(n2))*AA_vx(15,n1)*2*DT
	   B3= (meshes_gij_faces(1,2,1,n2)-meshes_gij_faces(1,2,2,n1))*(faces_w(n1_faces(4))-faces_w(n1_faces(3)) + &
										                            faces_w(n2_faces(4))-faces_w(n2_faces(3)))/4.0
	   B4= (meshes_gij_faces(1,3,1,n2)-meshes_gij_faces(1,3,2,n1))*(faces_w(n1_faces(6))-faces_w(n1_faces(5)) + & 
										                            faces_w(n2_faces(6))-faces_w(n2_faces(5)))/4.0


       if(nj1==0) then
         w_y=(-3*faces_w(meshes_faces(2,n1))+4*faces_w(meshes_faces(2,nj2))-faces_w(meshes_faces(2,nj3)))/4.
       elseif(nj2==0) then
         w_y=(3*faces_w(meshes_faces(2,n1))-4*faces_w(meshes_faces(2,nj1))+faces_w(meshes_faces(2,nj3)))/4.
       else
         w_y=(faces_w(meshes_faces(2,nj2))-faces_w(meshes_faces(2,nj1)))/4.
       endif

       if(nk1==0) then
         w_z=(-3*faces_w(meshes_faces(2,n1))+4*faces_w(meshes_faces(2,nk2))-faces_w(meshes_faces(2,nk3)))/4.
       elseif(nk2==0) then
         w_z=(3*faces_w(meshes_faces(2,n1))-4*faces_w(meshes_faces(2,nk1))+faces_w(meshes_faces(2,nk3)))/4.
       else
         w_z=(faces_w(meshes_faces(2,nk2))-faces_w(meshes_faces(2,nk1)))/4.
       endif


       B3= (meshes_gij_faces(1,2,1,n2)-meshes_gij_faces(1,2,2,n1))*w_y
       B4= (meshes_gij_faces(1,3,1,n2)-meshes_gij_faces(1,3,2,n1))*w_z

	   faces_w(n1_faces(2))=(1-w_w)*faces_w(n1_faces(2)) + w_w/(AA_vx(2,n1)) * ( &
											-AA_vx(1 ,n1)*faces_w(n1_faces(1)) &
											-AA_vx(3 ,n1)*faces_w(n2_faces(2)) &
											-AA_vx(4 ,n1)*faces_w(n1_faces(3)) &
											-AA_vx(5 ,n1)*faces_w(n1_faces(4)) &
											-AA_vx(6 ,n1)*faces_w(n2_faces(3)) &
											-AA_vx(7 ,n1)*faces_w(n2_faces(4)) &
											-AA_vx(8 ,n1)*faces_w(n1_faces(5)) &
											-AA_vx(9 ,n1)*faces_w(n1_faces(6)) &
											-AA_vx(10,n1)*faces_w(n2_faces(5)) &
											-AA_vx(11,n1)*faces_w(n2_faces(6)) &
											-AA_vx(12,n1)*(meshes_wxyz(n1)-meshes_wxyz_old(n1)) &
											-AA_vx(14,n1)*(meshes_wxyz(n2)-meshes_wxyz_old(n2)) &
											+B1+B2+B3+B4)
	endif
	!                               ======================
	! ============================== ||F||A||C||E|| ||4|| =============================
	!                               ======================
	if (faces_BC(n1_faces(4)) == BC_id(1)) THEN
	   n2=meshes_neighbors(4,n1)
	   n2_faces(:)=meshes_faces(:,n2)
       ni1=meshes_neighbors(1,n1)
       ni2=meshes_neighbors(2,n1)
       nk1=meshes_neighbors(5,n1)
       nk2=meshes_neighbors(6,n1)


       if(ni1==0) ni3=meshes_neighbors(2,ni2)
       if(ni2==0) ni3=meshes_neighbors(1,ni1)
       if(nk1==0) nk3=meshes_neighbors(6,nk2)
       if(nk2==0) nk3=meshes_neighbors(5,nk1)




! =========================== INTEGRATE U ============================================
	   B1=(meshes_R_C_u(n1)+meshes_R_D_u(n1)+meshes_R_B_u(n1))*AA_vy(13,n1)*2*DT
	   B2=(meshes_R_C_u(n2)+meshes_R_D_u(n2)+meshes_R_B_u(n2))*AA_vy(15,n1)*2*DT
	   B3= (meshes_gij_faces(2,1,3,n2)-meshes_gij_faces(2,1,4,n1))*(faces_u(n1_faces(2))-faces_u(n1_faces(1)) + &
										                            faces_u(n2_faces(2))-faces_u(n2_faces(1)))/4.0
	   B4= (meshes_gij_faces(2,3,3,n2)-meshes_gij_faces(2,3,4,n1))*(faces_u(n1_faces(6))-faces_u(n1_faces(5)) + &
										                            faces_u(n2_faces(6))-faces_u(n2_faces(5)))/4.0


       if(ni1==0) then
         u_x=(-3*faces_u(meshes_faces(4,n1))+4*faces_u(meshes_faces(4,ni2))-faces_u(meshes_faces(4,ni3)))/4.
       elseif(ni2==0) then
         u_x=(3*faces_u(meshes_faces(4,n1))-4*faces_u(meshes_faces(4,ni1))+faces_u(meshes_faces(4,ni3)))/4.
       else
         u_x=(faces_u(meshes_faces(4,ni2))-faces_u(meshes_faces(4,ni1)))/4.
       endif

       if(nk1==0) then
         u_z=(-3*faces_u(meshes_faces(4,n1))+4*faces_u(meshes_faces(4,nk2))-faces_u(meshes_faces(4,nk3)))/4.
       elseif(nk2==0) then
         u_z=(3*faces_u(meshes_faces(4,n1))-4*faces_u(meshes_faces(4,nk1))+faces_u(meshes_faces(4,nk3)))/4.
       else
         u_z=(faces_u(meshes_faces(4,nk2))-faces_u(meshes_faces(4,nk1)))/4.
       endif


       B3= (meshes_gij_faces(2,1,3,n2)-meshes_gij_faces(2,1,4,n1))*u_x
       B4= (meshes_gij_faces(2,3,3,n2)-meshes_gij_faces(2,3,4,n1))*u_z
	   faces_u(n1_faces(4))=(1-w_u)*faces_u(n1_faces(4)) + w_u/(AA_vy(6,n1)) * ( &
											-AA_vy(1 ,n1)*faces_u(n1_faces(1)) &
											-AA_vy(2 ,n1)*faces_u(n1_faces(2)) &
											-AA_vy(3 ,n1)*faces_u(n2_faces(1)) &
											-AA_vy(4 ,n1)*faces_u(n2_faces(2)) &
											-AA_vy(5 ,n1)*faces_u(n1_faces(3)) &
											-AA_vy(7 ,n1)*faces_u(n2_faces(4)) &
											-AA_vy(8 ,n1)*faces_u(n1_faces(5)) &
											-AA_vy(9 ,n1)*faces_u(n1_faces(6)) &
											-AA_vy(10,n1)*faces_u(n2_faces(5)) &
											-AA_vy(11,n1)*faces_u(n2_faces(6)) &
											-AA_vy(12,n1)*(meshes_uxyz(n1)-meshes_uxyz_old(n1)) &
											-AA_vy(14,n1)*(meshes_uxyz(n2)-meshes_uxyz_old(n2)) &
											+B1+B2+B3+B4)
											
! =========================== INTEGRATE V ============================================
	   B1=(meshes_R_C_v(n1)+meshes_R_D_v(n1)+meshes_R_B_v(n1))*AA_vy(13,n1)*2*DT
	   B2=(meshes_R_C_v(n2)+meshes_R_D_v(n2)+meshes_R_B_v(n2))*AA_vy(15,n1)*2*DT
	   B3= (meshes_gij_faces(2,1,3,n2)-meshes_gij_faces(2,1,4,n1))*(faces_v(n1_faces(2))-faces_v(n1_faces(1)) + &
										                            faces_v(n2_faces(2))-faces_v(n2_faces(1)))/4.0
	   B4= (meshes_gij_faces(2,3,3,n2)-meshes_gij_faces(2,3,4,n1))*(faces_v(n1_faces(6))-faces_v(n1_faces(5)) + &
										                            faces_v(n2_faces(6))-faces_v(n2_faces(5)))/4.0


       if(ni1==0) then
         v_x=(-3*faces_v(meshes_faces(4,n1))+4*faces_v(meshes_faces(4,ni2))-faces_v(meshes_faces(4,ni3)))/4.
       elseif(ni2==0) then
         v_x=(3*faces_v(meshes_faces(4,n1))-4*faces_v(meshes_faces(4,ni1))+faces_v(meshes_faces(4,ni3)))/4.
       else
         v_x=(faces_v(meshes_faces(4,ni2))-faces_v(meshes_faces(4,ni1)))/4.
       endif

       if(nk1==0) then
         v_z=(-3*faces_v(meshes_faces(4,n1))+4*faces_v(meshes_faces(4,nk2))-faces_v(meshes_faces(4,nk3)))/4.
       elseif(nk2==0) then
         v_z=(3*faces_v(meshes_faces(4,n1))-4*faces_v(meshes_faces(4,nk1))+faces_v(meshes_faces(4,nk3)))/4.
       else
         v_z=(faces_v(meshes_faces(4,nk2))-faces_v(meshes_faces(4,nk1)))/4.
       endif


       B3= (meshes_gij_faces(2,1,3,n2)-meshes_gij_faces(2,1,4,n1))*v_x
       B4= (meshes_gij_faces(2,3,3,n2)-meshes_gij_faces(2,3,4,n1))*v_z
	   faces_v(n1_faces(4))=(1-w_v)*faces_v(n1_faces(4)) + w_v/(AA_vy(6,n1)) * ( &
											-AA_vy(1 ,n1)*faces_v(n1_faces(1)) &
											-AA_vy(2 ,n1)*faces_v(n1_faces(2)) &
											-AA_vy(3 ,n1)*faces_v(n2_faces(1)) &
											-AA_vy(4 ,n1)*faces_v(n2_faces(2)) &
											-AA_vy(5 ,n1)*faces_v(n1_faces(3)) &
											-AA_vy(7 ,n1)*faces_v(n2_faces(4)) &
											-AA_vy(8 ,n1)*faces_v(n1_faces(5)) &
											-AA_vy(9 ,n1)*faces_v(n1_faces(6)) &
											-AA_vy(10,n1)*faces_v(n2_faces(5)) &
											-AA_vy(11,n1)*faces_v(n2_faces(6)) &
											-AA_vy(12,n1)*(meshes_vxyz(n1)-meshes_vxyz_old(n1)) &
											-AA_vy(14,n1)*(meshes_vxyz(n2)-meshes_vxyz_old(n2)) &
											+B1+B2+B3+B4)
											
! =========================== INTEGRATE W ============================================
	   B1=(meshes_R_C_w(n1)+meshes_R_D_w(n1)+meshes_R_B_w(n1))*AA_vy(13,n1)*2*DT
	   B2=(meshes_R_C_w(n2)+meshes_R_D_w(n2)+meshes_R_B_w(n2))*AA_vy(15,n1)*2*DT
	   B3= (meshes_gij_faces(2,1,3,n2)-meshes_gij_faces(2,1,4,n1))*(faces_w(n1_faces(2))-faces_w(n1_faces(1)) + &
										                            faces_w(n2_faces(2))-faces_w(n2_faces(1)))/4.0
	   B4= (meshes_gij_faces(2,3,3,n2)-meshes_gij_faces(2,3,4,n1))*(faces_w(n1_faces(6))-faces_w(n1_faces(5)) + &
										                            faces_w(n2_faces(6))-faces_w(n2_faces(5)))/4.0



       if(ni1==0) then
         w_x=(-3*faces_w(meshes_faces(4,n1))+4*faces_w(meshes_faces(4,ni2))-faces_w(meshes_faces(4,ni3)))/4.
       elseif(ni2==0) then
         w_x=(3*faces_w(meshes_faces(4,n1))-4*faces_w(meshes_faces(4,ni1))+faces_w(meshes_faces(4,ni3)))/4.
       else
         w_x=(faces_w(meshes_faces(4,ni2))-faces_w(meshes_faces(4,ni1)))/4.
       endif

       if(nk1==0) then
         w_z=(-3*faces_w(meshes_faces(4,n1))+4*faces_w(meshes_faces(4,nk2))-faces_w(meshes_faces(4,nk3)))/4.
       elseif(nk2==0) then
         w_z=(3*faces_w(meshes_faces(4,n1))-4*faces_w(meshes_faces(4,nk1))+faces_w(meshes_faces(4,nk3)))/4.
       else
         w_z=(faces_w(meshes_faces(4,nk2))-faces_w(meshes_faces(4,nk1)))/4.
       endif


       B3= (meshes_gij_faces(2,1,3,n2)-meshes_gij_faces(2,1,4,n1))*w_x
       B4= (meshes_gij_faces(2,3,3,n2)-meshes_gij_faces(2,3,4,n1))*w_z

	   faces_w(n1_faces(4))=(1-w_w)*faces_w(n1_faces(4)) + w_w/(AA_vy(6,n1)) * ( &
											-AA_vy(1 ,n1)*faces_w(n1_faces(1)) &
											-AA_vy(2 ,n1)*faces_w(n1_faces(2)) &
											-AA_vy(3 ,n1)*faces_w(n2_faces(1)) &
											-AA_vy(4 ,n1)*faces_w(n2_faces(2)) &
											-AA_vy(5 ,n1)*faces_w(n1_faces(3)) &
											-AA_vy(7 ,n1)*faces_w(n2_faces(4)) &
											-AA_vy(8 ,n1)*faces_w(n1_faces(5)) &
											-AA_vy(9 ,n1)*faces_w(n1_faces(6)) &
											-AA_vy(10,n1)*faces_w(n2_faces(5)) &
											-AA_vy(11,n1)*faces_w(n2_faces(6)) &
											-AA_vy(12,n1)*(meshes_wxyz(n1)-meshes_wxyz_old(n1)) &
											-AA_vy(14,n1)*(meshes_wxyz(n2)-meshes_wxyz_old(n2)) &
											+B1+B2+B3+B4)
	endif
	!                               ======================
	! ============================== ||F||A||C||E|| ||6|| =============================
	!                               ======================
	if (faces_BC(n1_faces(6)) == BC_id(1)) THEN
	   n2=meshes_neighbors(6,n1)
	   n2_faces(:)=meshes_faces(:,n2)
       ni1=meshes_neighbors(1,n1)
       ni2=meshes_neighbors(2,n1)
       nj1=meshes_neighbors(3,n1)
       nj2=meshes_neighbors(4,n1)


       if(ni1==0) ni3=meshes_neighbors(2,ni2)
       if(ni2==0) ni3=meshes_neighbors(1,ni1)
       if(nj1==0) nj3=meshes_neighbors(4,nj2)
       if(nj2==0) nj3=meshes_neighbors(3,nj1)


	   
	   
! =========================== INTEGRATE U ============================================
	   B1=(meshes_R_C_u(n1)+meshes_R_D_u(n1)+meshes_R_B_u(n1))*AA_vz(13,n1)*2*DT
	   B2=(meshes_R_C_u(n2)+meshes_R_D_u(n2)+meshes_R_B_u(n2))*AA_vz(15,n1)*2*DT
	   B3= (meshes_gij_faces(3,1,5,n2)-meshes_gij_faces(3,1,6,n1))*(faces_u(n1_faces(2))-faces_u(n1_faces(1)) + &
										                            faces_u(n2_faces(2))-faces_u(n2_faces(1)))/4.0
	   B4= (meshes_gij_faces(3,2,5,n2)-meshes_gij_faces(3,2,6,n1))*(faces_u(n1_faces(4))-faces_u(n1_faces(3)) + &
										                            faces_u(n2_faces(4))-faces_u(n2_faces(3)))/4.0

       if(ni1==0) then
         u_x=(-3*faces_u(meshes_faces(6,n1))+4*faces_u(meshes_faces(6,ni2))-faces_u(meshes_faces(6,ni3)))/4.
       elseif(ni2==0) then
         u_x=(3*faces_u(meshes_faces(6,n1))-4*faces_u(meshes_faces(6,ni1))+faces_u(meshes_faces(6,ni3)))/4.
       else
         u_x=(faces_u(meshes_faces(6,ni2))-faces_u(meshes_faces(6,ni1)))/4.
       endif

       if(nj1==0) then
         u_y=(-3*faces_u(meshes_faces(6,n1))+4*faces_u(meshes_faces(6,nj2))-faces_u(meshes_faces(6,nj3)))/4.
       elseif(nj2==0) then
         u_y=(3*faces_u(meshes_faces(6,n1))-4*faces_u(meshes_faces(6,nj1))+faces_u(meshes_faces(6,nj3)))/4.
       else
         u_y=(faces_u(meshes_faces(6,nj2))-faces_u(meshes_faces(6,nj1)))/4.
       endif


       B3= (meshes_gij_faces(3,1,5,n2)-meshes_gij_faces(3,1,6,n1))*u_x
       B4= (meshes_gij_faces(3,2,5,n2)-meshes_gij_faces(3,2,6,n1))*u_y

	   faces_u(n1_faces(6))=(1-w_u)*faces_u(n1_faces(6)) + w_u/(AA_vz(10,n1)) * ( &
											-AA_vz(1 ,n1)*faces_u(n1_faces(1)) &
											-AA_vz(2 ,n1)*faces_u(n1_faces(2)) &
											-AA_vz(3 ,n1)*faces_u(n2_faces(1)) &
											-AA_vz(4 ,n1)*faces_u(n2_faces(2)) &
											-AA_vz(5 ,n1)*faces_u(n1_faces(3)) &
											-AA_vz(6 ,n1)*faces_u(n1_faces(4)) &
											-AA_vz(7 ,n1)*faces_u(n2_faces(3)) &
											-AA_vz(8 ,n1)*faces_u(n2_faces(4)) &
											-AA_vz(9 ,n1)*faces_u(n1_faces(5)) &
											-AA_vz(11,n1)*faces_u(n2_faces(6)) &
											-AA_vz(12,n1)*(meshes_uxyz(n1)-meshes_uxyz_old(n1)) &
											-AA_vz(14,n1)*(meshes_uxyz(n2)-meshes_uxyz_old(n2)) &
											+B1+B2+B3+B4)
											
! =========================== INTEGRATE V ============================================
	   B1=(meshes_R_C_v(n1)+meshes_R_D_v(n1)+meshes_R_B_v(n1))*AA_vz(13,n1)*2*DT
	   B2=(meshes_R_C_v(n2)+meshes_R_D_v(n2)+meshes_R_B_v(n2))*AA_vz(15,n1)*2*DT
	   B3= (meshes_gij_faces(3,1,5,n2)-meshes_gij_faces(3,1,6,n1))*(faces_v(n1_faces(2))-faces_v(n1_faces(1)) + &
										                            faces_v(n2_faces(2))-faces_v(n2_faces(1)))/4.0
	   B4= (meshes_gij_faces(3,2,5,n2)-meshes_gij_faces(3,2,6,n1))*(faces_v(n1_faces(4))-faces_v(n1_faces(3)) + &
										                            faces_v(n2_faces(4))-faces_v(n2_faces(3)))/4.0

       if(ni1==0) then
         v_x=(-3*faces_v(meshes_faces(6,n1))+4*faces_v(meshes_faces(6,ni2))-faces_v(meshes_faces(6,ni3)))/4.
       elseif(ni2==0) then
         v_x=(3*faces_v(meshes_faces(6,n1))-4*faces_v(meshes_faces(6,ni1))+faces_v(meshes_faces(6,ni3)))/4.
       else
         v_x=(faces_v(meshes_faces(6,ni2))-faces_v(meshes_faces(6,ni1)))/4.
       endif

       if(nj1==0) then
         v_y=(-3*faces_v(meshes_faces(6,n1))+4*faces_v(meshes_faces(6,nj2))-faces_v(meshes_faces(6,nj3)))/4.
       elseif(nj2==0) then
         v_y=(3*faces_v(meshes_faces(6,n1))-4*faces_v(meshes_faces(6,nj1))+faces_v(meshes_faces(6,nj3)))/4.
       else
         v_y=(faces_v(meshes_faces(6,nj2))-faces_v(meshes_faces(6,nj1)))/4.
       endif


       B3= (meshes_gij_faces(3,1,5,n2)-meshes_gij_faces(3,1,6,n1))*v_x
       B4= (meshes_gij_faces(3,2,5,n2)-meshes_gij_faces(3,2,6,n1))*v_y

	   faces_v(n1_faces(6))=(1-w_v)*faces_v(n1_faces(6)) + w_v/(AA_vz(10,n1)) * ( &
											-AA_vz(1 ,n1)*faces_v(n1_faces(1)) &
											-AA_vz(2 ,n1)*faces_v(n1_faces(2)) &
											-AA_vz(3 ,n1)*faces_v(n2_faces(1)) &
											-AA_vz(4 ,n1)*faces_v(n2_faces(2)) &
											-AA_vz(5 ,n1)*faces_v(n1_faces(3)) &
											-AA_vz(6 ,n1)*faces_v(n1_faces(4)) &
											-AA_vz(7 ,n1)*faces_v(n2_faces(3)) &
											-AA_vz(8 ,n1)*faces_v(n2_faces(4)) &
											-AA_vz(9 ,n1)*faces_v(n1_faces(5)) &
											-AA_vz(11,n1)*faces_v(n2_faces(6)) &
											-AA_vz(12,n1)*(meshes_vxyz(n1)-meshes_vxyz_old(n1)) &
											-AA_vz(14,n1)*(meshes_vxyz(n2)-meshes_vxyz_old(n2)) &
											+B1+B2+B3+B4)
											
! =========================== INTEGRATE W ============================================
	   B1=(meshes_R_C_w(n1)+meshes_R_D_w(n1)+meshes_R_B_w(n1))*AA_vz(13,n1)*2*DT
	   B2=(meshes_R_C_w(n2)+meshes_R_D_w(n2)+meshes_R_B_w(n2))*AA_vz(15,n1)*2*DT
	   B3= (meshes_gij_faces(3,1,5,n2)-meshes_gij_faces(3,1,6,n1))*(faces_w(n1_faces(2))-faces_w(n1_faces(1)) + &
										                            faces_w(n2_faces(2))-faces_w(n2_faces(1)))/4.0
	   B4= (meshes_gij_faces(3,2,5,n2)-meshes_gij_faces(3,2,6,n1))*(faces_w(n1_faces(4))-faces_w(n1_faces(3)) + &
										                            faces_w(n2_faces(4))-faces_w(n2_faces(3)))/4.0


       if(ni1==0) then
         w_x=(-3*faces_w(meshes_faces(6,n1))+4*faces_w(meshes_faces(6,ni2))-faces_w(meshes_faces(6,ni3)))/4.
       elseif(ni2==0) then
         w_x=(3*faces_w(meshes_faces(6,n1))-4*faces_w(meshes_faces(6,ni1))+faces_w(meshes_faces(6,ni3)))/4.
       else
         w_x=(faces_w(meshes_faces(6,ni2))-faces_w(meshes_faces(6,ni1)))/4.
       endif

       if(nj1==0) then
         w_y=(-3*faces_w(meshes_faces(6,n1))+4*faces_w(meshes_faces(6,nj2))-faces_w(meshes_faces(6,nj3)))/4.
       elseif(nj2==0) then
         w_y=(3*faces_w(meshes_faces(6,n1))-4*faces_w(meshes_faces(6,nj1))+faces_w(meshes_faces(6,nj3)))/4.
       else
         w_y=(faces_w(meshes_faces(6,nj2))-faces_w(meshes_faces(6,nj1)))/4.
       endif


       B3= (meshes_gij_faces(3,1,5,n2)-meshes_gij_faces(3,1,6,n1))*w_x
       B4= (meshes_gij_faces(3,2,5,n2)-meshes_gij_faces(3,2,6,n1))*w_y

	   faces_w(n1_faces(6))=(1-w_w)*faces_w(n1_faces(6)) + w_w/(AA_vz(10,n1)) * ( &
											-AA_vz(1 ,n1)*faces_w(n1_faces(1)) &
											-AA_vz(2 ,n1)*faces_w(n1_faces(2)) &
											-AA_vz(3 ,n1)*faces_w(n2_faces(1)) &
											-AA_vz(4 ,n1)*faces_w(n2_faces(2)) &
											-AA_vz(5 ,n1)*faces_w(n1_faces(3)) &
											-AA_vz(6 ,n1)*faces_w(n1_faces(4)) &
											-AA_vz(7 ,n1)*faces_w(n2_faces(3)) &
											-AA_vz(8 ,n1)*faces_w(n2_faces(4)) &
											-AA_vz(9 ,n1)*faces_w(n1_faces(5)) &
											-AA_vz(11,n1)*faces_w(n2_faces(6)) &
											-AA_vz(12,n1)*(meshes_wxyz(n1)-meshes_wxyz_old(n1)) &
											-AA_vz(14,n1)*(meshes_wxyz(n2)-meshes_wxyz_old(n2)) &
											+B1+B2+B3+B4)
	endif
end subroutine integrate_velocity_field_i


subroutine integrate_velocity_field_formulation2_i(i)
	integer, intent(in)                  :: i
	real                                 :: B1,B2,B3,B4,B_ij
	integer                              :: fnn, i2
	integer, dimension(1:6)              :: i_faces, i2_faces
	i_faces(:)=meshes_faces(:,i)

! =========================== INTEGRATE U ============================================
	B_ij=(meshes_R_C_u(i)+meshes_R_D_u(i)+meshes_R_B_u(i))
	meshes_uxyz(i)=(1-w_u)*meshes_uxyz(i)+(w_u/M_F4(1,i))*( &
					   +M_F4(2,i)*meshes_uxyz_old(i) &
					   +M_F4(3,i)*faces_u(i_faces(1)) &
					   +M_F4(4,i)*faces_u(i_faces(2)) &
					   +M_F4(5,i)*faces_u(i_faces(3)) &
					   +M_F4(6,i)*faces_u(i_faces(4)) &
					   +M_F4(7,i)*faces_u(i_faces(5)) &
					   +M_F4(8,i)*faces_u(i_faces(6)) &
					   +B_ij)
! =========================== INTEGRATE V ============================================
								 
	B_ij=(meshes_R_C_v(i)+meshes_R_D_v(i)+meshes_R_B_v(i))
	meshes_vxyz(i)=(1-w_v)*meshes_vxyz(i)+(w_v/M_F4(1,i))*( &
					   +M_F4(2,i)*meshes_vxyz_old(i) &
					   +M_F4(3,i)*faces_v(i_faces(1)) &
					   +M_F4(4,i)*faces_v(i_faces(2)) &
					   +M_F4(5,i)*faces_v(i_faces(3)) &
					   +M_F4(6,i)*faces_v(i_faces(4)) &
					   +M_F4(7,i)*faces_v(i_faces(5)) &
					   +M_F4(8,i)*faces_v(i_faces(6)) &
					   +B_ij)
! =========================== INTEGRATE W ============================================
								 
	B_ij=(meshes_R_C_w(i)+meshes_R_D_w(i)+meshes_R_B_w(i))
	meshes_wxyz(i)=(1-w_w)*meshes_wxyz(i)+(w_w/M_F4(1,i))*( &
					   +M_F4(2,i)*meshes_wxyz_old(i) &
					   +M_F4(3,i)*faces_w(i_faces(1)) &
					   +M_F4(4,i)*faces_w(i_faces(2)) &
					   +M_F4(5,i)*faces_w(i_faces(3)) &
					   +M_F4(6,i)*faces_w(i_faces(4)) &
					   +M_F4(7,i)*faces_w(i_faces(5)) &
					   +M_F4(8,i)*faces_w(i_faces(6)) &
					   +B_ij)

	!                               ======================
	! ============================== ||F||A||C||E|| ||2|| =============================
	!                               ======================

	if (meshes_neighbors(2,i)/=0) THEN
		i2=meshes_neighbors(2,i)
	   i2_faces(:)=meshes_faces(:,i2)
! =========================== INTEGRATE U ============================================
	   B3= (meshes_gij_faces(1,2,1,i2)-meshes_gij_faces(1,2,2,i))*(faces_u(i_faces(4))-faces_u(i_faces(3)) + &
										                            faces_u(i2_faces(4))-faces_u(i2_faces(3)))/4.0
	   B4= (meshes_gij_faces(1,3,1,i2)-meshes_gij_faces(1,3,2,i))*(faces_u(i_faces(6))-faces_u(i_faces(5)) + & 
									                        	    faces_u(i2_faces(6))-faces_u(i2_faces(5)))/4.0
																	
		faces_u(i_faces(2))=(1-w_u)*faces_u(i_faces(2))+(w_u/M_F1(1,i)) * ( &
								 +M_F1(2,i) *  faces_u(i_faces(1)) &
								 +M_F1(3,i) *  faces_u(i2_faces(2)) &
								 +M_F1(4,i) * (meshes_uxyz(i)+meshes_uxyz_old(i))     &
								 +M_F1(5,i) * (meshes_uxyz(i2)+meshes_uxyz_old(i2)) &
								 -B3-B4)
! =========================== INTEGRATE V ============================================

	   B3= (meshes_gij_faces(1,2,1,i2)-meshes_gij_faces(1,2,2,i))*(faces_v(i_faces(4))-faces_v(i_faces(3)) + &
										                            faces_v(i2_faces(4))-faces_v(i2_faces(3)))/4.0
	   B4= (meshes_gij_faces(1,3,1,i2)-meshes_gij_faces(1,3,2,i))*(faces_v(i_faces(6))-faces_v(i_faces(5)) + & 
									                        	    faces_v(i2_faces(6))-faces_v(i2_faces(5)))/4.0
																	
		faces_v(i_faces(2))=(1-w_v)*faces_v(i_faces(2))+(w_v/M_F1(1,i)) * ( &
								 +M_F1(2,i) *  faces_v(i_faces(1)) &
								 +M_F1(3,i) *  faces_v(i2_faces(2)) &
								 +M_F1(4,i) * (meshes_vxyz(i)+meshes_vxyz_old(i))     &
								 +M_F1(5,i) * (meshes_vxyz(i2)+meshes_vxyz_old(i2)) &
								 -B3-B4)
! =========================== INTEGRATE W ============================================
								 
	   B3= (meshes_gij_faces(1,2,1,i2)-meshes_gij_faces(1,2,2,i))*(faces_w(i_faces(4))-faces_w(i_faces(3)) + &
										                            faces_w(i2_faces(4))-faces_w(i2_faces(3)))/4.0
	   B4= (meshes_gij_faces(1,3,1,i2)-meshes_gij_faces(1,3,2,i))*(faces_w(i_faces(6))-faces_w(i_faces(5)) + & 
									                        	    faces_w(i2_faces(6))-faces_w(i2_faces(5)))/4.0
																	
		faces_w(i_faces(2))=(1-w_w)*faces_w(i_faces(2))+(w_w/M_F1(1,i)) * ( &
								 +M_F1(2,i) *  faces_w(i_faces(1)) &
								 +M_F1(3,i) *  faces_w(i2_faces(2)) &
								 +M_F1(4,i) * (meshes_wxyz(i)+meshes_wxyz_old(i))     &
								 +M_F1(5,i) * (meshes_wxyz(i2)+meshes_wxyz_old(i2)) &
								 -B3-B4)
	endif
	!                               ======================
	! ============================== ||F||A||C||E|| ||4|| =============================
	!                               ======================
	if (meshes_neighbors(4,i)/=0) THEN
		i2=meshes_neighbors(4,i)
	   i2_faces(:)=meshes_faces(:,i2)
! =========================== INTEGRATE U ============================================
	   B3= (meshes_gij_faces(2,1,3,i2)-meshes_gij_faces(2,1,4,i))*(faces_u(i_faces(2))-faces_u(i_faces(1)) + &
										                            faces_u(i2_faces(2))-faces_u(i2_faces(1)))/4.0
	   B4= (meshes_gij_faces(2,3,3,i2)-meshes_gij_faces(2,3,4,i))*(faces_u(i_faces(6))-faces_u(i_faces(5)) + &
										                            faces_u(i2_faces(6))-faces_u(i2_faces(5)))/4.0
		faces_u(i_faces(4))=(1-w_u)*faces_u(i_faces(4))+(w_u/M_F2(1,i)) * ( &
								 +M_F2(2,i) *  faces_u(i_faces(3)) &
								 +M_F2(3,i) *  faces_u(i2_faces(4)) &
								 +M_F2(4,i) * (meshes_uxyz(i)+meshes_uxyz_old(i))     &
								 +M_F2(5,i) * (meshes_uxyz(i2)+meshes_uxyz_old(i2)) &
								 -B3-B4)

! =========================== INTEGRATE V ============================================
	   B3= (meshes_gij_faces(2,1,3,i2)-meshes_gij_faces(2,1,4,i))*(faces_v(i_faces(2))-faces_v(i_faces(1)) + &
										                            faces_v(i2_faces(2))-faces_v(i2_faces(1)))/4.0
	   B4= (meshes_gij_faces(2,3,3,i2)-meshes_gij_faces(2,3,4,i))*(faces_v(i_faces(6))-faces_v(i_faces(5)) + &
										                            faces_v(i2_faces(6))-faces_v(i2_faces(5)))/4.0
		faces_v(i_faces(4))=(1-w_v)*faces_v(i_faces(4))+(w_v/M_F2(1,i)) * ( &
								 +M_F2(2,i) *  faces_v(i_faces(3)) &
								 +M_F2(3,i) *  faces_v(i2_faces(4)) &
								 +M_F2(4,i) * (meshes_vxyz(i)+meshes_vxyz_old(i))     &
								 +M_F2(5,i) * (meshes_vxyz(i2)+meshes_vxyz_old(i2)) &
								 -B3-B4)

! =========================== INTEGRATE W ============================================
	   B3= (meshes_gij_faces(2,1,3,i2)-meshes_gij_faces(2,1,4,i))*(faces_w(i_faces(2))-faces_w(i_faces(1)) + &
										                            faces_w(i2_faces(2))-faces_w(i2_faces(1)))/4.0
	   B4= (meshes_gij_faces(2,3,3,i2)-meshes_gij_faces(2,3,4,i))*(faces_w(i_faces(6))-faces_w(i_faces(5)) + &
										                            faces_w(i2_faces(6))-faces_w(i2_faces(5)))/4.0
		faces_w(i_faces(4))=(1-w_w)*faces_w(i_faces(4))+(w_w/M_F2(1,i)) * ( &
								 +M_F2(2,i) *  faces_w(i_faces(3)) &
								 +M_F2(3,i) *  faces_w(i2_faces(4)) &
								 +M_F2(4,i) * (meshes_wxyz(i)+meshes_wxyz_old(i))     &
								 +M_F2(5,i) * (meshes_wxyz(i2)+meshes_wxyz_old(i2)) &
								 -B3-B4)

	endif
	!                               ======================
	! ============================== ||F||A||C||E|| ||6|| =============================
	!                               ======================
	
	if (meshes_neighbors(6,i)/=0) THEN
		i2=meshes_neighbors(6,i)
	   i2_faces(:)=meshes_faces(:,i2)
! =========================== INTEGRATE U ============================================	   
	   B3= (meshes_gij_faces(3,1,5,i2)-meshes_gij_faces(3,1,6,i))*(faces_u(i_faces(2))-faces_u(i_faces(1)) + &
										                            faces_u(i2_faces(2))-faces_u(i2_faces(1)))/4.0
	   B4= (meshes_gij_faces(3,2,5,i2)-meshes_gij_faces(3,2,6,i))*(faces_u(i_faces(4))-faces_u(i_faces(3)) + &
										                            faces_u(i2_faces(4))-faces_u(i2_faces(3)))/4.0
		faces_u(i_faces(6))=(1-w_u)*faces_u(i_faces(6))+(w_u/M_F3(1,i)) * ( &
								 +M_F3(2,i) *  faces_u(i_faces(5)) &
								 +M_F3(3,i) *  faces_u(i2_faces(6)) &
								 +M_F3(4,i) * (meshes_uxyz(i)+meshes_uxyz_old(i))     &
								 +M_F3(5,i) * (meshes_uxyz(i2)+meshes_uxyz_old(i2)) &
								 -B3-B4)	   
	   
! =========================== INTEGRATE V ============================================	   
	   B3= (meshes_gij_faces(3,1,5,i2)-meshes_gij_faces(3,1,6,i))*(faces_v(i_faces(2))-faces_v(i_faces(1)) + &
										                            faces_v(i2_faces(2))-faces_v(i2_faces(1)))/4.0
	   B4= (meshes_gij_faces(3,2,5,i2)-meshes_gij_faces(3,2,6,i))*(faces_v(i_faces(4))-faces_v(i_faces(3)) + &
										                            faces_v(i2_faces(4))-faces_v(i2_faces(3)))/4.0
		faces_v(i_faces(6))=(1-w_v)*faces_v(i_faces(6))+(w_v/M_F3(1,i)) * ( &
								 +M_F3(2,i) *  faces_v(i_faces(5)) &
								 +M_F3(3,i) *  faces_v(i2_faces(6)) &
								 +M_F3(4,i) * (meshes_vxyz(i)+meshes_vxyz_old(i))     &
								 +M_F3(5,i) * (meshes_vxyz(i2)+meshes_vxyz_old(i2)) &
								 -B3-B4)	   
! =========================== INTEGRATE W ============================================	   
	   B3= (meshes_gij_faces(3,1,5,i2)-meshes_gij_faces(3,1,6,i))*(faces_w(i_faces(2))-faces_w(i_faces(1)) + &
										                            faces_w(i2_faces(2))-faces_w(i2_faces(1)))/4.0
	   B4= (meshes_gij_faces(3,2,5,i2)-meshes_gij_faces(3,2,6,i))*(faces_w(i_faces(4))-faces_w(i_faces(3)) + &
										                            faces_w(i2_faces(4))-faces_w(i2_faces(3)))/4.0
		faces_w(i_faces(6))=(1-w_w)*faces_w(i_faces(6))+(w_w/M_F3(1,i)) * ( &
								 +M_F3(2,i) *  faces_w(i_faces(5)) &
								 +M_F3(3,i) *  faces_w(i2_faces(6)) &
								 +M_F3(4,i) * (meshes_wxyz(i)+meshes_wxyz_old(i))     &
								 +M_F3(5,i) * (meshes_wxyz(i2)+meshes_wxyz_old(i2)) &
								 -B3-B4)		   	   
	   
	endif	
end subroutine integrate_velocity_field_formulation2_i


subroutine integrate_temperature_field_formulation2_i(i)
	integer, intent(in)                  :: i
	real                                 :: B1,B2,B3,B4,B_ij
	integer                              :: fnn, i2
	integer, dimension(1:6)              :: i_faces, i2_faces
	i_faces(:)=meshes_faces(:,i)

! =========================== INTEGRATE T ============================================
	B_ij=(meshes_R_D_T(i)+meshes_R_B_T(i))
	meshes_Txyz(i)=(1-w_T)*meshes_Txyz(i)+(w_T/MT_F4(1,i))*( &
					   +MT_F4(2,i)*meshes_Txyz_old(i) &
					   +MT_F4(3,i)*faces_T(i_faces(1)) &
					   +MT_F4(4,i)*faces_T(i_faces(2)) &
					   +MT_F4(5,i)*faces_T(i_faces(3)) &
					   +MT_F4(6,i)*faces_T(i_faces(4)) &
					   +MT_F4(7,i)*faces_T(i_faces(5)) &
					   +MT_F4(8,i)*faces_T(i_faces(6)) &
					   +B_ij)


	!                               ======================
	! ============================== ||F||A||C||E|| ||2|| =============================
	!                               ======================

	if (meshes_neighbors(2,i)/=0) THEN
		i2=meshes_neighbors(2,i)
	   i2_faces(:)=meshes_faces(:,i2)
! =========================== INTEGRATE U ============================================
	   B3= (meshes_gij_faces(1,2,1,i2)-meshes_gij_faces(1,2,2,i))*(faces_T(i_faces(4))-faces_T(i_faces(3)) + &
										                            faces_T(i2_faces(4))-faces_T(i2_faces(3)))/4.0
	   B4= (meshes_gij_faces(1,3,1,i2)-meshes_gij_faces(1,3,2,i))*(faces_T(i_faces(6))-faces_T(i_faces(5)) + & 
									                        	    faces_T(i2_faces(6))-faces_T(i2_faces(5)))/4.0
																	
		faces_T(i_faces(2))=(1-w_T)*faces_T(i_faces(2))+(w_T/MT_F1(1,i)) * ( &
								 +MT_F1(2,i) *  faces_T(i_faces(1)) &
								 +MT_F1(3,i) *  faces_T(i2_faces(2)) &
								 +MT_F1(4,i) * (meshes_Txyz(i)+meshes_Txyz_old(i))     &
								 +MT_F1(5,i) * (meshes_Txyz(i2)+meshes_Txyz_old(i2)) &
								 -B3-B4)
	endif
	!                               ======================
	! ============================== ||F||A||C||E|| ||4|| =============================
	!                               ======================
	if (meshes_neighbors(4,i)/=0) THEN
		i2=meshes_neighbors(4,i)
	   i2_faces(:)=meshes_faces(:,i2)
! =========================== INTEGRATE U ============================================
	   B3= (meshes_gij_faces(2,1,3,i2)-meshes_gij_faces(2,1,4,i))*(faces_T(i_faces(2))-faces_T(i_faces(1)) + &
										                            faces_T(i2_faces(2))-faces_T(i2_faces(1)))/4.0
	   B4= (meshes_gij_faces(2,3,3,i2)-meshes_gij_faces(2,3,4,i))*(faces_T(i_faces(6))-faces_T(i_faces(5)) + &
										                            faces_T(i2_faces(6))-faces_T(i2_faces(5)))/4.0
		faces_T(i_faces(4))=(1-w_T)*faces_T(i_faces(4))+(w_T/MT_F2(1,i)) * ( &
								 +MT_F2(2,i) *  faces_T(i_faces(3)) &
								 +MT_F2(3,i) *  faces_T(i2_faces(4)) &
								 +MT_F2(4,i) * (meshes_Txyz(i)+meshes_Txyz_old(i))     &
								 +MT_F2(5,i) * (meshes_Txyz(i2)+meshes_Txyz_old(i2)) &
								 -B3-B4)


	endif
	!                               ======================
	! ============================== ||F||A||C||E|| ||6|| =============================
	!                               ======================
	
	if (meshes_neighbors(6,i)/=0) THEN
		i2=meshes_neighbors(6,i)
	   i2_faces(:)=meshes_faces(:,i2)
! =========================== INTEGRATE U ============================================	   
	   B3= (meshes_gij_faces(3,1,5,i2)-meshes_gij_faces(3,1,6,i))*(faces_T(i_faces(2))-faces_T(i_faces(1)) + &
										                            faces_T(i2_faces(2))-faces_T(i2_faces(1)))/4.0
	   B4= (meshes_gij_faces(3,2,5,i2)-meshes_gij_faces(3,2,6,i))*(faces_T(i_faces(4))-faces_T(i_faces(3)) + &
										                            faces_T(i2_faces(4))-faces_T(i2_faces(3)))/4.0
		faces_T(i_faces(6))=(1-w_T)*faces_T(i_faces(6))+(w_T/MT_F3(1,i)) * ( &
								 +MT_F3(2,i) *  faces_T(i_faces(5)) &
								 +MT_F3(3,i) *  faces_T(i2_faces(6)) &
								 +MT_F3(4,i) * (meshes_Txyz(i)+meshes_Txyz_old(i))     &
								 +MT_F3(5,i) * (meshes_Txyz(i2)+meshes_Txyz_old(i2)) &
								 -B3-B4)	   
	   
	endif	
end subroutine integrate_temperature_field_formulation2_i

subroutine integrate_pressure_field_i(n)
	integer, intent(in)                  :: n
	real                                 :: B1,B2,B3,B4,B_ij
	integer                              :: fnn,n1,n2
	integer, dimension(1:6)              :: n1_faces, n2_faces
    integer                              :: ni1,ni2,ni3,nj1,nj2,nj3,nk1,nk2,nk3
    real                                 :: p_x,p_y,p_z
	n1=n

	n1_faces(:)=meshes_faces(:,n1)

	!                               ======================
	! ============================== ||F||A||C||E|| ||2|| =============================
	!                               ======================
	if (faces_BC(n1_faces(2)) == BC_id(1)) THEN
	   n2=meshes_neighbors(2,n1)
	   n2_faces(:)=meshes_faces(:,n2)

       nj1=meshes_neighbors(3,n1)
       nj2=meshes_neighbors(4,n1)
       nk1=meshes_neighbors(5,n1)
       nk2=meshes_neighbors(6,n1)


       if(nj1==0) nj3=meshes_neighbors(4,nj2)
       if(nj2==0) nj3=meshes_neighbors(3,nj1)
       if(nk1==0) nk3=meshes_neighbors(6,nk2)
       if(nk2==0) nk3=meshes_neighbors(5,nk1)



	   B1=(meshes_R_D_p(n1)+meshes_R_B_p(n1))*AA_px(13,n1)
	   B2=(meshes_R_D_p(n2)+meshes_R_B_p(n2))*AA_px(15,n1)
	   B3= (meshes_gij_faces(1,2,1,n2)-meshes_gij_faces(1,2,2,n1))*(faces_p(n1_faces(4))-faces_p(n1_faces(3)) + &
										                            faces_p(n2_faces(4))-faces_p(n2_faces(3)))/4.0
	   B4= (meshes_gij_faces(1,3,1,n2)-meshes_gij_faces(1,3,2,n1))*(faces_p(n1_faces(6))-faces_p(n1_faces(5)) + & 
									                        	    faces_p(n2_faces(6))-faces_p(n2_faces(5)))/4.0


       if(nj1==0) then
         p_y=(-3*faces_p(meshes_faces(2,n1))+4*faces_p(meshes_faces(2,nj2))-faces_p(meshes_faces(2,nj3)))/4.
       elseif(nj2==0) then
         p_y=(3*faces_p(meshes_faces(2,n1))-4*faces_p(meshes_faces(2,nj1))+faces_p(meshes_faces(2,nj3)))/4.
       else
         p_y=(faces_p(meshes_faces(2,nj2))-faces_p(meshes_faces(2,nj1)))/4.
       endif

       if(nk1==0) then
         p_z=(-3*faces_p(meshes_faces(2,n1))+4*faces_p(meshes_faces(2,nk2))-faces_p(meshes_faces(2,nk3)))/4.
       elseif(nk2==0) then
         p_z=(3*faces_p(meshes_faces(2,n1))-4*faces_p(meshes_faces(2,nk1))+faces_p(meshes_faces(2,nk3)))/4.
       else
         p_z=(faces_p(meshes_faces(2,nk2))-faces_p(meshes_faces(2,nk1)))/4.
       endif


       B3= (meshes_gij_faces(1,2,1,n2)-meshes_gij_faces(1,2,2,n1))*p_y
       B4= (meshes_gij_faces(1,3,1,n2)-meshes_gij_faces(1,3,2,n1))*p_z

	   faces_p(n1_faces(2))=(1-w_p)*faces_p(n1_faces(2)) + w_p/(AA_px(2,n1)) * ( &
											-AA_px(1 ,n1)*faces_p(n1_faces(1)) &
											-AA_px(3 ,n1)*faces_p(n2_faces(2)) &
											-AA_px(4 ,n1)*faces_p(n1_faces(3)) &
											-AA_px(5 ,n1)*faces_p(n1_faces(4)) &
											-AA_px(6 ,n1)*faces_p(n2_faces(3)) &
											-AA_px(7 ,n1)*faces_p(n2_faces(4)) &
											-AA_px(8 ,n1)*faces_p(n1_faces(5)) &
											-AA_px(9 ,n1)*faces_p(n1_faces(6)) &
											-AA_px(10,n1)*faces_p(n2_faces(5)) &
											-AA_px(11,n1)*faces_p(n2_faces(6)) &
											+B1+B2+B3+B4)
!		faces_p(n1_faces(2))=faces_pT(n1_faces(2))

	endif
	!                               ======================
	! ============================== ||F||A||C||E|| ||4|| =============================
	!                               ======================
	if (faces_BC(n1_faces(4)) == BC_id(1)) THEN
	   n2=meshes_neighbors(4,n1)
	   n2_faces(:)=meshes_faces(:,n2)
       ni1=meshes_neighbors(1,n1)
       ni2=meshes_neighbors(2,n1)
       nk1=meshes_neighbors(5,n1)
       nk2=meshes_neighbors(6,n1)


       if(ni1==0) ni3=meshes_neighbors(2,ni2)
       if(ni2==0) ni3=meshes_neighbors(1,ni1)
       if(nk1==0) nk3=meshes_neighbors(6,nk2)
       if(nk2==0) nk3=meshes_neighbors(5,nk1)




	   B1=(meshes_R_D_p(n1)+meshes_R_B_p(n1))*AA_py(13,n1)
	   B2=(meshes_R_D_p(n2)+meshes_R_B_p(n2))*AA_py(15,n1)
	   B3= (meshes_gij_faces(2,1,3,n2)-meshes_gij_faces(2,1,4,n1))*(faces_p(n1_faces(2))-faces_p(n1_faces(1)) + &
										                            faces_p(n2_faces(2))-faces_p(n2_faces(1)))/4.0
	   B4= (meshes_gij_faces(2,3,3,n2)-meshes_gij_faces(2,3,4,n1))*(faces_p(n1_faces(6))-faces_p(n1_faces(5)) + &
										                            faces_p(n2_faces(6))-faces_p(n2_faces(5)))/4.0

       if(ni1==0) then
         p_x=(-3*faces_p(meshes_faces(4,n1))+4*faces_p(meshes_faces(4,ni2))-faces_p(meshes_faces(4,ni3)))/4.
       elseif(ni2==0) then
         p_x=(3*faces_p(meshes_faces(4,n1))-4*faces_p(meshes_faces(4,ni1))+faces_p(meshes_faces(4,ni3)))/4.
       else
         p_x=(faces_p(meshes_faces(4,ni2))-faces_p(meshes_faces(4,ni1)))/4.
       endif

       if(nk1==0) then
         p_z=(-3*faces_p(meshes_faces(4,n1))+4*faces_p(meshes_faces(4,nk2))-faces_p(meshes_faces(4,nk3)))/4.
       elseif(nk2==0) then
         p_z=(3*faces_p(meshes_faces(4,n1))-4*faces_p(meshes_faces(4,nk1))+faces_p(meshes_faces(4,nk3)))/4.
       else
         p_z=(faces_p(meshes_faces(4,nk2))-faces_p(meshes_faces(4,nk1)))/4.
       endif

       B3= (meshes_gij_faces(2,1,3,n2)-meshes_gij_faces(2,1,4,n1))*p_x
       B4= (meshes_gij_faces(2,3,3,n2)-meshes_gij_faces(2,3,4,n1))*p_z

	   faces_p(n1_faces(4))=(1-w_p)*faces_p(n1_faces(4)) + w_p/(AA_py(6,n1)) * ( &
											-AA_py(1 ,n1)*faces_p(n1_faces(1)) &
											-AA_py(2 ,n1)*faces_p(n1_faces(2)) &
											-AA_py(3 ,n1)*faces_p(n2_faces(1)) &
											-AA_py(4 ,n1)*faces_p(n2_faces(2)) &
											-AA_py(5 ,n1)*faces_p(n1_faces(3)) &
											-AA_py(7 ,n1)*faces_p(n2_faces(4)) &
											-AA_py(8 ,n1)*faces_p(n1_faces(5)) &
											-AA_py(9 ,n1)*faces_p(n1_faces(6)) &
											-AA_py(10,n1)*faces_p(n2_faces(5)) &
											-AA_py(11,n1)*faces_p(n2_faces(6)) &
											+B1+B2+B3+B4)

!		faces_p(n1_faces(4))=faces_pT(n1_faces(4))
	endif
	!                               ======================
	! ============================== ||F||A||C||E|| ||6|| =============================
	!                               ======================
	if (faces_BC(n1_faces(6)) == BC_id(1)) THEN
	   n2=meshes_neighbors(6,n1)
	   n2_faces(:)=meshes_faces(:,n2)
       ni1=meshes_neighbors(1,n1)
       ni2=meshes_neighbors(2,n1)
       nj1=meshes_neighbors(3,n1)
       nj2=meshes_neighbors(4,n1)


       if(ni1==0) ni3=meshes_neighbors(2,ni2)
       if(ni2==0) ni3=meshes_neighbors(1,ni1)
       if(nj1==0) nj3=meshes_neighbors(4,nj2)
       if(nj2==0) nj3=meshes_neighbors(3,nj1)


	   
	   B1=(meshes_R_D_p(n1)+meshes_R_B_p(n1))*AA_pz(13,n1)
	   B2=(meshes_R_D_p(n2)+meshes_R_B_p(n2))*AA_pz(15,n1)
	   B3= (meshes_gij_faces(3,1,5,n2)-meshes_gij_faces(3,1,6,n1))*(faces_p(n1_faces(2))-faces_p(n1_faces(1)) + &
										                            faces_p(n2_faces(2))-faces_p(n2_faces(1)))/4.0
	   B4= (meshes_gij_faces(3,2,5,n2)-meshes_gij_faces(3,2,6,n1))*(faces_p(n1_faces(4))-faces_p(n1_faces(3)) + &
										                            faces_p(n2_faces(4))-faces_p(n2_faces(3)))/4.0

	   
       if(ni1==0) then
         p_x=(-3*faces_p(meshes_faces(6,n1))+4*faces_p(meshes_faces(6,ni2))-faces_p(meshes_faces(6,ni3)))/4.
       elseif(ni2==0) then
         p_x=(3*faces_p(meshes_faces(6,n1))-4*faces_p(meshes_faces(6,ni1))+faces_p(meshes_faces(6,ni3)))/4.
       else
         p_x=(faces_p(meshes_faces(6,ni2))-faces_p(meshes_faces(6,ni1)))/4.
       endif

       if(nj1==0) then
         p_y=(-3*faces_p(meshes_faces(6,n1))+4*faces_p(meshes_faces(6,nj2))-faces_p(meshes_faces(6,nj3)))/4.
       elseif(nj2==0) then
         p_y=(3*faces_p(meshes_faces(6,n1))-4*faces_p(meshes_faces(6,nj1))+faces_p(meshes_faces(6,nj3)))/4.
       else
         p_y=(faces_p(meshes_faces(6,nj2))-faces_p(meshes_faces(6,nj1)))/4.
       endif


       B3= (meshes_gij_faces(3,1,5,n2)-meshes_gij_faces(3,1,6,n1))*p_x
       B4= (meshes_gij_faces(3,2,5,n2)-meshes_gij_faces(3,2,6,n1))*p_y


	   faces_p(n1_faces(6))=(1-w_p)*faces_p(n1_faces(6)) + w_p/(AA_pz(10,n1)) * ( &
											-AA_pz(1 ,n1)*faces_p(n1_faces(1)) &
											-AA_pz(2 ,n1)*faces_p(n1_faces(2)) &
											-AA_pz(3 ,n1)*faces_p(n2_faces(1)) &
											-AA_pz(4 ,n1)*faces_p(n2_faces(2)) &
											-AA_pz(5 ,n1)*faces_p(n1_faces(3)) &
											-AA_pz(6 ,n1)*faces_p(n1_faces(4)) &
											-AA_pz(7 ,n1)*faces_p(n2_faces(3)) &
											-AA_pz(8 ,n1)*faces_p(n2_faces(4)) &
											-AA_pz(9 ,n1)*faces_p(n1_faces(5)) &
											-AA_pz(11,n1)*faces_p(n2_faces(6)) &
											+B1+B2+B3+B4)

!		faces_p(n1_faces(6))=faces_pT(n1_faces(6))
	endif
end subroutine integrate_pressure_field_i


subroutine update_time_step
	integer i
! UPDATE TIME STEP
	STEP=STEP+1
    TI=TI+2*DT
!$OMP     PARALLEL DO PRIVATE (i)
	do i=1,NMESH
	meshes_uxyzt(i)=(meshes_uxyz_old(i)+meshes_uxyz(i))/2.
	meshes_uxyz_old(i)=meshes_uxyz(i)
	meshes_vxyzt(i)=(meshes_vxyz_old(i)+meshes_vxyz(i))/2.
	meshes_vxyz_old(i)=meshes_vxyz(i)
	meshes_wxyzt(i)=(meshes_wxyz_old(i)+meshes_wxyz(i))/2.
	meshes_wxyz_old(i)=meshes_wxyz(i)
	meshes_Txyzt(i)=(meshes_Txyz_old(i)+meshes_Txyz(i))/2.
	meshes_Txyz_old(i)=meshes_Txyz(i)
	meshes_p(i)=(faces_p(meshes_faces(1,i))+faces_p(meshes_faces(2,i))+ faces_p(meshes_faces(3,i)) + &
			     faces_p(meshes_faces(4,i))+faces_p(meshes_faces(5,i))+ faces_p(meshes_faces(6,i)))/6.0
	meshes_dila_old(i)=meshes_dila(i)
	enddo
!$OMP     END PARALLEL DO

end subroutine update_time_step


subroutine check_steady_state
	integer                     :: i
	real                        :: max_dudt, max_dvdt, max_dwdt, max_dTdt
	max_dudt=0.
	max_dvdt=0.
	max_dwdt=0.
	!$OMP     PARALLEL DO  PRIVATE (i) REDUCTION(max: max_dudt,max_dvdt,max_dwdt, max_dTdt)
	do i=1,NMESH
		max_dudt=MAX(max_dudt,abs(meshes_uxyz(i)-meshes_uxyz_old(i))/(2*DT))
		max_dvdt=MAX(max_dvdt,abs(meshes_vxyz(i)-meshes_vxyz_old(i))/(2*DT))
		max_dwdt=MAX(max_dwdt,abs(meshes_wxyz(i)-meshes_wxyz_old(i))/(2*DT))
		max_dTdt=MAX(max_dTdt,abs(meshes_Txyz(i)-meshes_Txyz_old(i))/(2*DT))
	enddo
	!$OMP     END PARALLEL DO
	dudt_max=MAX(max_dwdt,MAX(max_dudt,max_dvdt))
	dTdt_max=max_dTdt
end subroutine check_steady_state


subroutine calculate_R_c(n)
	integer, intent(in)         :: n
	real                        :: du_dxi1,du_dxi2,du_dxi3,dv_dxi1,dv_dxi2,dv_dxi3,dw_dxi1,dw_dxi2,dw_dxi3
    integer, dimension(1:6)     :: n_faces
    real, dimension(1:6)        :: uu_faces,uv_faces,uw_faces,vv_faces,vw_faces,ww_faces
	integer                     :: k
	
	
	n_faces(:)=meshes_faces(:,n)
	do k=1,6
		uu_faces(k)=faces_u(n_faces(k))*faces_u(n_faces(k))
		uv_faces(k)=faces_u(n_faces(k))*faces_v(n_faces(k))
		uw_faces(k)=faces_u(n_faces(k))*faces_w(n_faces(k))
		vv_faces(k)=faces_v(n_faces(k))*faces_v(n_faces(k))
		vw_faces(k)=faces_v(n_faces(k))*faces_w(n_faces(k))
		ww_faces(k)=faces_w(n_faces(k))*faces_w(n_faces(k))
	enddo 		
	
	
	n_faces(:)=meshes_faces(:,n)
! ==============================================================================================
    du_dxi1=(faces_u(n_faces(2))-faces_u(n_faces(1)))/2.
    du_dxi2=(faces_u(n_faces(4))-faces_u(n_faces(3)))/2.
    du_dxi3=(faces_u(n_faces(6))-faces_u(n_faces(5)))/2.

    dv_dxi1=(faces_v(n_faces(2))-faces_v(n_faces(1)))/2.
    dv_dxi2=(faces_v(n_faces(4))-faces_v(n_faces(3)))/2.
    dv_dxi3=(faces_v(n_faces(6))-faces_v(n_faces(5)))/2.
	
    dw_dxi1=(faces_w(n_faces(2))-faces_w(n_faces(1)))/2.
    dw_dxi2=(faces_w(n_faces(4))-faces_w(n_faces(3)))/2.
    dw_dxi3=(faces_w(n_faces(6))-faces_w(n_faces(5)))/2.

! ==============================================================================================

! ==============================================================================================


!	meshes_R_C_u(n)=-(meshes_u(n)*du_dxi1 +meshes_v(n)*du_dxi2 +meshes_w(n)*du_dxi3 )+ &
!					 (dphi_dx2(uu_faces,n_faces,n)+ &
!					  dphi_dy2(uv_faces,n_faces,n)+ &
!					  dphi_dz2(uw_faces,n_faces,n))

	meshes_R_C_u(n)=(dphi_dx2(uu_faces,n_faces,n)-meshes_u(n)*du_dxi1) &
				   +(dphi_dy2(uv_faces,n_faces,n)-meshes_v(n)*du_dxi2) &
				   +(dphi_dz2(uw_faces,n_faces,n)-meshes_w(n)*du_dxi3)
				   
	meshes_R_C_v(n)=(dphi_dx2(uv_faces,n_faces,n)-meshes_u(n)*dv_dxi1) & 
				   +(dphi_dy2(vv_faces,n_faces,n)-meshes_v(n)*dv_dxi2) &
				   +(dphi_dz2(vw_faces,n_faces,n)-meshes_w(n)*dv_dxi3)

	meshes_R_C_w(n)=(dphi_dx2(uw_faces,n_faces,n)-meshes_u(n)*dw_dxi1) &
				   +(dphi_dy2(vw_faces,n_faces,n)-meshes_v(n)*dw_dxi2) &
				   +(dphi_dz2(ww_faces,n_faces,n)-meshes_w(n)*dw_dxi3)

end subroutine calculate_R_c


subroutine calculate_R_c_new(n)
	integer, intent(in)         :: n
	real                        :: du_dxi1,du_dxi2,du_dxi3,dv_dxi1,dv_dxi2,dv_dxi3,dw_dxi1,dw_dxi2,dw_dxi3, uc,vc,wc
    integer, dimension(1:6)     :: n_faces
    real, dimension(1:6)        :: u_faces,v_faces,w_faces
	integer                     :: k
	
	
	n_faces(:)=meshes_faces(:,n)
	do k=1,6
		u_faces(k)=faces_u(n_faces(k))
		v_faces(k)=faces_v(n_faces(k))
		w_faces(k)=faces_w(n_faces(k))
	enddo 		
	
	
	n_faces(:)=meshes_faces(:,n)
! ==============================================================================================
    du_dxi1=(faces_u(n_faces(2))-faces_u(n_faces(1)))/2.
    du_dxi2=(faces_u(n_faces(4))-faces_u(n_faces(3)))/2.
    du_dxi3=(faces_u(n_faces(6))-faces_u(n_faces(5)))/2.

    dv_dxi1=(faces_v(n_faces(2))-faces_v(n_faces(1)))/2.
    dv_dxi2=(faces_v(n_faces(4))-faces_v(n_faces(3)))/2.
    dv_dxi3=(faces_v(n_faces(6))-faces_v(n_faces(5)))/2.
	
    dw_dxi1=(faces_w(n_faces(2))-faces_w(n_faces(1)))/2.
    dw_dxi2=(faces_w(n_faces(4))-faces_w(n_faces(3)))/2.
    dw_dxi3=(faces_w(n_faces(6))-faces_w(n_faces(5)))/2.

! ==============================================================================================

! ==============================================================================================
		uc=meshes_uxyzt(n)*meshes_xi_x(1,1,n)+meshes_vxyzt(n)*meshes_xi_x(1,2,n)+meshes_wxyzt(n)*meshes_xi_x(1,3,n)
		vc=meshes_uxyzt(n)*meshes_xi_x(2,1,n)+meshes_vxyzt(n)*meshes_xi_x(2,2,n)+meshes_wxyzt(n)*meshes_xi_x(2,3,n)
		wc=meshes_uxyzt(n)*meshes_xi_x(3,1,n)+meshes_vxyzt(n)*meshes_xi_x(3,2,n)+meshes_wxyzt(n)*meshes_xi_x(3,3,n)

	meshes_R_C_u(n)=-((meshes_u(n)-uc)*du_dxi1 +(meshes_v(n)-vc)*du_dxi2+(meshes_w(n)-wc)*du_dxi3)
	meshes_R_C_v(n)=-((meshes_u(n)-uc)*dv_dxi1 +(meshes_v(n)-vc)*dv_dxi2+(meshes_w(n)-wc)*dv_dxi3)
	meshes_R_C_w(n)=-((meshes_u(n)-uc)*dw_dxi1 +(meshes_v(n)-vc)*dw_dxi2+(meshes_w(n)-wc)*dw_dxi3)

end subroutine calculate_R_c_new



subroutine calculate_R_c2(n)
	integer, intent(in)         :: n
	real                        :: du_dxi1,du_dxi2,du_dxi3,dv_dxi1,dv_dxi2,dv_dxi3,dw_dxi1,dw_dxi2,dw_dxi3, u1,v1,w1
    integer, dimension(1:6)     :: n_faces
    real, dimension(1:6)        :: uu_faces,uv_faces,uw_faces,vv_faces,vw_faces,ww_faces
	integer                     :: k
	

	n_faces(:)=meshes_faces(:,n)
! ==============================================================================================
    du_dxi1=(faces_u(n_faces(2))-faces_u(n_faces(1)))/2.
    du_dxi2=(faces_u(n_faces(4))-faces_u(n_faces(3)))/2.
    du_dxi3=(faces_u(n_faces(6))-faces_u(n_faces(5)))/2.

    dv_dxi1=(faces_v(n_faces(2))-faces_v(n_faces(1)))/2.
    dv_dxi2=(faces_v(n_faces(4))-faces_v(n_faces(3)))/2.
    dv_dxi3=(faces_v(n_faces(6))-faces_v(n_faces(5)))/2.
	
    dw_dxi1=(faces_w(n_faces(2))-faces_w(n_faces(1)))/2.
    dw_dxi2=(faces_w(n_faces(4))-faces_w(n_faces(3)))/2.
    dw_dxi3=(faces_w(n_faces(6))-faces_w(n_faces(5)))/2.

! ==============================================================================================

! ==============================================================================================
    u1=1./6.*(faces_u(n_faces(1))+faces_u(n_faces(2))+faces_u(n_faces(3))+faces_u(n_faces(4))+faces_u(n_faces(5))+faces_u(n_faces(6)))*meshes_xi_x(1,1,n)
    v1=1./6.*(faces_v(n_faces(1))+faces_v(n_faces(2))+faces_v(n_faces(3))+faces_v(n_faces(4))+faces_v(n_faces(5))+faces_v(n_faces(6)))*meshes_xi_x(2,2,n)
    w1=1./6.*(faces_w(n_faces(1))+faces_w(n_faces(2))+faces_w(n_faces(3))+faces_w(n_faces(4))+faces_w(n_faces(5))+faces_w(n_faces(6)))*meshes_xi_x(3,3,n)

	meshes_R_C_u(n)=-((meshes_u(n)-u1)*du_dxi1 +(meshes_v(n)-v1)*du_dxi2 +(meshes_w(n)-w1)*du_dxi3 )

	meshes_R_C_v(n)=-((meshes_u(n)-u1)*dv_dxi1 +(meshes_v(n)-v1)*dv_dxi2 +(meshes_w(n)-w1)*dv_dxi3 )
	meshes_R_C_w(n)=-((meshes_u(n)-u1)*dw_dxi1 +(meshes_v(n)-v1)*dw_dxi2 +(meshes_w(n)-w1)*dw_dxi3 )

end subroutine calculate_R_c2

subroutine calculate_R_c3(n)
	integer, intent(in)         :: n
	real                        :: du_dxi1,du_dxi2,du_dxi3,dv_dxi1,dv_dxi2,dv_dxi3,dw_dxi1,dw_dxi2,dw_dxi3
    real                        :: duu_dxi1,duu_dxi2,duu_dxi3,duv_dxi1,duv_dxi2,duv_dxi3,duw_dxi1,duw_dxi2,duw_dxi3
	real                        :: dvv_dxi1,dvv_dxi2,dvv_dxi3,dvw_dxi1,dvw_dxi2,dvw_dxi3,dww_dxi1,dww_dxi2,dww_dxi3
    integer, dimension(1:6)     :: n_faces
	
	
	n_faces(:)=meshes_faces(:,n)
! ==============================================================================================
    du_dxi1=(faces_u(n_faces(2))-faces_u(n_faces(1)))/2.
    du_dxi2=(faces_u(n_faces(4))-faces_u(n_faces(3)))/2.
    du_dxi3=(faces_u(n_faces(6))-faces_u(n_faces(5)))/2.

    dv_dxi1=(faces_v(n_faces(2))-faces_v(n_faces(1)))/2.
    dv_dxi2=(faces_v(n_faces(4))-faces_v(n_faces(3)))/2.
    dv_dxi3=(faces_v(n_faces(6))-faces_v(n_faces(5)))/2.
	
    dw_dxi1=(faces_w(n_faces(2))-faces_w(n_faces(1)))/2.
    dw_dxi2=(faces_w(n_faces(4))-faces_w(n_faces(3)))/2.
    dw_dxi3=(faces_w(n_faces(6))-faces_w(n_faces(5)))/2.

! ==============================================================================================

	duu_dxi1=(faces_u(n_faces(2))*faces_u(n_faces(2))- &
			  faces_u(n_faces(1))*faces_u(n_faces(1)))/2.
	duu_dxi2=(faces_u(n_faces(4))*faces_u(n_faces(4))- &
			  faces_u(n_faces(3))*faces_u(n_faces(3)))/2.			  
	duu_dxi3=(faces_u(n_faces(6))*faces_u(n_faces(6))- &
			  faces_u(n_faces(5))*faces_u(n_faces(5)))/2.			  

	duv_dxi1=(faces_u(n_faces(2))*faces_v(n_faces(2))- &
			  faces_u(n_faces(1))*faces_v(n_faces(1)))/2.
	duv_dxi2=(faces_u(n_faces(4))*faces_v(n_faces(4))- &
			  faces_u(n_faces(3))*faces_v(n_faces(3)))/2.			  
	duv_dxi3=(faces_u(n_faces(6))*faces_v(n_faces(6))- &
			  faces_u(n_faces(5))*faces_v(n_faces(5)))/2.	

	duw_dxi1=(faces_u(n_faces(2))*faces_w(n_faces(2))- &
			  faces_u(n_faces(1))*faces_w(n_faces(1)))/2.
	duw_dxi2=(faces_u(n_faces(4))*faces_w(n_faces(4))- &
			  faces_u(n_faces(3))*faces_w(n_faces(3)))/2.			  
	duw_dxi3=(faces_u(n_faces(6))*faces_w(n_faces(6))- &
			  faces_u(n_faces(5))*faces_w(n_faces(5)))/2.	

	dvv_dxi1=(faces_v(n_faces(2))*faces_v(n_faces(2))- &
			  faces_v(n_faces(1))*faces_v(n_faces(1)))/2.
	dvv_dxi2=(faces_v(n_faces(4))*faces_v(n_faces(4))- &
			  faces_v(n_faces(3))*faces_v(n_faces(3)))/2.			  
	dvv_dxi3=(faces_v(n_faces(6))*faces_v(n_faces(6))- &
			  faces_v(n_faces(5))*faces_v(n_faces(5)))/2.	

	dvw_dxi1=(faces_v(n_faces(2))*faces_w(n_faces(2))- &
			  faces_v(n_faces(1))*faces_w(n_faces(1)))/2.
	dvw_dxi2=(faces_v(n_faces(4))*faces_w(n_faces(4))- &
			  faces_v(n_faces(3))*faces_w(n_faces(3)))/2.			  
	dvw_dxi3=(faces_v(n_faces(6))*faces_w(n_faces(6))- &
			  faces_v(n_faces(5))*faces_w(n_faces(5)))/2.	
	
	dww_dxi1=(faces_w(n_faces(2))*faces_w(n_faces(2))- &
			  faces_w(n_faces(1))*faces_w(n_faces(1)))/2.
	dww_dxi2=(faces_w(n_faces(4))*faces_w(n_faces(4))- &
			  faces_w(n_faces(3))*faces_w(n_faces(3)))/2.			  
	dww_dxi3=(faces_w(n_faces(6))*faces_w(n_faces(6))- &
			  faces_w(n_faces(5))*faces_w(n_faces(5)))/2.	

! ==============================================================================================


	meshes_R_C_u(n)=-(meshes_u(n)*du_dxi1 +meshes_v(n)*du_dxi2 +meshes_w(n)*du_dxi3 )+ &
					 dphi_dx(duu_dxi1,duu_dxi2,duu_dxi3,n)+ &
					 dphi_dy(duv_dxi1,duv_dxi2,duv_dxi3,n)+ &
					 dphi_dz(duw_dxi1,duw_dxi2,duw_dxi3,n)

	meshes_R_C_v(n)=-(meshes_u(n)*dv_dxi1 +meshes_v(n)*dv_dxi2 +meshes_w(n)*dv_dxi3 )+ &
					 dphi_dx(duv_dxi1,duv_dxi2,duv_dxi3,n)+ &
					 dphi_dy(dvv_dxi1,dvv_dxi2,dvv_dxi3,n)+ &
					 dphi_dz(dvw_dxi1,dvw_dxi2,dvw_dxi3,n)

	meshes_R_C_w(n)=-(meshes_u(n)*dw_dxi1 +meshes_v(n)*dw_dxi2 +meshes_w(n)*dw_dxi3 )+ &
					 dphi_dx(duw_dxi1,duw_dxi2,duw_dxi3,n)+ &
					 dphi_dy(dvw_dxi1,dvw_dxi2,dvw_dxi3,n)+ &
					 dphi_dz(dww_dxi1,dww_dxi2,dww_dxi3,n)


end subroutine calculate_R_c3


subroutine calculate_R_b(n)
	integer, intent(in)         :: n
	real                        :: dp_dx, dp_dy, dp_dz
    integer, dimension(1:6)     :: n_faces
    real, dimension(1:6)        :: p_faces
	integer                     :: k
	
	
	n_faces(:)=meshes_faces(:,n)
	do k=1,6
		p_faces(k)=faces_p(n_faces(k))
	enddo 		

! ==============================================================================================
    ! dp_dxi1=(faces_p(n_faces(2))-faces_p(n_faces(1)))/2.
    ! dp_dxi2=(faces_p(n_faces(4))-faces_p(n_faces(3)))/2.
    ! dp_dxi3=(faces_p(n_faces(6))-faces_p(n_faces(5)))/2.

    dp_dx= dphi_dx2(p_faces,n_faces,n)
    dp_dy= dphi_dy2(p_faces,n_faces,n)
	dp_dz= dphi_dz2(p_faces,n_faces,n)

! ==============================================================================================


	meshes_R_B_u(n)=meshes_bx(n)+1/rho*dp_dx
	meshes_R_B_v(n)=meshes_by(n)+1/rho*dp_dy
	meshes_R_B_w(n)=meshes_bz(n)+1/rho*dp_dz


end subroutine calculate_R_b


subroutine calculate_R_b2(n)
	integer, intent(in)         :: n
	real                        :: dp_dx, dp_dy, dp_dz,dp_dxi1,dp_dxi2,dp_dxi3
    integer, dimension(1:6)     :: n_faces
    real, dimension(1:6)        :: p_faces
	integer                     :: k
	
	
	n_faces(:)=meshes_faces(:,n)
	do k=1,6
		p_faces(k)=faces_p(n_faces(k))
	enddo 		

! ==============================================================================================
    dp_dxi1=(p_faces(2)-p_faces(1))/2.
    dp_dxi2=(p_faces(4)-p_faces(3))/2.
    dp_dxi3=(p_faces(6)-p_faces(5))/2.

    dp_dx= dphi_dx(dp_dxi1,dp_dxi2,dp_dxi3,n) 
    dp_dy= dphi_dy(dp_dxi1,dp_dxi2,dp_dxi3,n) 
	dp_dz= dphi_dz(dp_dxi1,dp_dxi2,dp_dxi3,n) 

! ==============================================================================================


	meshes_R_B_u(n)=meshes_bx(n)+1/rho*dp_dx
	meshes_R_B_v(n)=meshes_by(n)+1/rho*dp_dy
	meshes_R_B_w(n)=meshes_bz(n)+1/rho*dp_dz


end subroutine calculate_R_b2


subroutine calculate_R_D2(n)
	integer, intent(in)         :: n
    real                        :: uxy,uxz,uyz,ux,uy,uz,ux1,ux2,uy1,uy2,uz1,uz2,uyx,uzx,uzy	
    real                        :: vxy,vxz,vyz,vx,vy,vz,vx1,vx2,vy1,vy2,vz1,vz2,vyx,vzx,vzy	
    real                        :: wxy,wxz,wyz,wx,wy,wz,wx1,wx2,wy1,wy2,wz1,wz2,wyx,wzx,wzy	
    integer                     :: i,nj1,nj2,nj3,ni1,ni2,ni3,nk1,nk2,nk3,n1
	
    nj1=meshes_neighbors(3,n)
    nj2=meshes_neighbors(4,n)

    if(nj1==0) nj3=meshes_neighbors(4,nj2)
    if(nj2==0) nj3=meshes_neighbors(3,nj1)
	
    if(nj1==0) then
      uy1=(-3*faces_u(meshes_faces(1,n))+4*faces_u(meshes_faces(1,nj2))-faces_u(meshes_faces(1,nj3)))/4.
	  vy1=(-3*faces_v(meshes_faces(1,n))+4*faces_v(meshes_faces(1,nj2))-faces_v(meshes_faces(1,nj3)))/4.
	  wy1=(-3*faces_w(meshes_faces(1,n))+4*faces_w(meshes_faces(1,nj2))-faces_w(meshes_faces(1,nj3)))/4.
    elseif(nj2==0) then
      uy1=(3*faces_u(meshes_faces(1,n))-4*faces_u(meshes_faces(1,nj1))+faces_u(meshes_faces(1,nj3)))/4.
      vy1=(3*faces_v(meshes_faces(1,n))-4*faces_v(meshes_faces(1,nj1))+faces_v(meshes_faces(1,nj3)))/4.
      wy1=(3*faces_w(meshes_faces(1,n))-4*faces_w(meshes_faces(1,nj1))+faces_w(meshes_faces(1,nj3)))/4.
    else
      uy1=(faces_u(meshes_faces(1,nj2))-faces_u(meshes_faces(1,nj1)))/4.
      vy1=(faces_v(meshes_faces(1,nj2))-faces_v(meshes_faces(1,nj1)))/4.
      wy1=(faces_w(meshes_faces(1,nj2))-faces_w(meshes_faces(1,nj1)))/4.
    endif
    if(nj1==0) then
      uy2=(-3*faces_u(meshes_faces(2,n))+4*faces_u(meshes_faces(2,nj2))-faces_u(meshes_faces(2,nj3)))/4.
	  vy2=(-3*faces_v(meshes_faces(2,n))+4*faces_v(meshes_faces(2,nj2))-faces_v(meshes_faces(2,nj3)))/4.
	  wy2=(-3*faces_w(meshes_faces(2,n))+4*faces_w(meshes_faces(2,nj2))-faces_w(meshes_faces(2,nj3)))/4.
    elseif(nj2==0) then
      uy2=(3*faces_u(meshes_faces(2,n))-4*faces_u(meshes_faces(2,nj1))+faces_u(meshes_faces(2,nj3)))/4.
      vy2=(3*faces_v(meshes_faces(2,n))-4*faces_v(meshes_faces(2,nj1))+faces_v(meshes_faces(2,nj3)))/4.
      wy2=(3*faces_w(meshes_faces(2,n))-4*faces_w(meshes_faces(2,nj1))+faces_w(meshes_faces(2,nj3)))/4.
    else
      uy2=(faces_u(meshes_faces(2,nj2))-faces_u(meshes_faces(2,nj1)))/4.
      vy2=(faces_v(meshes_faces(2,nj2))-faces_v(meshes_faces(2,nj1)))/4.
      wy2=(faces_w(meshes_faces(2,nj2))-faces_w(meshes_faces(2,nj1)))/4.
    endif

    uxy=0.5*(uy2-uy1)
	vxy=0.5*(vy2-vy1)
	wxy=0.5*(wy2-wy1)
	
	
    if(nj1==0) then
      uy1=(-3*faces_u(meshes_faces(5,n))+4*faces_u(meshes_faces(5,nj2))-faces_u(meshes_faces(5,nj3)))/4.
	  vy1=(-3*faces_v(meshes_faces(5,n))+4*faces_v(meshes_faces(5,nj2))-faces_v(meshes_faces(5,nj3)))/4.
	  wy1=(-3*faces_w(meshes_faces(5,n))+4*faces_w(meshes_faces(5,nj2))-faces_w(meshes_faces(5,nj3)))/4.
    elseif(nj2==0) then
      uy1=(3*faces_u(meshes_faces(5,n))-4*faces_u(meshes_faces(5,nj1))+faces_u(meshes_faces(5,nj3)))/4.
      vy1=(3*faces_v(meshes_faces(5,n))-4*faces_v(meshes_faces(5,nj1))+faces_v(meshes_faces(5,nj3)))/4.
      wy1=(3*faces_w(meshes_faces(5,n))-4*faces_w(meshes_faces(5,nj1))+faces_w(meshes_faces(5,nj3)))/4.
    else
      uy1=(faces_u(meshes_faces(5,nj2))-faces_u(meshes_faces(5,nj1)))/4.
      vy1=(faces_v(meshes_faces(5,nj2))-faces_v(meshes_faces(5,nj1)))/4.
      wy1=(faces_w(meshes_faces(5,nj2))-faces_w(meshes_faces(5,nj1)))/4.
    endif
    if(nj1==0) then
      uy2=(-3*faces_u(meshes_faces(6,n))+4*faces_u(meshes_faces(6,nj2))-faces_u(meshes_faces(6,nj3)))/4.
	  vy2=(-3*faces_v(meshes_faces(6,n))+4*faces_v(meshes_faces(6,nj2))-faces_v(meshes_faces(6,nj3)))/4.
	  wy2=(-3*faces_w(meshes_faces(6,n))+4*faces_w(meshes_faces(6,nj2))-faces_w(meshes_faces(6,nj3)))/4.
    elseif(nj2==0) then
      uy2=(3*faces_u(meshes_faces(6,n))-4*faces_u(meshes_faces(6,nj1))+faces_u(meshes_faces(6,nj3)))/4.
      vy2=(3*faces_v(meshes_faces(6,n))-4*faces_v(meshes_faces(6,nj1))+faces_v(meshes_faces(6,nj3)))/4.
      wy2=(3*faces_w(meshes_faces(6,n))-4*faces_w(meshes_faces(6,nj1))+faces_w(meshes_faces(6,nj3)))/4.
    else
      uy2=(faces_u(meshes_faces(6,nj2))-faces_u(meshes_faces(6,nj1)))/4.
      vy2=(faces_v(meshes_faces(6,nj2))-faces_v(meshes_faces(6,nj1)))/4.
      wy2=(faces_w(meshes_faces(6,nj2))-faces_w(meshes_faces(6,nj1)))/4.
    endif
    uzy=0.5*(uy2-uy1)	
	vzy=0.5*(vy2-vy1)
	wzy=0.5*(wy2-wy1)
	

    ni1=meshes_neighbors(1,n)
    ni2=meshes_neighbors(2,n)

    if(ni1==0) ni3=meshes_neighbors(2,ni2)
    if(ni2==0) ni3=meshes_neighbors(1,ni1)


    if(ni1==0) then
      ux1=(-3*faces_u(meshes_faces(3,n))+4*faces_u(meshes_faces(3,ni2))-faces_u(meshes_faces(3,ni3)))/4.
	  vx1=(-3*faces_v(meshes_faces(3,n))+4*faces_v(meshes_faces(3,ni2))-faces_v(meshes_faces(3,ni3)))/4.
	  wx1=(-3*faces_w(meshes_faces(3,n))+4*faces_w(meshes_faces(3,ni2))-faces_w(meshes_faces(3,ni3)))/4.
    elseif(ni2==0) then
      ux1=(3*faces_u(meshes_faces(3,n))-4*faces_u(meshes_faces(3,ni1))+faces_u(meshes_faces(3,ni3)))/4.
      vx1=(3*faces_v(meshes_faces(3,n))-4*faces_v(meshes_faces(3,ni1))+faces_v(meshes_faces(3,ni3)))/4.
      wx1=(3*faces_w(meshes_faces(3,n))-4*faces_w(meshes_faces(3,ni1))+faces_w(meshes_faces(3,ni3)))/4.
    else
      ux1=(faces_u(meshes_faces(3,ni2))-faces_u(meshes_faces(3,ni1)))/4.
      vx1=(faces_v(meshes_faces(3,ni2))-faces_v(meshes_faces(3,ni1)))/4.
      wx1=(faces_w(meshes_faces(3,ni2))-faces_w(meshes_faces(3,ni1)))/4.
    endif
    if(ni1==0) then
      ux2=(-3*faces_u(meshes_faces(4,n))+4*faces_u(meshes_faces(4,ni2))-faces_u(meshes_faces(4,ni3)))/4.
	  vx2=(-3*faces_v(meshes_faces(4,n))+4*faces_v(meshes_faces(4,ni2))-faces_v(meshes_faces(4,ni3)))/4.
	  wx2=(-3*faces_w(meshes_faces(4,n))+4*faces_w(meshes_faces(4,ni2))-faces_w(meshes_faces(4,ni3)))/4.
    elseif(ni2==0) then
      ux2=(3*faces_u(meshes_faces(4,n))-4*faces_u(meshes_faces(4,ni1))+faces_u(meshes_faces(4,ni3)))/4.
      vx2=(3*faces_v(meshes_faces(4,n))-4*faces_v(meshes_faces(4,ni1))+faces_v(meshes_faces(4,ni3)))/4.
      wx2=(3*faces_w(meshes_faces(4,n))-4*faces_w(meshes_faces(4,ni1))+faces_w(meshes_faces(4,ni3)))/4.
    else
      ux2=(faces_u(meshes_faces(4,ni2))-faces_u(meshes_faces(4,ni1)))/4.
      vx2=(faces_v(meshes_faces(4,ni2))-faces_v(meshes_faces(4,ni1)))/4.
      wx2=(faces_w(meshes_faces(4,ni2))-faces_w(meshes_faces(4,ni1)))/4.
    endif

    uyx=0.5*(ux2-ux1)
	vyx=0.5*(vx2-vx1)
	wyx=0.5*(wx2-wx1)


    if(ni1==0) then
      ux1=(-3*faces_u(meshes_faces(5,n))+4*faces_u(meshes_faces(5,ni2))-faces_u(meshes_faces(5,ni3)))/4.
	  vx1=(-3*faces_v(meshes_faces(5,n))+4*faces_v(meshes_faces(5,ni2))-faces_v(meshes_faces(5,ni3)))/4.
	  wx1=(-3*faces_w(meshes_faces(5,n))+4*faces_w(meshes_faces(5,ni2))-faces_w(meshes_faces(5,ni3)))/4.
    elseif(ni2==0) then
      ux1=(3*faces_u(meshes_faces(5,n))-4*faces_u(meshes_faces(5,ni1))+faces_u(meshes_faces(5,ni3)))/4.
      vx1=(3*faces_v(meshes_faces(5,n))-4*faces_v(meshes_faces(5,ni1))+faces_v(meshes_faces(5,ni3)))/4.
      wx1=(3*faces_w(meshes_faces(5,n))-4*faces_w(meshes_faces(5,ni1))+faces_w(meshes_faces(5,ni3)))/4.
    else
      ux1=(faces_u(meshes_faces(5,ni2))-faces_u(meshes_faces(5,ni1)))/4.
      vx1=(faces_v(meshes_faces(5,ni2))-faces_v(meshes_faces(5,ni1)))/4.
      wx1=(faces_w(meshes_faces(5,ni2))-faces_w(meshes_faces(5,ni1)))/4.
    endif
    if(ni1==0) then
      ux2=(-3*faces_u(meshes_faces(6,n))+4*faces_u(meshes_faces(6,ni2))-faces_u(meshes_faces(6,ni3)))/4.
	  vx2=(-3*faces_v(meshes_faces(6,n))+4*faces_v(meshes_faces(6,ni2))-faces_v(meshes_faces(6,ni3)))/4.
	  wx2=(-3*faces_w(meshes_faces(6,n))+4*faces_w(meshes_faces(6,ni2))-faces_w(meshes_faces(6,ni3)))/4.
    elseif(ni2==0) then
      ux2=(3*faces_u(meshes_faces(6,n))-4*faces_u(meshes_faces(6,ni1))+faces_u(meshes_faces(6,ni3)))/4.
      vx2=(3*faces_v(meshes_faces(6,n))-4*faces_v(meshes_faces(6,ni1))+faces_v(meshes_faces(6,ni3)))/4.
      wx2=(3*faces_w(meshes_faces(6,n))-4*faces_w(meshes_faces(6,ni1))+faces_w(meshes_faces(6,ni3)))/4.
    else
      ux2=(faces_u(meshes_faces(6,ni2))-faces_u(meshes_faces(6,ni1)))/4.
      vx2=(faces_v(meshes_faces(6,ni2))-faces_v(meshes_faces(6,ni1)))/4.
      wx2=(faces_w(meshes_faces(6,ni2))-faces_w(meshes_faces(6,ni1)))/4.
    endif
    uzx=0.5*(ux2-ux1)
	vzx=0.5*(vx2-vx1)
	wzx=0.5*(wx2-wx1)
	
	
    nk1=meshes_neighbors(5,n)
    nk2=meshes_neighbors(6,n)

    if(nk1==0) nk3=meshes_neighbors(6,nk2)
    if(nk2==0) nk3=meshes_neighbors(5,nk1)

    if(nk1==0) then
      uz1=(-3*faces_u(meshes_faces(1,n))+4*faces_u(meshes_faces(1,nk2))-faces_u(meshes_faces(1,nk3)))/4.
	  vz1=(-3*faces_v(meshes_faces(1,n))+4*faces_v(meshes_faces(1,nk2))-faces_v(meshes_faces(1,nk3)))/4.
	  wz1=(-3*faces_w(meshes_faces(1,n))+4*faces_w(meshes_faces(1,nk2))-faces_w(meshes_faces(1,nk3)))/4.
    elseif(nk2==0) then
      uz1=(3*faces_u(meshes_faces(1,n))-4*faces_u(meshes_faces(1,nk1))+faces_u(meshes_faces(1,nk3)))/4.
      vz1=(3*faces_v(meshes_faces(1,n))-4*faces_v(meshes_faces(1,nk1))+faces_v(meshes_faces(1,nk3)))/4.
      wz1=(3*faces_w(meshes_faces(1,n))-4*faces_w(meshes_faces(1,nk1))+faces_w(meshes_faces(1,nk3)))/4.
    else
      uz1=(faces_u(meshes_faces(1,nk2))-faces_u(meshes_faces(1,nk1)))/4.
      vz1=(faces_v(meshes_faces(1,nk2))-faces_v(meshes_faces(1,nk1)))/4.
      wz1=(faces_w(meshes_faces(1,nk2))-faces_w(meshes_faces(1,nk1)))/4.
    endif
    if(nk1==0) then
      uz2=(-3*faces_u(meshes_faces(2,n))+4*faces_u(meshes_faces(2,nk2))-faces_u(meshes_faces(2,nk3)))/4.
	  vz2=(-3*faces_v(meshes_faces(2,n))+4*faces_v(meshes_faces(2,nk2))-faces_v(meshes_faces(2,nk3)))/4.
	  wz2=(-3*faces_w(meshes_faces(2,n))+4*faces_w(meshes_faces(2,nk2))-faces_w(meshes_faces(2,nk3)))/4.
    elseif(nk2==0) then
      uz2=(3*faces_u(meshes_faces(2,n))-4*faces_u(meshes_faces(2,nk1))+faces_u(meshes_faces(2,nk3)))/4.
      vz2=(3*faces_v(meshes_faces(2,n))-4*faces_v(meshes_faces(2,nk1))+faces_v(meshes_faces(2,nk3)))/4.
      wz2=(3*faces_w(meshes_faces(2,n))-4*faces_w(meshes_faces(2,nk1))+faces_w(meshes_faces(2,nk3)))/4.
    else
      uz2=(faces_u(meshes_faces(2,nk2))-faces_u(meshes_faces(2,nk1)))/4.
      vz2=(faces_v(meshes_faces(2,nk2))-faces_v(meshes_faces(2,nk1)))/4.
      wz2=(faces_w(meshes_faces(2,nk2))-faces_w(meshes_faces(2,nk1)))/4.
    endif
    uxz=0.5*(uz2-uz1)
	vxz=0.5*(vz2-vz1)	
	wxz=0.5*(wz2-wz1)	

    if(nk1==0) then
      uz1=(-3*faces_u(meshes_faces(3,n))+4*faces_u(meshes_faces(3,nk2))-faces_u(meshes_faces(3,nk3)))/4.
	  vz1=(-3*faces_v(meshes_faces(3,n))+4*faces_v(meshes_faces(3,nk2))-faces_v(meshes_faces(3,nk3)))/4.
	  wz1=(-3*faces_w(meshes_faces(3,n))+4*faces_w(meshes_faces(3,nk2))-faces_w(meshes_faces(3,nk3)))/4.
    elseif(nk2==0) then
      uz1=(3*faces_u(meshes_faces(3,n))-4*faces_u(meshes_faces(3,nk1))+faces_u(meshes_faces(3,nk3)))/4.
      vz1=(3*faces_v(meshes_faces(3,n))-4*faces_v(meshes_faces(3,nk1))+faces_v(meshes_faces(3,nk3)))/4.
      wz1=(3*faces_w(meshes_faces(3,n))-4*faces_w(meshes_faces(3,nk1))+faces_w(meshes_faces(3,nk3)))/4.
    else
      uz1=(faces_u(meshes_faces(3,nk2))-faces_u(meshes_faces(3,nk1)))/4.
      vz1=(faces_v(meshes_faces(3,nk2))-faces_v(meshes_faces(3,nk1)))/4.
      wz1=(faces_w(meshes_faces(3,nk2))-faces_w(meshes_faces(3,nk1)))/4.
    endif
    if(nk1==0) then
      uz2=(-3*faces_u(meshes_faces(4,n))+4*faces_u(meshes_faces(4,nk2))-faces_u(meshes_faces(4,nk3)))/4.
	  vz2=(-3*faces_v(meshes_faces(4,n))+4*faces_v(meshes_faces(4,nk2))-faces_v(meshes_faces(4,nk3)))/4.
	  wz2=(-3*faces_w(meshes_faces(4,n))+4*faces_w(meshes_faces(4,nk2))-faces_w(meshes_faces(4,nk3)))/4.
    elseif(nk2==0) then
      uz2=(3*faces_u(meshes_faces(4,n))-4*faces_u(meshes_faces(4,nk1))+faces_u(meshes_faces(4,nk3)))/4.
      vz2=(3*faces_v(meshes_faces(4,n))-4*faces_v(meshes_faces(4,nk1))+faces_v(meshes_faces(4,nk3)))/4.
      wz2=(3*faces_w(meshes_faces(4,n))-4*faces_w(meshes_faces(4,nk1))+faces_w(meshes_faces(4,nk3)))/4.
    else
      uz2=(faces_u(meshes_faces(4,nk2))-faces_u(meshes_faces(4,nk1)))/4.
      vz2=(faces_v(meshes_faces(4,nk2))-faces_v(meshes_faces(4,nk1)))/4.
      wz2=(faces_w(meshes_faces(4,nk2))-faces_w(meshes_faces(4,nk1)))/4.
    endif
    uyz=0.5*(uz2-uz1)
	vyz=0.5*(vz2-vz1)
	wyz=0.5*(wz2-wz1)


    ux=(faces_u(meshes_faces(2,n))-faces_u(meshes_faces(1,n)))/2.
    uy=(faces_u(meshes_faces(4,n))-faces_u(meshes_faces(3,n)))/2.
    uz=(faces_u(meshes_faces(6,n))-faces_u(meshes_faces(5,n)))/2.
    vx=(faces_v(meshes_faces(2,n))-faces_v(meshes_faces(1,n)))/2.
    vy=(faces_v(meshes_faces(4,n))-faces_v(meshes_faces(3,n)))/2.
    vz=(faces_v(meshes_faces(6,n))-faces_v(meshes_faces(5,n)))/2.
    wx=(faces_w(meshes_faces(2,n))-faces_w(meshes_faces(1,n)))/2.
    wy=(faces_w(meshes_faces(4,n))-faces_w(meshes_faces(3,n)))/2.
    wz=(faces_w(meshes_faces(6,n))-faces_w(meshes_faces(5,n)))/2.


    meshes_R_D_u(n)=1/Re*(meshes_G1(n)*ux+meshes_G2(n)*uy+meshes_G3(n)*uz-meshes_gij(1,2,n)*uxy-meshes_gij(2,1,n)*uyx &
					     -meshes_gij(1,3,n)*uxz-meshes_gij(3,1,n)*uzx-meshes_gij(2,3,n)*uyz-meshes_gij(3,2,n)*uzy)
						 
    meshes_R_D_v(n)=1/Re*(meshes_G1(n)*vx+meshes_G2(n)*vy+meshes_G3(n)*vz-meshes_gij(1,2,n)*vxy-meshes_gij(2,1,n)*vyx &
					     -meshes_gij(1,3,n)*vxz-meshes_gij(3,1,n)*vzx-meshes_gij(2,3,n)*vyz-meshes_gij(3,2,n)*vzy)	
						 
    meshes_R_D_w(n)=1/Re*(meshes_G1(n)*wx+meshes_G2(n)*wy+meshes_G3(n)*wz-meshes_gij(1,2,n)*wxy-meshes_gij(2,1,n)*wyx &
					     -meshes_gij(1,3,n)*wxz-meshes_gij(3,1,n)*wzx-meshes_gij(2,3,n)*wyz-meshes_gij(3,2,n)*wzy)
end subroutine calculate_R_D2


subroutine calculate_R_D_p2(n)
	integer, intent(in)         :: n
    real                        :: pxy,pxz,pyz,px,py,pz,px1,px2,py1,py2,pz1,pz2,pyx,pzx,pzy	
    integer                     :: i,nj1,nj2,nj3,ni1,ni2,ni3,nk1,nk2,nk3,n1
	
    nj1=meshes_neighbors(3,n)
    nj2=meshes_neighbors(4,n)

    if(nj1==0) nj3=meshes_neighbors(4,nj2)
    if(nj2==0) nj3=meshes_neighbors(3,nj1)
	
    if(nj1==0) then
      py1=(-3*faces_p(meshes_faces(1,n))+4*faces_p(meshes_faces(1,nj2))-faces_p(meshes_faces(1,nj3)))/4.
    elseif(nj2==0) then
      py1=(3*faces_p(meshes_faces(1,n))-4*faces_p(meshes_faces(1,nj1))+faces_p(meshes_faces(1,nj3)))/4.
    else
      py1=(faces_p(meshes_faces(1,nj2))-faces_p(meshes_faces(1,nj1)))/4.
    endif
    if(nj1==0) then
      py2=(-3*faces_p(meshes_faces(2,n))+4*faces_p(meshes_faces(2,nj2))-faces_p(meshes_faces(2,nj3)))/4.
    elseif(nj2==0) then
      py2=(3*faces_p(meshes_faces(2,n))-4*faces_p(meshes_faces(2,nj1))+faces_p(meshes_faces(2,nj3)))/4.
    else
      py2=(faces_p(meshes_faces(2,nj2))-faces_p(meshes_faces(2,nj1)))/4.
    endif

    pxy=0.5*(py2-py1)
	
	
    if(nj1==0) then
      py1=(-3*faces_p(meshes_faces(5,n))+4*faces_p(meshes_faces(5,nj2))-faces_p(meshes_faces(5,nj3)))/4.
    elseif(nj2==0) then
      py1=(3*faces_p(meshes_faces(5,n))-4*faces_p(meshes_faces(5,nj1))+faces_p(meshes_faces(5,nj3)))/4.
    else
      py1=(faces_p(meshes_faces(5,nj2))-faces_p(meshes_faces(5,nj1)))/4.
    endif
    if(nj1==0) then
      py2=(-3*faces_p(meshes_faces(6,n))+4*faces_p(meshes_faces(6,nj2))-faces_p(meshes_faces(6,nj3)))/4.
    elseif(nj2==0) then
      py2=(3*faces_p(meshes_faces(6,n))-4*faces_p(meshes_faces(6,nj1))+faces_p(meshes_faces(6,nj3)))/4.
    else
      py2=(faces_p(meshes_faces(6,nj2))-faces_p(meshes_faces(6,nj1)))/4.
    endif
    pzy=0.5*(py2-py1)	
	

    ni1=meshes_neighbors(1,n)
    ni2=meshes_neighbors(2,n)

    if(ni1==0) ni3=meshes_neighbors(2,ni2)
    if(ni2==0) ni3=meshes_neighbors(1,ni1)


    if(ni1==0) then
      px1=(-3*faces_p(meshes_faces(3,n))+4*faces_p(meshes_faces(3,ni2))-faces_p(meshes_faces(3,ni3)))/4.
    elseif(ni2==0) then
      px1=(3*faces_p(meshes_faces(3,n))-4*faces_p(meshes_faces(3,ni1))+faces_p(meshes_faces(3,ni3)))/4.
    else
      px1=(faces_p(meshes_faces(3,ni2))-faces_p(meshes_faces(3,ni1)))/4.
    endif
    if(ni1==0) then
      px2=(-3*faces_p(meshes_faces(4,n))+4*faces_p(meshes_faces(4,ni2))-faces_p(meshes_faces(4,ni3)))/4.
    elseif(ni2==0) then
      px2=(3*faces_p(meshes_faces(4,n))-4*faces_p(meshes_faces(4,ni1))+faces_p(meshes_faces(4,ni3)))/4.
    else
      px2=(faces_p(meshes_faces(4,ni2))-faces_p(meshes_faces(4,ni1)))/4.
    endif

    pyx=0.5*(px2-px1)


    if(ni1==0) then
      px1=(-3*faces_p(meshes_faces(5,n))+4*faces_p(meshes_faces(5,ni2))-faces_p(meshes_faces(5,ni3)))/4.
    elseif(ni2==0) then
      px1=(3*faces_p(meshes_faces(5,n))-4*faces_p(meshes_faces(5,ni1))+faces_p(meshes_faces(5,ni3)))/4.
    else
      px1=(faces_p(meshes_faces(5,ni2))-faces_p(meshes_faces(5,ni1)))/4.
    endif
    if(ni1==0) then
      px2=(-3*faces_p(meshes_faces(6,n))+4*faces_p(meshes_faces(6,ni2))-faces_p(meshes_faces(6,ni3)))/4.
    elseif(ni2==0) then
      px2=(3*faces_p(meshes_faces(6,n))-4*faces_p(meshes_faces(6,ni1))+faces_p(meshes_faces(6,ni3)))/4.
    else
      px2=(faces_p(meshes_faces(6,ni2))-faces_p(meshes_faces(6,ni1)))/4.
    endif
    pzx=0.5*(px2-px1)
	
	
    nk1=meshes_neighbors(5,n)
    nk2=meshes_neighbors(6,n)

    if(nk1==0) nk3=meshes_neighbors(6,nk2)
    if(nk2==0) nk3=meshes_neighbors(5,nk1)

    if(nk1==0) then
      pz1=(-3*faces_p(meshes_faces(1,n))+4*faces_p(meshes_faces(1,nk2))-faces_p(meshes_faces(1,nk3)))/4.
    elseif(nk2==0) then
      pz1=(3*faces_p(meshes_faces(1,n))-4*faces_p(meshes_faces(1,nk1))+faces_p(meshes_faces(1,nk3)))/4.
    else
      pz1=(faces_p(meshes_faces(1,nk2))-faces_p(meshes_faces(1,nk1)))/4.
    endif
    if(nk1==0) then
      pz2=(-3*faces_p(meshes_faces(2,n))+4*faces_p(meshes_faces(2,nk2))-faces_p(meshes_faces(2,nk3)))/4.
    elseif(nk2==0) then
      pz2=(3*faces_p(meshes_faces(2,n))-4*faces_p(meshes_faces(2,nk1))+faces_p(meshes_faces(2,nk3)))/4.
    else
      pz2=(faces_p(meshes_faces(2,nk2))-faces_p(meshes_faces(2,nk1)))/4.
    endif
    pxz=0.5*(pz2-pz1)

    if(nk1==0) then
      pz1=(-3*faces_p(meshes_faces(3,n))+4*faces_p(meshes_faces(3,nk2))-faces_p(meshes_faces(3,nk3)))/4.
    elseif(nk2==0) then
      pz1=(3*faces_p(meshes_faces(3,n))-4*faces_p(meshes_faces(3,nk1))+faces_p(meshes_faces(3,nk3)))/4.
    else
      pz1=(faces_p(meshes_faces(3,nk2))-faces_p(meshes_faces(3,nk1)))/4.
    endif
    if(nk1==0) then
      pz2=(-3*faces_p(meshes_faces(4,n))+4*faces_p(meshes_faces(4,nk2))-faces_p(meshes_faces(4,nk3)))/4.
    elseif(nk2==0) then
      pz2=(3*faces_p(meshes_faces(4,n))-4*faces_p(meshes_faces(4,nk1))+faces_p(meshes_faces(4,nk3)))/4.
    else
      pz2=(faces_p(meshes_faces(4,nk2))-faces_p(meshes_faces(4,nk1)))/4.
    endif
    pyz=0.5*(pz2-pz1)


    px=(faces_p(meshes_faces(2,n))-faces_p(meshes_faces(1,n)))/2.
    py=(faces_p(meshes_faces(4,n))-faces_p(meshes_faces(3,n)))/2.
    pz=(faces_p(meshes_faces(6,n))-faces_p(meshes_faces(5,n)))/2.


    meshes_R_D_p(n)=1/rho*(meshes_G1(n)*px+meshes_G2(n)*py+meshes_G3(n)*pz-meshes_gij(1,2,n)*pxy-meshes_gij(2,1,n)*pyx &
					      -meshes_gij(1,3,n)*pxz-meshes_gij(3,1,n)*pzx-meshes_gij(2,3,n)*pyz-meshes_gij(3,2,n)*pzy)
end subroutine calculate_R_D_p2


subroutine calculate_R_D_T(n)
	integer, intent(in)         :: n
    real                        :: Txy,Txz,Tyz,Tx,Ty,Tz,Tx1,Tx2,Ty1,Ty2,Tz1,Tz2,Tyx,Tzx,Tzy	
    integer                     :: i,nj1,nj2,nj3,ni1,ni2,ni3,nk1,nk2,nk3,n1
	
    nj1=meshes_neighbors(3,n)
    nj2=meshes_neighbors(4,n)

    if(nj1==0) nj3=meshes_neighbors(4,nj2)
    if(nj2==0) nj3=meshes_neighbors(3,nj1)
	
    if(nj1==0) then
      Ty1=(-3*faces_T(meshes_faces(1,n))+4*faces_T(meshes_faces(1,nj2))-faces_T(meshes_faces(1,nj3)))/4.
    elseif(nj2==0) then
      Ty1=(3*faces_T(meshes_faces(1,n))-4*faces_T(meshes_faces(1,nj1))+faces_T(meshes_faces(1,nj3)))/4.
    else
      Ty1=(faces_T(meshes_faces(1,nj2))-faces_T(meshes_faces(1,nj1)))/4.
    endif
    if(nj1==0) then
      Ty2=(-3*faces_T(meshes_faces(2,n))+4*faces_T(meshes_faces(2,nj2))-faces_T(meshes_faces(2,nj3)))/4.
    elseif(nj2==0) then
      Ty2=(3*faces_T(meshes_faces(2,n))-4*faces_T(meshes_faces(2,nj1))+faces_T(meshes_faces(2,nj3)))/4.
    else
      Ty2=(faces_T(meshes_faces(2,nj2))-faces_T(meshes_faces(2,nj1)))/4.
    endif

    Txy=0.5*(Ty2-Ty1)
	
	
    if(nj1==0) then
      Ty1=(-3*faces_T(meshes_faces(5,n))+4*faces_T(meshes_faces(5,nj2))-faces_T(meshes_faces(5,nj3)))/4.
    elseif(nj2==0) then
      Ty1=(3*faces_T(meshes_faces(5,n))-4*faces_T(meshes_faces(5,nj1))+faces_T(meshes_faces(5,nj3)))/4.
    else
      Ty1=(faces_T(meshes_faces(5,nj2))-faces_T(meshes_faces(5,nj1)))/4.
    endif
    if(nj1==0) then
      Ty2=(-3*faces_T(meshes_faces(6,n))+4*faces_T(meshes_faces(6,nj2))-faces_T(meshes_faces(6,nj3)))/4.
    elseif(nj2==0) then
      Ty2=(3*faces_T(meshes_faces(6,n))-4*faces_T(meshes_faces(6,nj1))+faces_T(meshes_faces(6,nj3)))/4.
    else
      Ty2=(faces_T(meshes_faces(6,nj2))-faces_T(meshes_faces(6,nj1)))/4.
    endif
    Tzy=0.5*(Ty2-Ty1)	
	

    ni1=meshes_neighbors(1,n)
    ni2=meshes_neighbors(2,n)

    if(ni1==0) ni3=meshes_neighbors(2,ni2)
    if(ni2==0) ni3=meshes_neighbors(1,ni1)


    if(ni1==0) then
      Tx1=(-3*faces_T(meshes_faces(3,n))+4*faces_T(meshes_faces(3,ni2))-faces_T(meshes_faces(3,ni3)))/4.
    elseif(ni2==0) then
      Tx1=(3*faces_T(meshes_faces(3,n))-4*faces_T(meshes_faces(3,ni1))+faces_T(meshes_faces(3,ni3)))/4.
    else
      Tx1=(faces_T(meshes_faces(3,ni2))-faces_T(meshes_faces(3,ni1)))/4.
    endif
    if(ni1==0) then
      Tx2=(-3*faces_T(meshes_faces(4,n))+4*faces_T(meshes_faces(4,ni2))-faces_T(meshes_faces(4,ni3)))/4.
    elseif(ni2==0) then
      Tx2=(3*faces_T(meshes_faces(4,n))-4*faces_T(meshes_faces(4,ni1))+faces_T(meshes_faces(4,ni3)))/4.
    else
      Tx2=(faces_T(meshes_faces(4,ni2))-faces_T(meshes_faces(4,ni1)))/4.
    endif

    Tyx=0.5*(Tx2-Tx1)


    if(ni1==0) then
      Tx1=(-3*faces_T(meshes_faces(5,n))+4*faces_T(meshes_faces(5,ni2))-faces_T(meshes_faces(5,ni3)))/4.
    elseif(ni2==0) then
      Tx1=(3*faces_T(meshes_faces(5,n))-4*faces_T(meshes_faces(5,ni1))+faces_T(meshes_faces(5,ni3)))/4.
    else
      Tx1=(faces_T(meshes_faces(5,ni2))-faces_T(meshes_faces(5,ni1)))/4.
    endif
    if(ni1==0) then
      Tx2=(-3*faces_T(meshes_faces(6,n))+4*faces_T(meshes_faces(6,ni2))-faces_T(meshes_faces(6,ni3)))/4.
    elseif(ni2==0) then
      Tx2=(3*faces_T(meshes_faces(6,n))-4*faces_T(meshes_faces(6,ni1))+faces_T(meshes_faces(6,ni3)))/4.
    else
      Tx2=(faces_T(meshes_faces(6,ni2))-faces_T(meshes_faces(6,ni1)))/4.
    endif
    Tzx=0.5*(Tx2-Tx1)
	
	
    nk1=meshes_neighbors(5,n)
    nk2=meshes_neighbors(6,n)

    if(nk1==0) nk3=meshes_neighbors(6,nk2)
    if(nk2==0) nk3=meshes_neighbors(5,nk1)

    if(nk1==0) then
      Tz1=(-3*faces_T(meshes_faces(1,n))+4*faces_T(meshes_faces(1,nk2))-faces_T(meshes_faces(1,nk3)))/4.
    elseif(nk2==0) then
      Tz1=(3*faces_T(meshes_faces(1,n))-4*faces_T(meshes_faces(1,nk1))+faces_T(meshes_faces(1,nk3)))/4.
    else
      Tz1=(faces_T(meshes_faces(1,nk2))-faces_T(meshes_faces(1,nk1)))/4.
    endif
    if(nk1==0) then
      Tz2=(-3*faces_T(meshes_faces(2,n))+4*faces_T(meshes_faces(2,nk2))-faces_T(meshes_faces(2,nk3)))/4.
    elseif(nk2==0) then
      Tz2=(3*faces_T(meshes_faces(2,n))-4*faces_T(meshes_faces(2,nk1))+faces_T(meshes_faces(2,nk3)))/4.
    else
      Tz2=(faces_T(meshes_faces(2,nk2))-faces_T(meshes_faces(2,nk1)))/4.
    endif
    Txz=0.5*(Tz2-Tz1)

    if(nk1==0) then
      Tz1=(-3*faces_T(meshes_faces(3,n))+4*faces_T(meshes_faces(3,nk2))-faces_T(meshes_faces(3,nk3)))/4.
    elseif(nk2==0) then
      Tz1=(3*faces_T(meshes_faces(3,n))-4*faces_T(meshes_faces(3,nk1))+faces_T(meshes_faces(3,nk3)))/4.
    else
      Tz1=(faces_T(meshes_faces(3,nk2))-faces_T(meshes_faces(3,nk1)))/4.
    endif
    if(nk1==0) then
      Tz2=(-3*faces_T(meshes_faces(4,n))+4*faces_T(meshes_faces(4,nk2))-faces_T(meshes_faces(4,nk3)))/4.
    elseif(nk2==0) then
      Tz2=(3*faces_T(meshes_faces(4,n))-4*faces_T(meshes_faces(4,nk1))+faces_T(meshes_faces(4,nk3)))/4.
    else
      Tz2=(faces_T(meshes_faces(4,nk2))-faces_T(meshes_faces(4,nk1)))/4.
    endif
    Tyz=0.5*(Tz2-Tz1)


    Tx=(faces_T(meshes_faces(2,n))-faces_T(meshes_faces(1,n)))/2.
    Ty=(faces_T(meshes_faces(4,n))-faces_T(meshes_faces(3,n)))/2.
    Tz=(faces_T(meshes_faces(6,n))-faces_T(meshes_faces(5,n)))/2.


    meshes_R_D_T(n)=alpha*(meshes_G1(n)*Tx+meshes_G2(n)*Ty+meshes_G3(n)*Tz-meshes_gij(1,2,n)*Txy-meshes_gij(2,1,n)*Tyx &
					      -meshes_gij(1,3,n)*Txz-meshes_gij(3,1,n)*Tzx-meshes_gij(2,3,n)*Tyz-meshes_gij(3,2,n)*Tzy)
end subroutine calculate_R_D_T

  subroutine calculate_R_D (n)
    integer, intent(in)                  :: n
    integer                              :: i
    real,dimension(1:4)                  :: uxy_i,uxz_i,uyz_i        ! The derivatives
    real,dimension(1:4)                  :: vxy_i,vxz_i,vyz_i        ! The derivatives
    real,dimension(1:4)                  :: wxy_i,wxz_i,wyz_i        ! The derivatives
    integer,dimension(1:3)               :: counter 
    real, dimension(1:3)                 :: sums_u,sums_v,sums_w
    real                                 :: uxy,uxz,uyz,ux,uy,uz,uyx,uzy,uzx
    real                                 :: vxy,vxz,vyz,vx,vy,vz,vyx,vzy,vzx
    real                                 :: wxy,wxz,wyz,wx,wy,wz,wyx,wzy,wzx
    integer, dimension(1:6)              :: n_faces,n_faces1,n_faces2,n_faces3,n_faces4,n_faces5,n_faces6
! i, -x mesh, ii +x mesh

    n_faces(:)=meshes_faces(:,n)
    counter=(/0,0,0/)
    uxy_i  =(/0.,0.,0.,0./)
    uxz_i  =(/0.,0.,0.,0./)
    uyz_i  =(/0.,0.,0.,0./)
    sums_u   =(/0.,0.,0./)
    vxy_i  =(/0.,0.,0.,0./)
    vxz_i  =(/0.,0.,0.,0./)
    vyz_i  =(/0.,0.,0.,0./)
    sums_v   =(/0.,0.,0./)
    wxy_i  =(/0.,0.,0.,0./)
    wxz_i  =(/0.,0.,0.,0./)
    wyz_i  =(/0.,0.,0.,0./)
    sums_w   =(/0.,0.,0./)
!         =  FACE 1 =
    if (faces_BC(n_faces(1)) == BC_id(1)) Then
       n_faces1(:)=meshes_faces(:,meshes_neighbors(1,n))
       uxy_i(1)=(faces_u(n_faces(4))-faces_u(n_faces(3))-faces_u(n_faces1(4))+faces_u(n_faces1(3)))/4.
       vxy_i(1)=(faces_v(n_faces(4))-faces_v(n_faces(3))-faces_v(n_faces1(4))+faces_v(n_faces1(3)))/4.
       wxy_i(1)=(faces_w(n_faces(4))-faces_w(n_faces(3))-faces_w(n_faces1(4))+faces_w(n_faces1(3)))/4.

       counter(1)=counter(1)+1
       uxz_i(1)=(faces_u(n_faces(6))-faces_u(n_faces(5))-faces_u(n_faces1(6))+faces_u(n_faces1(5)))/4.
       vxz_i(1)=(faces_v(n_faces(6))-faces_v(n_faces(5))-faces_v(n_faces1(6))+faces_v(n_faces1(5)))/4.
       wxz_i(1)=(faces_w(n_faces(6))-faces_w(n_faces(5))-faces_w(n_faces1(6))+faces_w(n_faces1(5)))/4.

       counter(2)=counter(2)+1
    endif 
    if (faces_BC(n_faces(2)) == BC_id(1)) Then
       n_faces2(:)=meshes_faces(:,meshes_neighbors(2,n))
       uxy_i(2)=-(faces_u(n_faces(4))-faces_u(n_faces(3))-faces_u(n_faces2(4))+faces_u(n_faces2(3)))/4
       vxy_i(2)=-(faces_v(n_faces(4))-faces_v(n_faces(3))-faces_v(n_faces2(4))+faces_v(n_faces2(3)))/4
       wxy_i(2)=-(faces_w(n_faces(4))-faces_w(n_faces(3))-faces_w(n_faces2(4))+faces_w(n_faces2(3)))/4
       counter(1)=counter(1)+1
       uxz_i(2)=-(faces_u(n_faces(6))-faces_u(n_faces(5))-faces_u(n_faces2(6))+faces_u(n_faces2(5)))/4
       vxz_i(2)=-(faces_v(n_faces(6))-faces_v(n_faces(5))-faces_v(n_faces2(6))+faces_v(n_faces2(5)))/4
       wxz_i(2)=-(faces_w(n_faces(6))-faces_w(n_faces(5))-faces_w(n_faces2(6))+faces_w(n_faces2(5)))/4
       counter(2)=counter(2)+1
    endif
    if (faces_BC(n_faces(3)) == BC_id(1)) Then
       n_faces3(:)=meshes_faces(:,meshes_neighbors(3,n))
       uxy_i(3)=(faces_u(n_faces(2))-faces_u(n_faces(1))-faces_u(n_faces3(2))+faces_u(n_faces3(1)))/4.
       vxy_i(3)=(faces_v(n_faces(2))-faces_v(n_faces(1))-faces_v(n_faces3(2))+faces_v(n_faces3(1)))/4.
       wxy_i(3)=(faces_w(n_faces(2))-faces_w(n_faces(1))-faces_w(n_faces3(2))+faces_w(n_faces3(1)))/4.
       counter(1)=counter(1)+1
       uyz_i(1)=(faces_u(n_faces(6))-faces_u(n_faces(5))-faces_u(n_faces3(6))+faces_u(n_faces3(5)))/4.
       vyz_i(1)=(faces_v(n_faces(6))-faces_v(n_faces(5))-faces_v(n_faces3(6))+faces_v(n_faces3(5)))/4.
       wyz_i(1)=(faces_w(n_faces(6))-faces_w(n_faces(5))-faces_w(n_faces3(6))+faces_w(n_faces3(5)))/4.
       counter(3)=counter(3)+1
    endif
    if (faces_BC(n_faces(4)) == BC_id(1)) Then
       n_faces4(:)=meshes_faces(:,meshes_neighbors(4,n))
       uxy_i(4)=-(faces_u(n_faces(2))-faces_u(n_faces(1))-faces_u(n_faces4(2))+faces_u(n_faces4(1)))/4
       vxy_i(4)=-(faces_v(n_faces(2))-faces_v(n_faces(1))-faces_v(n_faces4(2))+faces_v(n_faces4(1)))/4
       wxy_i(4)=-(faces_w(n_faces(2))-faces_w(n_faces(1))-faces_w(n_faces4(2))+faces_w(n_faces4(1)))/4
       counter(1)=counter(1)+1
       uyz_i(2)=-(faces_u(n_faces(6))-faces_u(n_faces(5))-faces_u(n_faces4(6))+faces_u(n_faces4(5)))/4
       vyz_i(2)=-(faces_v(n_faces(6))-faces_v(n_faces(5))-faces_v(n_faces4(6))+faces_v(n_faces4(5)))/4
       wyz_i(2)=-(faces_w(n_faces(6))-faces_w(n_faces(5))-faces_w(n_faces4(6))+faces_w(n_faces4(5)))/4
       counter(3)=counter(3)+1
    endif
    if (faces_BC(n_faces(5)) == BC_id(1)) Then
       n_faces5(:)=meshes_faces(:,meshes_neighbors(5,n))
       uxz_i(3)=(faces_u(n_faces(2))-faces_u(n_faces(1))-faces_u(n_faces5(2))+faces_u(n_faces5(1)))/4.
       vxz_i(3)=(faces_v(n_faces(2))-faces_v(n_faces(1))-faces_v(n_faces5(2))+faces_v(n_faces5(1)))/4.
       wxz_i(3)=(faces_w(n_faces(2))-faces_w(n_faces(1))-faces_w(n_faces5(2))+faces_w(n_faces5(1)))/4.
       counter(2)=counter(2)+1
       uyz_i(3)=(faces_u(n_faces(4))-faces_u(n_faces(3))-faces_u(n_faces5(4))+faces_u(n_faces5(3)))/4.
       vyz_i(3)=(faces_v(n_faces(4))-faces_v(n_faces(3))-faces_v(n_faces5(4))+faces_v(n_faces5(3)))/4.
       wyz_i(3)=(faces_w(n_faces(4))-faces_w(n_faces(3))-faces_w(n_faces5(4))+faces_w(n_faces5(3)))/4.
       counter(3)=counter(3)+1
    endif
    if (faces_BC(n_faces(6)) == BC_id(1)) Then
       n_faces6(:)=meshes_faces(:,meshes_neighbors(6,n))
       uxz_i(4)=-(faces_u(n_faces(2))-faces_u(n_faces(1))-faces_u(n_faces6(2))+faces_u(n_faces6(1)))/4
       vxz_i(4)=-(faces_v(n_faces(2))-faces_v(n_faces(1))-faces_v(n_faces6(2))+faces_v(n_faces6(1)))/4
       wxz_i(4)=-(faces_w(n_faces(2))-faces_w(n_faces(1))-faces_w(n_faces6(2))+faces_w(n_faces6(1)))/4
       counter(2)=counter(2)+1
       uyz_i(4)=-(faces_u(n_faces(4))-faces_u(n_faces(3))-faces_u(n_faces6(4))+faces_u(n_faces6(3)))/4
       vyz_i(4)=-(faces_v(n_faces(4))-faces_v(n_faces(3))-faces_v(n_faces6(4))+faces_v(n_faces6(3)))/4
       wyz_i(4)=-(faces_w(n_faces(4))-faces_w(n_faces(3))-faces_w(n_faces6(4))+faces_w(n_faces6(3)))/4
       counter(3)=counter(3)+1
    endif

    do i=1,4
      sums_u(1)=sums_u(1)+uxy_i(i)
      sums_u(2)=sums_u(2)+uxz_i(i)
      sums_u(3)=sums_u(3)+uyz_i(i)
	  
      sums_v(1)=sums_v(1)+vxy_i(i)
      sums_v(2)=sums_v(2)+vxz_i(i)
      sums_v(3)=sums_v(3)+vyz_i(i)
	  
      sums_w(1)=sums_w(1)+wxy_i(i)
      sums_w(2)=sums_w(2)+wxz_i(i)
      sums_w(3)=sums_w(3)+wyz_i(i)	  
	  
    enddo
    
    uxy=sums_u(1)/counter(1)
    uxz=sums_u(2)/counter(2)
    uyz=sums_u(3)/counter(3)

    vxy=sums_v(1)/counter(1)
    vxz=sums_v(2)/counter(2)
    vyz=sums_v(3)/counter(3)

    wxy=sums_w(1)/counter(1)
    wxz=sums_w(2)/counter(2)
    wyz=sums_w(3)/counter(3)

    ux=(faces_u(meshes_faces(2,n))-faces_u(meshes_faces(1,n)))/2.
    uy=(faces_u(meshes_faces(4,n))-faces_u(meshes_faces(3,n)))/2.
    uz=(faces_u(meshes_faces(6,n))-faces_u(meshes_faces(5,n)))/2.
    vx=(faces_v(meshes_faces(2,n))-faces_v(meshes_faces(1,n)))/2.
    vy=(faces_v(meshes_faces(4,n))-faces_v(meshes_faces(3,n)))/2.
    vz=(faces_v(meshes_faces(6,n))-faces_v(meshes_faces(5,n)))/2.
    wx=(faces_w(meshes_faces(2,n))-faces_w(meshes_faces(1,n)))/2.
    wy=(faces_w(meshes_faces(4,n))-faces_w(meshes_faces(3,n)))/2.
    wz=(faces_w(meshes_faces(6,n))-faces_w(meshes_faces(5,n)))/2.


    meshes_R_D_u(n)=1/Re*(meshes_G1(n)*ux+meshes_G2(n)*uy+meshes_G3(n)*uz-meshes_gij(1,2,n)*uxy-meshes_gij(2,1,n)*uyx &
					      -meshes_gij(1,3,n)*uxz-meshes_gij(3,1,n)*uzx-meshes_gij(2,3,n)*uyz-meshes_gij(3,2,n)*uzy)
						  
    meshes_R_D_v(n)=1/Re*(meshes_G1(n)*vx+meshes_G2(n)*vy+meshes_G3(n)*vz-meshes_gij(1,2,n)*vxy-meshes_gij(2,1,n)*vyx &
					      -meshes_gij(1,3,n)*vxz-meshes_gij(3,1,n)*vzx-meshes_gij(2,3,n)*vyz-meshes_gij(3,2,n)*vzy)

    meshes_R_D_w(n)=1/Re*(meshes_G1(n)*wx+meshes_G2(n)*wy+meshes_G3(n)*wz-meshes_gij(1,2,n)*wxy-meshes_gij(2,1,n)*wyx &
					      -meshes_gij(1,3,n)*wxz-meshes_gij(3,1,n)*wzx-meshes_gij(2,3,n)*wyz-meshes_gij(3,2,n)*wzy)						  

  end subroutine calculate_R_D


  subroutine calculate_R_D_p (n)
    integer, intent(in)                  :: n
    integer                              :: i
    real,dimension(1:4)                  :: pxy_i,pxz_i,pyz_i        ! The derivatives
    integer,dimension(1:3)               :: counter 
    real, dimension(1:3)                 :: sums
    real                                 :: pxy,pxz,pyz,px,py,pz,pyx,pzy,pzx
    integer, dimension(1:6)              :: n_faces,n_faces1,n_faces2,n_faces3,n_faces4,n_faces5,n_faces6
! i, -x mesh, ii +x mesh

    n_faces(:)=meshes_faces(:,n)
    counter=(/0,0,0/)
    pxy_i  =(/0.,0.,0.,0./)
    pxz_i  =(/0.,0.,0.,0./)
    pyz_i  =(/0.,0.,0.,0./)
    sums   =(/0.,0.,0./)
!         =  FACE 1 =
    if (faces_BC(n_faces(1)) == BC_id(1)) Then
       n_faces1(:)=meshes_faces(:,meshes_neighbors(1,n))
       pxy_i(1)=(faces_p(n_faces(4))-faces_p(n_faces(3))-faces_p(n_faces1(4))+faces_p(n_faces1(3)))/4.
       counter(1)=counter(1)+1
       pxz_i(1)=(faces_p(n_faces(6))-faces_p(n_faces(5))-faces_p(n_faces1(6))+faces_p(n_faces1(5)))/4.
       counter(2)=counter(2)+1
    endif 
    if (faces_BC(n_faces(2)) == BC_id(1)) Then
       n_faces2(:)=meshes_faces(:,meshes_neighbors(2,n))
       pxy_i(2)=-(faces_p(n_faces(4))-faces_p(n_faces(3))-faces_p(n_faces2(4))+faces_p(n_faces2(3)))/4
       counter(1)=counter(1)+1
       pxz_i(2)=-(faces_p(n_faces(6))-faces_p(n_faces(5))-faces_p(n_faces2(6))+faces_p(n_faces2(5)))/4
       counter(2)=counter(2)+1
    endif
    if (faces_BC(n_faces(3)) == BC_id(1)) Then
       n_faces3(:)=meshes_faces(:,meshes_neighbors(3,n))
       pxy_i(3)=(faces_p(n_faces(2))-faces_p(n_faces(1))-faces_p(n_faces3(2))+faces_p(n_faces3(1)))/4.
       counter(1)=counter(1)+1
       pyz_i(1)=(faces_p(n_faces(6))-faces_p(n_faces(5))-faces_p(n_faces3(6))+faces_p(n_faces3(5)))/4.
       counter(3)=counter(3)+1
    endif
    if (faces_BC(n_faces(4)) == BC_id(1)) Then
       n_faces4(:)=meshes_faces(:,meshes_neighbors(4,n))
       pxy_i(4)=-(faces_p(n_faces(2))-faces_p(n_faces(1))-faces_p(n_faces4(2))+faces_p(n_faces4(1)))/4
       counter(1)=counter(1)+1
       pyz_i(2)=-(faces_p(n_faces(6))-faces_p(n_faces(5))-faces_p(n_faces4(6))+faces_p(n_faces4(5)))/4
       counter(3)=counter(3)+1
    endif
    if (faces_BC(n_faces(5)) == BC_id(1)) Then
       n_faces5(:)=meshes_faces(:,meshes_neighbors(5,n))
       pxz_i(3)=(faces_p(n_faces(2))-faces_p(n_faces(1))-faces_p(n_faces5(2))+faces_p(n_faces5(1)))/4.
       counter(2)=counter(2)+1
       pyz_i(3)=(faces_p(n_faces(4))-faces_p(n_faces(3))-faces_p(n_faces5(4))+faces_p(n_faces5(3)))/4.
       counter(3)=counter(3)+1
    endif
    if (faces_BC(n_faces(6)) == BC_id(1)) Then
       n_faces6(:)=meshes_faces(:,meshes_neighbors(6,n))
       pxz_i(4)=-(faces_p(n_faces(2))-faces_p(n_faces(1))-faces_p(n_faces6(2))+faces_p(n_faces6(1)))/4
       counter(2)=counter(2)+1
       pyz_i(4)=-(faces_p(n_faces(4))-faces_p(n_faces(3))-faces_p(n_faces6(4))+faces_p(n_faces6(3)))/4
       counter(3)=counter(3)+1
    endif

    do i=1,4
      sums(1)=sums(1)+pxy_i(i)
      sums(2)=sums(2)+pxz_i(i)
      sums(3)=sums(3)+pyz_i(i)
    enddo
    
    pxy=sums(1)/counter(1)
    pxz=sums(2)/counter(2)
    pyz=sums(3)/counter(3)

    px=(faces_p(meshes_faces(2,n))-faces_p(meshes_faces(1,n)))/2.
    py=(faces_p(meshes_faces(4,n))-faces_p(meshes_faces(3,n)))/2.
    pz=(faces_p(meshes_faces(6,n))-faces_p(meshes_faces(5,n)))/2.

    meshes_R_D_p(n)=1/rho*(meshes_G1(n)*px+meshes_G2(n)*py+meshes_G3(n)*pz-meshes_gij(1,2,n)*pxy-meshes_gij(2,1,n)*pyx &
					      -meshes_gij(1,3,n)*pxz-meshes_gij(3,1,n)*pzx-meshes_gij(2,3,n)*pyz-meshes_gij(3,2,n)*pzy)
  end subroutine calculate_R_D_p

subroutine calculate_R_b_p(n)
	integer, intent(in)         :: n
	real                        :: du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz
	real                        :: dbx_dx, dby_dy, dbz_dz, Dila
    integer, dimension(1:6)     :: n_faces
	integer                     :: k
    real, dimension(1:6)        :: u_faces,v_faces,w_faces,bx_faces,by_faces,bz_faces
	
	
	n_faces(:)=meshes_faces(:,n)
	do k=1,6
		u_faces(k)=faces_u(n_faces(k))
		v_faces(k)=faces_v(n_faces(k))
		w_faces(k)=faces_w(n_faces(k))
		bx_faces(k)=faces_bx(n_faces(k))
		by_faces(k)=faces_by(n_faces(k))
		bz_faces(k)=faces_bz(n_faces(k))
	enddo 		


! ==============================================================================================
! Derivatives in Cartesian system
    du_dx= dphi_dx2(u_faces,n_faces,n)
    du_dy= dphi_dy2(u_faces,n_faces,n)
	du_dz= dphi_dz2(u_faces,n_faces,n)

    dv_dx= dphi_dx2(v_faces,n_faces,n)
    dv_dy= dphi_dy2(v_faces,n_faces,n)
	dv_dz= dphi_dz2(v_faces,n_faces,n)
	
    dw_dx= dphi_dx2(w_faces,n_faces,n)
    dw_dy= dphi_dy2(w_faces,n_faces,n)
	dw_dz= dphi_dz2(w_faces,n_faces,n)

! bx, by, bz
	dbx_dx= dphi_dx2(bx_faces,n_faces,n)
	dby_dy= dphi_dy2(by_faces,n_faces,n)
	dbz_dz= dphi_dz2(bz_faces,n_faces,n)
	
! ==============================================================================================
	Dila=du_dx+dv_dy+dw_dz
	meshes_dila(n)=Dila

	meshes_R_b_p(n)=-(du_dx**2. + dv_dy**2. + dw_dz**2.)-2.0*(du_dy*dv_dx+dv_dz*dw_dy+dw_dx*du_dz)-dbx_dx-dby_dy-dbz_dz+Dila/(2*DT)

end subroutine calculate_R_b_p

subroutine calculate_R_b_T(n)
	integer, intent(in)         :: n

	meshes_R_b_T(n)=meshes_q(n)

end subroutine calculate_R_b_T


subroutine calculate_R_b_p2(n)
	integer, intent(in)         :: n
	real                        :: du_dxi1,du_dxi2,du_dxi3,dv_dxi1,dv_dxi2,dv_dxi3,dw_dxi1,dw_dxi2,dw_dxi3
	real                        :: du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz
	real                        :: dbx_dxi1,dbx_dxi2,dbx_dxi3,dby_dxi1,dby_dxi2,dby_dxi3,dbz_dxi1,dbz_dxi2,dbz_dxi3
	real                        :: dbx_dx, dby_dy, dbz_dz, Dila
    integer, dimension(1:6)     :: n_faces
	
	
	n_faces(:)=meshes_faces(:,n)
! ==============================================================================================
! Derivatives in curvilinear coordinates
    du_dxi1=(faces_u(n_faces(2))-faces_u(n_faces(1)))/2.
    du_dxi2=(faces_u(n_faces(4))-faces_u(n_faces(3)))/2.
    du_dxi3=(faces_u(n_faces(6))-faces_u(n_faces(5)))/2.

    dv_dxi1=(faces_v(n_faces(2))-faces_v(n_faces(1)))/2.
    dv_dxi2=(faces_v(n_faces(4))-faces_v(n_faces(3)))/2.
    dv_dxi3=(faces_v(n_faces(6))-faces_v(n_faces(5)))/2.
	
    dw_dxi1=(faces_w(n_faces(2))-faces_w(n_faces(1)))/2.
    dw_dxi2=(faces_w(n_faces(4))-faces_w(n_faces(3)))/2.
    dw_dxi3=(faces_w(n_faces(6))-faces_w(n_faces(5)))/2.

! bx, by, bz

    dbx_dxi1=(faces_bx(n_faces(2))-faces_bx(n_faces(1)))/2.
    dbx_dxi2=(faces_bx(n_faces(4))-faces_bx(n_faces(3)))/2.
    dbx_dxi3=(faces_bx(n_faces(6))-faces_bx(n_faces(5)))/2.

    dby_dxi1=(faces_by(n_faces(2))-faces_by(n_faces(1)))/2.
    dby_dxi2=(faces_by(n_faces(4))-faces_by(n_faces(3)))/2.
    dby_dxi3=(faces_by(n_faces(6))-faces_by(n_faces(5)))/2.

    dbz_dxi1=(faces_bz(n_faces(2))-faces_bz(n_faces(1)))/2.
    dbz_dxi2=(faces_bz(n_faces(4))-faces_bz(n_faces(3)))/2.
    dbz_dxi3=(faces_bz(n_faces(6))-faces_bz(n_faces(5)))/2.


! ==============================================================================================
! Derivatives in Cartesian system
    du_dx= dphi_dx(du_dxi1,du_dxi2,du_dxi3,n)
    du_dy= dphi_dy(du_dxi1,du_dxi2,du_dxi3,n)
	du_dz= dphi_dz(du_dxi1,du_dxi2,du_dxi3,n)

    dv_dx= dphi_dx(dv_dxi1,dv_dxi2,dv_dxi3,n)
    dv_dy= dphi_dy(dv_dxi1,dv_dxi2,dv_dxi3,n)
	dv_dz= dphi_dz(dv_dxi1,dv_dxi2,dv_dxi3,n)
	
    dw_dx= dphi_dx(dw_dxi1,dw_dxi2,dw_dxi3,n)
    dw_dy= dphi_dy(dw_dxi1,dw_dxi2,dw_dxi3,n)
	dw_dz= dphi_dz(dw_dxi1,dw_dxi2,dw_dxi3,n)

! bx, by, bz
	dbx_dx= dphi_dx(dbx_dxi1,dbx_dxi2,dbx_dxi3,n)
	dby_dy= dphi_dy(dby_dxi1,dby_dxi2,dby_dxi3,n)
	dbz_dz= dphi_dz(dbz_dxi1,dbz_dxi2,dbz_dxi3,n)
	
! ==============================================================================================
!	Dila=du_dx+dv_dy+dw_dz
!	meshes_dila(n)=Dila

	meshes_R_b_p(n)=-(du_dx**2. + dv_dy**2. + dw_dz**2.)-2.0*(du_dy*dv_dx+dv_dz*dw_dy+dw_dx*du_dz)-dbx_dx-dby_dy-dbz_dz+meshes_dila(n)/(2*DT)

end subroutine calculate_R_b_p2


subroutine calculate_gforces
	integer i
	!$OMP     PARALLEL DO  PRIVATE (i)
	do i=1,NMESH
		meshes_bx(i)=0.
		meshes_by(i)=0.
		meshes_bz(i)=-g_beta*0.5*(meshes_Txyz(i)+meshes_Txyz_old(i))
	enddo
     !$OMP     END PARALLEL DO

	!$OMP     PARALLEL DO  PRIVATE (i)
	do i=1,NFACE
		faces_bx(i)=0.
		faces_by(i)=0.
		faces_bz(i)=-g_beta*(faces_T(i))
	enddo
	!$OMP     END PARALLEL DO

end subroutine calculate_gforces

subroutine evaluate_D
	integer                     :: n,k
	real                        :: du_dx, dv_dy, dw_dz, Dila
    integer, dimension(1:6)     :: n_faces
	real                        :: sum1
    real, dimension(1:6)        :: u_faces,v_faces,w_faces
	
	sum1=0.0
	!$OMP PARALLEL DO PRIVATE (n,k,n_faces,u_faces,v_faces,w_faces,du_dx, dv_dy, dw_dz,Dila) REDUCTION(+: sum1)
	do n=1,NMESH
		n_faces(:)=meshes_faces(:,n)
		do k=1,6
			u_faces(k)=faces_u(n_faces(k))
			v_faces(k)=faces_v(n_faces(k))
			w_faces(k)=faces_w(n_faces(k))
		enddo 	
		! ==============================================================================================
		! Derivatives in Cartesian system
		du_dx= dphi_dx2(u_faces,n_faces,n)

		dv_dy= dphi_dy2(v_faces,n_faces,n)

		dw_dz= dphi_dz2(w_faces,n_faces,n)

		! ==============================================================================================
		Dila=du_dx+dv_dy+dw_dz
		meshes_dila(n)=Dila
		
		sum1=sum1+Dila**2

    enddo
    !$OMP     END PARALLEL DO
	Dilatation=sqrt(sum1/NMESH)
end subroutine evaluate_D

subroutine evaluate_D2
	integer                     :: n
	real                        :: du_dxi1,du_dxi2,du_dxi3,dv_dxi1,dv_dxi2,dv_dxi3,dw_dxi1,dw_dxi2,dw_dxi3
	real                        :: du_dx, dv_dy, dw_dz, Dila
    integer, dimension(1:6)     :: n_faces
	real                        :: sum1
	
	sum1=0.0
	!$OMP PARALLEL DO PRIVATE (n,n_faces,du_dxi1,du_dxi2,du_dxi3,dv_dxi1,dv_dxi2,dv_dxi3,dw_dxi1,dw_dxi2,dw_dxi3,du_dx, dv_dy, dw_dz,Dila) REDUCTION(+: sum1)
	do n=1,NMESH
		n_faces(:)=meshes_faces(:,n)
		! ==============================================================================================
		! Derivatives in curvilinear coordinates
		du_dxi1=(faces_u(n_faces(2))-faces_u(n_faces(1)))/2.
		du_dxi2=(faces_u(n_faces(4))-faces_u(n_faces(3)))/2.
		du_dxi3=(faces_u(n_faces(6))-faces_u(n_faces(5)))/2.

		dv_dxi1=(faces_v(n_faces(2))-faces_v(n_faces(1)))/2.
		dv_dxi2=(faces_v(n_faces(4))-faces_v(n_faces(3)))/2.
		dv_dxi3=(faces_v(n_faces(6))-faces_v(n_faces(5)))/2.

		dw_dxi1=(faces_w(n_faces(2))-faces_w(n_faces(1)))/2.
		dw_dxi2=(faces_w(n_faces(4))-faces_w(n_faces(3)))/2.
		dw_dxi3=(faces_w(n_faces(6))-faces_w(n_faces(5)))/2.

		! ==============================================================================================
		! Derivatives in Cartesian system
		du_dx= dphi_dx(du_dxi1,du_dxi2,du_dxi3,n)

		dv_dy= dphi_dy(dv_dxi1,dv_dxi2,dv_dxi3,n)

		dw_dz= dphi_dz(dw_dxi1,dw_dxi2,dw_dxi3,n)

		! ==============================================================================================
		Dila=du_dx+dv_dy+dw_dz
		meshes_dila(n)=Dila
		
		sum1=sum1+Dila**2

    enddo
    !$OMP     END PARALLEL DO
	Dilatation=sqrt(sum1/NMESH)
end subroutine evaluate_D2


subroutine find_Nu
        integer        :: i,j,M, fn1,fn0
        real           :: sum1, sum2,sumA1, sumA2, Nu0, Nu1, dTdx, max_nu0,max_nu1,B3,B4
        integer        ::  f1,f2,f3,f4
! Nu0 : Nu at x=0
! Nu1 : Nu at x=1
        integer, dimension(1:6)              :: n1_faces
        real,dimension(1:4)   :: h

        sum1=0.
        sum2=0.
        sumA1=0.
        sumA2=0.
        max_nu0=0.
        max_nu1=0.
        fn0=1
        fn1=1
        ! Right
    do i=1,tot_faces_bcs
! ==================== |F|A|C|E| |2| =====================
                M=BCs_faces_mesh_num(i)
                n1_faces=meshes_faces(:,M)
                if (BCs_faces_position(i)==1) THEN
                          f1=BCs_faces(i)
                          f2=meshes_faces(2,M)
                          f3=meshes_faces(2,meshes_neighbors(2,M))
                          f4=meshes_faces(2,meshes_neighbors(2,meshes_neighbors(2,M)))

           B3= (-meshes_gij_faces(1,2,1,M))*(faces_T(meshes_faces(4,M))-faces_T(meshes_faces(3,M)))/2.0
           B4= (-meshes_gij_faces(1,3,1,M))*(faces_T(meshes_faces(6,M))-faces_T(meshes_faces(5,M)))/2.0


!                          h(1)=sqrt((faces_mid_coordinates(1,f1)-faces_mid_coordinates(1,f2))**2+(faces_mid_coordinates(2,f1)-faces_mid_coordinates(2,f2))**2+(faces_mid_coordinates(3,f1)-faces_mid_coordinates(3,f2))**2)
!                          h(2)=sqrt((faces_mid_coordinates(1,f2)-faces_mid_coordinates(1,f3))**2+(faces_mid_coordinates(2,f2)-faces_mid_coordinates(2,f3))**2+(faces_mid_coordinates(3,f2)-faces_mid_coordinates(3,f3))**2)
!                          h(3)=sqrt((faces_mid_coordinates(1,f3)-faces_mid_coordinates(1,f4))**2+(faces_mid_coordinates(2,f3)-faces_mid_coordinates(2,f4))**2+(faces_mid_coordinates(3,f3)-faces_mid_coordinates(3,f4))**2)
                          ! h(1)=-11./12.*meshes_xi_x(1,1,Mi,Mj)!/4.
                          ! h(2)=+18./12.*meshes_xi_x(1,1,Mi,Mj)!/4.
                          ! h(3)=-9. /12.*meshes_xi_x(1,1,Mi,Mj)!/4.
                          ! h(4)= 2. /12.*meshes_xi_x(1,1,Mi,Mj)!/4.
!                          dTdx=-(-(1./h(1)+1./(h(1)+h(2))+1./(h(1)+h(2)+h(3)))*faces_T(f1) &
!                                        +((h(1)+h(2))*(h(1)+h(2)+h(3)))/(h(1)*h(2)*(h(2)+h(3)))*faces_T(f2) &
!                                        -(h(1)*(h(1)+h(2)+h(3)))/(h(2)*h(3)*(h(1)+h(2)))*faces_T(f3) &
!                                        +(h(1)*(h(1)+h(2)))/(h(3)*(h(2)+h(3))*(h(1)+h(2)+h(3)))*faces_T(f4))
                          dTdx=BCs_faces_aT(1,i)*faces_T(f1)-BCs_faces_aT(2,i)*faces_T(f2)-BCs_faces_aT(3,i)*(meshes_Txyz(M)+meshes_Txyz_old(M))+B3+B4

                        if(dTdx>max_nu0) then
                                max_nu0=dTdx
                                fn0=i
                        endif
!                       max_nu0=MAX(max_nu0,dTdx)
                        sum2=sum2+faces_area(n1_faces(1))*dTdx
                        sumA2=sumA2+faces_area(n1_faces(1))
                        ! Nu_local(Mj)=dTdx
                        ! Nu_coord(Mj)=faces_mid_coordinates(2,meshes_faces(1,M))


                elseif(BCs_faces_position(i)==2) THEN
                          ! h(1)=11./12.*meshes_xi_x(1,1,Mi,Mj)!/4.
                          ! h(2)=-18./12.*meshes_xi_x(1,1,Mi,Mj)!/4.
                          ! h(3)=+9. /12.*meshes_xi_x(1,1,Mi,Mj)!/4.
                          ! h(4)=-2. /12.*meshes_xi_x(1,1,Mi,Mj)!/4.
!                       max_nu1=MAX(max_nu1,dTdx)
                          f1=BCs_faces(i)
                          f2=meshes_faces(1,M)
                          f3=meshes_faces(1,meshes_neighbors(1,M))
                          f4=meshes_faces(1,meshes_neighbors(1,meshes_neighbors(1,M)))

           B3= (-meshes_gij_faces(1,2,2,M))*(faces_T(meshes_faces(4,M))-faces_T(meshes_faces(3,M)))/2.0
           B4= (-meshes_gij_faces(1,3,2,M))*(faces_T(meshes_faces(6,M))-faces_T(meshes_faces(5,M)))/2.0
!                          h(1)=sqrt((faces_mid_coordinates(1,f1)-faces_mid_coordinates(1,f2))**2+(faces_mid_coordinates(2,f1)-faces_mid_coordinates(2,f2))**2+(faces_mid_coordinates(3,f1)-faces_mid_coordinates(3,f2))**2)
!                          h(2)=sqrt((faces_mid_coordinates(1,f2)-faces_mid_coordinates(1,f3))**2+(faces_mid_coordinates(2,f2)-faces_mid_coordinates(2,f3))**2+(faces_mid_coordinates(3,f2)-faces_mid_coordinates(3,f3))**2)
!                          h(3)=sqrt((faces_mid_coordinates(1,f3)-faces_mid_coordinates(1,f4))**2+(faces_mid_coordinates(2,f3)-faces_mid_coordinates(2,f4))**2+(faces_mid_coordinates(3,f3)-faces_mid_coordinates(3,f4))**2)
!                          dTdx=(-(1./h(1)+1./(h(1)+h(2))+1./(h(1)+h(2)+h(3)))*faces_T(f1) &
!                                        +((h(1)+h(2))*(h(1)+h(2)+h(3)))/(h(1)*h(2)*(h(2)+h(3)))*faces_T(f2) &
!                                        -(h(1)*(h(1)+h(2)+h(3)))/(h(2)*h(3)*(h(1)+h(2)))*faces_T(f3) &
!                                        +(h(1)*(h(1)+h(2)))/(h(3)*(h(2)+h(3))*(h(1)+h(2)+h(3)))*faces_T(f4))
                          dTdx=BCs_faces_aT(1,i)*faces_T(f1)-BCs_faces_aT(2,i)*faces_T(f2)-BCs_faces_aT(3,i)*(meshes_Txyz(M)+meshes_Txyz_old(M))+B3+B4

                        if(dTdx>max_nu1) then
                                max_nu1=dTdx
                                fn1=i
                        endif
                        sum1=sum1+faces_area(n1_faces(2))*dTdx
                        sumA1=sumA1+faces_area(n1_faces(2))
                        !Nu_local(Mj)=dTdx
                        !Nu_coord(Mj)=1-faces_mid_coordinates(2,meshes_edges(1,M))*sqrt(2.)

!                       print*, Nu_coord(Mj)!,Nu_local(Mj)

                elseif(BCs_faces_position(i)==3) THEN

                elseif(BCs_faces_position(i)==4) THEN

                endif
    enddo

        ! do j=1,NMESH_j
                ! n1_edges(:)=meshes_edges(:,i,j)
                ! dTdx=(MT_F12(i,j)+2*MT_F14(i,j))*edges_T(n1_edges(2))-MT_F12(i,j)*edges_T(n1_edges(1))-MT_F14(i,j)*(meshes_Txy(i,j)+meshes_Txy_old(i,j))
                ! PRINT*, (MT_F12(i,j)),(2*MT_F14(i,j)),edges_T(n1_edges(2)),edges_T(n1_edges(1)), dTdx
                ! sum1=sum1+edges_length(n1_edges(2))*dTdx
        ! enddo



        print*, "================================================================================="
        write (*,'(A, 1p,E11.4, A, 1p,E11.4, A, 1p,E11.4, A, 1p,E11.4, A, 1p,E11.4)') &
        'Nu_H=', -sum1/sumA1*R_o ,' |  Nu_C=', -sum2/sumA2*R_i, '   |    Average = ',  -((sum1/sumA1*R_o)+(sum2/sumA2*R_i))/2., '   |    Area1 = ',sumA1, '   |    Area2 = ',sumA2
!       'Nu_H=', -sum1 ,' |  Nu_C=', -sum2, '   |    Average = ',  -((sum1/sumA1*R_o)+(sum2/sumA2*R_i))/2., '   |    Area1 = ',sumA1, '   |    Area2 = ',sumA2
        print*, "================================================================================="

end subroutine find_Nu



subroutine find_Nu2
	integer        :: i,j,M, fn1,fn0
	real           :: sum1, sum2,sumA1, sumA2, Nu0, Nu1, dTdx, max_nu0,max_nu1
	integer        ::  f1,f2,f3,f4
! Nu0 : Nu at x=0
! Nu1 : Nu at x=1
	integer, dimension(1:6)              :: n1_faces
	real,dimension(1:4)   :: h

	sum1=0.
	sum2=0.
	sumA1=0.
	sumA2=0.
	max_nu0=0.
	max_nu1=0.
	fn0=1
	fn1=1
	! Right
    do i=1,tot_faces_bcs
! ==================== |F|A|C|E| |2| =====================
		M=BCs_faces_mesh_num(i)
		n1_faces=meshes_faces(:,M)
		if (BCs_faces_position(i)==1) THEN
			  f1=BCs_faces(i)
			  f2=meshes_faces(2,M)
			  f3=meshes_faces(2,meshes_neighbors(2,M))
			  f4=meshes_faces(2,meshes_neighbors(2,meshes_neighbors(2,M)))
			  h(1)=sqrt((faces_mid_coordinates(1,f1)-faces_mid_coordinates(1,f2))**2+(faces_mid_coordinates(2,f1)-faces_mid_coordinates(2,f2))**2+(faces_mid_coordinates(3,f1)-faces_mid_coordinates(3,f2))**2)
			  h(2)=sqrt((faces_mid_coordinates(1,f2)-faces_mid_coordinates(1,f3))**2+(faces_mid_coordinates(2,f2)-faces_mid_coordinates(2,f3))**2+(faces_mid_coordinates(3,f2)-faces_mid_coordinates(3,f3))**2)
			  h(3)=sqrt((faces_mid_coordinates(1,f3)-faces_mid_coordinates(1,f4))**2+(faces_mid_coordinates(2,f3)-faces_mid_coordinates(2,f4))**2+(faces_mid_coordinates(3,f3)-faces_mid_coordinates(3,f4))**2)
			  ! h(1)=-11./12.*meshes_xi_x(1,1,Mi,Mj)!/4.
			  ! h(2)=+18./12.*meshes_xi_x(1,1,Mi,Mj)!/4.
			  ! h(3)=-9. /12.*meshes_xi_x(1,1,Mi,Mj)!/4.
			  ! h(4)= 2. /12.*meshes_xi_x(1,1,Mi,Mj)!/4.
			  dTdx=-(-(1./h(1)+1./(h(1)+h(2))+1./(h(1)+h(2)+h(3)))*faces_T(f1) &
					+((h(1)+h(2))*(h(1)+h(2)+h(3)))/(h(1)*h(2)*(h(2)+h(3)))*faces_T(f2) &
					-(h(1)*(h(1)+h(2)+h(3)))/(h(2)*h(3)*(h(1)+h(2)))*faces_T(f3) &
					+(h(1)*(h(1)+h(2)))/(h(3)*(h(2)+h(3))*(h(1)+h(2)+h(3)))*faces_T(f4))
			if(dTdx>max_nu0) then
				max_nu0=dTdx
				fn0=i
			endif
!			max_nu0=MAX(max_nu0,dTdx)
			sum2=sum2+faces_area(n1_faces(1))*dTdx
			sumA2=sumA2+faces_area(n1_faces(1))
			! Nu_local(Mj)=dTdx
			! Nu_coord(Mj)=faces_mid_coordinates(2,meshes_faces(1,M))

			
		elseif(BCs_faces_position(i)==2) THEN
			  ! h(1)=11./12.*meshes_xi_x(1,1,Mi,Mj)!/4.
			  ! h(2)=-18./12.*meshes_xi_x(1,1,Mi,Mj)!/4.
			  ! h(3)=+9. /12.*meshes_xi_x(1,1,Mi,Mj)!/4.
			  ! h(4)=-2. /12.*meshes_xi_x(1,1,Mi,Mj)!/4.
!			max_nu1=MAX(max_nu1,dTdx)
			  f1=BCs_faces(i)
			  f2=meshes_faces(1,M)
			  f3=meshes_faces(1,meshes_neighbors(1,M))
			  f4=meshes_faces(1,meshes_neighbors(1,meshes_neighbors(1,M)))
			  h(1)=sqrt((faces_mid_coordinates(1,f1)-faces_mid_coordinates(1,f2))**2+(faces_mid_coordinates(2,f1)-faces_mid_coordinates(2,f2))**2+(faces_mid_coordinates(3,f1)-faces_mid_coordinates(3,f2))**2)
			  h(2)=sqrt((faces_mid_coordinates(1,f2)-faces_mid_coordinates(1,f3))**2+(faces_mid_coordinates(2,f2)-faces_mid_coordinates(2,f3))**2+(faces_mid_coordinates(3,f2)-faces_mid_coordinates(3,f3))**2)
			  h(3)=sqrt((faces_mid_coordinates(1,f3)-faces_mid_coordinates(1,f4))**2+(faces_mid_coordinates(2,f3)-faces_mid_coordinates(2,f4))**2+(faces_mid_coordinates(3,f3)-faces_mid_coordinates(3,f4))**2)
			  dTdx=(-(1./h(1)+1./(h(1)+h(2))+1./(h(1)+h(2)+h(3)))*faces_T(f1) &
					+((h(1)+h(2))*(h(1)+h(2)+h(3)))/(h(1)*h(2)*(h(2)+h(3)))*faces_T(f2) &
					-(h(1)*(h(1)+h(2)+h(3)))/(h(2)*h(3)*(h(1)+h(2)))*faces_T(f3) &
					+(h(1)*(h(1)+h(2)))/(h(3)*(h(2)+h(3))*(h(1)+h(2)+h(3)))*faces_T(f4))
			if(dTdx>max_nu1) then
				max_nu1=dTdx
				fn1=i
			endif
			sum1=sum1+faces_area(n1_faces(2))*dTdx
			sumA1=sumA1+faces_area(n1_faces(2))
			!Nu_local(Mj)=dTdx
			!Nu_coord(Mj)=1-faces_mid_coordinates(2,meshes_edges(1,M))*sqrt(2.)

!			print*, Nu_coord(Mj)!,Nu_local(Mj)
			
		elseif(BCs_faces_position(i)==3) THEN
		
		elseif(BCs_faces_position(i)==4) THEN

		endif
    enddo

	! do j=1,NMESH_j
		! n1_edges(:)=meshes_edges(:,i,j)
		! dTdx=(MT_F12(i,j)+2*MT_F14(i,j))*edges_T(n1_edges(2))-MT_F12(i,j)*edges_T(n1_edges(1))-MT_F14(i,j)*(meshes_Txy(i,j)+meshes_Txy_old(i,j))
		! PRINT*, (MT_F12(i,j)),(2*MT_F14(i,j)),edges_T(n1_edges(2)),edges_T(n1_edges(1)), dTdx
		! sum1=sum1+edges_length(n1_edges(2))*dTdx
	! enddo



	print*, "================================================================================="
	write (*,'(A, 1p,E11.4, A, 1p,E11.4, A, 1p,E11.4, A, 1p,E11.4, A, 1p,E11.4)') & 
	'Nu_H=', -sum1/sumA1*R_o ,' |  Nu_C=', -sum2/sumA2*R_i, '   |    Average = ',  -((sum1/sumA1*R_o)+(sum2/sumA2*R_i))/2., '   |    Area1 = ',sumA1, '   |    Area2 = ',sumA2
!	'Nu_H=', -sum1 ,' |  Nu_C=', -sum2, '   |    Average = ',  -((sum1/sumA1*R_o)+(sum2/sumA2*R_i))/2., '   |    Area1 = ',sumA1, '   |    Area2 = ',sumA2
	print*, "================================================================================="

end subroutine find_Nu2


subroutine find_rms
	integer        :: i
	real           :: sum1, sum2, sum3
	
	sum1=0.
	sum2=0.
	sum3=0.
	!$OMP     PARALLEL DO  PRIVATE (i) REDUCTION(+: sum1,sum2,sum3)

	do i=1,NMESH
		sum1=sum1+abs(meshes_uT(i)-meshes_uxyz(i))**2.
		sum2=sum2+abs(meshes_vT(i)-meshes_vxyz(i))**2.
		sum3=sum3+abs(meshes_wT(i)-meshes_wxyz(i))**2.
	enddo

    !$OMP     END PARALLEL DO

	RMS_meshes_u=sqrt(sum1/(NMESH))
	RMS_meshes_v=sqrt(sum2/(NMESH))
	RMS_meshes_w=sqrt(sum3/(NMESH))
	sum1=0.
	sum2=0.
	sum3=0.
	!$OMP     PARALLEL DO  PRIVATE (i) REDUCTION(+: sum1,sum2,sum3)
	do i=1,NFACE
		sum1=sum1+abs(faces_uT(i)-faces_u(i))**2.
		sum2=sum2+abs(faces_vT(i)-faces_v(i))**2.
		sum3=sum3+abs(faces_wT(i)-faces_w(i))**2.
	enddo
    !$OMP     END PARALLEL DO
	RMS_faces_u=sqrt(sum1/NFACE)
	RMS_faces_v=sqrt(sum2/NFACE)
	RMS_faces_w=sqrt(sum3/NFACE)
	print*, "================================================================================="
	write (*,'(A, 1p,E11.4, A, 1p,E11.4, A, 1p,E11.4, A, 1p,E11.4, A, 1p,E11.4,A, 1p,E11.4)') & 
	'RMS : uxy=', RMS_meshes_u , ' | vxy=', RMS_meshes_v, ' | wxy=', RMS_meshes_w, ' | faces u=', RMS_faces_u, ' | faces v=', RMS_faces_v, ' | faces w=', RMS_faces_w
	print*, "================================================================================="

end subroutine find_rms


! This function is used to find dphi_dx in cartesian system from the derivatives in curvilinear coordinates
function dphi_dx(dphi_dxi1,dphi_dxi2,dphi_dxi3,n) result(T)
	real, intent(in)    :: dphi_dxi1,dphi_dxi2,dphi_dxi3
	integer, intent(in) :: n 							! mesh umber
	real                :: T
    T= meshes_xi_x(1,1,n)*dphi_dxi1+meshes_xi_x(2,1,n)*dphi_dxi2+meshes_xi_x(3,1,n)*dphi_dxi3
end function dphi_dx

function dphi_dy(dphi_dxi1,dphi_dxi2,dphi_dxi3,n) result(T)
	real, intent(in)    :: dphi_dxi1,dphi_dxi2,dphi_dxi3
	integer, intent(in) :: n 							! mesh umber
	real                :: T
    T= meshes_xi_x(1,2,n)*dphi_dxi1+meshes_xi_x(2,2,n)*dphi_dxi2+meshes_xi_x(3,2,n)*dphi_dxi3
end function dphi_dy

function dphi_dz(dphi_dxi1,dphi_dxi2,dphi_dxi3,n) result(T)
	real, intent(in)    :: dphi_dxi1,dphi_dxi2,dphi_dxi3
	integer, intent(in) :: n 							! mesh umber
	real                :: T
    T= meshes_xi_x(1,3,n)*dphi_dxi1+meshes_xi_x(2,3,n)*dphi_dxi2+meshes_xi_x(3,3,n)*dphi_dxi3
end function dphi_dz


! To find first derivative with respect to x in cartesian coordinates using the corner points of a mesh
function dphi_dx2(phi,n_faces,i) result(T)
	real,dimension(1:6),intent(in)   :: phi
	integer,dimension(1:6),intent(in):: n_faces
	integer, intent(in)           :: i 							! mesh umber
	real                          :: T
    T= 1./meshes_volume(i)*(phi(2)*faces_area(n_faces(2))*(faces_normal(1,n_faces(2))) &
	                       -phi(1)*faces_area(n_faces(1))*(faces_normal(1,n_faces(1))) &
	                       +phi(4)*faces_area(n_faces(4))*(faces_normal(1,n_faces(4))) &
	                       -phi(3)*faces_area(n_faces(3))*(faces_normal(1,n_faces(3))) &
	                       +phi(6)*faces_area(n_faces(6))*(faces_normal(1,n_faces(6))) &
	                       -phi(5)*faces_area(n_faces(5))*(faces_normal(1,n_faces(5))))
end function dphi_dx2

function dphi_dy2(phi,n_faces,i) result(T)
	real,dimension(1:6),intent(in)   :: phi
	integer,dimension(1:6),intent(in):: n_faces
	integer, intent(in)           :: i 							! mesh umber
	real                          :: T
    T= 1./meshes_volume(i)*(phi(2)*faces_area(n_faces(2))*(faces_normal(2,n_faces(2))) &
	                       -phi(1)*faces_area(n_faces(1))*(faces_normal(2,n_faces(1))) &
	                       +phi(4)*faces_area(n_faces(4))*(faces_normal(2,n_faces(4))) &
	                       -phi(3)*faces_area(n_faces(3))*(faces_normal(2,n_faces(3))) &
	                       +phi(6)*faces_area(n_faces(6))*(faces_normal(2,n_faces(6))) &
	                       -phi(5)*faces_area(n_faces(5))*(faces_normal(2,n_faces(5))))
end function dphi_dy2

function dphi_dz2(phi,n_faces,i) result(T)
	real,dimension(1:6),intent(in)   :: phi
	integer,dimension(1:6),intent(in):: n_faces
	integer, intent(in)           :: i 							! mesh umber
	real                          :: T
    T= 1./meshes_volume(i)*(phi(2)*faces_area(n_faces(2))*(faces_normal(3,n_faces(2))) &
	                       -phi(1)*faces_area(n_faces(1))*(faces_normal(3,n_faces(1))) &
	                       +phi(4)*faces_area(n_faces(4))*(faces_normal(3,n_faces(4))) &
	                       -phi(3)*faces_area(n_faces(3))*(faces_normal(3,n_faces(3))) &
	                       +phi(6)*faces_area(n_faces(6))*(faces_normal(3,n_faces(6))) &
	                       -phi(5)*faces_area(n_faces(5))*(faces_normal(3,n_faces(5))))
end function dphi_dz2
end module class_nse

