module class_general
  use class_mesh
  use class_parameters
!  implicit none

contains

  function analytical_xy (x,y,z,t,id) result (TT)
    real, intent (in)            :: x,y,z,t
	integer, intent(in)          :: id
	real                         :: TT
	TT=0.

  end function analytical_xy

  function initial_condition_xyz (x,y,z,id) result (TT)
    real, intent (in)            :: x,y,z
	integer, intent(in)          :: id
    real                         :: TT
    if (id==velx1) then
       TT=0.
    elseif(id==velx2) then
       TT=0.
    elseif(id==velx3) then
       TT=0.
    elseif(id==PPE) then
       TT=0.
    elseif(id==CDE) then
       TT=0.
    endif
  end function initial_condition_xyz

  function boundary_condition_1 (x,y,z,id) result (TT)
    real, intent (in)            :: x,y,z
	integer, intent(in)          :: id
    real                         :: TT
    if (id==velx1) then
       TT=0.
    elseif(id==velx2) then
       TT=0.
    elseif(id==velx3) then
       TT=0.
    elseif(id==PPE) then
       TT=0.
    elseif(id==CDE) then
       TT=0.
    endif
  end function boundary_condition_1

  function boundary_condition_2 (x,y,z,id) result (TT)
    real, intent (in)            :: x,y,z
	integer, intent(in)          :: id
    real                         :: TT
    if (id==velx1) then
       TT=0.
    elseif(id==velx2) then
       TT=0.
    elseif(id==velx3) then
       TT=0.
    elseif(id==PPE) then
       TT=0.
    elseif(id==CDE) then
       TT=0.
    endif
  end function boundary_condition_2

  function boundary_condition_3 (x,y,z,id) result (TT)
    real, intent (in)            :: x,y,z
	integer, intent(in)          :: id
    real                         :: TT
    if (id==velx1) then
       TT=0.
    elseif(id==velx2) then
       TT=0.
    elseif(id==velx3) then
       TT=0.
    elseif(id==PPE) then
       TT=0.
    elseif(id==CDE) then
       TT=0.
    endif
  end function boundary_condition_3


  function boundary_condition_4 (x,y,z,id) result (TT)
    real, intent (in)            :: x,y,z
	integer, intent(in)          :: id
    real                         :: TT
    if (id==velx1) then
       TT=0.
    elseif(id==velx2) then
       TT=0.
    elseif(id==velx3) then
       TT=0.
    elseif(id==PPE) then
       TT=0.
    elseif(id==CDE) then
       TT=1.
    endif
  end function boundary_condition_4

  function boundary_condition_5 (x,y,z,id) result (TT)
    real, intent (in)            :: x,y,z
	integer, intent(in)          :: id
    real                         :: TT
    if (id==velx1) then
       TT=0.
    elseif(id==velx2) then
       TT=0.
    elseif(id==velx3) then
       TT=0.
    elseif(id==PPE) then
       TT=0.
    elseif(id==CDE) then
       TT=0.
    endif
  end function boundary_condition_5


  function boundary_condition_6 (x,y,z,id) result (TT)
    real, intent (in)            :: x,y,z
	integer, intent(in)          :: id
    real                         :: TT
    if (id==velx1) then
       TT=0.
    elseif(id==velx2) then
       TT=0.
    elseif(id==velx3) then
       TT=0.
    elseif(id==PPE) then
       TT=0.
    elseif(id==CDE) then
       TT=0.
    endif
  end function boundary_condition_6

  
  function boundary_condition_7 (x,y,z,id) result (TT)
    real, intent (in)            :: x,y,z
	integer, intent(in)          :: id
    real                         :: TT
    if (id==velx1) then
       TT=0.
    elseif(id==velx2) then
       TT=0.
    elseif(id==velx3) then
       TT=0.
    elseif(id==PPE) then
       TT=0.
    elseif(id==CDE) then
       TT=0.
    endif
  end function boundary_condition_7


  function gforce_x (x,y,z,t) result (S)
    real, intent (in)            :: x,y,z,t
    real                         :: S
    S=0.

  end function gforce_x
  
  function gforce_y (x,y,z,t) result (S)
    real, intent (in)            :: x,y,z,t
    real                         :: S
    S=0.

  end function gforce_y
  
  function gforce_z (x,y,z,t) result (S)
    real, intent (in)            :: x,y,z,t
    real                         :: S
    S=0.

  end function gforce_z
  
  function gforce_q (x,y,z,t) result (S)
    real, intent (in)            :: x,y,z,t
    real                         :: S
    S=0.

  end function gforce_q
  

  function Diffusion_coefficient (x,y,z,t) result (Diff)
    real, intent (in)            :: x,y,z
    real, optional               :: t
    real                         :: Diff
    Diff=1.0/Re
  end function Diffusion_coefficient



  function functions_id (x,y,z,t,id) result (fun)
    real, intent (in)            :: x,y,z
    real, optional               :: t
	integer, intent(in)          :: id
    real                         :: fun

    SELECT CASE (id)
      CASE (analytical_velx1)
           fun=analytical_xy (x,y,z,t,velx1)
      CASE (analytical_velx2)
           fun=analytical_xy (x,y,z,t,velx2)
      CASE (analytical_velx3)
           fun=analytical_xy (x,y,z,t,velx3)
      CASE (analytical_p)
           fun=analytical_xy (x,y,z,t,PPE)
      CASE (analytical_CDE)
           fun=analytical_xy (x,y,z,t,CDE)
      CASE (IC_velx1)
           fun=initial_condition_xyz(x,y,z,velx1)
      CASE (IC_velx2)
           fun=initial_condition_xyz(x,y,z,velx2)
      CASE (IC_velx3)
           fun=initial_condition_xyz(x,y,z,velx3)
      CASE (IC_PPE)
           fun=initial_condition_xyz(x,y,z,PPE)
      CASE (IC_CDE)
           fun=initial_condition_xyz(x,y,z,CDE)
	  CASE (bx_fun)
           fun=gforce_x(x,y,z,t)
      CASE (by_fun)
           fun=gforce_y(x,y,z,t)
      CASE (bz_fun)
           fun=gforce_z(x,y,z,t)
      CASE (q_fun)
           fun=gforce_q(x,y,z,t)
      CASE DEFAULT
           print*, "ERROR : Unknown function call with id : ", id
    END SELECT

  end function functions_id


  function BCs_id (x,y,z,field_id,bc_id) result (fun)
    real, intent (in)            :: x,y,z
	integer, intent(in)          :: field_id, bc_id
    real                         :: fun

    SELECT CASE (bc_id)
      CASE (1)
           fun=boundary_condition_1(x,y,z,field_id)
      CASE (2)
           fun=boundary_condition_2(x,y,z,field_id)
      CASE (3)
           fun=boundary_condition_3(x,y,z,field_id)
      CASE (4)
           fun=boundary_condition_4(x,y,z,field_id)
      CASE (5)
           fun=boundary_condition_5(x,y,z,field_id)
      CASE (6)
           fun=boundary_condition_6(x,y,z,field_id)
      CASE (7)
           fun=boundary_condition_7(x,y,z,field_id)
      CASE DEFAULT
           print*, "ERROR : Unknown function call with id : ", bc_id
    END SELECT

  end function BCs_id

  function BCs_derivatives (x,y,z,bc_id,dvdt_id) result (fun)
    real, intent (in)            :: x,y,z
	integer, intent(in)          :: dvdt_id, bc_id
    real                         :: fun

    SELECT CASE (bc_id)
      CASE (1)
           fun=BCs_derivatives_1(x,y,z,dvdt_id)
      CASE (2)
           fun=BCs_derivatives_2(x,y,z,dvdt_id)
      CASE (3)
           fun=BCs_derivatives_3(x,y,z,dvdt_id)
      CASE (4)
           fun=BCs_derivatives_4(x,y,z,dvdt_id)
      CASE (5)
           fun=BCs_derivatives_5(x,y,z,dvdt_id)
      CASE (6)
           fun=BCs_derivatives_6(x,y,z,dvdt_id)
      CASE (7)
           fun=BCs_derivatives_7(x,y,z,dvdt_id)
      CASE DEFAULT
           print*, "ERROR : Unknown function call (BCs_derivatives) with id : ", bc_id
    END SELECT

  end function BCs_derivatives


  function BCs_derivatives_1 (x,y,z,id) result (TT)
    real, intent (in)            :: x,y,z
	integer, intent(in)          :: id
    real                         :: TT
    if (id==bc_dvdn) then
       TT=0.
    elseif(id==bc_dvdnt) then
       TT=0.
    elseif(id==bc_dvdns) then
       TT=0.
    elseif(id==bc_dvdntt) then
       TT=0.
    elseif(id==bc_dvdnss) then
       TT=0.
    elseif(id==bc_dvdnts) then
       TT=0.
    endif
  end function BCs_derivatives_1

  function BCs_derivatives_2 (x,y,z,id) result (TT)
    real, intent (in)            :: x,y,z
	integer, intent(in)          :: id
    real                         :: TT
    if (id==bc_dvdn) then
       TT=0.
    elseif(id==bc_dvdnt) then
       TT=0.
    elseif(id==bc_dvdns) then
       TT=0.
    elseif(id==bc_dvdntt) then
       TT=0.
    elseif(id==bc_dvdnss) then
       TT=0.
    elseif(id==bc_dvdnts) then
       TT=0.
    endif
  end function BCs_derivatives_2

  function BCs_derivatives_3 (x,y,z,id) result (TT)
    real, intent (in)            :: x,y,z
	integer, intent(in)          :: id
    real                         :: TT
    if (id==bc_dvdn) then
       TT=0.
    elseif(id==bc_dvdnt) then
       TT=0.
    elseif(id==bc_dvdns) then
       TT=0.
    elseif(id==bc_dvdntt) then
       TT=0.
    elseif(id==bc_dvdnss) then
       TT=0.
    elseif(id==bc_dvdnts) then
       TT=0.
    endif
  end function BCs_derivatives_3

  function BCs_derivatives_4 (x,y,z,id) result (TT)
    real, intent (in)            :: x,y,z
	integer, intent(in)          :: id
    real                         :: TT
    if (id==bc_dvdn) then
       TT=0.
    elseif(id==bc_dvdnt) then
       TT=0.
    elseif(id==bc_dvdns) then
       TT=0.
    elseif(id==bc_dvdntt) then
       TT=0.
    elseif(id==bc_dvdnss) then
       TT=0.
    elseif(id==bc_dvdnts) then
       TT=0.
    endif
  end function BCs_derivatives_4

  function BCs_derivatives_5 (x,y,z,id) result (TT)
    real, intent (in)            :: x,y,z
	integer, intent(in)          :: id
    real                         :: TT
    if (id==bc_dvdn) then
       TT=0.
    elseif(id==bc_dvdnt) then
       TT=0.
    elseif(id==bc_dvdns) then
       TT=0.
    elseif(id==bc_dvdntt) then
       TT=0.
    elseif(id==bc_dvdnss) then
       TT=0.
    elseif(id==bc_dvdnts) then
       TT=0.
    endif
  end function BCs_derivatives_5
  
  function BCs_derivatives_6 (x,y,z,id) result (TT)
    real, intent (in)            :: x,y,z
	integer, intent(in)          :: id
    real                         :: TT
    if (id==bc_dvdn) then
       TT=0.
    elseif(id==bc_dvdnt) then
       TT=0.
    elseif(id==bc_dvdns) then
       TT=0.
    elseif(id==bc_dvdntt) then
       TT=0.
    elseif(id==bc_dvdnss) then
       TT=0.
    elseif(id==bc_dvdnts) then
       TT=0.
    endif
  end function BCs_derivatives_6

  function BCs_derivatives_7 (x,y,z,id) result (TT)
    real, intent (in)            :: x,y,z
	integer, intent(in)          :: id
    real                         :: TT
    if (id==bc_dvdn) then
       TT=0.
    elseif(id==bc_dvdnt) then
       TT=0.
    elseif(id==bc_dvdns) then
       TT=0.
    elseif(id==bc_dvdntt) then
       TT=0.
    elseif(id==bc_dvdnss) then
       TT=0.
    elseif(id==bc_dvdnts) then
       TT=0.
    endif
  end function BCs_derivatives_7
 
! =================================================================================
! =================================================================================
! ============================ NODE DEFORMATION FUNCTIONS =========================
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 


subroutine deform_nodes
  implicit none
   call refine_walls_radial_order(1.5)

!   call create_square_from_annulus_order_n_flipped_three_regions(1.25,1.75,-1.)
   call create_hexagon_from_annulus_order_n_flipped_three_regions(1.25,1.75,-1.)
!   call create_hexagon_from_annulus_order_n_flipped_three_regions(1.1,1.9,0.25)





!  call geometry_dilation(10.,10.,10.)
!  call move_geometry(-.5,-.5,-.5)


!   call create_square_from_annulus_order_n_flipped_mine(R_i,R_o,0.25)
!   call create_square_from_annulus_order_n_flipped_three_regions(R_i,R_o,0.25)
!   call create_hexagon_from_annulus_order_n_flipped_mine(R_i,R_o,0.25)



!   call create_square_from_annulus_order_n_flipped(R_i,R_o,0.25)


!  call move_geometry(0.0,0.0,0.1)
!  call h2_perturb
!  call cube_to_circle
!  call geometry_dilation(.14142,.14142,.1)
!  call geometry_dilation(1.,1.,1.)
  
!  call create_cylinder_from_cube
end subroutine deform_nodes


subroutine refine_walls_radial_order(b1)
! The cube should be 1x1x1 centered at 0,0,0
  real, intent(in)                           :: b1
  real                                       :: x,y,z,r,theta,Roi
  integer i,j,k
  do i=1,NNODE
         x=NODES_ARRAY(i)%x
         y=NODES_ARRAY(i)%y
                 r=sqrt(x**2+y**2)-R_i
                Roi=R_o-R_i
                theta=atan(y/(x))
                if(x<0) theta=theta+MATH_PI
                r=r/Roi
         if(r<=0.5 .and. r>0.) r=((2*r)**b1)/2.
         if(r>0.5  .and. r<1.)  r=1.-((-2.*(r-1))**b1)/2.
                r=(r*Roi+R_i)
                NODES_ARRAY(i)%x=r*cos(theta)
                NODES_ARRAY(i)%y=r*sin(theta)

  enddo

end subroutine refine_walls_radial_order


subroutine refine_walls_radial_order2(b1)
! The cube should be 1x1x1 centered at 0,0,0
  real, intent(in)                           :: b1
  real                                       :: x,y,z,r,theta
  integer i,j,k
  do i=1,NNODE
         x=NODES_ARRAY(i)%x
         y=NODES_ARRAY(i)%y
		 r=sqrt(x**2+y**2)-R_i
		theta=atan(y/(x))
		if(x<0) theta=theta+MATH_PI
         if(r<=0.5 .and. r>0.) r=((2*r)**b1)/2.
         if(r>0.5 .and. r<1.)  r=1.-((-2.*(r-1))**b1)/2.
		r=r+R_i
		NODES_ARRAY(i)%x=r*cos(theta)
		NODES_ARRAY(i)%y=r*sin(theta)

  enddo

end subroutine refine_walls_radial_order2



subroutine create_cylinder_from_cube
  implicit none

  call geometry_dilation(20.,20.,20.)
  call move_geometry(-1.0,-1.0,-1.0)
!  call geometry_dilation(sqrt(2.),sqrt(2.),1.)
!  call geometry_dilation(1.414,1.414,1.)


  call cube_to_circle2(sqrt(2.))
!  call geometry_dilation(sqrt(2.),sqrt(2.),1.)

!  call geometry_dilation(1.,1.,1.)


end subroutine create_cylinder_from_cube


subroutine geometry_dilation(xx,xy,xz)
  implicit none
  real, intent(in) :: xx,xy,xz
  integer          :: i
  do i=1,NNODE
    NODES_ARRAY(i)%x=NODES_ARRAY(i)%x*xx
    NODES_ARRAY(i)%y=NODES_ARRAY(i)%y*xy
    NODES_ARRAY(i)%z=NODES_ARRAY(i)%z*xz
  enddo
end subroutine geometry_dilation


subroutine cube_to_circle2(R)
  implicit none
  real, intent(in) :: R
  real             :: x,y,z
  integer          :: i
  do i=1,NNODE
    x=NODES_ARRAY(i)%x
    y=NODES_ARRAY(i)%y
    z=NODES_ARRAY(i)%z

     NODES_ARRAY(i)%x=R*x*sqrt(1-(y**2/2))
     NODES_ARRAY(i)%y=R*y*sqrt(1-(x**2/2))
!     NODES_ARRAY(i)%z=z+0.05*sin(2*PI*x)*sin(2*PI*y)*sin(2*PI*z)
  enddo
end subroutine cube_to_circle2

subroutine create_square_from_annulus_order_n(Ri,Ro,n)
  implicit none
  real, intent(in) :: Ri,Ro, n
  integer          :: i
  real             :: R,theta,x,y, x1,y1,M, p1,p2
  do i=1,NNODE
    x=NODES_ARRAY(i)%x
    y=NODES_ARRAY(i)%y
	theta=atan(y/(x))
	if(x<0) theta=theta+MATH_PI
	R=sqrt(x**2+y**2)
	x1=cos(theta)
	y1=sin(theta)
	M=max(abs(x1),abs(y1))
	p1=R*x1/M
	p2=R*y1/M
	
	
     NODES_ARRAY(i)%x=(x-p1)/(Ri**n-Ro**n)*R**n-(x*Ro**n-p1*Ri**n)/(Ri**n-Ro**n)
     NODES_ARRAY(i)%y=(y-p2)/(Ri**n-Ro**n)*R**n-(y*Ro**n-p2*Ri**n)/(Ri**n-Ro**n)
!     NODES_ARRAY(i)%z=z+0.05*sin(2*PI*x)*sin(2*PI*y)*sin(2*PI*z)
  enddo
end subroutine create_square_from_annulus_order_n


subroutine create_square_from_annulus_order_n_flipped(Ri,Ro,n)
  implicit none
  real, intent(in) :: Ri,Ro, n
  integer          :: i
  real             :: R,theta,x,y, x1,y1,M, p1,p2
  do i=1,NNODE
    x=NODES_ARRAY(i)%x
    y=NODES_ARRAY(i)%y
	theta=atan(y/(x))
	if(x<0) theta=theta+MATH_PI
	R=sqrt(x**2+y**2)
	x1=cos(theta)
	y1=sin(theta)
	M=max(abs(x1),abs(y1))
	p1=R*x1/M
	p2=R*y1/M
	
	
     NODES_ARRAY(i)%x=-(x-p1)/(Ri**n-Ro**n)*R**n-(p1*Ro**n-x*Ri**n)/(Ri**n-Ro**n)
     NODES_ARRAY(i)%y=-(y-p2)/(Ri**n-Ro**n)*R**n-(p2*Ro**n-y*Ri**n)/(Ri**n-Ro**n)
!     NODES_ARRAY(i)%z=z+0.05*sin(2*PI*x)*sin(2*PI*y)*sin(2*PI*z)
  enddo
end subroutine create_square_from_annulus_order_n_flipped


subroutine create_square_from_annulus_order_n_flipped_three_regions(Ri,Ro,n)
  implicit none
  real, intent(in) :: Ri,Ro, n
  integer          :: i
  real             :: R,theta,x,y, x1,y1,M, p1,p2, Reg1,Reg2
  do i=1,NNODE
    x=NODES_ARRAY(i)%x
    y=NODES_ARRAY(i)%y
	theta=atan(y/(x))
	if(x<0) theta=theta+MATH_PI
	R=sqrt(x**2+y**2)
	call convert_circle_to_square_mine(x,y,p1,p2)
	if(R<=Ri) then
		NODES_ARRAY(i)%x=p1
		NODES_ARRAY(i)%y=p2
	elseif(R>Ri .and. R<Ro) then
     NODES_ARRAY(i)%x=-(x-p1)/(Ri**n-Ro**n)*R**n-(p1*Ro**n-x*Ri**n)/(Ri**n-Ro**n)
     NODES_ARRAY(i)%y=-(y-p2)/(Ri**n-Ro**n)*R**n-(p2*Ro**n-y*Ri**n)/(Ri**n-Ro**n)
	endif

!     NODES_ARRAY(i)%z=z+0.05*sin(2*PI*x)*sin(2*PI*y)*sin(2*PI*z)
  enddo
end subroutine create_square_from_annulus_order_n_flipped_three_regions


subroutine create_square_from_annulus_order_n_flipped_mine(Ri,Ro,n)
  implicit none
  real, intent(in) :: Ri,Ro, n
  integer          :: i
  real             :: R,theta,x,y, x1,y1,M, p1,p2
  do i=1,NNODE
    x=NODES_ARRAY(i)%x
    y=NODES_ARRAY(i)%y
	theta=atan(y/(x))
	if(x<0) theta=theta+MATH_PI
	R=sqrt(x**2+y**2)
	call convert_circle_to_square_mine(x,y,p1,p2)
	
     NODES_ARRAY(i)%x=-(x-p1)/(Ri**n-Ro**n)*R**n-(p1*Ro**n-x*Ri**n)/(Ri**n-Ro**n)
     NODES_ARRAY(i)%y=-(y-p2)/(Ri**n-Ro**n)*R**n-(p2*Ro**n-y*Ri**n)/(Ri**n-Ro**n)
!     NODES_ARRAY(i)%z=z+0.05*sin(2*PI*x)*sin(2*PI*y)*sin(2*PI*z)
  enddo
end subroutine create_square_from_annulus_order_n_flipped_mine


! x1 and y1 are coordinates in a circle
! x2 and y2 are coordinates in a square
subroutine convert_circle_to_square_mine(x1,y1,x2,y2)
	real, intent(in)      :: x1,y1
	real, intent(inout)   :: x2,y2
	real                  :: theta, R
	theta=atan(y1/(x1))
	if(x1<0) theta=theta+MATH_PI
	R=sqrt(x1**2+y1**2)
	x2=x1
	y2=y1
	if(theta>MATH_PI/4. .and. theta<3*MATH_PI/4.) then
		y2=R*sin(MATH_PI/4.)
		x2=x1*y2/y1
	elseif(theta>3*MATH_PI/4. .and. theta<5*MATH_PI/4.) then
		x2=R*cos(3*MATH_PI/4.)
		y2=y1*x2/x1
	elseif((theta>5*MATH_PI/4. .and. theta<7*MATH_PI/4.) .or. (theta>-2*MATH_PI/4. .and. theta<-MATH_PI/4.)) then
		y2=R*sin(5*MATH_PI/4.)
		x2=x1*y2/y1
	elseif(theta>-MATH_PI/4. .and. theta<MATH_PI/4.) then
		x2=R*cos(MATH_PI/4.)
		y2=y1*x2/x1
	endif
	
end subroutine convert_circle_to_square_mine


subroutine create_hexagon_from_annulus_order_n_flipped_mine(Ri,Ro,n)
  implicit none
  real, intent(in) :: Ri,Ro, n
  integer          :: i
  real             :: R,theta,x,y, x1,y1,M, p1,p2
  do i=1,NNODE
    x=NODES_ARRAY(i)%x
    y=NODES_ARRAY(i)%y
	theta=atan(y/(x))
	if(x<0) theta=theta+MATH_PI
	R=sqrt(x**2+y**2)
	call convert_circle_to_hexagon_mine(x,y,p1,p2)
	
     NODES_ARRAY(i)%x=-(x-p1)/(Ri**n-Ro**n)*R**n-(p1*Ro**n-x*Ri**n)/(Ri**n-Ro**n)
     NODES_ARRAY(i)%y=-(y-p2)/(Ri**n-Ro**n)*R**n-(p2*Ro**n-y*Ri**n)/(Ri**n-Ro**n)	
!     NODES_ARRAY(i)%x=p1
!     NODES_ARRAY(i)%y=p2
!     NODES_ARRAY(i)%z=z+0.05*sin(2*PI*x)*sin(2*PI*y)*sin(2*PI*z)
  enddo
end subroutine create_hexagon_from_annulus_order_n_flipped_mine

subroutine create_hexagon_from_annulus_order_n_flipped_three_regions(Ri,Ro,n)
  implicit none
  real, intent(in) :: Ri,Ro, n
  integer          :: i
  real             :: R,theta,x,y, x1,y1,M, p1,p2, Reg1,Reg2
  do i=1,NNODE
    x=NODES_ARRAY(i)%x
    y=NODES_ARRAY(i)%y
	theta=atan(y/(x))
	if(x<0) theta=theta+MATH_PI
	R=sqrt(x**2+y**2)
	call convert_circle_to_hexagon_mine(x,y,p1,p2)
	if(R<=Ri) then
		NODES_ARRAY(i)%x=p1
		NODES_ARRAY(i)%y=p2
	elseif(R>Ri .and. R<Ro) then
     NODES_ARRAY(i)%x=-(x-p1)/(Ri**n-Ro**n)*R**n-(p1*Ro**n-x*Ri**n)/(Ri**n-Ro**n)
     NODES_ARRAY(i)%y=-(y-p2)/(Ri**n-Ro**n)*R**n-(p2*Ro**n-y*Ri**n)/(Ri**n-Ro**n)
	endif

!     NODES_ARRAY(i)%z=z+0.05*sin(2*PI*x)*sin(2*PI*y)*sin(2*PI*z)
  enddo
end subroutine create_hexagon_from_annulus_order_n_flipped_three_regions

! x1 and y1 are coordinates in a circle
! x2 and y2 are coordinates in a square
subroutine convert_circle_to_hexagon_mine(x1,y1,x2,y2)
	real, intent(in)      :: x1,y1
	real, intent(inout)   :: x2,y2
	real                  :: theta, R
	real                  :: a1,a2,b1,b2,c1,c2,pp1,pp2
	theta=atan(y1/(x1))
	if(x1<0) theta=theta+MATH_PI
	R=sqrt(x1**2+y1**2)
	x2=x1
	y2=y1
	if(theta>0. .and. theta<MATH_PI/3.) then
		pp1=R
		pp2=0.
		a1=(sin(0.)-sin(MATH_PI/3.))/(cos(0.)-cos(MATH_PI/3.))
		b1=-1.
		c1=pp2-a1*pp1
		a2=y1/x1
		b2=-1.
		c2=0.
		x2=(b1*c2-b2*c1)/(a1*b2-a2*b1)
		y2=(c1*a2-c2*a1)/(a1*b2-a2*b1)
	elseif(theta>MATH_PI/3. .and. theta<2*MATH_PI/3.) then
		y2=R*sin(MATH_PI/3.)
		x2=x1*y2/y1
	elseif(theta>2*MATH_PI/3. .and. theta<3*MATH_PI/3.) then
		pp1=-R
		pp2=0.
		a1=(sin(2*MATH_PI/3.)-sin(3*MATH_PI/3.))/(cos(2*MATH_PI/3.)-cos(3*MATH_PI/3.))
		b1=-1.
		c1=pp2-a1*pp1
		a2=y1/x1
		b2=-1.
		c2=0.
		x2=(b1*c2-b2*c1)/(a1*b2-a2*b1)
		y2=(c1*a2-c2*a1)/(a1*b2-a2*b1)
	elseif(theta>3*MATH_PI/3. .and. theta<4*MATH_PI/3.) then
		pp1=-R
		pp2=0.
		a1=(sin(3*MATH_PI/3.)-sin(4*MATH_PI/3.))/(cos(3*MATH_PI/3.)-cos(4*MATH_PI/3.))
		b1=-1.
		c1=pp2-a1*pp1
		a2=y1/x1
		b2=-1.
		c2=0.
		x2=(b1*c2-b2*c1)/(a1*b2-a2*b1)
		y2=(c1*a2-c2*a1)/(a1*b2-a2*b1)
	elseif((theta>4*MATH_PI/3. .and. theta<5*MATH_PI/3.) .or. (theta>-2*MATH_PI/3. .and. theta<-MATH_PI/3.)) then
		y2=R*sin(4*MATH_PI/3.)
		x2=x1*y2/y1
	elseif(theta>-MATH_PI/3. .and. theta<0.) then
		pp1=R
		pp2=0.
		a1=(sin(-MATH_PI/3.)-sin(0.))/(cos(-MATH_PI/3.)-cos(0.))
		b1=-1.
		c1=pp2-a1*pp1
		a2=y1/x1
		b2=-1.
		c2=0.
		x2=(b1*c2-b2*c1)/(a1*b2-a2*b1)
		y2=(c1*a2-c2*a1)/(a1*b2-a2*b1)
	endif
	
end subroutine convert_circle_to_hexagon_mine


subroutine move_geometry(xx,xy,xz)
  implicit none
  real, intent(in) :: xx,xy,xz
  integer          :: i
  do i=1,NNODE
    NODES_ARRAY(i)%x=NODES_ARRAY(i)%x+xx
    NODES_ARRAY(i)%y=NODES_ARRAY(i)%y+xy
    NODES_ARRAY(i)%z=NODES_ARRAY(i)%z+xz
  enddo
end subroutine move_geometry


subroutine h2_perturb
  real    :: x,y,z,PI,a
  integer :: i
  do i=1,NNODE
    PI=4.D0*DATAN(1.D0)
    x=NODES_ARRAY(i)%x
    y=NODES_ARRAY(i)%y
    z=NODES_ARRAY(i)%z
!    a=0.125!+0.075
    a=0.
    a=0.1
!    a=0.225
!     NODES_ARRAY(i)%x=x+0.06*sin(2*PI*x)*sin(2*PI*y)*sin(2*PI*z)
!     NODES_ARRAY(i)%y=y+0.06*sin(2*PI*x)*sin(2*PI*y)*sin(2*PI*z)
!    NODES_ARRAY(i)%z=z+0.06*sin(2*PI*x)*sin(2*PI*y)*sin(2*PI*z)

     NODES_ARRAY(i)%x=x+a*sin(2*MATH_PI*x)*sin(2*MATH_PI*y)*sin(2*MATH_PI*z)
     NODES_ARRAY(i)%y=y+a*sin(2*MATH_PI*x)*sin(2*MATH_PI*y)*sin(2*MATH_PI*z)
     NODES_ARRAY(i)%z=z+a*sin(2*MATH_PI*x)*sin(2*MATH_PI*y)*sin(2*MATH_PI*z)
	 
 !    NODES_ARRAY(i)%x=x+a*sin(2*MATH_PI*x)*sin(2*MATH_PI*z)
 !    NODES_ARRAY(i)%y=y-a*sin(2*MATH_PI*x)*sin(2*MATH_PI*y)
 !    NODES_ARRAY(i)%z=z-a*sin(2*MATH_PI*x)*sin(2*MATH_PI*z)

!     NODES_ARRAY(i)%z=z+a*sin(2*MATH_PI*z)

     ! NODES_ARRAY(i)%x=x+a*sin(2*PI*x)*sin(2*PI*y)
     ! NODES_ARRAY(i)%y=y-a*sin(2*PI*x)*sin(2*PI*y)

  enddo
end subroutine h2_perturb





! *************** Determinant of 3x3 matrix ****************************
  function det(a) result (x)
    implicit none
    real,dimension(0:2,0:2), intent (in)  :: a         ! if omitted ==> 2D
    real                                  :: x
    x=a(0,0)*a(1,1)*a(2,2) + a(0,1)*a(1,2)*a(2,0) + a(0,2)*a(1,0)*a(2,1) &
    - a(0,2)*a(1,1)*a(2,0) - a(0,1)*a(1,0)*a(2,2) - a(0,0)*a(1,2)*a(2,1)
  end function det

! *************** inverse of 3x3 matrix ****************************
  function inv3x3(a) result (x)
    implicit none
    real,dimension(0:2,0:2), intent (in)  :: a
    real,dimension(0:2,0:2)               :: x
    x = (1/det(a))*transpose(reshape((/ &
         a(1,1)*a(2,2)-a(2,1)*a(1,2), a(0,2)*a(2,1)-a(0,1)*a(2,2), a(0,1)*a(1,2)-a(0,2)*a(1,1), &
         a(1,2)*a(2,0)-a(1,0)*a(2,2), a(0,0)*a(2,2)-a(0,2)*a(2,0), a(0,2)*a(1,0)-a(0,0)*a(1,2), &
         a(1,0)*a(2,1)-a(1,1)*a(2,0), a(0,1)*a(2,0)-a(0,0)*a(2,1), a(0,0)*a(1,1)-a(0,1)*a(1,0) /), shape(x)))
  end function inv3x3


! *************** unit normal of three vectors ****************************
  function unit_normal(a,b,c) result (f)
    implicit none
    real,dimension(0:2), intent (in)      :: a,b,c
    real,dimension(0:2)                   :: f
    real                                  :: x,y,z,temp
    real,dimension(0:2,0:2)               :: x1, x2, x3
    x1 = transpose(reshape((/ 1.0, a(1), a(2), &
                              1.0, b(1), b(2), &
                              1.0, c(1), c(2) /), shape(x1)))
    
    x2 = transpose(reshape((/ a(0), 1.0, a(2), &
                              b(0), 1.0, b(2), &
                              c(0), 1.0, c(2) /), shape(x2)))
    x3 = transpose(reshape((/ a(0), a(1), 1.0, &
                              b(0), b(1), 1.0, &
                              c(0), c(1), 1.0 /), shape(x3)))
    x=det(x1)
    y=det(x2)
    z=det(x3)
    temp=(x**2 + y**2 + z**2)**.5
    f(0)=x/temp
    f(1)=y/temp
    f(2)=z/temp
  end function unit_normal
! *************** dot product of vectors ****************************
  function dot(a,b) result (x)
    implicit none
    real,dimension(1:3), intent (in)  :: a,b
    real                              :: x
    x=a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
  end function dot


! *************** cross product of vectors ****************************
  function cross(a,b) result (x)
    real,dimension(0:2), intent (in)  :: a,b 
    real,dimension(0:2)               :: x
    x(0)=a(1) * b(2) - a(2) * b(1)
    x(1)=a(2) * b(0) - a(0) * b(2)
    x(2)=a(0) * b(1) - a(1) * b(0)
  end function cross

! ********************** finding numbe rof lines in text file **********
  function number_of_lines(infile) result(n)
    implicit none
    CHARACTER(len=*), intent(in) :: infile
    integer :: n, i,ios
    CHARACTER(len=200)  :: line
    print*, infile
    open(unit=3, file=infile, iostat=ios)
    if ( ios /= 0 ) stop "Err100 error opening file input"

    n = 0
    do
        read(3, '(A)', iostat=ios) line
        if (ios /= 0) exit
        n = n + 1
    end do
    print*, 'number of lines in file ',infile,' is : ',n
    close(3)
  end function number_of_lines



  subroutine find_mesh_number(n_vertices,n_faces,n_elements)
    implicit none
    integer, intent(inout) :: n_vertices,n_faces,n_elements
    integer :: n, i,ios
    CHARACTER(len=200)  :: line
    
    open(unit=1, file='mesh//mesh.msh', iostat=ios)
    if ( ios /= 0 ) stop "Err100 error opening file input"

    n = 0    
    do
        read(1, '(A)', iostat=ios) line
        if (ios /= 0) exit
        n = n + 1
    end do
    print*, 'number of lines: ',n
    close(1)
    
  end subroutine find_mesh_number


! This function is to check if x element is existed in the A array or not
! Output is true or false
! input: array (A), array length (n), element to check (x)
! The array is integer array
  function is_exist (A,n,x) result (isexist)
     integer, intent (in)                :: n,x
     integer, dimension(1:n), intent(in) :: A
     logical                             :: isexist
     integer                             :: i
     isexist=.FALSE.
     do i=1,n
        if (A(i)==x) isexist=.TRUE.
     enddo
!     if ( ANY( A==x ) ) isexist=.TRUE.
  end function is_exist


! Position of element in an array
  function position_in_array (A,n,x) result (i)
     integer, intent (in)                :: n,x
     integer, dimension(1:n), intent(in) :: A
     integer                             :: i
     do i=1,n
        if (A(i)==x) THEN
            isexist=i
            return
        endif
     enddo
!     if ( ANY( A==x ) ) isexist=.TRUE.
  end function position_in_array

! ======================= ARRAY SORTING ==============================
! Sort an array in ascending order 
  subroutine array_ascending (A,n)
     integer, intent (in)                :: n
     real, dimension(1:n), intent(inout) :: A
     integer                             :: i,j
     real                                :: t
     do i=1,n
        do j=i,n
           if(A(i)>A(j)) then
              t=a(j)
              a(j)=a(i)
              a(i)=t
           endif
        enddo
     enddo
  end subroutine array_ascending
! ======================= ARRAY SORTING ==============================
! Sort an array in ascending order 
  subroutine array_ascending2 (A,n,m,col)
! A is two dimensional array to be sorted
! n is the number of rows
! m is the number of columns
! col is the column to be sorted based on
     integer, intent (in)                    :: n,m,col
     real, dimension(1:n,1:m), intent(inout) :: A
     integer                                 :: i,j
     real,dimension(1:m)                     :: t
     do i=1,n
        do j=i,n
           if(A(i,col)>A(j,col)) then
              t(1:m)=A(j,1:m)
              A(j,1:m)=A(i,1:m)
              A(i,1:m)=t
           endif
        enddo
     enddo
  end subroutine array_ascending2


! ======================= ARRAY SORTING ==============================
! Sort an array in ascending order 
  subroutine array_descending (A,n)
     integer, intent (in)                :: n
     real, dimension(1:n), intent(inout) :: A
     integer                             :: i,j
     real                                :: t
     do i=1,n
        do j=i,n
           if(A(i)<A(j)) then
              t=a(j)
              a(j)=a(i)
              a(i)=t
           endif
        enddo
     enddo
  end subroutine array_descending



! **************** Convert string ---> Integer *****************
  function str2int(str)  result(int)
    implicit none
    ! Arguments
    character(len=*),intent(in) :: str
    integer         :: int
    integer         :: stat

    read(str,*,iostat=stat)  int
    if ( stat /= 0 ) then
      print *,'ERR101: Conversion of string failed... FUNCTION (str2int-class_GeneralFunctions)'
    endif
  end function str2int




  subroutine map_to_xyz(x2,y2,z2,x,y,z,C,D,E)
! Mapping from x2,y2,z2 ----------> x,y,z (good for triple intergrals)
    real, dimension(0:7), intent(in) :: C,D,E     ! Coefficient arrays (mapping coeficcient for each element)
    real, intent(in)                 :: x2,y2,z2
    real, intent(inout)              :: x,y,z
!    real                             :: x
! -----------------------------------
     
    x=C(0)+C(1)*x2+C(2)*y2+C(3)*z2+C(4)*x2*y2+C(5)*x2*z2+C(6)*y2*z2+C(7)*x2*y2*z2
    y=D(0)+D(1)*x2+D(2)*y2+D(3)*z2+D(4)*x2*y2+D(5)*x2*z2+D(6)*y2*z2+D(7)*x2*y2*z2
    z=E(0)+E(1)*x2+E(2)*y2+E(3)*z2+E(4)*x2*y2+E(5)*x2*z2+E(6)*y2*z2+E(7)*x2*y2*z2
  end subroutine map_to_xyz


  subroutine map_ti(time,t_i)
! This function is to map the time variable to the -1 -> 1 domain
    real, intent(inout)              :: t_i
    real, intent(in)                 :: time
    real                             :: a, b
    a = TI-DT
    b = TI+DT
    t_i= (b-a)/2. *time +(b+a)/2.
  end subroutine map_ti


! The following function perform gauss integration in 3d
  function Gauss_integration (A) result (Int)
    real, intent (in), dimension(1:Gauss_int_p,1:Gauss_int_p,1:Gauss_int_p) :: A
    real                                      :: Int
    integer :: i,j,k
    Int=0.
    do i=1,Gauss_int_p
       do j=1,Gauss_int_p
          do k=1,Gauss_int_p
             Int=Int+Gauss_w_i(i)*Gauss_w_i(j)*Gauss_w_i(k)*A(i,j,k)
          enddo
       enddo
    enddo
!  print*, Int
  end function Gauss_integration

! The following function perform gauss integration in 3d with time
! Time variables should be fixed to make the integration limits: -1 -> 1
  function Gauss_integration_time (A) result (Int)
    real, intent (in), dimension(1:Gauss_int_p,1:Gauss_int_p,1:Gauss_int_p,1:Gauss_int_p) :: A
    real                                      :: Int
    integer :: i,j,k,l
    Int=0.
    do i=1,Gauss_int_p
       do j=1,Gauss_int_p
          do k=1,Gauss_int_p
             do l=1,Gauss_int_p
                Int=Int+Gauss_w_i(i)*Gauss_w_i(j)*Gauss_w_i(k)*Gauss_w_i(l)*A(i,j,k,l)
             enddo
          enddo
       enddo
    enddo
!  print*, Int
  end function Gauss_integration_time

  function Gauss_integration2d (A) result (Int)
    real, intent (in), dimension(1:Gauss_int_p,1:Gauss_int_p) :: A
    real                                      :: Int
    integer :: i,j
    Int=0.
    do i=1,Gauss_int_p
       do j=1,Gauss_int_p
             Int=Int+Gauss_w_i(i)*Gauss_w_i(j)*A(i,j)
       enddo
    enddo
!  print*, Int
  end function Gauss_integration2d




  function Gauss_integration2d_face (A,fn) result (Int)
    real, intent (in), dimension(1:Gauss_int_p,1:Gauss_int_p,1:Gauss_int_p) :: A
    integer, intent(in)                         :: fn
    real                                        :: Int
    real, dimension(1:Gauss_int_p,1:Gauss_int_p) :: A2d


    if (fn == 1 .OR. fn==2) then
       A2d=A(1,:,:)
    else if (fn == 3 .OR. fn==4) then
       A2d=A(:,1,:)
    else if (fn == 5 .OR. fn==6) then
       A2d=A(:,:,1)
    endif
    INT=Gauss_integration2d(A2d)
  end function Gauss_integration2d_face

  function Gauss_integration2d_face_time (A,fn) result (Int)
    real, intent (in), dimension(1:Gauss_int_p,1:Gauss_int_p,1:Gauss_int_p,1:Gauss_int_p) :: A
    integer, intent(in)                         :: fn
    real                                        :: Int
    real, dimension(1:Gauss_int_p,1:Gauss_int_p,1:Gauss_int_p) :: A2d


    if (fn == 1 .OR. fn==2) then
       A2d=A(1,:,:,:)
    else if (fn == 3 .OR. fn==4) then
       A2d=A(:,1,:,:)
    else if (fn == 5 .OR. fn==6) then
       A2d=A(:,:,1,:)
    endif
    INT=Gauss_integration(A2d)
  end function Gauss_integration2d_face_time


recursive subroutine quicksort(a, first, last)
  implicit none
  integer  a(*), x, t
  integer first, last
  integer i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort(a, first, i-1)
  if (j+1 < last)  call quicksort(a, j+1, last)
end subroutine quicksort


end module class_general



