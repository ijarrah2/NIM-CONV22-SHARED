! This class contains all the constants associated with the convection-diffusion equation
module class_nim
use class_general
use class_parameters
use class_mesh
use class_curvilinear
use class_create_mesh
use class_prepare_mesh
use class_geometry
implicit none

contains

! ++++++++++++++++++++++++++++ Velocity +++++++++++++++++++++++++++++
! *****************  Calculating average constants ******************

  subroutine calculate_velocity_nim_constants
    implicit none
    integer :: i,i2
	real    :: B_x,B_y,B_z,a,b,c,Re_x,Re_y,Re_z,K11,K22,K33
	real,dimension(1:NMESH)    :: Det
	real,dimension(1:3,1:NMESH)    :: CC,EE,FF,GG,HH,II,JJ

    a=1.
    b=1.
    c=1.	

	!$OMP     PARALLEL DO PRIVATE (i,Re_x,Re_y,Re_z,B_x,B_y,B_z,K11,K22,K33)
    do i=1,NMESH
	
		meshes_u(i)=meshes_uxyzt(i)*meshes_xi_x(1,1,i)+meshes_vxyzt(i)*meshes_xi_x(1,2,i)+meshes_wxyzt(i)*meshes_xi_x(1,3,i)
		meshes_v(i)=meshes_uxyzt(i)*meshes_xi_x(2,1,i)+meshes_vxyzt(i)*meshes_xi_x(2,2,i)+meshes_wxyzt(i)*meshes_xi_x(2,3,i)
		meshes_w(i)=meshes_uxyzt(i)*meshes_xi_x(3,1,i)+meshes_vxyzt(i)*meshes_xi_x(3,2,i)+meshes_wxyzt(i)*meshes_xi_x(3,3,i)
                K11=meshes_gij(1,1,i)*nu
		K22=meshes_gij(2,2,i)*nu
		K33=meshes_gij(3,3,i)*nu
		
		Re_x=2*meshes_u(i)/K11
		Re_y=2*meshes_v(i)/K22
		Re_z=2*meshes_w(i)/K33

		B_x=B_nim(meshes_u(i),K11, Re_x)
		B_y=B_nim(meshes_v(i),K22, Re_y)
		B_z=B_nim(meshes_w(i),K33, Re_z)

		CC(1,i)=C_nim(meshes_u(i),K11,Re_x,B_x,a)
		CC(2,i)=C_nim(meshes_v(i),K22,Re_y,B_y,b)
		CC(3,i)=C_nim(meshes_w(i),K33,Re_z,B_z,c)

		Det(i)=1./(CC(1,i)*CC(2,i)+CC(2,i)*CC(3,i)+CC(1,i)*CC(3,i))

		EE(1,i)=E_nim(meshes_u(i),K11,Re_x,B_x,Det(i),a)
		EE(2,i)=E_nim(meshes_v(i),K22,Re_y,B_y,Det(i),b)
		EE(3,i)=E_nim(meshes_w(i),K33,Re_z,B_z,Det(i),c)

		FF(1,i)=F_nim(meshes_u(i),K11,Re_x,B_x,Det(i),a)
		FF(2,i)=F_nim(meshes_v(i),K22,Re_y,B_y,Det(i),b)
		FF(3,i)=F_nim(meshes_w(i),K33,Re_z,B_z,Det(i),c)

		GG(1,i)=G_nim(meshes_u(i),K11,Re_x,a)
		GG(2,i)=G_nim(meshes_v(i),K22,Re_y,b)
		GG(3,i)=G_nim(meshes_w(i),K33,Re_z,c)

		HH(1,i)=H_nim(meshes_u(i),K11,Re_x,a)
		HH(2,i)=H_nim(meshes_v(i),K22,Re_y,b)
		HH(3,i)=H_nim(meshes_w(i),K33,Re_z,c)

		II(1,i)=I_nim(meshes_u(i),K11,Re_x,a)
		II(2,i)=I_nim(meshes_v(i),K22,Re_y,b)
		II(3,i)=I_nim(meshes_w(i),K33,Re_z,c)

		JJ(1,i)=J_nim(meshes_u(i),K11,Re_x,a)
		JJ(2,i)=J_nim(meshes_v(i),K22,Re_y,b)
		JJ(3,i)=J_nim(meshes_w(i),K33,Re_z,c)

    enddo
  !$OMP     END PARALLEL DO 

!$OMP     PARALLEL DO PRIVATE (i,i2)
    do i=1,NMESH
! ==================== |F|A|C|E| |2| =====================
    if (meshes_neighbors(2,i)/=0) THEN
       i2=meshes_neighbors(2,i)
!       M2=MESHES_ARRAY(MESHES_ARRAY(i)%neighbors(2))%const_V
       AA_vx(1,i)=-(GG(1,i)+(CC(2,i)+CC(3,i))*EE(1,i)*HH(1,i))*meshes_gij_faces(1,1,2,i)
       AA_vx(2,i)= (GG(1,i)-(CC(2,i)+CC(3,i))*FF(1,i)*HH(1,i))*meshes_gij_faces(1,1,2,i) &
                                    +(II(1,i2)+(CC(2,i2)+CC(3,i2))*EE(1,i2)*JJ(1,i2))*meshes_gij_faces(1,1,1,i2)
       AA_vx(3,i)=-(II(1,i2)-(CC(2,i2)+CC(3,i2))*FF(1,i2)*JJ(1,i2))*meshes_gij_faces(1,1,1,i2)
       AA_vx(4,i)=  CC(3,i)*EE(2,i)*HH(1,i)*meshes_gij_faces(1,1,2,i)
       AA_vx(5,i)=  CC(3,i)*FF(2,i)*HH(1,i)*meshes_gij_faces(1,1,2,i)
       AA_vx(6,i)= -CC(3,i2)*EE(2,i2)*JJ(1,i2)*meshes_gij_faces(1,1,1,i2)
       AA_vx(7,i)= -CC(3,i2)*FF(2,i2)*JJ(1,i2)*meshes_gij_faces(1,1,1,i2)
       AA_vx(8,i)=  CC(2,i)*EE(3,i)*HH(1,i)*meshes_gij_faces(1,1,2,i)        
       AA_vx(9,i)=  CC(2,i)*FF(3,i)*HH(1,i)*meshes_gij_faces(1,1,2,i)
       AA_vx(10,i)=-CC(2,i2)*EE(3,i2)*JJ(1,i2)*meshes_gij_faces(1,1,1,i2)
       AA_vx(11,i)=-CC(2,i2)*FF(3,i2)*JJ(1,i2)*meshes_gij_faces(1,1,1,i2)
       AA_vx(12,i)= CC(2,i)*CC(3,i)*Det(i)*HH(1,i)*meshes_gij_faces(1,1,2,i)/(2*DT)
       AA_vx(13,i)=-CC(2,i)*CC(3,i)*Det(i)*HH(1,i)*meshes_gij_faces(1,1,2,i)/(2*DT)
       AA_vx(14,i)=-CC(2,i2)*CC(3,i2)*Det(i2)*JJ(1,i2)*meshes_gij_faces(1,1,1,i2)/(2*DT)
       AA_vx(15,i)= CC(2,i2)*CC(3,i2)*Det(i2)*JJ(1,i2)*meshes_gij_faces(1,1,1,i2)/(2*DT)
    endif

! ==================== |F|A|C|E| |4| =====================
    if (meshes_neighbors(4,i)/=0) THEN
       i2=meshes_neighbors(4,i)
       AA_vy(1,i)=  CC(3,i)*EE(1,i)*HH(2,i)*meshes_gij_faces(2,2,4,i)
       AA_vy(2,i)=  CC(3,i)*FF(1,i)*HH(2,i)*meshes_gij_faces(2,2,4,i)
       AA_vy(3,i)= -CC(3,i2)*EE(1,i2)*JJ(2,i2)*meshes_gij_faces(2,2,3,i2)
       AA_vy(4,i)= -CC(3,i2)*FF(1,i2)*JJ(2,i2)*meshes_gij_faces(2,2,3,i2)
       AA_vy(5,i)=-(GG(2,i)+(CC(1,i)+CC(3,i))*EE(2,i)*HH(2,i))*meshes_gij_faces(2,2,4,i)
       AA_vy(6,i)= (GG(2,i)-(CC(1,i)+CC(3,i))*FF(2,i)*HH(2,i))*meshes_gij_faces(2,2,4,i) &
                                    +(II(2,i2)+(CC(1,i2)+CC(3,i2))*EE(2,i2)*JJ(2,i2))*meshes_gij_faces(2,2,3,i2)
       AA_vy(7,i)=-(II(2,i2)-(CC(1,i2)+CC(3,i2))*FF(2,i2)*JJ(2,i2))*meshes_gij_faces(2,2,3,i2)
       AA_vy(8,i)=  CC(1,i)*EE(3,i)*HH(2,i)*meshes_gij_faces(2,2,4,i)
       AA_vy(9,i)=  CC(1,i)*FF(3,i)*HH(2,i)*meshes_gij_faces(2,2,4,i)
       AA_vy(10,i)=-CC(1,i2)*EE(3,i2)*JJ(2,i2)*meshes_gij_faces(2,2,3,i2)
       AA_vy(11,i)=-CC(1,i2)*FF(3,i2)*JJ(2,i2)*meshes_gij_faces(2,2,3,i2)
       AA_vy(12,i)= CC(1,i)*CC(3,i)*Det(i)*HH(2,i)*meshes_gij_faces(2,2,4,i)/(2*DT)
       AA_vy(13,i)=-CC(1,i)*CC(3,i)*Det(i)*HH(2,i)*meshes_gij_faces(2,2,4,i)/(2*DT)
       AA_vy(14,i)=-CC(1,i2)*CC(3,i2)*Det(i2)*JJ(2,i2)*meshes_gij_faces(2,2,3,i2)/(2*DT)
       AA_vy(15,i)= CC(1,i2)*CC(3,i2)*Det(i2)*JJ(2,i2)*meshes_gij_faces(2,2,3,i2)/(2*DT)
    endif
! ==================== |F|A|C|E| |6| =====================
    if (meshes_neighbors(6,i)/=0) THEN
       i2=meshes_neighbors(6,i)
       AA_vz(1,i) =  CC(2,i)*EE(1,i)*HH(3,i)*meshes_gij_faces(3,3,6,i)
       AA_vz(2,i) =  CC(2,i)*FF(1,i)*HH(3,i)*meshes_gij_faces(3,3,6,i)
       AA_vz(3,i) = -CC(2,i2)*EE(1,i2)*JJ(3,i2)*meshes_gij_faces(3,3,5,i2)
       AA_vz(4,i) = -CC(2,i2)*FF(1,i2)*JJ(3,i2)*meshes_gij_faces(3,3,5,i2)
       AA_vz(5,i) =  CC(1,i)*EE(2,i)*HH(3,i)*meshes_gij_faces(3,3,6,i)
       AA_vz(6,i) =  CC(1,i)*FF(2,i)*HH(3,i)*meshes_gij_faces(3,3,6,i)
       AA_vz(7,i) = -CC(1,i2)*EE(2,i2)*JJ(3,i2)*meshes_gij_faces(3,3,5,i2)
       AA_vz(8,i) = -CC(1,i2)*FF(2,i2)*JJ(3,i2)*meshes_gij_faces(3,3,5,i2)
       AA_vz(9,i) =-(GG(3,i)+(CC(1,i)+CC(2,i))*EE(3,i)*HH(3,i))*meshes_gij_faces(3,3,6,i)
       AA_vz(10,i) =(GG(3,i)-(CC(1,i)+CC(2,i))*FF(3,i)*HH(3,i))*meshes_gij_faces(3,3,6,i) &
                                     +(II(3,i2)+(CC(1,i2)+CC(2,i2))*EE(3,i2)*JJ(3,i2))*meshes_gij_faces(3,3,5,i2)
       AA_vz(11,i)=-(II(3,i2)-(CC(1,i2)+CC(2,i2))*FF(3,i2)*JJ(3,i2))*meshes_gij_faces(3,3,5,i2)
       AA_vz(12,i)= CC(1,i)*CC(2,i)*Det(i)*HH(3,i)*meshes_gij_faces(3,3,6,i)/(2*DT)
       AA_vz(13,i)=-CC(1,i)*CC(2,i)*Det(i)*HH(3,i)*meshes_gij_faces(3,3,6,i)/(2*DT)
       AA_vz(14,i)=-CC(1,i2)*CC(2,i2)*Det(i2)*JJ(3,i2)*meshes_gij_faces(3,3,5,i2)/(2*DT)
       AA_vz(15,i)= CC(1,i2)*CC(2,i2)*Det(i2)*JJ(3,i2)*meshes_gij_faces(3,3,5,i2)/(2*DT)

    endif
       a_t(1,i)=2*EE(1,i)/Det(i) - 2*(CC(2,i)+CC(3,i))*CC(1,i)*EE(1,i)
       a_t(2,i)=2*FF(1,i)/Det(i) - 2*(CC(2,i)+CC(3,i))*CC(1,i)*FF(1,i)
       a_t(3,i)=2*EE(2,i)*CC(1,i)*CC(3,i)
       a_t(4,i)=2*FF(2,i)*CC(1,i)*CC(3,i)
       a_t(5,i)=2*EE(3,i)*CC(1,i)*CC(2,i)
       a_t(6,i)=2*FF(3,i)*CC(1,i)*CC(2,i)
       a_t(7,i)=2*CC(1,i)*CC(2,i)*CC(3,i)*Det(i)
       a_t(8,i)=-a_t(7,i)/(2*DT)-1
       a_t(9,i)= a_t(7,i)/(2*DT)-1

    enddo
  !$OMP     END PARALLEL DO 

  end subroutine calculate_velocity_nim_constants


  subroutine calculate_velocity_nim_constants_formulation2
! -----------------------------------
    integer                         :: i,i2
	real							:: Re_x,Re_y,Re_z,u,v,w,a,b,c,k11,k22,k33
	!$OMP     PARALLEL DO PRIVATE (i,Re_x,Re_y,Re_z,u,v,w,a,b,c,k11,k22,k33)
	do i=1,NMESH
		a=1.
		b=1.
		c=1.
		meshes_u(i)=meshes_uxyzt(i)*meshes_xi_x(1,1,i)+meshes_vxyzt(i)*meshes_xi_x(1,2,i)+meshes_wxyzt(i)*meshes_xi_x(1,3,i)
		meshes_v(i)=meshes_uxyzt(i)*meshes_xi_x(2,1,i)+meshes_vxyzt(i)*meshes_xi_x(2,2,i)+meshes_wxyzt(i)*meshes_xi_x(2,3,i)
		meshes_w(i)=meshes_uxyzt(i)*meshes_xi_x(3,1,i)+meshes_vxyzt(i)*meshes_xi_x(3,2,i)+meshes_wxyzt(i)*meshes_xi_x(3,3,i)

		u=meshes_u(i)
		v=meshes_v(i)
		w=meshes_w(i)
		
		k11=meshes_gij(1,1,i)*nu
		k22=meshes_gij(2,2,i)*nu
		k33=meshes_gij(3,3,i)*nu
		
		Re_x=2*a*u/k11
		Re_y=2*b*v/k22
		Re_z=2*c*w/k33
		if(abs(Re_x)>v_limit) then
		  M_A(1,1,i)=(2*a/(k11*(1./exp(Re_x)-1.))+1./u)*meshes_gij_faces(1,1,2,i)
		  M_A(1,2,i)=(-2*a/(k11*(1.-exp(Re_x)))-1./u)*meshes_gij_faces(1,1,1,i)
		  M_A(1,3,i)=(u/(k11*(1./exp(Re_x)-1.)))*meshes_gij_faces(1,1,2,i)
		  M_A(1,4,i)=(u/(k11*(1.-exp(Re_x))))*meshes_gij_faces(1,1,1,i)
		  
		  if(Re_x>0) M_A(1,5,i)=(-k11+(a*u*(1./exp(Re_x)+1.)/(-1./exp(Re_x)+1.)))/(u**2)
		  if(Re_x<0) M_A(1,5,i)=(-k11+(a*u*(1.+exp(Re_x))/(-1.+exp(Re_x))))/(u**2)

		  if(Re_x>0) M_A(1,6,i)=1./(-1./exp(Re_x)+1.)-1./(Re_x)
		  if(Re_x<0) M_A(1,6,i)=exp(Re_x)/(-1.+exp(Re_x))-1./(Re_x)
		else

		  M_A(1,1,i)= (-(a/k11) - (a**2*u)/(3.*k11**2) + (a**4*u**3)/(45.*k11**4) - (2*a**6*u**5)/(945.*k11**6) + &
		 -  (a**8*u**7)/(4725.*k11**8) - (2*a**10*u**9)/(93555.*k11**10))*meshes_gij_faces(1,1,2,i)

		  M_A(1,2,i)=(-(a/k11) + (a**2*u)/(3.*k11**2) - (a**4*u**3)/(45.*k11**4) + (2*a**6*u**5)/(945.*k11**6) - &
		 -  (a**8*u**7)/(4725.*k11**8) + (2*a**10*u**9)/(93555.*k11**10))*meshes_gij_faces(1,1,1,i)

		  M_A(1,3,i)= (-1/(2.*a) - u/(2.*k11) - (a*u**2)/(6.*k11**2) + (a**3*u**4)/(90.*k11**4) - &
		 -  (a**5*u**6)/(945.*k11**6) + (a**7*u**8)/(9450.*k11**8) - (a**9*u**10)/(93555.*k11**10))*meshes_gij_faces(1,1,2,i)

		  M_A(1,4,i)=(-1/(2.*a) + u/(2.*k11) - (a*u**2)/(6.*k11**2) + (a**3*u**4)/(90.*k11**4) - &
		 -  (a**5*u**6)/(945.*k11**6) + (a**7*u**8)/(9450.*k11**8) - (a**9*u**10)/(93555.*k11**10))*meshes_gij_faces(1,1,1,i)

		  M_A(1,5,i)= a**2/(3.*k11) - (a**4*u**2)/(45.*k11**3) + (2*a**6*u**4)/(945.*k11**5) - &
		 -  (a**8*u**6)/(4725.*k11**7) + (2*a**10*u**8)/(93555.*k11**9) - &
		 -  (1382*a**12*u**10)/(6.38512875e8*k11**11)

		  M_A(1,6,i)= 0.5 + (a*u)/(6.*k11) - (a**3*u**3)/(90.*k11**3) + (a**5*u**5)/(945.*k11**5) - &
		 -  (a**7*u**7)/(9450.*k11**7) + (a**9*u**9)/(93555.*k11**9)

		endif

		if(abs(Re_y)>v_limit) then
		  M_A(2,1,i)=(2*b/(k22*(1./exp(Re_y)-1.))+1./v)*meshes_gij_faces(2,2,4,i)
		  M_A(2,2,i)=(-2*b/(k22*(1-exp(Re_y)))-1./v)*meshes_gij_faces(2,2,3,i)
		  M_A(2,3,i)=(v/(k22*(1./exp(Re_y)-1.)))*meshes_gij_faces(2,2,4,i)
		  M_A(2,4,i)=(v/(k22*(1-exp(Re_y))))*meshes_gij_faces(2,2,3,i)

		  if(Re_y>0) M_A(2,5,i)=(-k22+(b*v*(1./exp(Re_y)+1.)/(-1./exp(Re_y)+1.)))/(v**2)
		  if(Re_y<0) M_A(2,5,i)=(-k22+(b*v*(1+exp(Re_y))/(-1.+exp(Re_y))))/(v**2)

		  if(Re_y>0) M_A(2,6,i)=1./(-1./exp(Re_y)+1.)-1./(Re_y)
		  if(Re_y<0) M_A(2,6,i)=exp(Re_y)/(-1+exp(Re_y))-1./(Re_y)
		!      print*, M_i,M_j,M_B1,M_B2,M_B3,M_B4,M_D1,Re_y


		else
		  M_A(2,1,i)= (-(b/k22) - (b**2*v)/(3.*k22**2) + (b**4*v**3)/(45.*k22**4) - (2*b**6*v**5)/(945.*k22**6) + &
		 -  (b**8*v**7)/(4725.*k22**8) - (2*b**10*v**9)/(93555.*k22**10))*meshes_gij_faces(2,2,4,i)

		  M_A(2,2,i)=(-(b/k22) + (b**2*v)/(3.*k22**2) - (b**4*v**3)/(45.*k22**4) + (2*b**6*v**5)/(945.*k22**6) - &
		 -  (b**8*v**7)/(4725.*k22**8) + (2*b**10*v**9)/(93555.*k22**10))*meshes_gij_faces(2,2,3,i)

		  M_A(2,3,i)= (-1/(2.*b) - v/(2.*k22) - (b*v**2)/(6.*k22**2) + (b**3*v**4)/(90.*k22**4) - &
		 -  (b**5*v**6)/(945.*k22**6) + (b**7*v**8)/(9450.*k22**8) - (b**9*v**10)/(93555.*k22**10))*meshes_gij_faces(2,2,4,i)

		  M_A(2,4,i)=(-1/(2.*b) + v/(2.*k22) - (b*v**2)/(6.*k22**2) + (b**3*v**4)/(90.*k22**4) - &
		 -  (b**5*v**6)/(945.*k22**6) + (b**7*v**8)/(9450.*k22**8) - (b**9*v**10)/(93555.*k22**10))*meshes_gij_faces(2,2,3,i)

		  M_A(2,5,i)= b**2/(3.*k22) - (b**4*v**2)/(45.*k22**3) + (2*b**6*v**4)/(945.*k22**5) - &
		 -  (b**8*v**6)/(4725.*k22**7) + (2*b**10*v**8)/(93555.*k22**9) - &
		 -  (1382*b**12*v**10)/(6.38512875e8*k22**11)

		  M_A(2,6,i)= 0.5 + (b*v)/(6.*k22) - (b**3*v**3)/(90.*k22**3) + (b**5*v**5)/(945.*k22**5) - &
		 -  (b**7*v**7)/(9450.*k22**7) + (b**9*v**9)/(93555.*k22**9)

		endif
		if(abs(Re_z)>v_limit) then
		  M_A(3,1,i)=(2*c/(k33*(1./exp(Re_z)-1.))+1./w)*meshes_gij_faces(3,3,6,i)
		  M_A(3,2,i)=(-2*c/(k33*(1-exp(Re_z)))-1./w)*meshes_gij_faces(3,3,5,i)
		  M_A(3,3,i)=(w/(k33*(1./exp(Re_z)-1.)))*meshes_gij_faces(3,3,6,i)
		  M_A(3,4,i)=(w/(k33*(1-exp(Re_z))))*meshes_gij_faces(3,3,5,i)

		  if(Re_z>0) M_A(3,5,i)=(-k33+(c*w*(1./exp(Re_z)+1.)/(-1./exp(Re_z)+1.)))/(w**2)
		  if(Re_z<0) M_A(3,5,i)=(-k33+(c*w*(1+exp(Re_z))/(-1.+exp(Re_z))))/(w**2)

		  if(Re_z>0) M_A(3,6,i)=1./(-1./exp(Re_z)+1.)-1./(Re_z)
		  if(Re_z<0) M_A(3,6,i)=exp(Re_z)/(-1+exp(Re_z))-1./(Re_z)
		!      print*, M_i,M_j,M_B1,M_B2,M_B3,M_B4,M_D1,Re_z


		else
		  M_A(3,1,i)= (-(c/k33) - (c**2*w)/(3.*k33**2) + (c**4*w**3)/(45.*k33**4) - (2*c**6*w**5)/(945.*k33**6) + &
		 -  (c**8*w**7)/(4725.*k33**8) - (2*c**10*w**9)/(93555.*k33**10))*meshes_gij_faces(3,3,6,i)

		  M_A(3,2,i)=(-(c/k33) + (c**2*w)/(3.*k33**2) - (c**4*w**3)/(45.*k33**4) + (2*c**6*w**5)/(945.*k33**6) - &
		 -  (c**8*w**7)/(4725.*k33**8) + (2*c**10*w**9)/(93555.*k33**10))*meshes_gij_faces(3,3,5,i)
		 
		  M_A(3,3,i)= (-1/(2.*c) - w/(2.*k33) - (c*w**2)/(6.*k33**2) + (c**3*w**4)/(90.*k33**4) - &
		 -  (c**5*w**6)/(945.*k33**6) + (c**7*w**8)/(9450.*k33**8) - (c**9*w**10)/(93555.*k33**10))*meshes_gij_faces(3,3,6,i)

		  M_A(3,4,i)=(-1/(2.*c) + w/(2.*k33) - (c*w**2)/(6.*k33**2) + (c**3*w**4)/(90.*k33**4) - &
		 -  (c**5*w**6)/(945.*k33**6) + (c**7*w**8)/(9450.*k33**8) - (c**9*w**10)/(93555.*k33**10))*meshes_gij_faces(3,3,5,i)

		  M_A(3,5,i)= c**2/(3.*k33) - (c**4*w**2)/(45.*k33**3) + (2*c**6*w**4)/(945.*k33**5) - &
		 -  (c**8*w**6)/(4725.*k33**7) + (2*c**10*w**8)/(93555.*k33**9) - &
		 -  (1382*c**12*w**10)/(6.38512875e8*k33**11)

		  M_A(3,6,i)= 0.5 + (c*w)/(6.*k33) - (c**3*w**3)/(90.*k33**3) + (c**5*w**5)/(945.*k33**5) - &
		 -  (c**7*w**7)/(9450.*k33**7) + (c**9*w**9)/(93555.*k33**9)

		endif
		
	enddo
  !$OMP     END PARALLEL DO 

	!$OMP     PARALLEL DO PRIVATE (i,i2)
	do i=1,NMESH
		M_F4(2,i)=0.5*(1./M_A(1,5,i)+1./M_A(2,5,i)+1./M_A(3,5,i)-1./DT)
		M_F4(3,i)=-M_A(1,6,i)/M_A(1,5,i)
		M_F4(4,i)=(M_A(1,6,i)-1.)/M_A(1,5,i)
		M_F4(5,i)=-M_A(2,6,i)/M_A(2,5,i)
		M_F4(6,i)=(M_A(2,6,i)-1.)/M_A(2,5,i)
		M_F4(7,i)=-M_A(3,6,i)/M_A(3,5,i)
		M_F4(8,i)=(M_A(3,6,i)-1.)/M_A(3,5,i)
		M_F4(1,i)=M_F4(2,i)+M_F4(3,i)+M_F4(4,i)+M_F4(5,i)+M_F4(6,i)+M_F4(7,i)+M_F4(8,i)


		if (meshes_neighbors(2,i)/=0) THEN
			i2=meshes_neighbors(2,i)
			M_F1(2,i)=M_A(1,3,i)-(M_A(1,1,i)*M_A(1,6,i)/M_A(1,5,i))
			M_F1(3,i)=M_A(1,4,i2)-M_A(1,2,i2)/M_A(1,5,i2)+M_A(1,2,i2)*M_A(1,6,i2)/M_A(1,5,i2)
			M_F1(4,i)=M_A(1,1,i)/(2.*M_A(1,5,i))
			M_F1(5,i)=M_A(1,2,i2)/(2.*M_A(1,5,i2))
			M_F1(1,i)=M_F1(2,i)+M_F1(3,i)+2.*M_F1(4,i)+2.*M_F1(5,i)
		endif
		if (meshes_neighbors(4,i)/=0) THEN
			i2=meshes_neighbors(4,i)
			M_F2(2,i)=M_A(2,3,i)-(M_A(2,1,i)*M_A(2,6,i)/M_A(2,5,i))
			M_F2(3,i)=M_A(2,4,i2)-M_A(2,2,i2)/M_A(2,5,i2)+M_A(2,2,i2)*M_A(2,6,i2)/M_A(2,5,i2)
			M_F2(4,i)=M_A(2,1,i)/(2.*M_A(2,5,i))
			M_F2(5,i)=M_A(2,2,i2)/(2.*M_A(2,5,i2))
			M_F2(1,i)=M_F2(2,i)+M_F2(3,i)+2.*M_F2(4,i)+2.*M_F2(5,i)
		endif
		if (meshes_neighbors(6,i)/=0) THEN
			i2=meshes_neighbors(6,i)
			M_F3(2,i)=M_A(3,3,i)-(M_A(3,1,i)*M_A(3,6,i)/M_A(3,5,i))
			M_F3(3,i)=M_A(3,4,i2)-M_A(3,2,i2)/M_A(3,5,i2)+M_A(3,2,i2)*M_A(3,6,i2)/M_A(3,5,i2)
			M_F3(4,i)=M_A(3,1,i)/(2.*M_A(3,5,i))
			M_F3(5,i)=M_A(3,2,i2)/(2.*M_A(3,5,i2))
			M_F3(1,i)=M_F3(2,i)+M_F3(3,i)+2.*M_F3(4,i)+2.*M_F3(5,i)
		endif
	enddo
  !$OMP     END PARALLEL DO 
  end subroutine calculate_velocity_nim_constants_formulation2



  subroutine calculate_temperature_nim_constants_formulation2
! -----------------------------------
    integer                         :: i,i2,M
	real							:: Pe_x,Pe_y,Pe_z,u,v,w,ut,vt,wt,a,b,c,k11,k22,k33
	!$OMP     PARALLEL DO PRIVATE (i,Pe_x,Pe_y,Pe_z,u,v,w,ut,vt,wt,a,b,c,k11,k22,k33)
	do i=1,NMESH
		a=1.
		b=1.
		c=1.

		ut=0.5*(meshes_uxyz(i)+meshes_uxyz_old(i))
		vt=0.5*(meshes_vxyz(i)+meshes_vxyz_old(i))
		wt=0.5*(meshes_wxyz(i)+meshes_wxyz_old(i))
		u=ut*meshes_xi_x(1,1,i)+vt*meshes_xi_x(1,2,i)+wt*meshes_xi_x(1,3,i)
		v=ut*meshes_xi_x(2,1,i)+vt*meshes_xi_x(2,2,i)+wt*meshes_xi_x(2,3,i)
		w=ut*meshes_xi_x(3,1,i)+vt*meshes_xi_x(3,2,i)+wt*meshes_xi_x(3,3,i)
		k11=meshes_gij(1,1,i)*alpha
		k22=meshes_gij(2,2,i)*alpha
		k33=meshes_gij(3,3,i)*alpha
		
		Pe_x=2*a*u/k11
		Pe_y=2*b*v/k22
		Pe_z=2*c*w/k33
		if(abs(Pe_x)>v_limit) then
		  MT_A(1,1,i)=(2*a/(k11*(1./exp(Pe_x)-1.))+1./u)*meshes_gij_faces(1,1,2,i)
		  MT_A(1,2,i)=(-2*a/(k11*(1.-exp(Pe_x)))-1./u)*meshes_gij_faces(1,1,1,i)
		  MT_A(1,3,i)=(u/(k11*(1./exp(Pe_x)-1.)))*meshes_gij_faces(1,1,2,i)
		  MT_A(1,4,i)=(u/(k11*(1.-exp(Pe_x))))*meshes_gij_faces(1,1,1,i)
		  
		  if(Pe_x>0) MT_A(1,5,i)=(-k11+(a*u*(1./exp(Pe_x)+1.)/(-1./exp(Pe_x)+1.)))/(u**2)
		  if(Pe_x<0) MT_A(1,5,i)=(-k11+(a*u*(1.+exp(Pe_x))/(-1.+exp(Pe_x))))/(u**2)

		  if(Pe_x>0) MT_A(1,6,i)=1./(-1./exp(Pe_x)+1.)-1./(Pe_x)
		  if(Pe_x<0) MT_A(1,6,i)=exp(Pe_x)/(-1.+exp(Pe_x))-1./(Pe_x)
		else

		  MT_A(1,1,i)= (-(a/k11) - (a**2*u)/(3.*k11**2) + (a**4*u**3)/(45.*k11**4) - (2*a**6*u**5)/(945.*k11**6) + &
		 -  (a**8*u**7)/(4725.*k11**8) - (2*a**10*u**9)/(93555.*k11**10))*meshes_gij_faces(1,1,2,i)

		  MT_A(1,2,i)=(-(a/k11) + (a**2*u)/(3.*k11**2) - (a**4*u**3)/(45.*k11**4) + (2*a**6*u**5)/(945.*k11**6) - &
		 -  (a**8*u**7)/(4725.*k11**8) + (2*a**10*u**9)/(93555.*k11**10))*meshes_gij_faces(1,1,1,i)

		  MT_A(1,3,i)= (-1/(2.*a) - u/(2.*k11) - (a*u**2)/(6.*k11**2) + (a**3*u**4)/(90.*k11**4) - &
		 -  (a**5*u**6)/(945.*k11**6) + (a**7*u**8)/(9450.*k11**8) - (a**9*u**10)/(93555.*k11**10))*meshes_gij_faces(1,1,2,i)

		  MT_A(1,4,i)=(-1/(2.*a) + u/(2.*k11) - (a*u**2)/(6.*k11**2) + (a**3*u**4)/(90.*k11**4) - &
		 -  (a**5*u**6)/(945.*k11**6) + (a**7*u**8)/(9450.*k11**8) - (a**9*u**10)/(93555.*k11**10))*meshes_gij_faces(1,1,1,i)

		  MT_A(1,5,i)= a**2/(3.*k11) - (a**4*u**2)/(45.*k11**3) + (2*a**6*u**4)/(945.*k11**5) - &
		 -  (a**8*u**6)/(4725.*k11**7) + (2*a**10*u**8)/(93555.*k11**9) - &
		 -  (1382*a**12*u**10)/(6.38512875e8*k11**11)

		  MT_A(1,6,i)= 0.5 + (a*u)/(6.*k11) - (a**3*u**3)/(90.*k11**3) + (a**5*u**5)/(945.*k11**5) - &
		 -  (a**7*u**7)/(9450.*k11**7) + (a**9*u**9)/(93555.*k11**9)

		endif

		if(abs(Pe_y)>v_limit) then
		  MT_A(2,1,i)=(2*b/(k22*(1./exp(Pe_y)-1.))+1./v)*meshes_gij_faces(2,2,4,i)
		  MT_A(2,2,i)=(-2*b/(k22*(1-exp(Pe_y)))-1./v)*meshes_gij_faces(2,2,3,i)
		  MT_A(2,3,i)=(v/(k22*(1./exp(Pe_y)-1.)))*meshes_gij_faces(2,2,4,i)
		  MT_A(2,4,i)=(v/(k22*(1-exp(Pe_y))))*meshes_gij_faces(2,2,3,i)

		  if(Pe_y>0) MT_A(2,5,i)=(-k22+(b*v*(1./exp(Pe_y)+1.)/(-1./exp(Pe_y)+1.)))/(v**2)
		  if(Pe_y<0) MT_A(2,5,i)=(-k22+(b*v*(1+exp(Pe_y))/(-1.+exp(Pe_y))))/(v**2)

		  if(Pe_y>0) MT_A(2,6,i)=1./(-1./exp(Pe_y)+1.)-1./(Pe_y)
		  if(Pe_y<0) MT_A(2,6,i)=exp(Pe_y)/(-1+exp(Pe_y))-1./(Pe_y)
		!      print*, MT_i,MT_j,MT_B1,MT_B2,MT_B3,MT_B4,MT_D1,Pe_y


		else
		  MT_A(2,1,i)= (-(b/k22) - (b**2*v)/(3.*k22**2) + (b**4*v**3)/(45.*k22**4) - (2*b**6*v**5)/(945.*k22**6) + &
		 -  (b**8*v**7)/(4725.*k22**8) - (2*b**10*v**9)/(93555.*k22**10))*meshes_gij_faces(2,2,4,i)

		  MT_A(2,2,i)=(-(b/k22) + (b**2*v)/(3.*k22**2) - (b**4*v**3)/(45.*k22**4) + (2*b**6*v**5)/(945.*k22**6) - &
		 -  (b**8*v**7)/(4725.*k22**8) + (2*b**10*v**9)/(93555.*k22**10))*meshes_gij_faces(2,2,3,i)

		  MT_A(2,3,i)= (-1/(2.*b) - v/(2.*k22) - (b*v**2)/(6.*k22**2) + (b**3*v**4)/(90.*k22**4) - &
		 -  (b**5*v**6)/(945.*k22**6) + (b**7*v**8)/(9450.*k22**8) - (b**9*v**10)/(93555.*k22**10))*meshes_gij_faces(2,2,4,i)

		  MT_A(2,4,i)=(-1/(2.*b) + v/(2.*k22) - (b*v**2)/(6.*k22**2) + (b**3*v**4)/(90.*k22**4) - &
		 -  (b**5*v**6)/(945.*k22**6) + (b**7*v**8)/(9450.*k22**8) - (b**9*v**10)/(93555.*k22**10))*meshes_gij_faces(2,2,3,i)

		  MT_A(2,5,i)= b**2/(3.*k22) - (b**4*v**2)/(45.*k22**3) + (2*b**6*v**4)/(945.*k22**5) - &
		 -  (b**8*v**6)/(4725.*k22**7) + (2*b**10*v**8)/(93555.*k22**9) - &
		 -  (1382*b**12*v**10)/(6.38512875e8*k22**11)

		  MT_A(2,6,i)= 0.5 + (b*v)/(6.*k22) - (b**3*v**3)/(90.*k22**3) + (b**5*v**5)/(945.*k22**5) - &
		 -  (b**7*v**7)/(9450.*k22**7) + (b**9*v**9)/(93555.*k22**9)

		endif
		if(abs(Pe_z)>v_limit) then
		  MT_A(3,1,i)=(2*c/(k33*(1./exp(Pe_z)-1.))+1./w)*meshes_gij_faces(3,3,6,i)
		  MT_A(3,2,i)=(-2*c/(k33*(1-exp(Pe_z)))-1./w)*meshes_gij_faces(3,3,5,i)
		  MT_A(3,3,i)=(w/(k33*(1./exp(Pe_z)-1.)))*meshes_gij_faces(3,3,6,i)
		  MT_A(3,4,i)=(w/(k33*(1-exp(Pe_z))))*meshes_gij_faces(3,3,5,i)

		  if(Pe_z>0) MT_A(3,5,i)=(-k33+(c*w*(1./exp(Pe_z)+1.)/(-1./exp(Pe_z)+1.)))/(w**2)
		  if(Pe_z<0) MT_A(3,5,i)=(-k33+(c*w*(1+exp(Pe_z))/(-1.+exp(Pe_z))))/(w**2)

		  if(Pe_z>0) MT_A(3,6,i)=1./(-1./exp(Pe_z)+1.)-1./(Pe_z)
		  if(Pe_z<0) MT_A(3,6,i)=exp(Pe_z)/(-1+exp(Pe_z))-1./(Pe_z)
		!      print*, MT_i,MT_j,MT_B1,MT_B2,MT_B3,MT_B4,MT_D1,Pe_z


		else
		  MT_A(3,1,i)= (-(c/k33) - (c**2*w)/(3.*k33**2) + (c**4*w**3)/(45.*k33**4) - (2*c**6*w**5)/(945.*k33**6) + &
		 -  (c**8*w**7)/(4725.*k33**8) - (2*c**10*w**9)/(93555.*k33**10))*meshes_gij_faces(3,3,6,i)

		  MT_A(3,2,i)=(-(c/k33) + (c**2*w)/(3.*k33**2) - (c**4*w**3)/(45.*k33**4) + (2*c**6*w**5)/(945.*k33**6) - &
		 -  (c**8*w**7)/(4725.*k33**8) + (2*c**10*w**9)/(93555.*k33**10))*meshes_gij_faces(3,3,5,i)
		 
		  MT_A(3,3,i)= (-1/(2.*c) - w/(2.*k33) - (c*w**2)/(6.*k33**2) + (c**3*w**4)/(90.*k33**4) - &
		 -  (c**5*w**6)/(945.*k33**6) + (c**7*w**8)/(9450.*k33**8) - (c**9*w**10)/(93555.*k33**10))*meshes_gij_faces(3,3,6,i)

		  MT_A(3,4,i)=(-1/(2.*c) + w/(2.*k33) - (c*w**2)/(6.*k33**2) + (c**3*w**4)/(90.*k33**4) - &
		 -  (c**5*w**6)/(945.*k33**6) + (c**7*w**8)/(9450.*k33**8) - (c**9*w**10)/(93555.*k33**10))*meshes_gij_faces(3,3,5,i)

		  MT_A(3,5,i)= c**2/(3.*k33) - (c**4*w**2)/(45.*k33**3) + (2*c**6*w**4)/(945.*k33**5) - &
		 -  (c**8*w**6)/(4725.*k33**7) + (2*c**10*w**8)/(93555.*k33**9) - &
		 -  (1382*c**12*w**10)/(6.38512875e8*k33**11)

		  MT_A(3,6,i)= 0.5 + (c*w)/(6.*k33) - (c**3*w**3)/(90.*k33**3) + (c**5*w**5)/(945.*k33**5) - &
		 -  (c**7*w**7)/(9450.*k33**7) + (c**9*w**9)/(93555.*k33**9)

		endif
		
	enddo
  !$OMP     END PARALLEL DO 

	!$OMP     PARALLEL DO PRIVATE (i,i2)
	do i=1,NMESH
		MT_F4(2,i)=0.5*(1./MT_A(1,5,i)+1./MT_A(2,5,i)+1./MT_A(3,5,i)-1./DT)
		MT_F4(3,i)=-MT_A(1,6,i)/MT_A(1,5,i)
		MT_F4(4,i)=(MT_A(1,6,i)-1.)/MT_A(1,5,i)
		MT_F4(5,i)=-MT_A(2,6,i)/MT_A(2,5,i)
		MT_F4(6,i)=(MT_A(2,6,i)-1.)/MT_A(2,5,i)
		MT_F4(7,i)=-MT_A(3,6,i)/MT_A(3,5,i)
		MT_F4(8,i)=(MT_A(3,6,i)-1.)/MT_A(3,5,i)
		MT_F4(1,i)=MT_F4(2,i)+MT_F4(3,i)+MT_F4(4,i)+MT_F4(5,i)+MT_F4(6,i)+MT_F4(7,i)+MT_F4(8,i)


		if (meshes_neighbors(2,i)/=0) THEN
			i2=meshes_neighbors(2,i)
			MT_F1(2,i)=MT_A(1,3,i)-(MT_A(1,1,i)*MT_A(1,6,i)/MT_A(1,5,i))
			MT_F1(3,i)=MT_A(1,4,i2)-MT_A(1,2,i2)/MT_A(1,5,i2)+MT_A(1,2,i2)*MT_A(1,6,i2)/MT_A(1,5,i2)
			MT_F1(4,i)=MT_A(1,1,i)/(2.*MT_A(1,5,i))
			MT_F1(5,i)=MT_A(1,2,i2)/(2.*MT_A(1,5,i2))
			MT_F1(1,i)=MT_F1(2,i)+MT_F1(3,i)+2.*MT_F1(4,i)+2.*MT_F1(5,i)
		endif
		if (meshes_neighbors(4,i)/=0) THEN
			i2=meshes_neighbors(4,i)
			MT_F2(2,i)=MT_A(2,3,i)-(MT_A(2,1,i)*MT_A(2,6,i)/MT_A(2,5,i))
			MT_F2(3,i)=MT_A(2,4,i2)-MT_A(2,2,i2)/MT_A(2,5,i2)+MT_A(2,2,i2)*MT_A(2,6,i2)/MT_A(2,5,i2)
			MT_F2(4,i)=MT_A(2,1,i)/(2.*MT_A(2,5,i))
			MT_F2(5,i)=MT_A(2,2,i2)/(2.*MT_A(2,5,i2))
			MT_F2(1,i)=MT_F2(2,i)+MT_F2(3,i)+2.*MT_F2(4,i)+2.*MT_F2(5,i)
		endif
		if (meshes_neighbors(6,i)/=0) THEN
			i2=meshes_neighbors(6,i)
			MT_F3(2,i)=MT_A(3,3,i)-(MT_A(3,1,i)*MT_A(3,6,i)/MT_A(3,5,i))
			MT_F3(3,i)=MT_A(3,4,i2)-MT_A(3,2,i2)/MT_A(3,5,i2)+MT_A(3,2,i2)*MT_A(3,6,i2)/MT_A(3,5,i2)
			MT_F3(4,i)=MT_A(3,1,i)/(2.*MT_A(3,5,i))
			MT_F3(5,i)=MT_A(3,2,i2)/(2.*MT_A(3,5,i2))
			MT_F3(1,i)=MT_F3(2,i)+MT_F3(3,i)+2.*MT_F3(4,i)+2.*MT_F3(5,i)
		endif
	enddo
  !$OMP     END PARALLEL DO 






!$OMP     PARALLEL DO PRIVATE (i,M)
    do i=1,tot_faces_bcs
! ==================== |F|A|C|E| |2| =====================
        M=BCs_faces_mesh_num(i)

!                if (BCs_faces_position(i)==1) THEN
!                        BCs_faces_aT(2,i)=MT_A(1,3,M)-MT_A(1,4,M)/MT_A(1,5,M)+MT_A(1,4,M)*MT_A(1,6,M)/MT_A(1,5,M)
!                        BCs_faces_aT(3,i)=MT_A(1,4,M)/(2.*MT_A(1,5,M))
!                        BCs_faces_aT(1,i)=BCs_faces_aT(2,i)+2*BCs_faces_aT(3,i)
!                elseif(BCs_faces_position(i)==2) THEN
!                        BCs_faces_aT(2,i)=MT_A(1,1,M)-(MT_A(1,2,M)*MT_A(1,6,M)/MT_A(1,5,M))
!                        BCs_faces_aT(3,i)=MT_A(1,2,M)/(2.*MT_A(1,5,M))
!                        BCs_faces_aT(1,i)=BCs_faces_aT(2,i)+2*BCs_faces_aT(3,i)
                if (BCs_faces_position(i)==1) THEN
                        BCs_faces_aT(2,i)=MT_A(1,4,M)-MT_A(1,2,M)/MT_A(1,5,M)+MT_A(1,2,M)*MT_A(1,6,M)/MT_A(1,5,M)
                        BCs_faces_aT(3,i)=MT_A(1,2,M)/(2.*MT_A(1,5,M))
                        BCs_faces_aT(1,i)=BCs_faces_aT(2,i)+2*BCs_faces_aT(3,i)
                elseif(BCs_faces_position(i)==2) THEN
                        BCs_faces_aT(2,i)=MT_A(1,3,M)-(MT_A(1,1,M)*MT_A(1,6,M)/MT_A(1,5,M))
                        BCs_faces_aT(3,i)=MT_A(1,1,M)/(2.*MT_A(1,5,M))
                        BCs_faces_aT(1,i)=BCs_faces_aT(2,i)+2*BCs_faces_aT(3,i)

                elseif(BCs_faces_position(i)==3) THEN
                        BCs_faces_aT(2,i)=MT_A(2,3,M)-MT_A(2,4,M)/MT_A(2,5,M)+MT_A(2,4,M)*MT_A(2,6,M)/MT_A(2,5,M)
                        BCs_faces_aT(3,i)=MT_A(2,4,M)/(2.*MT_A(2,5,M))
                        BCs_faces_aT(1,i)=BCs_faces_aT(2,i)+2*BCs_faces_aT(3,i)
                elseif(BCs_faces_position(i)==4) THEN
                        BCs_faces_aT(2,i)=MT_A(2,1,M)-(MT_A(2,2,M)*MT_A(2,6,M)/MT_A(2,5,M))
                        BCs_faces_aT(3,i)=MT_A(2,2,M)/(2.*MT_A(2,5,M))
                        BCs_faces_aT(1,i)=BCs_faces_aT(2,i)+2*BCs_faces_aT(3,i)
                elseif(BCs_faces_position(i)==5) THEN
                        BCs_faces_aT(2,i)=MT_A(3,3,M)-MT_A(3,4,M)/MT_A(3,5,M)+MT_A(3,4,M)*MT_A(3,6,M)/MT_A(3,5,M)
                        BCs_faces_aT(3,i)=MT_A(3,4,M)/(2.*MT_A(3,5,M))
                        BCs_faces_aT(1,i)=BCs_faces_aT(2,i)+2*BCs_faces_aT(3,i)
                elseif(BCs_faces_position(i)==6) THEN
                        BCs_faces_aT(2,i)=MT_A(3,1,M)-(MT_A(3,2,M)*MT_A(3,6,M)/MT_A(3,5,M))
                        BCs_faces_aT(3,i)=MT_A(3,2,M)/(2.*MT_A(3,5,M))
                        BCs_faces_aT(1,i)=BCs_faces_aT(2,i)+2*BCs_faces_aT(3,i)
                endif
    enddo
  !$OMP     END PARALLEL DO



  end subroutine calculate_temperature_nim_constants_formulation2


  subroutine calculate_pressure_nim_constants
    implicit none
    integer :: i,i2,M
	real    :: B_x,B_y,B_z,a,b,c,Re_x,Re_y,Re_z,K11,K22,K33
	real,dimension(1:NMESH)    :: Det
	real,dimension(1:3,1:NMESH)    :: CC,EE,FF,GG,HH,II,JJ

    a=1.
    b=1.
    c=1.	

	!$OMP     PARALLEL DO PRIVATE (i,Re_x,Re_y,Re_z,B_x,B_y,B_z,K11,K22,K33)
    do i=1,NMESH
	
        K11=meshes_gij(1,1,i)/rho
		K22=meshes_gij(2,2,i)/rho
		K33=meshes_gij(3,3,i)/rho
		
		Re_x=0.
		Re_y=0.
		Re_z=0.

		B_x=B_nim(0.,K11, 0.)
		B_y=B_nim(0.,K22, 0.)
		B_z=B_nim(0.,K33, 0.)

		CC(1,i)=C_nim(0.,K11,0.,B_x,a)
		CC(2,i)=C_nim(0.,K22,0.,B_y,b)
		CC(3,i)=C_nim(0.,K33,0.,B_z,c)

		Det(i)=1./(CC(1,i)*CC(2,i)+CC(2,i)*CC(3,i)+CC(1,i)*CC(3,i))

		EE(1,i)=E_nim(0.,K11,0.,B_x,Det(i),a)
		EE(2,i)=E_nim(0.,K22,0.,B_y,Det(i),b)
		EE(3,i)=E_nim(0.,K33,0.,B_z,Det(i),c)

		FF(1,i)=F_nim(0.,K11,0.,B_x,Det(i),a)
		FF(2,i)=F_nim(0.,K22,0.,B_y,Det(i),b)
		FF(3,i)=F_nim(0.,K33,0.,B_z,Det(i),c)

		GG(1,i)=G_nim(0.,K11,0.,a)
		GG(2,i)=G_nim(0.,K22,0.,b)
		GG(3,i)=G_nim(0.,K33,0.,c)

		HH(1,i)=H_nim(0.,K11,0.,a)
		HH(2,i)=H_nim(0.,K22,0.,b)
		HH(3,i)=H_nim(0.,K33,0.,c)

		II(1,i)=I_nim(0.,K11,0.,a)
		II(2,i)=I_nim(0.,K22,0.,b)
		II(3,i)=I_nim(0.,K33,0.,c)

		JJ(1,i)=J_nim(0.,K11,0.,a)
		JJ(2,i)=J_nim(0.,K22,0.,b)
		JJ(3,i)=J_nim(0.,K33,0.,c)

    enddo
  !$OMP     END PARALLEL DO 

!$OMP     PARALLEL DO PRIVATE (i,i2)
    do i=1,NMESH
! ==================== |F|A|C|E| |2| =====================
    if (meshes_neighbors(2,i)/=0) THEN
       i2=meshes_neighbors(2,i)
!       M2=MESHES_ARRAY(MESHES_ARRAY(i)%neighbors(2))%const_V
       AA_px(1,i)=-(GG(1,i)+(CC(2,i)+CC(3,i))*EE(1,i)*HH(1,i))*meshes_gij_faces(1,1,2,i)
       AA_px(2,i)= (GG(1,i)-(CC(2,i)+CC(3,i))*FF(1,i)*HH(1,i))*meshes_gij_faces(1,1,2,i) &
                                    +(II(1,i2)+(CC(2,i2)+CC(3,i2))*EE(1,i2)*JJ(1,i2))*meshes_gij_faces(1,1,1,i2)
       AA_px(3,i)=-(II(1,i2)-(CC(2,i2)+CC(3,i2))*FF(1,i2)*JJ(1,i2))*meshes_gij_faces(1,1,1,i2)
       AA_px(4,i)=  CC(3,i)*EE(2,i)*HH(1,i)*meshes_gij_faces(1,1,2,i)
       AA_px(5,i)=  CC(3,i)*FF(2,i)*HH(1,i)*meshes_gij_faces(1,1,2,i)
       AA_px(6,i)= -CC(3,i2)*EE(2,i2)*JJ(1,i2)*meshes_gij_faces(1,1,1,i2)
       AA_px(7,i)= -CC(3,i2)*FF(2,i2)*JJ(1,i2)*meshes_gij_faces(1,1,1,i2)
       AA_px(8,i)=  CC(2,i)*EE(3,i)*HH(1,i)*meshes_gij_faces(1,1,2,i)        
       AA_px(9,i)=  CC(2,i)*FF(3,i)*HH(1,i)*meshes_gij_faces(1,1,2,i)
       AA_px(10,i)=-CC(2,i2)*EE(3,i2)*JJ(1,i2)*meshes_gij_faces(1,1,1,i2)
       AA_px(11,i)=-CC(2,i2)*FF(3,i2)*JJ(1,i2)*meshes_gij_faces(1,1,1,i2)
       AA_px(12,i)= CC(2,i)*CC(3,i)*Det(i)*HH(1,i)*meshes_gij_faces(1,1,2,i)
       AA_px(13,i)=-CC(2,i)*CC(3,i)*Det(i)*HH(1,i)*meshes_gij_faces(1,1,2,i)
       AA_px(14,i)=-CC(2,i2)*CC(3,i2)*Det(i2)*JJ(1,i2)*meshes_gij_faces(1,1,1,i2)
       AA_px(15,i)= CC(2,i2)*CC(3,i2)*Det(i2)*JJ(1,i2)*meshes_gij_faces(1,1,1,i2)
    endif

! ==================== |F|A|C|E| |4| =====================
    if (meshes_neighbors(4,i)/=0) THEN
       i2=meshes_neighbors(4,i)
       AA_py(1,i)=  CC(3,i)*EE(1,i)*HH(2,i)*meshes_gij_faces(2,2,4,i)
       AA_py(2,i)=  CC(3,i)*FF(1,i)*HH(2,i)*meshes_gij_faces(2,2,4,i)
       AA_py(3,i)= -CC(3,i2)*EE(1,i2)*JJ(2,i2)*meshes_gij_faces(2,2,3,i2)
       AA_py(4,i)= -CC(3,i2)*FF(1,i2)*JJ(2,i2)*meshes_gij_faces(2,2,3,i2)
       AA_py(5,i)=-(GG(2,i)+(CC(1,i)+CC(3,i))*EE(2,i)*HH(2,i))*meshes_gij_faces(2,2,4,i)
       AA_py(6,i)= (GG(2,i)-(CC(1,i)+CC(3,i))*FF(2,i)*HH(2,i))*meshes_gij_faces(2,2,4,i) &
                                    +(II(2,i2)+(CC(1,i2)+CC(3,i2))*EE(2,i2)*JJ(2,i2))*meshes_gij_faces(2,2,3,i2)
       AA_py(7,i)=-(II(2,i2)-(CC(1,i2)+CC(3,i2))*FF(2,i2)*JJ(2,i2))*meshes_gij_faces(2,2,3,i2)
       AA_py(8,i)=  CC(1,i)*EE(3,i)*HH(2,i)*meshes_gij_faces(2,2,4,i)
       AA_py(9,i)=  CC(1,i)*FF(3,i)*HH(2,i)*meshes_gij_faces(2,2,4,i)
       AA_py(10,i)=-CC(1,i2)*EE(3,i2)*JJ(2,i2)*meshes_gij_faces(2,2,3,i2)
       AA_py(11,i)=-CC(1,i2)*FF(3,i2)*JJ(2,i2)*meshes_gij_faces(2,2,3,i2)
       AA_py(12,i)= CC(1,i)*CC(3,i)*Det(i)*HH(2,i)*meshes_gij_faces(2,2,4,i)
       AA_py(13,i)=-CC(1,i)*CC(3,i)*Det(i)*HH(2,i)*meshes_gij_faces(2,2,4,i)
       AA_py(14,i)=-CC(1,i2)*CC(3,i2)*Det(i2)*JJ(2,i2)*meshes_gij_faces(2,2,3,i2)
       AA_py(15,i)= CC(1,i2)*CC(3,i2)*Det(i2)*JJ(2,i2)*meshes_gij_faces(2,2,3,i2)
    endif
! ==================== |F|A|C|E| |6| =====================
    if (meshes_neighbors(6,i)/=0) THEN
       i2=meshes_neighbors(6,i)
       AA_pz(1,i) =  CC(2,i)*EE(1,i)*HH(3,i)*meshes_gij_faces(3,3,6,i)
       AA_pz(2,i) =  CC(2,i)*FF(1,i)*HH(3,i)*meshes_gij_faces(3,3,6,i)
       AA_pz(3,i) = -CC(2,i2)*EE(1,i2)*JJ(3,i2)*meshes_gij_faces(3,3,5,i2)
       AA_pz(4,i) = -CC(2,i2)*FF(1,i2)*JJ(3,i2)*meshes_gij_faces(3,3,5,i2)
       AA_pz(5,i) =  CC(1,i)*EE(2,i)*HH(3,i)*meshes_gij_faces(3,3,6,i)
       AA_pz(6,i) =  CC(1,i)*FF(2,i)*HH(3,i)*meshes_gij_faces(3,3,6,i)
       AA_pz(7,i) = -CC(1,i2)*EE(2,i2)*JJ(3,i2)*meshes_gij_faces(3,3,5,i2)
       AA_pz(8,i) = -CC(1,i2)*FF(2,i2)*JJ(3,i2)*meshes_gij_faces(3,3,5,i2)
       AA_pz(9,i) =-(GG(3,i)+(CC(1,i)+CC(2,i))*EE(3,i)*HH(3,i))*meshes_gij_faces(3,3,6,i)
       AA_pz(10,i) =(GG(3,i)-(CC(1,i)+CC(2,i))*FF(3,i)*HH(3,i))*meshes_gij_faces(3,3,6,i) &
                                     +(II(3,i2)+(CC(1,i2)+CC(2,i2))*EE(3,i2)*JJ(3,i2))*meshes_gij_faces(3,3,5,i2)
       AA_pz(11,i)=-(II(3,i2)-(CC(1,i2)+CC(2,i2))*FF(3,i2)*JJ(3,i2))*meshes_gij_faces(3,3,5,i2)
       AA_pz(12,i)= CC(1,i)*CC(2,i)*Det(i)*HH(3,i)*meshes_gij_faces(3,3,6,i)
       AA_pz(13,i)=-CC(1,i)*CC(2,i)*Det(i)*HH(3,i)*meshes_gij_faces(3,3,6,i)
       AA_pz(14,i)=-CC(1,i2)*CC(2,i2)*Det(i2)*JJ(3,i2)*meshes_gij_faces(3,3,5,i2)
       AA_pz(15,i)= CC(1,i2)*CC(2,i2)*Det(i2)*JJ(3,i2)*meshes_gij_faces(3,3,5,i2)
    endif
    enddo
  !$OMP     END PARALLEL DO 

!$OMP     PARALLEL DO PRIVATE (i,M)
    do i=1,tot_faces_bcs
! ==================== |F|A|C|E| |2| =====================
	M=BCs_faces_mesh_num(i)
    if (BCs_faces_position(i)==1) THEN
		BCs_faces_ap(1,i)=-(II(1,M)+(CC(2,M)+CC(3,M))*EE(1,M)*JJ(1,M))*meshes_gij_faces(1,1,1,M)
		BCs_faces_ap(2,i)= (II(1,M)-(CC(2,M)+CC(3,M))*FF(1,M)*JJ(1,M))*meshes_gij_faces(1,1,1,M)
		BCs_faces_ap(3,i)= CC(3,M)*EE(2,M)*JJ(1,M)*meshes_gij_faces(1,1,1,M)
		BCs_faces_ap(4,i)= CC(3,M)*FF(2,M)*JJ(1,M)*meshes_gij_faces(1,1,1,M)
		BCs_faces_ap(5,i)= CC(2,M)*EE(3,M)*JJ(1,M)*meshes_gij_faces(1,1,1,M)
		BCs_faces_ap(6,i)= CC(2,M)*FF(3,M)*JJ(1,M)*meshes_gij_faces(1,1,1,M)
		BCs_faces_ap(7,i)= CC(2,M)*CC(3,M)*Det(M)*JJ(1,M)*meshes_gij_faces(1,1,1,M)
		BCs_faces_ap(8,i)=-CC(2,M)*CC(3,M)*Det(M)*JJ(1,M)*meshes_gij_faces(1,1,1,M)
	elseif(BCs_faces_position(i)==2) THEN
		BCs_faces_ap(1,i)=-(GG(1,M)+(CC(2,M)+CC(3,M))*EE(1,M)*HH(1,M))*meshes_gij_faces(1,1,2,M)
		BCs_faces_ap(2,i)= (GG(1,M)-(CC(2,M)+CC(3,M))*FF(1,M)*HH(1,M))*meshes_gij_faces(1,1,2,M)
		BCs_faces_ap(3,i)= CC(3,M)*EE(2,M)*HH(1,M)*meshes_gij_faces(1,1,2,M)
		BCs_faces_ap(4,i)= CC(3,M)*FF(2,M)*HH(1,M)*meshes_gij_faces(1,1,2,M)
		BCs_faces_ap(5,i)= CC(2,M)*EE(3,M)*HH(1,M)*meshes_gij_faces(1,1,2,M)  
		BCs_faces_ap(6,i)= CC(2,M)*FF(3,M)*HH(1,M)*meshes_gij_faces(1,1,2,M)
		BCs_faces_ap(7,i)= CC(2,M)*CC(3,M)*Det(M)*HH(1,M)*meshes_gij_faces(1,1,2,M)
		BCs_faces_ap(8,i)=-CC(2,M)*CC(3,M)*Det(M)*HH(1,M)*meshes_gij_faces(1,1,2,M)
	elseif(BCs_faces_position(i)==3) THEN
		BCs_faces_ap(1,i)= CC(3,M)*EE(1,M)*JJ(2,M)*meshes_gij_faces(2,2,3,M)
		BCs_faces_ap(2,i)= CC(3,M)*FF(1,M)*JJ(2,M)*meshes_gij_faces(2,2,3,M)
		BCs_faces_ap(3,i)=-(II(2,M)+(CC(1,M)+CC(3,M))*EE(2,M)*JJ(2,M))*meshes_gij_faces(2,2,3,M)
		BCs_faces_ap(4,i)= (II(2,M)-(CC(1,M)+CC(3,M))*FF(2,M)*JJ(2,M))*meshes_gij_faces(2,2,3,M)
		BCs_faces_ap(5,i)= CC(1,M)*EE(3,M)*JJ(2,M)*meshes_gij_faces(2,2,3,M)
		BCs_faces_ap(6,i)= CC(1,M)*FF(3,M)*JJ(2,M)*meshes_gij_faces(2,2,3,M)
		BCs_faces_ap(7,i)= CC(1,M)*CC(3,M)*Det(M)*JJ(2,M)*meshes_gij_faces(2,2,3,M)
		BCs_faces_ap(8,i)=-CC(1,M)*CC(3,M)*Det(M)*JJ(2,M)*meshes_gij_faces(2,2,3,M)
	elseif(BCs_faces_position(i)==4) THEN
		BCs_faces_ap(1,i)= CC(3,M)*EE(1,M)*HH(2,M)*meshes_gij_faces(2,2,4,M)
		BCs_faces_ap(2,i)= CC(3,M)*FF(1,M)*HH(2,M)*meshes_gij_faces(2,2,4,M)
		BCs_faces_ap(3,i)=-(GG(2,M)+(CC(1,M)+CC(3,M))*EE(2,M)*HH(2,M))*meshes_gij_faces(2,2,4,M)
		BCs_faces_ap(4,i)= (GG(2,M)-(CC(1,M)+CC(3,M))*FF(2,M)*HH(2,M))*meshes_gij_faces(2,2,4,M)
		BCs_faces_ap(5,i)= CC(1,M)*EE(3,M)*HH(2,M)*meshes_gij_faces(2,2,4,M)
		BCs_faces_ap(6,i)= CC(1,M)*FF(3,M)*HH(2,M)*meshes_gij_faces(2,2,4,M)
		BCs_faces_ap(7,i)= CC(1,M)*CC(3,M)*Det(M)*HH(2,M)*meshes_gij_faces(2,2,4,M)
		BCs_faces_ap(8,i)=-CC(1,M)*CC(3,M)*Det(M)*HH(2,M)*meshes_gij_faces(2,2,4,M)
	elseif(BCs_faces_position(i)==5) THEN
		BCs_faces_ap(1,i)= CC(2,M)*EE(1,M)*JJ(3,M)*meshes_gij_faces(3,3,5,M)
		BCs_faces_ap(2,i)= CC(2,M)*FF(1,M)*JJ(3,M)*meshes_gij_faces(3,3,5,M)
		BCs_faces_ap(3,i)= CC(1,M)*EE(2,M)*JJ(3,M)*meshes_gij_faces(3,3,5,M)
		BCs_faces_ap(4,i)= CC(1,M)*FF(2,M)*JJ(3,M)*meshes_gij_faces(3,3,5,M)
		BCs_faces_ap(5,i)=-(II(3,M)+(CC(1,M)+CC(2,M))*EE(3,M)*JJ(3,M))*meshes_gij_faces(3,3,5,M)
		BCs_faces_ap(6,i)= (II(3,M)-(CC(1,M)+CC(2,M))*FF(3,M)*JJ(3,M))*meshes_gij_faces(3,3,5,M)
		BCs_faces_ap(7,i)= CC(1,M)*CC(2,M)*Det(M)*JJ(3,M)*meshes_gij_faces(3,3,5,M)
		BCs_faces_ap(8,i)=-CC(1,M)*CC(2,M)*Det(M)*JJ(3,M)*meshes_gij_faces(3,3,5,M)
	elseif(BCs_faces_position(i)==6) THEN
		BCs_faces_ap(1,i)= CC(2,M)*EE(1,M)*HH(3,M)*meshes_gij_faces(3,3,6,M)
		BCs_faces_ap(2,i)= CC(2,M)*FF(1,M)*HH(3,M)*meshes_gij_faces(3,3,6,M)
		BCs_faces_ap(3,i)= CC(1,M)*EE(2,M)*HH(3,M)*meshes_gij_faces(3,3,6,M)
		BCs_faces_ap(4,i)= CC(1,M)*FF(2,M)*HH(3,M)*meshes_gij_faces(3,3,6,M)
		BCs_faces_ap(5,i)=-(GG(3,M)+(CC(1,M)+CC(2,M))*EE(3,M)*HH(3,M))*meshes_gij_faces(3,3,6,M)
		BCs_faces_ap(6,i)= (GG(3,M)-(CC(1,M)+CC(2,M))*FF(3,M)*HH(3,M))*meshes_gij_faces(3,3,6,M)
		BCs_faces_ap(7,i)= CC(1,M)*CC(2,M)*Det(M)*HH(3,M)*meshes_gij_faces(3,3,6,M)
		BCs_faces_ap(8,i)=-CC(1,M)*CC(2,M)*Det(M)*HH(3,M)*meshes_gij_faces(3,3,6,M)
    endif
    enddo
  !$OMP     END PARALLEL DO 
  end subroutine calculate_pressure_nim_constants



! ========================================= NIM constants =========================================
  function condition(u,F,Re) result(cond)
     real, intent(in)        :: Re,u,F
     integer                    :: cond
    if (abs(Re) <v_limit) then     ! 0.01
       cond=1              ! For this case, it is undefined.
    else
       cond=3
    endif
  end function condition

  function B_nim(u,F,Re)  result(P)
     real, intent(in)        :: Re,u,F
     real                    :: P
    select case (condition(u,F,Re))
       case (1)
            P=0         ! For this case, it is undefined.
       case (2)
            P=0
       case default
            P=1./u*(1./Re - 1./(exp(Re)-1))
 !           p=1./u*((exp(Re)-1.)/Re-1.)/(exp(Re)-1)
    end select
  end function B_nim

  function C_nim(u,F,Re,B,a)  result(P)
     real, intent(in)        :: u,F,Re,B,a
     real                    :: P
    select case (condition(u,F,Re))
       case (1)
            P=-a**2/(3.*F)+(a**4*u**2)/(45.*F**3)-(2*a**6*u**4)/(945.*F**5) + &
             -(a**8*u**6)/(4725.*F**7) - (2*a**10*u**8)/(93555.*F**9) + &
             -(1382*a**12*u**10)/(6.38512875e8*F**11)
       case (2)
            P=-1./u
       case default
            P=(-1./u + 2.*B)*a
    end select
  end function C_nim

  function E_nim(u,F,Re,B,D,a)  result(P)
     real, intent(in)        :: u,F,Re,B,D,a
     real                    :: P
    select case (condition(u,F,Re))
       case (1)
            P=D/2. + (a*D*u)/(6.*F) - (a**3*D*u**3)/(90.*F**3) + (a**5*D*u**5)/(945.*F**5) - &
            -(a**7*D*u**7)/(9450.*F**7) + (a**9*D*u**9)/(93555.*F**9)
       case (2)
            P=D
       case default
            P=D*(1-B*u)
!           P=D*(1.-(1./Re - 1./(exp(Re)-1.)))
    end select
  end function E_nim

  function F_nim(u,FF,Re,B,D,a)  result(P)
     real, intent(in)        :: u,FF,Re,B,D,a
     real                    :: P
    select case (condition(u,FF,Re))
       case (1)
            P=D/2. - (a*D*u)/(6.*FF) + (a**3*D*u**3)/(90.*FF**3) - (a**5*D*u**5)/(945.*FF**5) + &
            -  (a**7*D*u**7)/(9450.*FF**7) - (a**9*D*u**9)/(93555.*FF**9)
       case (2)
            P=0
       case default
!            P=D*(1./Re - 1./(exp(Re)-1))
            P=D*(B*u)
    end select
  end function F_nim


  function G_nim(u,F,Re,a)  result(P)
     real, intent(in)        :: u,F,Re,a
     real                    :: P
    select case (condition(u,F,Re))
       case (1)
            P=1/(2.*a) + u/(2.*F) + (a*u**2)/(6.*F**2) - (a**3*u**4)/(90.*F**4) + &
             -  (a**5*u**6)/(945.*F**6) - (a**7*u**8)/(9450.*F**8) + (a**9*u**10)/(93555.*F**10)
       case (2)
            P=0
       case default
            P=exp(Re)*u/(F*(exp(Re)-1.))
            p=u/F/(1.-1./exp(Re))
    end select
  end function G_nim

  function H_nim(u,F,Re,a)  result(P)
     real, intent(in)        :: u,F,Re,a
     real                    :: P
    select case (condition(u,F,Re))
       case (1)
            P=a/F+(a**2*u)/(3.*F**2)-(a**4*u**3)/(45.*F**4)+(2*a**6*u**5)/(945.*F**6) - &
            -(a**8*u**7)/(4725.*F**8) + (2*a**10*u**9)/(93555.*F**10)
       case (2)
            P=0
       case default
            P=(2*a*exp(Re))/(F*(exp(Re)-1))-1/u
            p=2.*a/F/(1-1/exp(Re))-1./u
    end select
  end function H_nim

  function I_nim(u,F,Re,a)  result(P)
     real, intent(in)        :: u,F,Re,a
     real                    :: P
    select case (condition(u,F,Re))
       case (1)
            P=1/(2.*a) - u/(2.*F) + (a*u**2)/(6.*F**2) - (a**3*u**4)/(90.*F**4) + &
             -(a**5*u**6)/(945.*F**6) - (a**7*u**8)/(9450.*F**8) + (a**9*u**10)/(93555.*F**10)
       case (2)
            P=0
       case default
            P=u/(F*(exp(Re)-1))
    end select
  end function I_nim

  function J_nim(u,F,Re,a)  result(P)
     real, intent(in)        :: u,F,Re,a
     real                    :: P
    select case (condition(u,F,Re))
       case (1)
            P=-(a/F) + (a**2*u)/(3.*F**2) - (a**4*u**3)/(45.*F**4) + &
              -(2*a**6*u**5)/(945.*F**6) - (a**8*u**7)/(4725.*F**8) + &
              -(2*a**10*u**9)/(93555.*F**10)
       case (2)
            P=0
       case default
            P=2.*a/(F*(exp(Re)-1.))-1./u
    end select
  end function J_nim

end module class_nim

