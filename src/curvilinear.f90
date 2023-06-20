module class_curvilinear
use class_mesh
use class_general
use class_parameters
implicit none

contains

subroutine evaluate_curvilinear(M,C,x,y,z)
  implicit none
  real, intent(in)        :: x,y,z
  type(mesh), intent(in)              :: M
  type(curvi3d),intent(inout)         :: C


  C%x=M%C(0)+M%C(1)*x+M%C(2)*y+M%C(3)*z+M%C(4)*x*y+M%C(5)*x*z+M%C(6)*y*z+M%C(7)*x*y*z
  C%y=M%D(0)+M%D(1)*x+M%D(2)*y+M%D(3)*z+M%D(4)*x*y+M%D(5)*x*z+M%D(6)*y*z+M%D(7)*x*y*z
  C%z=M%E(0)+M%E(1)*x+M%E(2)*y+M%E(3)*z+M%E(4)*x*y+M%E(5)*x*z+M%E(6)*y*z+M%E(7)*x*y*z
  

  ! Find g1, g2 vectors
  C%jacmat= jacobian_matrix(M,x,y,z)
  C%g_1=C%jacmat(1,:)
  C%g_2=C%jacmat(2,:)
  C%g_3=C%jacmat(3,:)

  ! Find g_1, g_2 vectors
  C%jacinv= inverse_jocobian_matrix(M,C%jacmat)
  C%g1=C%jacinv(:,1)
  C%g2=C%jacinv(:,2)
  C%g3=C%jacinv(:,3)

  ! Find metric matrices
  C%g_ij= g_ij_matrix(C%g_1,C%g_2,C%g_3)
  C%gij = gij_matrix (C%g1 ,C%g2 ,C%g3 )

  ! Determinants  
  C%gamma_tensor = find_gamma_tensor(M,C%jacinv,C%g1,C%g2,C%g3,x,y,z)
  C%detG = det(C%g_ij)
  C%detJ = det(C%jacmat)
end subroutine evaluate_curvilinear


function jacobian_matrix(M,x,y,z)  result(J)
  implicit none
  real, intent(in)        :: x,y,z
  type(mesh), intent(in)              :: M
  real, dimension(1:3,1:3):: J

  J(1,1)=M%C(1)+M%C(4)*y+M%C(5)*z+M%C(7)*y*z
  J(1,2)=M%D(1)+M%D(4)*y+M%D(5)*z+M%D(7)*y*z
  J(1,3)=M%E(1)+M%E(4)*y+M%E(5)*z+M%E(7)*y*z

  J(2,1)=M%C(2)+M%C(4)*x+M%C(6)*z+M%C(7)*x*z
  J(2,2)=M%D(2)+M%D(4)*x+M%D(6)*z+M%D(7)*x*z
  J(2,3)=M%E(2)+M%E(4)*x+M%E(6)*z+M%E(7)*x*z

  J(3,1)=M%C(3)+M%C(5)*x+M%C(6)*y+M%C(7)*x*y
  J(3,2)=M%D(3)+M%D(5)*x+M%D(6)*y+M%D(7)*x*y
  J(3,3)=M%E(3)+M%E(5)*x+M%E(6)*y+M%E(7)*x*y


end function jacobian_matrix


function g_ij_matrix(g1,g2,g3)  result(A)
  implicit none
  real, dimension(1:3),intent(in):: g1,g2,g3
  real, dimension(1:3,1:3)       :: A

  A(1,1)=DOT_PRODUCT(g1,g1)
  A(1,2)=DOT_PRODUCT(g1,g2)
  A(1,3)=DOT_PRODUCT(g1,g3)
  A(2,1)=DOT_PRODUCT(g2,g1)
  A(2,2)=DOT_PRODUCT(g2,g2)
  A(2,3)=DOT_PRODUCT(g2,g3)
  A(3,1)=DOT_PRODUCT(g3,g1)
  A(3,2)=DOT_PRODUCT(g3,g2)
  A(3,3)=DOT_PRODUCT(g3,g3)
end function g_ij_matrix


function inverse_jocobian_matrix(M,J)  result(J2)
  implicit none
  real, dimension(1:3,1:3),intent(in):: J
  type(mesh), intent(in)              :: M
  real, dimension(1:3,1:3):: J2

  J2=inv3x3(J)
end function inverse_jocobian_matrix

function gij_matrix(g1,g2,g3)  result(A)
  implicit none
  real, dimension(1:3),intent(in):: g1,g2,g3
  real, dimension(1:3,1:3)       :: A
  A(1,1)=DOT_PRODUCT(g1,g1)
  A(1,2)=DOT_PRODUCT(g1,g2)
  A(1,3)=DOT_PRODUCT(g1,g3)
  A(2,1)=DOT_PRODUCT(g2,g1)
  A(2,2)=DOT_PRODUCT(g2,g2)
  A(2,3)=DOT_PRODUCT(g2,g3)
  A(3,1)=DOT_PRODUCT(g3,g1)
  A(3,2)=DOT_PRODUCT(g3,g2)
  A(3,3)=DOT_PRODUCT(g3,g3)

end function gij_matrix


function find_gamma_tensor(M,J,g1,g2,g3,x,y,z)  result(A)
  implicit none
  real, dimension(1:3,1:3),intent(in):: J
  real, dimension(1:3),intent(in)    :: g1,g2,g3
  type(mesh), intent(in)                         :: M
  real, dimension(1:3,1:3,1:3)       :: A
  real,dimension(1:3)                :: g1_2,g1_3,g2_3
  real                               :: x,y,z
  
  g1_2(1)= M%C(4)+M%C(7)*z
  g1_2(2)= M%D(4)+M%D(7)*z
  g1_2(3)= M%E(4)+M%E(7)*z

  g1_3(1)= M%C(5)+M%C(7)*y
  g1_3(2)= M%D(5)+M%D(7)*y
  g1_3(3)= M%E(5)+M%E(7)*y

  g2_3(1)= M%C(6)+M%C(7)*x
  g2_3(2)= M%D(6)+M%D(7)*x
  g2_3(3)= M%E(6)+M%E(7)*x

! =======
  A(1,1,1) = 0.
  A(1,1,2) = DOT_PRODUCT(g1,g1_2)
  A(1,1,3) = DOT_PRODUCT(g1,g1_3)
  A(1,2,1) = A(1,1,2)
  A(1,2,2) = 0.
  A(1,2,3) = DOT_PRODUCT(g1,g2_3)
  A(1,3,1) = A(1,1,3)
  A(1,3,2) = A(1,2,3)
  A(1,3,3) = 0.
! ========
  A(2,1,1) = 0.
  A(2,1,2) = DOT_PRODUCT(g2,g1_2)
  A(2,1,3) = DOT_PRODUCT(g2,g1_3)
  A(2,2,1) = A(2,1,2)
  A(2,2,2) = 0.
  A(2,2,3) = DOT_PRODUCT(g2,g2_3)
  A(2,3,1) = A(2,1,3)
  A(2,3,2) = A(2,2,3)
  A(2,3,3) = 0.
! ========
  A(3,1,1) = 0.
  A(3,1,2) = DOT_PRODUCT(g3,g1_2)
  A(3,1,3) = DOT_PRODUCT(g3,g1_3)
  A(3,2,1) = A(3,1,2)
  A(3,2,2) = 0.
  A(3,2,3) = DOT_PRODUCT(g3,g2_3)
  A(3,3,1) = A(3,1,3)
  A(3,3,2) = A(3,2,3)
  A(3,3,3) = 0.
end function find_gamma_tensor

end module class_curvilinear
