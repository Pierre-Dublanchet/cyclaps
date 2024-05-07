!---------------------------------------------------------!
!                                                         !
!                       CYCLAPS                           !
!                                                         !
! Copyright 2024 (c) Pierre Dublanchet                    !
! Author: Pierre Dublanchet (Mines Paris PSL/Armines)     !
! Contact: pierre.dublanchet@minesparis.psl.eu            !
! License: CeCILL-B                                       !
!                                                         !   
!---------------------------------------------------------!


subroutine normalize_rns_p(a_proc,b_proc,dc_proc,s_proc,u_proc,phi_proc,th_proc,p_proc,cm,chy,cc,cr)

use mod_cst_fault_rns
use mod_cst_fault_rns_p
implicit none
type(pmeca) :: cm
type(phydro) :: chy
type(pcomput) :: cc
type(prec) :: cr
integer :: i
real*8, dimension(cc%nvalpr) :: a_proc,b_proc,dc_proc,s_proc,u_proc,phi_proc,th_proc,p_proc

do i=1,cc%nvalpr
	a_proc(i)=a_proc(i)/cm%a0
	b_proc(i)=b_proc(i)/cm%b0
	dc_proc(i)=dc_proc(i)/cm%dc0
	s_proc(i)=s_proc(i)/cm%sig0
	p_proc(i)=p_proc(i)/cm%sig0
	u_proc(i)=u_proc(i)/cm%dc0
	phi_proc(i)=phi_proc(i)-log(cm%vb0)
	th_proc(i)=th_proc(i)+log(cm%vb0/cm%dc0)
end do
 cm%mu=0.5*cm%Y/(1+cm%nu)       !--shear modulus
 cm%r=cm%a0/cm%b0
 cm%beta=0.5*sqrt(cm%mu*cm%rho)*cm%vb0/(cm%b0*cm%sig0)   !--normalized damping
 cm%H= cm%H*cm%b0*cm%sig0/(cm%mu*cm%dc0)
 cm%mu0=cm%mu0/cm%b0
 cm%alpha=cm%alpha/cm%b0
 cm%vsis=cm%vsis/cm%vb0
 cm%vb=cm%vb/cm%vb0

 cm%paramdsext(1)= cm%paramdsext(1)*cm%dc0/(cm%b0*cm%sig0*cm%vb0)  !--normalized shear stressing rate
 cm%paramdsext(2)= cm%paramdsext(2)/(cm%b0*cm%sig0)  !--normalized amplitude of transient shear stress perturbation
 cm%paramdsext(3)= cm%paramdsext(3)*cm%b0*cm%sig0/(cm%mu*cm%dc0)  !--normalized radius of stressed patch
 cm%paramdsext(4)= cm%paramdsext(4)*cm%vb0/cm%dc0          !--normalized duration of transient shear stress perturbation
 cm%paramdsext(5)= cm%paramdsext(5)*cm%b0*cm%sig0/(cm%mu*cm%dc0)  !--normalized y coordinate of the center of stressed patch
 cm%paramdsext(6)= cm%paramdsext(6)*cm%b0*cm%sig0/(cm%mu*cm%dc0)  !--normalized y coordinate of the center of stressed patch

  !--Normalize hydraulic parameters here--!
 chy%D=chy%D*((cm%b0**2)*(cm%sig0**2))/(cm%vb0*(cm%mu**2)*cm%dc0)
 chy%paraminj(1)=chy%paraminj(1)*cm%vb0/cm%dc0                         !--duration of injection
 chy%paraminj(2)=chy%paraminj(2)*cm%b0*cm%sig0/(cm%mu*cm%dc0)    !--x coordinate of injection point
 chy%paraminj(3)=chy%paraminj(3)*cm%b0*cm%sig0/(cm%mu*cm%dc0)  !--y coordinate of injection point
 !chy%paraminj(4)=chy%paraminj(4)*cm%dc0/(cm%sig0*cm%vb0)           !--rate of mean pore pressure increase
 chy%paraminj(4)=chy%paraminj(4)*cm%b0/(cm%mu*cm%vb0*chy%cf*chy%phi)           !--Darcy Velocity at injection well
  chy%paraminj(5)=chy%paraminj(5)/cm%sig0                           !--pore pressure at injection point
    chy%paraminj(6)=chy%paraminj(6)*cm%b0*cm%sig0/(cm%mu*cm%dc0)     !--borehole radius
 !chy%paraminj(6)=chy%paraminj(6)*cm%b0*cm%sig0/(cm%mu*cm%dc0)     !--feedzone thickness
  !-------------------------------------------------!
 cc%dx=cc%dx*cm%b0*cm%sig0/(cm%mu*cm%dc0)
 cc%dy=cc%dy*cm%b0*cm%sig0/(cm%mu*cm%dc0)
 cc%paramdt(4)=cc%paramdt(4)*cm%vb0/cm%dc0
 cc%paramdt(5)=cc%paramdt(5)*cm%vb0/cm%dc0

 do i=1,4
 	cr%dtfrec(i)=cr%dtfrec(i)*cm%vb0/cm%dc0
 end do
 
 if (rang .eq. 0)  then
 	if (cc%dx .gt. 0.25*1.33) then
 		print*,'The model is discrete to some extent, for dx'
 	elseif (cc%dy .gt. 0.25*1.33) then
 	       print*,'The model is discrete to some extent, for dy'
	else
	      print*,'The model is continuous'
	end if 
 end if


 

end subroutine normalize_rns_p
