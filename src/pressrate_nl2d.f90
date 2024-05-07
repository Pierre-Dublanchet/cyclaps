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


subroutine pressrate_nl2d(prate_proc,p_proc,p0,pe,diffu_proc,d0,de,cm,chy,cc)

use mod_cst_fault_rns
use mod_cst_fault_rns_p
implicit none
type(pmeca) :: cm
type(phydro) :: chy
type(pcomput) :: cc
integer :: i,j,ic,ir,il,iu,id
real*8, dimension(cc%ny) :: p0,pe,d0,de
real*8, dimension(cc%nvalpr) :: prate_proc,p_proc,diffu_proc

!------------------------!
!--Diffusion term--------!
!------------------------!
!--2<=i<=ny-1--!
do i=2,cc%ny-1
   !--2<=j<=nx/np-1--!
   do j=2,cc%nx/nprocs-1
      ic=(j-1)*cc%ny+i
      iu=ic-1
      id=ic+1
      ir=ic+cc%ny
      il=ic-cc%ny
      prate_proc(ic)=((diffu_proc(ic)+diffu_proc(ir))*(p_proc(ir)-p_proc(ic))&
			-(diffu_proc(il)+diffu_proc(ic))*(p_proc(ic)-p_proc(il)))/(2*(cc%dx**2))
      prate_proc(ic)=prate_proc(ic)+((diffu_proc(ic)+diffu_proc(id))*(p_proc(id)-p_proc(ic))&
			-(diffu_proc(iu)+diffu_proc(ic))*(p_proc(ic)-p_proc(iu)))/(2*(cc%dy**2))
   end do
   !--j=1----!
   ic=i
   iu=ic-1
   id=ic+1
   ir=ic+cc%ny
   il=i
   prate_proc(ic)=((diffu_proc(ic)+diffu_proc(ir))*(p_proc(ir)-p_proc(ic))&
			-(d0(il)+diffu_proc(ic))*(p_proc(ic)-p0(il)))/(2*(cc%dx**2))
   prate_proc(ic)=prate_proc(ic)+((diffu_proc(ic)+diffu_proc(id))*(p_proc(id)-p_proc(ic))&
			-(diffu_proc(iu)+diffu_proc(ic))*(p_proc(ic)-p_proc(iu)))/(2*(cc%dy**2))
   !--j=nx/np----!
   ic=(cc%nx/nprocs-1)*cc%ny+i
   iu=ic-1
   id=ic+1
   ir=i
   il=ic-cc%ny
   prate_proc(ic)=((diffu_proc(ic)+de(ir))*(pe(ir)-p_proc(ic))&
			-(diffu_proc(il)+diffu_proc(ic))*(p_proc(ic)-p_proc(il)))/(2*(cc%dx**2))
   prate_proc(ic)=prate_proc(ic)+((diffu_proc(ic)+diffu_proc(id))*(p_proc(id)-p_proc(ic))&
			-(diffu_proc(iu)+diffu_proc(ic))*(p_proc(ic)-p_proc(iu)))/(2*(cc%dy**2))
end do

!--i=1--!
!--2<=j<=nx/np-1--!
do j=2,cc%nx/nprocs-1
   ic=(j-1)*cc%ny+1
   iu=ic+1
   id=ic+1
   ir=ic+cc%ny
   il=ic-cc%ny
   prate_proc(ic)=((diffu_proc(ic)+diffu_proc(ir))*(p_proc(ir)-p_proc(ic))&
			-(diffu_proc(il)+diffu_proc(ic))*(p_proc(ic)-p_proc(il)))/(2*(cc%dx**2))
   prate_proc(ic)=prate_proc(ic)+((diffu_proc(ic)+diffu_proc(id))*(p_proc(id)-p_proc(ic))&
			-(diffu_proc(iu)+diffu_proc(ic))*(p_proc(ic)-p_proc(iu)))/(2*(cc%dy**2))
end do
!--j=1--!
ic=1
iu=2
id=2
ir=1+cc%ny
il=1
prate_proc(ic)=((diffu_proc(ic)+diffu_proc(ir))*(p_proc(ir)-p_proc(ic))&
			-(d0(il)+diffu_proc(ic))*(p_proc(ic)-p0(il)))/(2*(cc%dx**2))
prate_proc(ic)=prate_proc(ic)+((diffu_proc(ic)+diffu_proc(id))*(p_proc(id)-p_proc(ic))&
			-(diffu_proc(iu)+diffu_proc(ic))*(p_proc(ic)-p_proc(iu)))/(2*(cc%dy**2))
   
!--j=nx/np----!
ic=cc%nvalpr-cc%ny+1
iu=ic+1
id=ic+1
ir=1
il=ic-cc%ny
prate_proc(ic)=((diffu_proc(ic)+de(ir))*(pe(ir)-p_proc(ic))&
		-(diffu_proc(il)+diffu_proc(ic))*(p_proc(ic)-p_proc(il)))/(2*(cc%dx**2))
prate_proc(ic)=prate_proc(ic)+((diffu_proc(ic)+diffu_proc(id))*(p_proc(id)-p_proc(ic))&
			-(diffu_proc(iu)+diffu_proc(ic))*(p_proc(ic)-p_proc(iu)))/(2*(cc%dy**2))

!--i=ny--!
!--2<=j<=nx/np-1--!
do j=2,cc%nx/nprocs-1
   ic=j*cc%ny
   iu=ic-1
   id=ic-1
   ir=ic+cc%ny
   il=ic-cc%ny
   prate_proc(ic)=((diffu_proc(ic)+diffu_proc(ir))*(p_proc(ir)-p_proc(ic))&
			-(diffu_proc(il)+diffu_proc(ic))*(p_proc(ic)-p_proc(il)))/(2*(cc%dx**2))
   prate_proc(ic)=prate_proc(ic)+((diffu_proc(ic)+diffu_proc(id))*(p_proc(id)-p_proc(ic))&
			-(diffu_proc(iu)+diffu_proc(ic))*(p_proc(ic)-p_proc(iu)))/(2*(cc%dy**2))
end do
!--j=1--!
ic=cc%ny
iu=ic-1
id=ic-1
ir=2*cc%ny
il=cc%ny
prate_proc(ic)=((diffu_proc(ic)+diffu_proc(ir))*(p_proc(ir)-p_proc(ic))&
			-(d0(il)+diffu_proc(ic))*(p_proc(ic)-p0(il)))/(2*(cc%dx**2))
prate_proc(ic)=prate_proc(ic)+((diffu_proc(ic)+diffu_proc(id))*(p_proc(id)-p_proc(ic))&
			-(diffu_proc(iu)+diffu_proc(ic))*(p_proc(ic)-p_proc(iu)))/(2*(cc%dy**2))
   
!--j=nx/np----!
ic=cc%nvalpr
iu=ic-1
id=ic-1
ir=cc%ny
il=ic-cc%ny
prate_proc(ic)=((diffu_proc(ic)+de(ir))*(pe(ir)-p_proc(ic))&
			-(diffu_proc(il)+diffu_proc(ic))*(p_proc(ic)-p_proc(il)))/(2*(cc%dx**2))
prate_proc(ic)=prate_proc(ic)+((diffu_proc(ic)+diffu_proc(id))*(p_proc(id)-p_proc(ic))&
			-(diffu_proc(iu)+diffu_proc(ic))*(p_proc(ic)-p_proc(iu)))/(2*(cc%dy**2))


!------------------------!
!--Source term-----------!
!------------------------!
if (rang .eq. chy%rang_inj) then
	prate_proc(chy%i_inj_proc)=prate_proc(chy%i_inj_proc)+chy%paraminj(4)*2*pi*chy%paraminj(6)/(cc%dx*cc%dy)
end if

end subroutine pressrate_nl2d
