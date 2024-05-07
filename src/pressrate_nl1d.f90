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


subroutine pressrate_nl1d(prate_proc,p_proc,p0,pe,diffu_proc,d0,de,cm,chy,cc,t)

use mod_cst_fault_rns
use mod_cst_fault_rns_p
implicit none
type(pmeca) :: cm
type(phydro) :: chy
type(pcomput) :: cc
integer :: i
real*8 :: t
real*8 :: p0,pe,d0,de
real*8, dimension(cc%nvalpr) :: prate_proc,p_proc,diffu_proc

!--Diffusion term with zero flux at fault ends--------!
do i=2,cc%nvalpr-1 
	prate_proc(i)=((diffu_proc(i)+diffu_proc(i+1))*(p_proc(i+1)-p_proc(i))-(diffu_proc(i-1)+&
		diffu_proc(i))*(p_proc(i)-p_proc(i-1)))/(2*cc%dx**2)	
end do
prate_proc(1)=((diffu_proc(1)+diffu_proc(2))*(p_proc(2)-p_proc(1))-(d0+&
		diffu_proc(1))*(p_proc(1)-p0))/(2*cc%dx**2)
prate_proc(cc%nvalpr)=((diffu_proc(cc%nvalpr)+de)*(pe-p_proc(cc%nvalpr))-(diffu_proc(cc%nvalpr-1)+&
		diffu_proc(cc%nvalpr))*(p_proc(cc%nvalpr)-p_proc(cc%nvalpr-1)))/(2*cc%dx**2)

!--Source term-----------!
if ((rang .eq. chy%rang_inj) .and. (t .le. chy%paraminj(1))) then
    !--if chy%paraminj(4): mean pore pressure rate
	!prate_proc(chy%i_inj_proc)=prate_proc(chy%i_inj_proc)+cc%nx*cc%dx*chy%paraminj(4)
	!--for SEAS benchamrk BP6
	prate_proc(chy%i_inj_proc)=prate_proc(chy%i_inj_proc)+chy%paraminj(4)/(cc%dx)
end if

end subroutine pressrate_nl1d
