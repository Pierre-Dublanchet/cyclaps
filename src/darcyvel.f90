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


subroutine darcyvel(p_proc,p0,pe,qdarcy_proc,cc,chy,cm)

use mod_cst_fault_rns
use mod_cst_fault_rns_p
implicit none
type(pcomput) :: cc
type(pmeca) :: cm
type(phydro) :: chy
integer :: i
real*8 :: p0,pe
real*8, dimension(cc%nvalpr) :: p_proc,qdarcy_proc

do i=2,cc%nvalpr-1
	qdarcy_proc(i)=-chy%D*(p_proc(i+1)-p_proc(i-1))/(2*cc%dx)
end do

qdarcy_proc(1)=-chy%D*(p_proc(2)-p0)/(2*cc%dx)
qdarcy_proc(cc%nvalpr)=-chy%D*(pe-p_proc(cc%nvalpr))/(2*cc%dx)


end subroutine darcyvel