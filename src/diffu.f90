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


subroutine diffu(diffu_proc,phi_proc,th_proc,p_proc,s_proc,a_proc,b_proc,dc_proc,cm,chy,cc)

use mod_cst_fault_rns
use mod_cst_fault_rns_p
implicit none
type(pmeca) :: cm
type(phydro) :: chy
type(pcomput) :: cc
integer :: i
real*8, dimension(cc%nvalpr) :: diffu_proc,phi_proc,th_proc,p_proc,s_proc,a_proc,b_proc,dc_proc

do i=1,cc%nvalpr
	diffu_proc(i)=chy%D   !*exp(-th_proc(i)/100.0)
end do

end subroutine diffu
