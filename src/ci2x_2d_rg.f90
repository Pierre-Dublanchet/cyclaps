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


subroutine ci2x_2d_rg(i,cm,cc,x)

use mod_cst_fault_rns
implicit none
type(pmeca) :: cm
type(pcomput) :: cc
integer :: i
real*8 :: x

x=-0.5*cc%nx*cc%dx+(i+rang*cc%nvalpr-1)*cc%dx+0.5*cc%dx

end subroutine ci2x_2d_rg
