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


subroutine x2ci_2d_rg(i,cc,x)

use mod_cst_fault_rns
implicit none
type(pcomput) :: cc
integer :: i
real*8 :: x

i=floor(x/cc%dx+cc%nx/2+1/2)+1




end subroutine x2ci_2d_rg
