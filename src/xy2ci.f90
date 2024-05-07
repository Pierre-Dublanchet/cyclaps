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


subroutine xy2ci(i,cc,x,y)

use mod_cst_fault_rns
implicit none
type(pcomput) :: cc
integer :: i,k,l
real*8 :: x,y


l=floor(x/cc%dx+cc%nx/2+1)
k=floor(y/cc%dy+cc%ny/2+1)

i=(l-1)*cc%ny+k


end subroutine xy2ci
