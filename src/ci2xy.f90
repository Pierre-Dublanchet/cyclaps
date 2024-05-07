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


subroutine ci2xy(i,cc,x,y)

use mod_cst_fault_rns
implicit none
type(pcomput) :: cc
integer :: i
real*8 :: x,y

x=-0.5*cc%nx*cc%dx+floor(real((i+rang*cc%nvalpr-1)/cc%ny))*cc%dx+0.5*cc%dx
y=-0.5*cc%ny*cc%dy+(i+rang*cc%nvalpr-1-floor(real((i+rang*cc%nvalpr-1)/cc%ny))*cc%ny)*cc%dy+0.5*cc%dy

end subroutine ci2xy
