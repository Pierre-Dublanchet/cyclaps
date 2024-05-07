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


subroutine elastic_kernel_2d_cr(kappa_proc,cc,cm)

use mod_cst_fault_rns
implicit none
type(pmeca) :: cm
type(pcomput) :: cc
integer :: i
integer*8 :: k
real*8 :: ci,si
real*8, dimension(2*cc%nvalpr) :: kappa_proc

  do i=1,2*cc%nvalpr
      if (i+rang*2*cc%nvalpr .le. cc%n) then
         k=(i+rang*2*cc%nvalpr)
     else
         k=2*cc%n+2-i-rang*2*cc%nvalpr
      end if
      call cisia(abs(k-1)*pi,ci,si)
    kappa_proc(i)=-abs(k-1)*si/(cc%dx*cc%nx)
     
   end do
   
   if (cm%slipmode .eq. 2) then
   	kappa_proc(:)=kappa_proc(:)/(1-cm%nu)
   end if

end subroutine elastic_kernel_2d_cr
