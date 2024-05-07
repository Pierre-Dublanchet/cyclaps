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


subroutine boxmonitor(cc,cr,vboxm_proc,vwbox_proc,a_proc,b_proc)

use mod_cst_fault_rns
implicit none
type(pcomput) :: cc
type(prec) :: cr
integer :: i
integer, dimension(cc%nvalpr) :: vboxm_proc,vwbox_proc
real*8, dimension(cc%nvalpr) :: a_proc,b_proc
real*8 :: x,y,r


	do i=1,cc%nvalpr
	   call ci2xy(i,cc,x,y)
	   r=sqrt(x**2+y**2)
	   !if ((abs(x) .le. 0.5*cc%nx*cc%dx-cr%hboxm) .and. (abs(y) .le. 0.5*cc%ny*cc%dy-cr%hboxm)) then
	   if (r .le. 1.5*cr%hboxm) then
	   	  vboxm_proc(i)=1
	   else
	   	  vboxm_proc(i)=0
	   end if

	   if (r .le. cr%hboxm) then
	   	  vwbox_proc(i)=1
	   else
	   	  vwbox_proc(i)=0
	   end if
	end do




end subroutine boxmonitor