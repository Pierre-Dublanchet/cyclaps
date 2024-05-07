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


subroutine pressrate_1d(prate_proc,t,cm,chy,cc)

use mod_cst_fault_rns
use mod_cst_fault_rns_p
implicit none
type(pmeca) :: cm
type(phydro) :: chy
type(pcomput) :: cc
integer :: i
real*8 :: t,x
real*8, dimension(cc%nvalpr) :: prate_proc

do i=1,cc%nvalpr
	call ci2x_2d_rg(i,cm,cc,x)
	if (t .gt. 0.0) then
		prate_proc(i)=chy%paraminj(5)*(abs(x-chy%paraminj(2))/(2*sqrt(pi*chy%D*t)*t))*exp(-abs((x-chy%paraminj(2))**2)/(4.0*chy%D*t))
	else
		prate_proc(i)=0.0
	end if
	
end do

end subroutine pressrate_1d
