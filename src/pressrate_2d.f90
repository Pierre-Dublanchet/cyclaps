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


subroutine pressrate_2d(prate_proc,t,cm,chy,cc)

use mod_cst_fault_rns
use mod_cst_fault_rns_p
implicit none
type(pmeca) :: cm
type(phydro) :: chy
type(pcomput) :: cc
integer :: i
real*8 :: t,x,y,r,eta
real*8, dimension(cc%nvalpr) :: prate_proc

do i=1,cc%nvalpr
	if (t .gt. 0.0) then
	          call ci2xy(i,cc,x,y)
     	         r=sqrt((x-chy%paraminj(2))**2+(y-chy%paraminj(3))**2)
     	         eta=r/(2*sqrt(chy%D*t))
	       !--Mauvaise solution pour constant flux at injection point
		prate_proc(i)=(chy%paraminj(5)*eta/(sqrt(pi)*t))*exp(-eta**2)
	else
		prate_proc(i)=0.0
	end if
	
end do

end subroutine pressrate_2d
