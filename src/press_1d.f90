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


subroutine press_1d(p_proc,t,cm,chy,cc)

use mod_cst_fault_rns
use mod_cst_fault_rns_p
implicit none
type(pmeca) :: cm
type(phydro) :: chy
type(pcomput) :: cc
integer :: i
real*8 :: t,x
real*8, dimension(cc%nvalpr) :: p_proc



do i=1,cc%nvalpr
       call ci2x_2d_rg(i,cm,cc,x)
       if (t .gt. 0.0) then
		p_proc(i)=chy%paraminj(5)*erfc(abs(x-chy%paraminj(2))/(2*sqrt(chy%D*t)))
			!if (i+rang*cc%nvalpr .eq. cc%nx/2) then
			!print*,p_proc(i),chy%paraminj(1),chy%paraminj(2),chy%D,t
			!endif
	else
	     	p_proc(i)=0.0
	endif
end do


	

end subroutine press_1d
