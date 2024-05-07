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


subroutine elastic_kernel_2d_freesurface(Kstatic_proc,cm,cc)

use mod_cst_fault_rns
implicit none
type(pmeca) :: cm
type(pcomput) :: cc
integer :: i,j
real*8, dimension(cc%nvalpr,cc%n) :: Kstatic_proc

do i=1,cc%nvalpr
   do j=1,cc%n
      Kstatic_proc(i,j)=(i+rang*cc%nvalpr+j-1.5)*(i+rang*cc%nvalpr+j-0.5)+(i+rang*cc%nvalpr-j-0.5)*(i+rang*cc%nvalpr-j+0.5)   
      Kstatic_proc(i,j)=Kstatic_proc(i,j)/(((i+rang*cc%nvalpr-0.5)**2-(j-1)**2)*((i+rang*cc%nvalpr-0.5)**2-j**2))
      Kstatic_proc(i,j)=Kstatic_proc(i,j)*0.5/(pi*cc%dx)
     ! Kstatic_proc(i,j)=0.0
       !print*,Kstatic_proc(i,j)
   end do
end do

end subroutine elastic_kernel_2d_freesurface
