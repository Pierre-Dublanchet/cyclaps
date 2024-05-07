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


subroutine elastic_kernel_2d_infperiodic(kappa_proc,cc,cm)

use mod_cst_fault_rns
implicit none
type(pmeca) :: cm
type(pcomput) :: cc
integer :: i
integer*8 :: k
real*8, dimension(cc%nvalpr) :: kappa_proc

   do i=1,cc%nvalpr
      if (i+rang*cc%nvalpr .le. cc%n/2) then
         k=(i+rang*cc%nvalpr)
      else
         k=cc%n+2-i-rang*cc%nvalpr
      end if

    kappa_proc(i)=abs(k-1)*pi/(cc%nx*cc%dx)


      if (k .eq. 1) then
         kappa_proc(i)=0.0
      else  
         kappa_proc(i)=-(kappa_proc(i)*(1.0+exp(-4.0*kappa_proc(i)*cm%H))/(1.0-exp(-4.0*kappa_proc(i)*cm%H)))
      end if
      
      
   end do
   
   if (cm%slipmode .eq. 2) then
   	kappa_proc(:)=kappa_proc(:)/(1-cm%nu)
   end if

end subroutine elastic_kernel_2d_infperiodic
