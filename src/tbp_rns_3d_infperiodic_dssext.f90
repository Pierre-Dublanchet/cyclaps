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


subroutine tbp_rns_3d_infperiodic_dssext(tbp_proc,phi_proc,cm,cc,t)

use mod_cst_fault_rns
implicit none
type(pmeca) :: cm
type(pcomput) :: cc
integer :: i
real*8 :: vm,t,x,y,r
real*8, dimension(cc%nvalpr) :: tbp_proc,phi_proc
real*8, dimension(cc%n) :: phi

   call MPI_GATHER(phi_proc,cc%nvalpr,MPI_DOUBLE_PRECISION,phi,cc%nvalpr,MPI_DOUBLE_PRECISION,0,&
                   MPI_COMM_WORLD,code)

   if (rang .eq. 0) then
      vm=sum(exp(phi))/cc%n
   end if

   call MPI_BCAST(vm,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,code)

   do i=1,cc%nvalpr
      
      tbp_proc(i)=cm%paramdsext(1)-(0.5/cm%H)*(vm-cm%vb)
      if ((t .le. cm%paramdsext(4)) .and. (t .gt. 0.0)) then
         call ci2xy(i,cc,x,y)
         r=sqrt((x-cm%paramdsext(5))**2+(y-cm%paramdsext(6))**2)
         if (r .le. cm%paramdsext(3)) then
            tbp_proc(i)=tbp_proc(i)+cm%paramdsext(2)*exp(r**2/(r**2-cm%paramdsext(3)**2))&
                        *(2*(cm%paramdsext(4)**2)*(cm%paramdsext(4)-t))/((t**2)*((t-2*cm%paramdsext(4))**2))&
                           *exp((t-cm%paramdsext(4))**2/(t*(t-2*cm%paramdsext(4))))
         end if
      end if

   end do

end subroutine tbp_rns_3d_infperiodic_dssext
