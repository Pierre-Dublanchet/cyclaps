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


subroutine tbp_rns_2d_infperiodic(tbp_proc,phi_proc,cm,cc)

use mod_cst_fault_rns
implicit none
type(pmeca) :: cm
type(pcomput) :: cc
integer :: i
real*8 :: vm
real*8, dimension(cc%nvalpr) :: tbp_proc,phi_proc
real*8, dimension(cc%n) :: phi



   call MPI_GATHER(phi_proc,cc%nvalpr,MPI_DOUBLE_PRECISION,phi,cc%nvalpr,MPI_DOUBLE_PRECISION,0,&
                   MPI_COMM_WORLD,code)

   if (rang .eq. 0) then
      vm=sum(exp(phi))/cc%n
   end if

   call MPI_BCAST(vm,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,code)

   do i=1,cc%nvalpr
      tbp_proc(i)=cm%paramdsext(1)-(0.5/(cm%H*cm%sm_factor))*(vm-cm%vb)
      
   end do






end subroutine tbp_rns_2d_infperiodic
