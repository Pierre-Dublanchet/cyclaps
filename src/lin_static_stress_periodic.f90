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


subroutine lin_static_stress_periodic(sk_proc,phi_proc,kappa_proc,cc)

use mod_cst_fault_rns
implicit none
type(pcomput) :: cc
integer :: i
real*8, dimension(cc%nvalpr) :: phi_proc,sk_proc,kappa_proc
real*8, dimension(cc%nvalpr) :: v_proc
complex*16, dimension(cc%nvalpr) :: iv_proc,vc_proc,w_proc

   do i=1,cc%nvalpr
      v_proc(i)=exp(phi_proc(i))
   end do

   call pzfft1d(dcmplx(v_proc),iv_proc,w_proc,cc%n1,MPI_COMM_WORLD,rang,nprocs,0)
   call pzfft1d(dcmplx(v_proc),iv_proc,w_proc,cc%n1,MPI_COMM_WORLD,rang,nprocs,-1)


   do i=1,cc%nvalpr
      
      iv_proc(i)=iv_proc(i)*kappa_proc(i)
      
   end do
   
   call pzfft1d(iv_proc,vc_proc,w_proc,cc%n1,MPI_COMM_WORLD,rang,nprocs,1)

   sk_proc=real(vc_proc)

end subroutine lin_static_stress_periodic
