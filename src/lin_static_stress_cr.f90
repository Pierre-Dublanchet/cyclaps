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


subroutine lin_static_stress_cr(sk_proc,phi_proc,kappa_proc,v0_proc,rkr,rks,ad,cm,cc)

use mod_cst_fault_rns
implicit none
type(pmeca) :: cm
type(pcomput) :: cc
integer :: i,k
integer, parameter :: etiquette=11
integer, dimension(nprocs) :: rkr,rks,ad
real*8, dimension(cc%nvalpr) :: phi_proc,sk_proc
real*8, dimension(2*cc%nvalpr) :: v0_proc,kappa_proc
complex*16, dimension(2*cc%nvalpr) :: iv_proc,vc_proc,w_proc


 !-----------------------------------------------!
 !--Padding of v------------------------------!
 !-----------------------------------------------!
 do i=1,2*cc%nvalpr
  	v0_proc(i)=0.0
  end do


     do k=1,nprocs

        if (rks(k) .eq. rkr(k)) then
            if (rang .eq. rks(k)) then
              v0_proc(ad(k):ad(k)+cc%nvalpr-1)=exp(phi_proc(1:cc%nvalpr))-cm%vb
            end if
        else
           if (rang .eq. rks(k)) then
              call MPI_SEND(exp(phi_proc(1:cc%nvalpr))-cm%vb,cc%nvalpr,MPI_DOUBLE_PRECISION,rkr(k),etiquette,MPI_COMM_WORLD,code)
           elseif (rang .eq. rkr(k)) then
              call MPI_RECV(v0_proc(ad(k)),cc%nvalpr,MPI_DOUBLE_PRECISION,rks(k),etiquette,MPI_COMM_WORLD,statut,code)
           end if
        end if

     end do
   
  call MPI_BARRIER(MPI_COMM_WORLD,code)
 !-------------------------------------------------------!
  !--Compute Fourier transform of slip rate-----!
 !-------------------------------------------------------!
 
   call pzfft1d(dcmplx(v0_proc),iv_proc,w_proc,cc%n1,MPI_COMM_WORLD,rang,nprocs,0)
   call pzfft1d(dcmplx(v0_proc),iv_proc,w_proc,cc%n1,MPI_COMM_WORLD,rang,nprocs,-1)

 !-------------------------------------------------------!
  !--Perform convolution with elastic kernel-----!
 !-------------------------------------------------------!
   do i=1,2*cc%nvalpr
      iv_proc(i)=kappa_proc(i)*iv_proc(i)  
   end do


 !-----------------------------------------------------------------------!
  !--Compute inverse Fourier transform of stressing rate-----!
 !-----------------------------------------------------------------------!
  call pzfft1d(iv_proc,vc_proc,w_proc,cc%n1,MPI_COMM_WORLD,rang,nprocs,1)
  v0_proc=real(vc_proc)

call MPI_BARRIER(MPI_COMM_WORLD,code)

!----------------------------------------------------------------!
!--Rearrange stressing rate--------------------------------!
!----------------------------------------------------------------!
  do k=1,nprocs

     if (rks(k) .eq. rkr(k)) then
        if (rang .eq. rks(k)) then
           sk_proc(1:cc%nvalpr)=v0_proc(ad(k):ad(k)+cc%nvalpr-1)
        end if
     else
        if (rang .eq. rkr(k)) then
           call MPI_SEND(v0_proc(ad(k):ad(k)+cc%nvalpr-1),cc%nvalpr,MPI_DOUBLE_PRECISION,rks(k),etiquette,MPI_COMM_WORLD,code)
        elseif (rang .eq. rks(k)) then
           call MPI_RECV(sk_proc(1),cc%nvalpr,MPI_DOUBLE_PRECISION,rkr(k),etiquette,MPI_COMM_WORLD,statut,code)
        end if
     end if
           
  end do
 
end subroutine lin_static_stress_cr
