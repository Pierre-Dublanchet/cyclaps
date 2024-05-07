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


subroutine plan_static_stress_periodic(sk_proc,phi_proc,kappa_proc,cc)

use mod_cst_fault_rns
implicit none
type(pcomput) :: cc
integer :: i,j
real*8, dimension(cc%nvalpr) :: phi_proc,sk_proc
real*8, dimension(cc%ny,cc%nvalpr/cc%ny) :: v_proc, kappa_proc
complex*16, dimension(cc%ny,cc%nvalpr/cc%ny) :: iv_proc,vc_proc

do i=1,cc%ny
   do j=1,cc%nvalpr/cc%ny
      v_proc(i,j)=exp(phi_proc(i+(j-1)*cc%ny))
   end do
end do

call pzfft2d(dcmplx(v_proc),iv_proc,cc%ny,cc%nx,MPI_COMM_WORLD,nprocs,0)
call pzfft2d(dcmplx(v_proc),iv_proc,cc%ny,cc%nx,MPI_COMM_WORLD,nprocs,-1)

do i=1,cc%nvalpr/cc%ny
	do j=1,cc%ny
		iv_proc(j,i)=iv_proc(j,i)*kappa_proc(j,i)
	end do
end do

call pzfft2d(iv_proc,vc_proc,cc%ny,cc%nx,MPI_COMM_WORLD,nprocs,1)

do i=1,cc%ny
   do j=1,cc%nvalpr/cc%ny
      sk_proc(i+(j-1)*cc%ny)=real(vc_proc(i,j))
   end do
end do

end subroutine plan_static_stress_periodic
