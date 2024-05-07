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


subroutine lin_static_freesurface(sk_proc,phi_proc,Kstatic_proc,cm,cc)

use mod_cst_fault_rns
implicit none
type(pmeca) :: cm
type(pcomput) :: cc
integer :: i,j
real*8, dimension(cc%nvalpr) :: phi_proc, sk_proc
real*8, dimension(cc%n) :: phi
real*8, dimension(cc%nvalpr,cc%n) :: Kstatic_proc

 call MPI_BARRIER(MPI_COMM_WORLD,code)
call MPI_ALLGATHER(phi_proc,cc%nvalpr,MPI_DOUBLE_PRECISION,phi,cc%nvalpr,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,code)

  !call MPI_GATHER(phi_proc,cc%nvalpr,MPI_DOUBLE_PRECISION,phi,cc%nvalpr,MPI_DOUBLE_PRECISION,0,&
 !                  MPI_COMM_WORLD,code)
                   
 do i=1,cc%nvalpr
     sk_proc(i)=0.0
 	   do j=1,cc%n
 		    sk_proc(i) =sk_proc(i)+ Kstatic_proc(i,j)*(exp(phi(j))-cm%vb)
 	   end do
 end do
call MPI_BARRIER(MPI_COMM_WORLD,code)
 		


end subroutine lin_static_freesurface
