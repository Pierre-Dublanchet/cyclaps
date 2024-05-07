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


subroutine tbp_rns_2d_cr(tbp_proc,phi_proc,cm,cc)

use mod_cst_fault_rns
implicit none
type(pmeca) :: cm
type(pcomput) :: cc
integer :: i
real*8 :: phi1,phi0
real*8, dimension(cc%nvalpr) :: tbp_proc,phi_proc


  if (rang .eq. 0) then
       phi0=phi_proc(1)
  end if 
  if (rang .eq. nprocs-1) then
       phi1=phi_proc(cc%nvalpr)
  end if 
    
 call MPI_BCAST(phi0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,code) 
 call MPI_BCAST(phi1,1,MPI_DOUBLE_PRECISION,nprocs-1,MPI_COMM_WORLD,code)



 call MPI_BARRIER(MPI_COMM_WORLD,code)
 
     do i=1,cc%nvalpr
       tbp_proc(i)=((exp(phi1)-cm%vb)/(i+rang*cc%nvalpr-cc%n-0.5))-((exp(phi0)-cm%vb)/(i+rang*cc%nvalpr-0.5))
       tbp_proc(i)=cm%paramdsext(1)+tbp_proc(i)/(2*pi*cc%dx*cm%sm_factor)
       !tbp_proc(i)=tbp_proc(i)/(2*pi*cc%dx*cm%sm_factor)
    end do


end subroutine tbp_rns_2d_cr
