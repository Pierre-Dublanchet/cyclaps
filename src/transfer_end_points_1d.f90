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


subroutine transfer_end_points_1d(diffu_proc,p_proc,p0,pe,d0,de,cc)

use mod_cst_fault_rns
implicit none
type(pcomput) :: cc
integer :: i
integer, parameter :: etiquettedp=11, etiquettegp=12,etiquettedd=13, etiquettegd=14
real*8 :: p0,pe,d0,de
real*8, dimension(cc%nvalpr) :: p_proc,diffu_proc

if (nprocs .eq. 1) then

   p0=p_proc(2)
   pe=p_proc(cc%nvalpr-1)
   d0=diffu_proc(2)
   de=diffu_proc(cc%nvalpr-1)

else

 if (rang .eq. 0) then
 	call MPI_SEND(p_proc(cc%nvalpr),1,MPI_DOUBLE_PRECISION,rang+1,etiquettedp,MPI_COMM_WORLD,code)
 	call MPI_RECV(pe,1,MPI_DOUBLE_PRECISION,rang+1,etiquettegp,MPI_COMM_WORLD,statut,code)
 	p0=p_proc(2)
 	
 	call MPI_SEND(diffu_proc(cc%nvalpr),1,MPI_DOUBLE_PRECISION,rang+1,etiquettedd,MPI_COMM_WORLD,code)
 	call MPI_RECV(de,1,MPI_DOUBLE_PRECISION,rang+1,etiquettegd,MPI_COMM_WORLD,statut,code)
 	d0=diffu_proc(2)
 elseif (rang .eq. nprocs-1) then
 	call MPI_SEND(p_proc(1),1,MPI_DOUBLE_PRECISION,rang-1,etiquettegp,MPI_COMM_WORLD,code)
	call MPI_RECV(p0,1,MPI_DOUBLE_PRECISION,rang-1,etiquettedp,MPI_COMM_WORLD,statut,code)
	pe=p_proc(cc%nvalpr-1)
	
	call MPI_SEND(diffu_proc(1),1,MPI_DOUBLE_PRECISION,rang-1,etiquettegd,MPI_COMM_WORLD,code)
	call MPI_RECV(d0,1,MPI_DOUBLE_PRECISION,rang-1,etiquettedd,MPI_COMM_WORLD,statut,code)
	de=diffu_proc(cc%nvalpr-1)
else
 
	call MPI_SEND(p_proc(cc%nvalpr),1,MPI_DOUBLE_PRECISION,rang+1,etiquettedp,MPI_COMM_WORLD,code)
	call MPI_SEND(p_proc(1),1,MPI_DOUBLE_PRECISION,rang-1,etiquettegp,MPI_COMM_WORLD,code)
	call MPI_RECV(p0,1,MPI_DOUBLE_PRECISION,rang-1,etiquettedp,MPI_COMM_WORLD,statut,code)
	call MPI_RECV(pe,1,MPI_DOUBLE_PRECISION,rang+1,etiquettegp,MPI_COMM_WORLD,statut,code)

	call MPI_SEND(diffu_proc(cc%nvalpr),1,MPI_DOUBLE_PRECISION,rang+1,etiquettedd,MPI_COMM_WORLD,code)
	call MPI_SEND(diffu_proc(1),1,MPI_DOUBLE_PRECISION,rang-1,etiquettegd,MPI_COMM_WORLD,code)
	call MPI_RECV(d0,1,MPI_DOUBLE_PRECISION,rang-1,etiquettedd,MPI_COMM_WORLD,statut,code)
	call MPI_RECV(de,1,MPI_DOUBLE_PRECISION,rang+1,etiquettegd,MPI_COMM_WORLD,statut,code)
	
end if

end if

call MPI_BARRIER(MPI_COMM_WORLD,code)

end subroutine transfer_end_points_1d
