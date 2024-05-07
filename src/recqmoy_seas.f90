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


subroutine recqmoy_seas(recok,cm,cc,cr,phi,th,u,tau,irec,t,dtrec,vmax,vboxm,vwbox)

use mod_cst_fault_rns
use netcdf
implicit none
type(pcomput) :: cc
type(pmeca) :: cm
type(prec) :: cr
integer :: recok,irec
integer, dimension(cc%n) :: vboxm,vwbox
real*8 :: t,dtrec,vmax
real*8, dimension(cc%n) :: phi,th,u,tau

     irec=irec+1

     if (rang .eq. 0) then
            !write(10,rec=irec) real(t*cm%dc0/cm%vb0,8)
            !write(9,rec=irec) real(dtrec*cm%dc0/cm%vb0,8)
            !write(30,rec=irec) real(vmax*cm%vb0,8)
            !write(40,rec=irec) real(sum(u)*cm%dc0/cc%n,8)
            !write(50,rec=irec) real(sum(exp(th))*cm%dc0/(cm%vb0*cc%n),8)
            !write(60,rec=irec) real(sum(exp(phi)*vboxm(:))*cm%vb0/cc%n,8)
            !write(61,rec=irec) real(sum(exp(phi)*vwbox(:))*cm%vb0/cc%n,8)
            !write(70,rec=irec) real(sum(tau)*cm%b0*cm%sig0/cc%n,8)
            !call flush(10)
            !call flush(9)
            !call flush(20)
            !call flush(30)
            !call flush(40)
            !call flush(50)
            !call flush(60)
            !call flush(61)
            !call flush(70)


            statusnetcdf=nf90_put_var(cr%ncqmoyid, cr%qmoy_varid(1), t*cm%dc0/cm%vb0, start=(/irec/))
            statusnetcdf=nf90_put_var(cr%ncqmoyid, cr%qmoy_varid(2), dtrec*cm%dc0/cm%vb0, start=(/irec/))
            statusnetcdf=nf90_put_var(cr%ncqmoyid, cr%qmoy_varid(3), sum(exp(phi)*vboxm(:))*cm%vb0/cc%n, start=(/irec/))
            statusnetcdf=nf90_put_var(cr%ncqmoyid, cr%qmoy_varid(4), vmax*cm%vb0, start=(/irec/))
            statusnetcdf=nf90_put_var(cr%ncqmoyid, cr%qmoy_varid(5), sum(u)*cm%dc0/cc%n, start=(/irec/))
            statusnetcdf=nf90_put_var(cr%ncqmoyid, cr%qmoy_varid(6), sum(tau)*cm%b0*cm%sig0/cc%n, start=(/irec/))
            statusnetcdf=nf90_put_var(cr%ncqmoyid, cr%qmoy_varid(7), sum(exp(th))*cm%dc0/(cm%vb0*cc%n), start=(/irec/))
            statusnetcdf=nf90_put_var(cr%ncqmoyid, cr%qmoy_varid(8), sum(exp(phi)*vwbox(:))*cm%vb0/cc%n, start=(/irec/))

     end if
     
     recok=0
     dtrec=0.0



end subroutine recqmoy_seas
