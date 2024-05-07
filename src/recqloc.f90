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


subroutine recqloc(reclocok,vrecpoints,cm,cc,cr,phi,th,u,tau,irecloc,t,dtrecloc)

use mod_cst_fault_rns
use netcdf
implicit none
type(pcomput) :: cc
type(pmeca) :: cm
type(prec) :: cr
integer :: i,reclocok,irecloc
integer, dimension(cr%nploc) :: vrecpoints
real*8 :: t,dtrecloc
real*8, dimension(cc%n) :: phi,th,u,tau

     irecloc=irecloc+1

     if (rang .eq. 0) then
          do i=1,cr%nploc
            !write(200+i,rec=irecloc) real(u(vrecpoints(i))*cm%dc0,8)
            !write(300+i,rec=irecloc) real(exp(th(vrecpoints(i)))*cm%dc0/(cm%vb0),8)
            !write(400+i,rec=irecloc) real(exp(phi(vrecpoints(i)))*cm%vb0,8)
            !write(500+i,rec=irecloc) real(tau(vrecpoints(i))*cm%b0*cm%sig0,8)
            !call flush(200+i)
            !call flush(300+i)
            !call flush(400+i)
            !call flush(500+i)

            statusnetcdf=nf90_put_var(cr%ncqlocid(i), cr%qloc_varid(i,1), t*cm%dc0/cm%vb0, start=(/irecloc/))
            statusnetcdf=nf90_put_var(cr%ncqlocid(i), cr%qloc_varid(i,2), dtrecloc*cm%dc0/cm%vb0, start=(/irecloc/))
            statusnetcdf=nf90_put_var(cr%ncqlocid(i), cr%qloc_varid(i,3), u(vrecpoints(i))*cm%dc0, start=(/irecloc/))
            statusnetcdf=nf90_put_var(cr%ncqlocid(i), cr%qloc_varid(i,4), exp(th(vrecpoints(i)))*cm%dc0/(cm%vb0), start=(/irecloc/))
            statusnetcdf=nf90_put_var(cr%ncqlocid(i), cr%qloc_varid(i,5), exp(phi(vrecpoints(i)))*cm%vb0, start=(/irecloc/))
            statusnetcdf=nf90_put_var(cr%ncqlocid(i), cr%qloc_varid(i,6), tau(vrecpoints(i))*cm%b0*cm%sig0, start=(/irecloc/))
          end do
     end if
     
     reclocok=0
     dtrecloc=0.0



end subroutine recqloc
