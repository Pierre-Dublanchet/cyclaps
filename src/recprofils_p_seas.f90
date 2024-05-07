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


subroutine recprofils_p_seas(recpok,cm,chy,cr,cc,phi,th,u,tau,p,qdarcy,irecp,t,dtrecp)

use mod_cst_fault_rns
use mod_cst_fault_rns_p
use netcdf
implicit none
type(pcomput) :: cc
type(pmeca) :: cm
type(phydro) :: chy
type(prec) :: cr
integer :: i,recpok,irecp
integer :: ncid,x_dimid,y_dimid,t_dimid
integer, dimension(8) :: varid
real*8 :: t,dtrecp
real*8, dimension(cc%n) :: phi,th,u,tau,p,qdarcy
 character(len=10) :: numf

     irecp=irecp+1

     if (rang .eq. 0) then
            !write(80,rec=irecp) real(t*cm%dc0/cm%vb0,8)
            !write(8,rec=irecp) real(dtrecp*cm%dc0/cm%vb0,8)
            !call flush(80)
            !call flush(8)

            select case (irecp)
               case (1:9)
                 write(numf,'(I1)') irecp
               case (10:99)
                 write(numf,'(I2)') irecp
               case (100:999)
                 write(numf,'(I3)') irecp
               case (1000:9999)
                 write(numf,'(I4)') irecp
               case (10000:99999)
                 write(numf,'(I5)') irecp
            end select
            !open(11,file=trim(cr%path_res)//'v'//trim(numf)//'.data',form='unformatted',access='direct',recl=8)
            !open(12,file=trim(cr%path_res)//'th'//trim(numf)//'.data',form='unformatted',access='direct',recl=8)
            !open(13,file=trim(cr%path_res)//'u'//trim(numf)//'.data',form='unformatted',access='direct',recl=8)
            !open(14,file=trim(cr%path_res)//'tau'//trim(numf)//'.data',form='unformatted',access='direct',recl=8)
            !open(15,file=trim(cr%path_res)//'press'//trim(numf)//'.data',form='unformatted',access='direct',recl=8)
            !open(16,file=trim(cr%path_res)//'qdarcy'//trim(numf)//'.data',form='unformatted',access='direct',recl=8)
            !do i=1,cc%n
            !   write(11,rec=i) real(exp(phi(i))*cm%vb0,8)
            !   write(12,rec=i) real(exp(th(i))*cm%dc0/cm%vb0,8)
            !   write(13,rec=i) real(u(i)*cm%dc0,8)
            !   write(14,rec=i) real(tau(i)*cm%b0*cm%sig0,8)
            !   write(15,rec=i) real(p(i)*cm%sig0,8)
            !   write(16,rec=i) real(qdarcy(i)*chy%phi*chy%cf*cm%mu*cm%vb0/cm%b0,8)
            !end do
            !close(11)
            !close(12)
            !close(13)
            !close(14)
            !close(15)
            !close(16)

            statusnetcdf=nf90_create(trim(cr%path_res)//'maps'//trim(numf)//'.nc', NF90_CLASSIC_MODEL, ncid)
            
            statusnetcdf=nf90_def_dim(ncid, "along strike",cc%n,x_dimid)
            statusnetcdf=nf90_def_dim(ncid, "along depth",1,y_dimid)
            statusnetcdf=nf90_def_dim(ncid, "time",1,t_dimid)

            statusnetcdf=nf90_def_var(ncid, "slip rate", NF90_DOUBLE,[x_dimid,y_dimid],varid(1))
            statusnetcdf=nf90_def_var(ncid, "state", NF90_DOUBLE,[x_dimid,y_dimid],varid(2))
            statusnetcdf=nf90_def_var(ncid, "slip", NF90_DOUBLE,[x_dimid,y_dimid],varid(3))
            statusnetcdf=nf90_def_var(ncid, "shear stress", NF90_DOUBLE,[x_dimid,y_dimid],varid(4))
            statusnetcdf=nf90_def_var(ncid, "pore pressure", NF90_DOUBLE,[x_dimid,y_dimid],varid(5))
            statusnetcdf=nf90_def_var(ncid, "darcy velocity", NF90_DOUBLE,[x_dimid,y_dimid],varid(6))
            statusnetcdf=nf90_def_var(ncid, "time", NF90_DOUBLE,[t_dimid],varid(7))
            statusnetcdf=nf90_def_var(ncid, "time delai", NF90_DOUBLE,[t_dimid],varid(8))

            statusnetcdf=nf90_enddef(ncid)

            statusnetcdf=nf90_put_var(ncid, varid(1), exp(phi(:))*cm%vb0)
            statusnetcdf=nf90_put_var(ncid, varid(2), exp(th(:))*cm%dc0/cm%vb0)
            statusnetcdf=nf90_put_var(ncid, varid(3), u(:)*cm%dc0)
            statusnetcdf=nf90_put_var(ncid, varid(4), tau(:)*cm%b0*cm%sig0)
            statusnetcdf=nf90_put_var(ncid, varid(5), p(:)*cm%sig0)
            statusnetcdf=nf90_put_var(ncid, varid(6), qdarcy(:)*chy%phi*chy%cf*cm%mu*cm%vb0/cm%b0)
            statusnetcdf=nf90_put_var(ncid, varid(7), t*cm%dc0/cm%vb0)
            statusnetcdf=nf90_put_var(ncid, varid(8), dtrecp*cm%dc0/cm%vb0)

            statusnetcdf=nf90_close(ncid)
      end if
      
      recpok=0
      dtrecp=0.0

end subroutine recprofils_p_seas
