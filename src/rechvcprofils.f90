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


subroutine rechvcprofils(recphvcok,cm,cr,cc,phi,th,u,tau,irecphvc,t,dtrecphvc)

use mod_cst_fault_rns
use netcdf
implicit none
type(pcomput) :: cc
type(pmeca) :: cm
type(prec) :: cr
integer :: i,recphvcok,irecphvc,k,l
integer :: ncid,x_dimid,y_dimid,t_dimid
integer, dimension(6) :: varid
real*8 :: t,dtrecphvc
real*8, dimension(cc%n) :: phi,th,u,tau
 character(len=10) :: numf

     irecphvc=irecphvc+1

     if (rang .eq. 0) then


            select case (irecphvc)
               case (1:9)
                 write(numf,'(I1)') irecphvc
               case (10:99)
                 write(numf,'(I2)') irecphvc
               case (100:999)
                 write(numf,'(I3)') irecphvc
               case (1000:9999)
                 write(numf,'(I4)') irecphvc
               case (10000:99999)
                 write(numf,'(I5)') irecphvc
            end select
            !----------------------!
            !--profils along x-----!
            !----------------------!
            !open(801,file=trim(cr%path_res)//'vx'//trim(numf)//'.data',form='unformatted',access='direct',recl=8)
            !open(802,file=trim(cr%path_res)//'thx'//trim(numf)//'.data',form='unformatted',access='direct',recl=8)
            !open(803,file=trim(cr%path_res)//'ux'//trim(numf)//'.data',form='unformatted',access='direct',recl=8)
            !open(804,file=trim(cr%path_res)//'taux'//trim(numf)//'.data',form='unformatted',access='direct',recl=8)
            !do i=1,cc%nx
            !   k=(i-1)*cc%ny+cc%ny/2-1
            !   l=k+1
            !   write(801,rec=i) real(exp(0.5*(phi(k)+phi(l)))*cm%vb0,8)
            !   write(802,rec=i) real(exp(0.5*(th(k)+th(l)))*cm%dc0/cm%vb0,8)
            !   write(803,rec=i) real(0.5*(u(k)+u(l))*cm%dc0,8)
            !   write(804,rec=i) real(0.5*(tau(k)+tau(l))*cm%b0*cm%sig0,8)
            !end do
            !close(801)
            !close(802)
            !close(803)
            !close(804)

            statusnetcdf=nf90_create(trim(cr%path_res)//'profilsx'//trim(numf)//'.nc', NF90_CLASSIC_MODEL, ncid)
            
            statusnetcdf=nf90_def_dim(ncid, "along strike",cc%nx,x_dimid)
            statusnetcdf=nf90_def_dim(ncid, "along depth",1,y_dimid)
            statusnetcdf=nf90_def_dim(ncid, "time",1,t_dimid)

            statusnetcdf=nf90_def_var(ncid, "slip rate", NF90_DOUBLE,[x_dimid,y_dimid],varid(1))
            statusnetcdf=nf90_def_var(ncid, "state", NF90_DOUBLE,[x_dimid,y_dimid],varid(2))
            statusnetcdf=nf90_def_var(ncid, "slip", NF90_DOUBLE,[x_dimid,y_dimid],varid(3))
            statusnetcdf=nf90_def_var(ncid, "shear stress", NF90_DOUBLE,[x_dimid,y_dimid],varid(4))
            statusnetcdf=nf90_def_var(ncid, "time", NF90_DOUBLE,[t_dimid],varid(5))
            statusnetcdf=nf90_def_var(ncid, "time delai", NF90_DOUBLE,[t_dimid],varid(6))

            statusnetcdf=nf90_enddef(ncid)

            do i=1,cc%nx
                k=(i-1)*cc%ny+cc%ny/2-1
                l=k+1

                statusnetcdf=nf90_put_var(ncid, varid(1), exp(0.5*(phi(k)+phi(l)))*cm%vb0, start=(/i/))
                statusnetcdf=nf90_put_var(ncid, varid(2), exp(0.5*(th(k)+th(l)))*cm%dc0/cm%vb0, start=(/i/))
                statusnetcdf=nf90_put_var(ncid, varid(3), 0.5*(u(k)+u(l))*cm%dc0, start=(/i/))
                statusnetcdf=nf90_put_var(ncid, varid(4), 0.5*(tau(k)+tau(l))*cm%b0*cm%sig0, start=(/i/))

            end do
            statusnetcdf=nf90_put_var(ncid, varid(5), t*cm%dc0/cm%vb0)
            statusnetcdf=nf90_put_var(ncid, varid(6), dtrecphvc*cm%dc0/cm%vb0)

            statusnetcdf=nf90_close(ncid)

            !----------------------!
            !--profils along y-----!
            !----------------------!
            !open(801,file=trim(cr%path_res)//'vy'//trim(numf)//'.data',form='unformatted',access='direct',recl=8)
            !open(802,file=trim(cr%path_res)//'thy'//trim(numf)//'.data',form='unformatted',access='direct',recl=8)
            !open(803,file=trim(cr%path_res)//'uy'//trim(numf)//'.data',form='unformatted',access='direct',recl=8)
            !open(804,file=trim(cr%path_res)//'tauy'//trim(numf)//'.data',form='unformatted',access='direct',recl=8)
            !do i=1,cc%ny
            !   k=(cc%nx/2-2)*cc%ny+i
            !   l=(cc%nx/2-1)*cc%ny+i
            !   write(801,rec=i) real(exp(0.5*(phi(k)+phi(l)))*cm%vb0,8)
            !   write(802,rec=i) real(exp(0.5*(th(k)+th(l)))*cm%dc0/cm%vb0,8)
            !   write(803,rec=i) real(0.5*(u(k)+u(l))*cm%dc0,8)
            !   write(804,rec=i) real(0.5*(tau(k)+tau(l))*cm%b0*cm%sig0,8)
            !end do
            !close(801)
            !close(802)
            !close(803)
            !close(804)


            statusnetcdf=nf90_create(trim(cr%path_res)//'profilsy'//trim(numf)//'.nc', NF90_CLASSIC_MODEL, ncid)
            
            statusnetcdf=nf90_def_dim(ncid, "along strike",cc%ny,x_dimid)
            statusnetcdf=nf90_def_dim(ncid, "along depth",1,y_dimid)
            statusnetcdf=nf90_def_dim(ncid, "time",1,t_dimid)

            statusnetcdf=nf90_def_var(ncid, "slip rate", NF90_DOUBLE,[x_dimid,y_dimid],varid(1))
            statusnetcdf=nf90_def_var(ncid, "state", NF90_DOUBLE,[x_dimid,y_dimid],varid(2))
            statusnetcdf=nf90_def_var(ncid, "slip", NF90_DOUBLE,[x_dimid,y_dimid],varid(3))
            statusnetcdf=nf90_def_var(ncid, "shear stress", NF90_DOUBLE,[x_dimid,y_dimid],varid(4))
            statusnetcdf=nf90_def_var(ncid, "time", NF90_DOUBLE,[t_dimid],varid(5))
            statusnetcdf=nf90_def_var(ncid, "time delai", NF90_DOUBLE,[t_dimid],varid(6))

            statusnetcdf=nf90_enddef(ncid)

            do i=1,cc%ny
               k=(cc%nx/2-2)*cc%ny+i
               l=(cc%nx/2-1)*cc%ny+i

                statusnetcdf=nf90_put_var(ncid, varid(1), exp(0.5*(phi(k)+phi(l)))*cm%vb0, start=(/i/))
                statusnetcdf=nf90_put_var(ncid, varid(2), exp(0.5*(th(k)+th(l)))*cm%dc0/cm%vb0, start=(/i/))
                statusnetcdf=nf90_put_var(ncid, varid(3), 0.5*(u(k)+u(l))*cm%dc0, start=(/i/))
                statusnetcdf=nf90_put_var(ncid, varid(4), 0.5*(tau(k)+tau(l))*cm%b0*cm%sig0, start=(/i/))

            end do
            statusnetcdf=nf90_put_var(ncid, varid(5), t*cm%dc0/cm%vb0)
            statusnetcdf=nf90_put_var(ncid, varid(6), dtrecphvc*cm%dc0/cm%vb0)

            statusnetcdf=nf90_close(ncid)

      end if
      
      recphvcok=0
      dtrecphvc=0.0


end subroutine rechvcprofils
