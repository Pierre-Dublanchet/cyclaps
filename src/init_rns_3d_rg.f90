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


subroutine init_rns_3d_rg(a_proc,b_proc,dc_proc,s_proc,phi_proc,th_proc,u_proc,cm,cc,cr)

use mod_cst_fault_rns
use netcdf
implicit none
type(pmeca) :: cm
type(pcomput) :: cc
type(prec) :: cr
integer :: i,j
integer :: ncid,x_dimid,y_dimid
integer, dimension(7) :: varid
real*8 :: x,y,r
real*8, dimension(cc%nvalpr) :: a_proc,b_proc,dc_proc,s_proc
real*8, dimension(cc%nvalpr) :: phi_proc,th_proc,u_proc
real*8, dimension(:), allocatable :: a,b,dc,s,phi,th,u

!-----------------------------------------------!
!--Initial Conditions defined here--------------!
!-----------------------------------------------!
if (cc%methinit .eq. 0) then

do i=1,cc%nvalpr
    !-----------------------------------------------!
    !--Spatial coordinates x along the fault-(non-normalized)-------!
    !-----------------------------------------------!
    call ci2xy(i,cc,x,y)
    r=sqrt(x**2+y**2)   !--radial distance from the  fault center
    !-----------------------------------------------!
    !--Initial conditions on VS regions-------------!
    !-----------------------------------------------!
    a_proc(i)=cm%a_asp(1)
    b_proc(i)=cm%b_asp(1)
    s_proc(i)=cm%sig_asp(1)
    dc_proc(i)=cm%dc_asp(1)

    phi_proc(i)=log(cm%vi_asp(1))
    th_proc(i)=log(cm%thi_asp(1))
    u_proc(i)=cm%ui_asp(1)
    
    !-----------------------------------------------!
    !--Initial conditions on VW regions-------------!
    !-----------------------------------------------!
    do j=1,cm%n_asp
       if (r .le. cm%R_asp(j+1))  then
           a_proc(i)=cm%a_asp(j+1)
           b_proc(i)=cm%b_asp(j+1)
           s_proc(i)=cm%sig_asp(j+1)
           dc_proc(i)=cm%dc_asp(j+1)

           phi_proc(i)=log(cm%vi_asp(j+1))
           th_proc(i)=log(cm%thi_asp(j+1))
           u_proc(i)=cm%ui_asp(j+1)

      end if
    end do
    


 end do

 if (rang .eq. 0) then
   allocate(a(cc%n),b(cc%n),s(cc%n),dc(cc%n),phi(cc%n),th(cc%n),u(cc%n))
 end if

 call MPI_GATHER(a_proc,cc%nvalpr,MPI_DOUBLE_PRECISION,a,cc%nvalpr,MPI_DOUBLE_PRECISION,0,&
                  MPI_COMM_WORLD,code)
 call MPI_GATHER(b_proc,cc%nvalpr,MPI_DOUBLE_PRECISION,b,cc%nvalpr,MPI_DOUBLE_PRECISION,0,&
                  MPI_COMM_WORLD,code)
 call MPI_GATHER(s_proc,cc%nvalpr,MPI_DOUBLE_PRECISION,s,cc%nvalpr,MPI_DOUBLE_PRECISION,0,&
                  MPI_COMM_WORLD,code)
 call MPI_GATHER(dc_proc,cc%nvalpr,MPI_DOUBLE_PRECISION,dc,cc%nvalpr,MPI_DOUBLE_PRECISION,0,&
                  MPI_COMM_WORLD,code)
 call MPI_GATHER(phi_proc,cc%nvalpr,MPI_DOUBLE_PRECISION,phi,cc%nvalpr,MPI_DOUBLE_PRECISION,0,&
                  MPI_COMM_WORLD,code)
 call MPI_GATHER(th_proc,cc%nvalpr,MPI_DOUBLE_PRECISION,th,cc%nvalpr,MPI_DOUBLE_PRECISION,0,&
                  MPI_COMM_WORLD,code)
 call MPI_GATHER(u_proc,cc%nvalpr,MPI_DOUBLE_PRECISION,u,cc%nvalpr,MPI_DOUBLE_PRECISION,0,&
                  MPI_COMM_WORLD,code)

 !--Extract reference mechanical parameters---!
  if (rang .eq. 0) then
  
   cm%a0=maxval(a(:))
   cm%b0=maxval(b(:))
   cm%sig0=maxval(s(:))
   cm%dc0=minval(dc(:))
   
   end if
   
   call MPI_BCAST(cm%a0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,code)
   call MPI_BCAST(cm%b0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,code)
   call MPI_BCAST(cm%dc0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,code)
   call MPI_BCAST(cm%sig0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,code)
   
   if (rang .eq. 0) then

    call system('rm -f '//trim(cc%path_init)//'init.nc')
    statusnetcdf=nf90_create(trim(cr%path_res)//'init.nc', NF90_CLASSIC_MODEL, ncid)
    statusnetcdf=nf90_def_dim(ncid, "index",cc%n,x_dimid)
    statusnetcdf=nf90_def_var(ncid, "a", NF90_DOUBLE,[x_dimid],varid(1))
    statusnetcdf=nf90_def_var(ncid, "b", NF90_DOUBLE,[x_dimid],varid(2))
    statusnetcdf=nf90_def_var(ncid, "s", NF90_DOUBLE,[x_dimid],varid(3))
    statusnetcdf=nf90_def_var(ncid, "dc", NF90_DOUBLE,[x_dimid],varid(4))
    statusnetcdf=nf90_def_var(ncid, "vi", NF90_DOUBLE,[x_dimid],varid(5))
    statusnetcdf=nf90_def_var(ncid, "thi", NF90_DOUBLE,[x_dimid],varid(6))
    statusnetcdf=nf90_def_var(ncid, "ui", NF90_DOUBLE,[x_dimid],varid(7))
    statusnetcdf=nf90_enddef(ncid)
    statusnetcdf=nf90_put_var(ncid, varid(1), a(:))
    statusnetcdf=nf90_put_var(ncid, varid(2), b(:))
    statusnetcdf=nf90_put_var(ncid, varid(3), s(:))
    statusnetcdf=nf90_put_var(ncid, varid(4), dc(:))
    statusnetcdf=nf90_put_var(ncid, varid(5), exp(phi(:)))
    statusnetcdf=nf90_put_var(ncid, varid(6), exp(th(:)))
    statusnetcdf=nf90_put_var(ncid, varid(7), u(:))
    statusnetcdf=nf90_close(ncid)

    deallocate(a,b,s,dc,phi,th,u)
  end if


!-----------------------------------------------!
!--Initial Conditions read from file------------!
!-----------------------------------------------!
 else

 if (rang .eq. 0) then
    allocate(a(cc%n),b(cc%n),s(cc%n),dc(cc%n),phi(cc%n),th(cc%n),u(cc%n))
    statusnetcdf=nf90_open(trim(cc%path_init)//'init.nc',NF90_NOWRITE,ncid)
    statusnetcdf=nf90_inq_dimid(ncid,'x',x_dimid)
    statusnetcdf=nf90_inq_dimid(ncid,'y',y_dimid)
    statusnetcdf=nf90_inq_varid(ncid,'a',varid(1))
    statusnetcdf=nf90_inq_varid(ncid,'b',varid(2))
    statusnetcdf=nf90_inq_varid(ncid,'s',varid(3))
    statusnetcdf=nf90_inq_varid(ncid,'dc',varid(4))
    statusnetcdf=nf90_inq_varid(ncid,'vi',varid(5))
    statusnetcdf=nf90_inq_varid(ncid,'thi',varid(6))
    statusnetcdf=nf90_inq_varid(ncid,'ui',varid(7))
    statusnetcdf=nf90_get_var(ncid,varid(1),a(:))
    statusnetcdf=nf90_get_var(ncid,varid(2),b(:))
    statusnetcdf=nf90_get_var(ncid,varid(3),s(:))
    statusnetcdf=nf90_get_var(ncid,varid(4),dc(:))
    statusnetcdf=nf90_get_var(ncid,varid(5),phi(:))
    statusnetcdf=nf90_get_var(ncid,varid(6),th(:))
    statusnetcdf=nf90_get_var(ncid,varid(7),u(:))
    statusnetcdf=nf90_close(ncid)
 end if
 
  !--Extract reference mechanical parameters---!
   if (rang .eq. 0) then
  
   cm%a0=minval(a(:))
   cm%b0=minval(b(:))
   cm%sig0=maxval(s(:))
   cm%dc0=minval(dc(:))
   
   end if
   
    call MPI_BCAST(cm%a0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,code)
   call MPI_BCAST(cm%b0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,code)
   call MPI_BCAST(cm%dc0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,code)
   call MPI_BCAST(cm%sig0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,code)

    call MPI_SCATTER(a,cc%nvalpr,MPI_DOUBLE_PRECISION,a_proc,cc%nvalpr,MPI_DOUBLE_PRECISION,0,&
                     MPI_COMM_WORLD,code)
    call MPI_SCATTER(b,cc%nvalpr,MPI_DOUBLE_PRECISION,b_proc,cc%nvalpr,MPI_DOUBLE_PRECISION,0,&
                     MPI_COMM_WORLD,code)
    call MPI_SCATTER(s,cc%nvalpr,MPI_DOUBLE_PRECISION,s_proc,cc%nvalpr,MPI_DOUBLE_PRECISION,0,&
                     MPI_COMM_WORLD,code)
    call MPI_SCATTER(dc,cc%nvalpr,MPI_DOUBLE_PRECISION,dc_proc,cc%nvalpr,MPI_DOUBLE_PRECISION,0,&
                     MPI_COMM_WORLD,code)
    call MPI_SCATTER(phi,cc%nvalpr,MPI_DOUBLE_PRECISION,phi_proc,cc%nvalpr,MPI_DOUBLE_PRECISION,0,&
                     MPI_COMM_WORLD,code)
    call MPI_SCATTER(th,cc%nvalpr,MPI_DOUBLE_PRECISION,th_proc,cc%nvalpr,MPI_DOUBLE_PRECISION,0,&
                     MPI_COMM_WORLD,code)
    call MPI_SCATTER(u,cc%nvalpr,MPI_DOUBLE_PRECISION,u_proc,cc%nvalpr,MPI_DOUBLE_PRECISION,0,&
                     MPI_COMM_WORLD,code)

 if (rang .eq. 0) then
    deallocate(a,b,s,dc,phi,th,u)
 end if

 end if

end subroutine init_rns_3d_rg
