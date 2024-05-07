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


subroutine reccatal_3d_rg(cc,cm,cr,phi_proc,th_proc,tau_proc,u_proc,s_proc,a_proc,b_proc,irecevt,t,vmax,isevt,&
				 t0evt,dt0evt,dtevt,m0_proc,&
                    m_proc,xa_proc,ya_proc,u0_proc,tau0_proc,du_proc,dtau_proc,nelevt_proc,nelevt_proc0,&
                    x0,y0,xa,ya,nelevt,xi,yi,m,du,dtau,tev)

use mod_cst_fault_rns
use netcdf
implicit none
type(pcomput) :: cc
type(pmeca) :: cm
type(prec) :: cr
integer :: i,irecevt,isevt,nelevt_proc,nelevt_proc0
integer, dimension(nprocs) :: nelevt
real*8 :: t0evt,dtevt,dt0evt,tev,m0_proc,m_proc,du_proc,dtau_proc,taui_proc,tauf_proc,t,vmax
real*8 :: x0,y0,x,y,xi,yi
real*8, dimension(nprocs) :: m,du,dtau,taui,tauf,xa,ya
real*8, dimension(cc%nvalpr) :: phi_proc,th_proc,tau_proc,u_proc,s_proc,a_proc,b_proc,xa_proc,ya_proc,u0_proc,tau0_proc

     !----------------------------------------------------------------!
     !--Eqk initiation parameters-------------------------------------!
     !----------------------------------------------------------------!
     if ((vmax .gt. cm%vsis) .and. (isevt .eq. 0)) then
        isevt=1                                 !--earthquake on
        t0evt=t                                 !--initiation time
        dt0evt=dtevt                            !--
        dtevt=0.0                               !--eqk duration 
        !--Initial seismic moment per proc.---!
        m0_proc=sum(u_proc(:))*cc%dx*cc%dy


        call MPI_BARRIER(MPI_COMM_WORLD,code)

        !-Eqk initiation point+initial slip and stress distribution-------!
        do i=1,cc%nvalpr
           xa_proc(i)=0.0
           ya_proc(i)=0.0
           u0_proc(i)=-1.0
           tau0_proc(i)=0.0
        end do

        call MPI_BARRIER(MPI_COMM_WORLD,code)    

        nelevt_proc=0  !-number of elements participating in the eqk rupture
        do i=1,cc%nvalpr
           if (exp(phi_proc(i)) .gt. cm%vsis) then 
              nelevt_proc=nelevt_proc+1
              call ci2xy(i,cc,x,y)
              xa_proc(i)=x
              ya_proc(i)=y
              u0_proc(i)=u_proc(i)
              tau0_proc(i)=(cm%mu0+cm%r*a_proc(i)*phi_proc(i)+b_proc(i)*th_proc(i))*s_proc(i)
           end if
        end do

        nelevt_proc0=nelevt_proc
        x0=sum(xa_proc)
        y0=sum(ya_proc)
        call MPI_GATHER(x0,1,MPI_DOUBLE_PRECISION,xa,1,MPI_DOUBLE_PRECISION,0,&
                   MPI_COMM_WORLD,code)
        call MPI_GATHER(y0,1,MPI_DOUBLE_PRECISION,ya,1,MPI_DOUBLE_PRECISION,0,&
                   MPI_COMM_WORLD,code)
        call MPI_GATHER(nelevt_proc0,1,MPI_INTEGER,nelevt,1,MPI_INTEGER,0,&
                   MPI_COMM_WORLD,code)
        if (rang .eq. 0) then
           xi=sum(xa)/sum(nelevt)
           yi=sum(ya)/sum(nelevt)
           print*,' New event (number ',irecevt+1,')'
        end if
      end if

      !-------------------------------------------------------------------!
      !--Track eqk development--------------------------------------------!
      !-------------------------------------------------------------------!
      if (isevt .eq. 1) then
        do i=1,cc%nvalpr
           if ((exp(phi_proc(i)) .gt. cm%vsis) .and. (u0_proc(i) .lt. 0.0)) then 
              nelevt_proc=nelevt_proc+1
              call ci2xy(i,cc,x,y)
              xa_proc(i)=x
              ya_proc(i)=y
              u0_proc(i)=u_proc(i)
              tau0_proc(i)=(cm%mu0+cm%r*a_proc(i)*phi_proc(i)+b_proc(i)*th_proc(i))*s_proc(i)
           end if
        end do
      end if

      !-------------------------------------------------------------------!
      !--Eqk termination construct catalogue entry------------------------!
      !-------------------------------------------------------------------!
      if ((isevt .eq. 1) .and. (vmax .lt. cm%vsis)) then
         isevt=0                                          !--eqk off
         tev=dtevt                                        !--final eqk duration
         irecevt=irecevt+1                                !--line number in eqk catalogue
        !--Eqk seismic moment/proc----------------------------------------!
		m_proc=sum(u_proc(:))*cc%dx*cc%dy-m0_proc

        !--Eqk coseismic slip+static stress drop distributions-----------! 
        du_proc=0.0
        dtau_proc=0.0
        tauf_proc=0.0
        taui_proc=0.0
        do i=1,cc%nvalpr
           if (u0_proc(i) .gt. 0.0) then
              du_proc=du_proc+u_proc(i)-u0_proc(i)
              taui_proc=taui_proc+tau0_proc(i)
              tauf_proc=tauf_proc+(cm%mu0+cm%r*a_proc(i)*phi_proc(i)+b_proc(i)*th_proc(i))*s_proc(i)
              dtau_proc=dtau_proc+(cm%mu0+cm%r*a_proc(i)*phi_proc(i)+b_proc(i)*th_proc(i))*s_proc(i)-tau0_proc(i)
           end if
        end do
         !--gather final quantities-------------------------------------!
         x0=sum(xa_proc)
         y0=sum(ya_proc)
         nelevt_proc0=nelevt_proc
         call MPI_GATHER(x0,1,MPI_DOUBLE_PRECISION,xa,1,MPI_DOUBLE_PRECISION,0,&
                      MPI_COMM_WORLD,code)
         call MPI_GATHER(y0,1,MPI_DOUBLE_PRECISION,ya,1,MPI_DOUBLE_PRECISION,0,&
                      MPI_COMM_WORLD,code)
         call MPI_GATHER(nelevt_proc0,1,MPI_INTEGER,nelevt,1,MPI_INTEGER,0,&
                      MPI_COMM_WORLD,code)
         call MPI_GATHER(m_proc,1,MPI_DOUBLE_PRECISION,m,1,MPI_DOUBLE_PRECISION,0,&
                      MPI_COMM_WORLD,code)
         call MPI_GATHER(du_proc,1,MPI_DOUBLE_PRECISION,du,1,MPI_DOUBLE_PRECISION,0,&
                      MPI_COMM_WORLD,code)
         call MPI_GATHER(dtau_proc,1,MPI_DOUBLE_PRECISION,dtau,1,MPI_DOUBLE_PRECISION,0,&
                      MPI_COMM_WORLD,code)
         call MPI_GATHER(taui_proc,1,MPI_DOUBLE_PRECISION,taui,1,MPI_DOUBLE_PRECISION,0,&
                      MPI_COMM_WORLD,code)
         call MPI_GATHER(tauf_proc,1,MPI_DOUBLE_PRECISION,tauf,1,MPI_DOUBLE_PRECISION,0,&
                      MPI_COMM_WORLD,code)
         !--Write eqk catalogue line-------------------------------------!
         !--Initiation time-duration-time since last eqk-x initiation-y initiation-x baricenter-y baricenter-number of points involved-coseismic moment-average stress drop-average stress before rupture-average stress after rupture-average coseismic slip--! 
         if (rang .eq. 0) then
            !-Catalogue entries with physical units----!
           t0evt=t0evt*cm%dc0/cm%vb0
           dt0evt=dt0evt*cm%dc0/cm%vb0
           tev=tev*cm%dc0/cm%vb0
           xi=xi*cm%mu*cm%dc0/(cm%b0*cm%sig0)
           yi=yi*cm%mu*cm%dc0/(cm%b0*cm%sig0)
           do i=1,nprocs
           	xa(i)=xa(i)*cm%mu*cm%dc0/(cm%b0*cm%sig0)
           	ya(i)=ya(i)*cm%mu*cm%dc0/(cm%b0*cm%sig0)
           	m(i)=m(i)*(cm%mu**3)*(cm%dc0**3)/((cm%b0**2)*(cm%sig0**2))
           	dtau(i)=dtau(i)*cm%b0*cm%sig0
            taui(i)=taui(i)*cm%b0*cm%sig0
            tauf(i)=tauf(i)*cm%b0*cm%sig0
            du(i)=du(i)*cm%dc0
           end do
            !--Write catalogue entries
            !write(90,rec=irecevt) real(t0evt,8), real(dt0evt,8), real(tev,8), real(xi,8), real(yi,8),&
            !                  real(sum(xa)/sum(nelevt),8), real(sum(ya)/sum(nelevt),8),real(sum(nelevt),8),&
            !                  real(sum(m),8),real(sum(dtau)/sum(nelevt),8),real(sum(taui)/sum(nelevt),8),&
            !                  real(sum(tauf)/sum(nelevt),8),real(sum(du)/sum(nelevt),8)
            !call flush(90)

            statusnetcdf=nf90_put_var(cr%nceqkcatid, cr%eqkcat_varid(1), t0evt, start=(/irecevt/))
            statusnetcdf=nf90_put_var(cr%nceqkcatid, cr%eqkcat_varid(2), dt0evt, start=(/irecevt/))
            statusnetcdf=nf90_put_var(cr%nceqkcatid, cr%eqkcat_varid(3), tev, start=(/irecevt/))
            statusnetcdf=nf90_put_var(cr%nceqkcatid, cr%eqkcat_varid(4), xi, start=(/irecevt/))
            statusnetcdf=nf90_put_var(cr%nceqkcatid, cr%eqkcat_varid(5), yi, start=(/irecevt/))
            statusnetcdf=nf90_put_var(cr%nceqkcatid, cr%eqkcat_varid(6), sum(xa)/sum(nelevt), start=(/irecevt/))
            statusnetcdf=nf90_put_var(cr%nceqkcatid, cr%eqkcat_varid(7), sum(ya)/sum(nelevt), start=(/irecevt/))
            statusnetcdf=nf90_put_var(cr%nceqkcatid, cr%eqkcat_varid(8), sum(nelevt), start=(/irecevt/))
            statusnetcdf=nf90_put_var(cr%nceqkcatid, cr%eqkcat_varid(9), sum(m), start=(/irecevt/))
            statusnetcdf=nf90_put_var(cr%nceqkcatid, cr%eqkcat_varid(10), sum(dtau)/sum(nelevt), start=(/irecevt/))
            statusnetcdf=nf90_put_var(cr%nceqkcatid, cr%eqkcat_varid(11), sum(du)/sum(nelevt), start=(/irecevt/))
            print*,' End of event, duration (s): ',tev,' magnitude: ',(2.0/3.0)*log10(sum(m))-6.0
         end if
         
      end if



end subroutine reccatal_3d_rg
