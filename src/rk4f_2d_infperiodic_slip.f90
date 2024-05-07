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


subroutine rk4f_2d_infperiodic_slip(phi_proc,th_proc,kappa_proc,&
                  a_proc,b_proc,dc_proc,s_proc,cc,cm,dt,dt0)

use mod_cst_fault_rns
implicit none
integer :: i,passtep,iterrk
real*8 :: dt,dt0,frns_slip,grns_slip
real*8 :: emax_proc,error
type(pmeca) :: cm
type(pcomput) :: cc
real*8, dimension(nprocs) :: emax
real*8, dimension(2*cc%nvalpr) :: e_proc
real*8, dimension(cc%nvalpr) :: a_proc,b_proc,dc_proc,s_proc
real*8, dimension(cc%nvalpr) :: phi_proc,th_proc
real*8, dimension(cc%nvalpr) :: dphi_proc,dth_proc
real*8, dimension(cc%nvalpr) :: kappa_proc,sk_proc,tbp_proc
real*8, dimension(cc%nvalpr) :: l1_proc,l2_proc,l3_proc,l4_proc,l5_proc,l6_proc
real*8, dimension(cc%nvalpr) :: m1_proc,m2_proc,m3_proc,m4_proc,m5_proc,m6_proc

passtep=0
iterrk=0

do while ((passtep .eq. 0) .and. (iterrk .le. cc%paramdt(1)))

  iterrk=iterrk+1  
  
  call tbp_rns_2d_infperiodic(tbp_proc,phi_proc,cm,cc)
  call lin_static_stress_periodic(sk_proc,phi_proc,kappa_proc,cc)

  do i=1,cc%nvalpr
        l1_proc(i)=frns_slip(phi_proc(i),th_proc(i),sk_proc(i),tbp_proc(i),a_proc(i),b_proc(i),dc_proc(i),s_proc(i),cm)
        m1_proc(i)=grns_slip(phi_proc(i),th_proc(i),dc_proc(i))
  end do

  call MPI_BARRIER(MPI_COMM_WORLD,code)

dphi_proc=cc%rkc(2,1)*dt*l1_proc
dth_proc=cc%rkc(2,1)*dt*m1_proc

  call tbp_rns_2d_infperiodic(tbp_proc,phi_proc+dphi_proc,cm,cc)
 call  lin_static_stress_periodic(sk_proc,phi_proc+dphi_proc,kappa_proc,cc)


!print*,' hello 7.3 '
  do i=1,cc%nvalpr
       
        l2_proc(i)=frns_slip(phi_proc(i)+dphi_proc(i),th_proc(i)+dth_proc(i),sk_proc(i),tbp_proc(i),&
        							a_proc(i),b_proc(i),dc_proc(i),s_proc(i),cm)
        m2_proc(i)=grns_slip(phi_proc(i)+dphi_proc(i),th_proc(i)+dth_proc(i),dc_proc(i))
  
  end do
!print*,' hello 7.4 '
  call MPI_BARRIER(MPI_COMM_WORLD,code)

dphi_proc=cc%rkc(3,1)*dt*l1_proc+cc%rkc(3,2)*dt*l2_proc
dth_proc=cc%rkc(3,1)*dt*m1_proc+cc%rkc(3,2)*dt*m2_proc

  call tbp_rns_2d_infperiodic(tbp_proc,phi_proc+dphi_proc,cm,cc)
 call  lin_static_stress_periodic(sk_proc,phi_proc+dphi_proc,kappa_proc,cc)

!print*,' hello 7.5 '
  do i=1,cc%nvalpr
        l3_proc(i)=frns_slip(phi_proc(i)+dphi_proc(i),th_proc(i)+dth_proc(i),sk_proc(i),tbp_proc(i),&
        							a_proc(i),b_proc(i),dc_proc(i),s_proc(i),cm)
        m3_proc(i)=grns_slip(phi_proc(i)+dphi_proc(i),th_proc(i)+dth_proc(i),dc_proc(i))
  end do
!print*,' hello 7.6 '
  call MPI_BARRIER(MPI_COMM_WORLD,code)
  
dphi_proc=cc%rkc(4,1)*dt*l1_proc+cc%rkc(4,2)*dt*l2_proc+cc%rkc(4,3)*dt*l3_proc
dth_proc=cc%rkc(4,1)*dt*m1_proc+cc%rkc(4,2)*dt*m2_proc+cc%rkc(4,3)*dt*m3_proc

  call tbp_rns_2d_infperiodic(tbp_proc,phi_proc+dphi_proc,cm,cc)
 call  lin_static_stress_periodic(sk_proc,phi_proc+dphi_proc,kappa_proc,cc)

!print*,' hello 7.7 '
  do i=1,cc%nvalpr
        l4_proc(i)=frns_slip(phi_proc(i)+dphi_proc(i),th_proc(i)+dth_proc(i),sk_proc(i),tbp_proc(i),&
        							a_proc(i),b_proc(i),dc_proc(i),s_proc(i),cm)
        m4_proc(i)=grns_slip(phi_proc(i)+dphi_proc(i),th_proc(i)+dth_proc(i),dc_proc(i))
  end do
!print*,' hello 7.8 '
  call MPI_BARRIER(MPI_COMM_WORLD,code)

dphi_proc=cc%rkc(5,1)*dt*l1_proc+cc%rkc(5,2)*dt*l2_proc+cc%rkc(5,3)*dt*l3_proc+cc%rkc(5,4)*dt*l4_proc
dth_proc=cc%rkc(5,1)*dt*m1_proc+cc%rkc(5,2)*dt*m2_proc+cc%rkc(5,3)*dt*m3_proc+cc%rkc(5,4)*dt*m4_proc

  call tbp_rns_2d_infperiodic(tbp_proc,phi_proc+dphi_proc,cm,cc)
 call  lin_static_stress_periodic(sk_proc,phi_proc+dphi_proc,kappa_proc,cc)

!print*,' hello 7.9 '
  do i=1,cc%nvalpr
        l5_proc(i)=frns_slip(phi_proc(i)+dphi_proc(i),th_proc(i)+dth_proc(i),sk_proc(i),tbp_proc(i),&
        							a_proc(i),b_proc(i),dc_proc(i),s_proc(i),cm)
        m5_proc(i)=grns_slip(phi_proc(i)+dphi_proc(i),th_proc(i)+dth_proc(i),dc_proc(i))
  end do
!print*,' hello 7.10 '
  call MPI_BARRIER(MPI_COMM_WORLD,code)

dphi_proc=cc%rkc(6,1)*dt*l1_proc+cc%rkc(6,2)*dt*l2_proc+cc%rkc(6,3)*dt*l3_proc+cc%rkc(6,4)*dt*l4_proc+cc%rkc(6,5)*dt*l5_proc
dth_proc=cc%rkc(6,1)*dt*m1_proc+cc%rkc(6,2)*dt*m2_proc+cc%rkc(6,3)*dt*m3_proc+cc%rkc(6,4)*dt*m4_proc+cc%rkc(6,5)*dt*m5_proc

  call tbp_rns_2d_infperiodic(tbp_proc,phi_proc+dphi_proc,cm,cc)
 call  lin_static_stress_periodic(sk_proc,phi_proc+dphi_proc,kappa_proc,cc)

!print*,' hello 7.11 '
  do i=1,cc%nvalpr
        l6_proc(i)=frns_slip(phi_proc(i)+dphi_proc(i),th_proc(i)+dth_proc(i),sk_proc(i),tbp_proc(i),&
        							a_proc(i),b_proc(i),dc_proc(i),s_proc(i),cm)
        m6_proc(i)=grns_slip(phi_proc(i)+dphi_proc(i),th_proc(i)+dth_proc(i),dc_proc(i))
  end do
!print*,' hello 7.12 ',dt
  call MPI_BARRIER(MPI_COMM_WORLD,code)

  do i=1,cc%nvalpr
!print*,' hello 7.12.0 '
	

     e_proc(i)=dt*abs(cc%rkr(1)*l1_proc(i)+cc%rkr(3)*l3_proc(i)+cc%rkr(4)*l4_proc(i)+cc%rkr(5)*l5_proc(i)+cc%rkr(6)*l6_proc(i))
     e_proc(i+cc%nvalpr)=dt*abs(cc%rkr(1)*m1_proc(i)+cc%rkr(3)*m3_proc(i)+cc%rkr(4)*m4_proc(i)&
     										+cc%rkr(5)*m5_proc(i)+cc%rkr(6)*m6_proc(i))
  end do

  emax_proc=maxval(e_proc)


  call MPI_GATHER(emax_proc,1,MPI_DOUBLE_PRECISION,emax,1,MPI_DOUBLE_PRECISION,0,&
                   MPI_COMM_WORLD,code)
!print*,' hello 7.14 '
  if (rang .eq. 0) then
     error=maxval(emax)
       if (error .le. cc%paramdt(2)) then  !--step accepted, variables updated, next time step increased
       passtep=1
	  dt0=dt
       dt=min(cc%paramdt(3)*dt*((cc%paramdt(2)/error)**(0.25)),cc%paramdt(4))
      else  !--step rejected, time step decreased
       dt=max(cc%paramdt(3)*dt*((cc%paramdt(2)/error)**(0.25)),cc%paramdt(5))
     end if
  end if

  call MPI_BCAST(dt0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,code)
  call MPI_BCAST(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,code)
  call MPI_BCAST(passtep,1,MPI_INTEGER,0,MPI_COMM_WORLD,code)

  call MPI_BARRIER(MPI_COMM_WORLD,code)

end do

phi_proc=phi_proc+dt0*(cc%rkf(1)*l1_proc+cc%rkf(3)*l3_proc+cc%rkf(4)*l4_proc+cc%rkf(5)*l5_proc+cc%rkf(6)*l6_proc)
th_proc=th_proc+dt0*(cc%rkf(1)*m1_proc+cc%rkf(3)*m3_proc+cc%rkf(4)*m4_proc+cc%rkf(5)*m5_proc+cc%rkf(6)*m6_proc)

  call MPI_BARRIER(MPI_COMM_WORLD,code)

end subroutine rk4f_2d_infperiodic_slip
