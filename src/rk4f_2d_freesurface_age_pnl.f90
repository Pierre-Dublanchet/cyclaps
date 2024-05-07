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


subroutine rk4f_2d_freesurface_age_pnl(t,phi_proc,th_proc,tbp_proc,Kstatic_proc,&
                  a_proc,b_proc,dc_proc,s_proc,p_proc,cc,cm,chy,dt,dt0)

use mod_cst_fault_rns
use mod_cst_fault_rns_p
implicit none
integer :: i,passtep,iterrk
real*8 :: t,dt,ddt,dt0
real*8 :: frns_age_p,grns_age_p
real*8 :: p0,pe,d0,de
real*8 :: emax_proc,error
type(pmeca) :: cm
type(phydro) :: chy
type(pcomput) :: cc
real*8, dimension(nprocs) :: emax
real*8, dimension(3*cc%nvalpr) :: e_proc
real*8, dimension(cc%nvalpr) :: a_proc,b_proc,dc_proc,s_proc,p_proc,diffu_proc
real*8, dimension(cc%nvalpr) :: phi_proc,th_proc
real*8, dimension(cc%nvalpr) :: dphi_proc,dth_proc,dp_proc
real*8, dimension(cc%nvalpr) :: sk_proc,tbp_proc
real*8, dimension(cc%nvalpr) :: l1_proc,l2_proc,l3_proc,l4_proc,l5_proc,l6_proc
real*8, dimension(cc%nvalpr) :: m1_proc,m2_proc,m3_proc,m4_proc,m5_proc,m6_proc
real*8, dimension(cc%nvalpr) :: k1_proc,k2_proc,k3_proc,k4_proc,k5_proc,k6_proc
real*8, dimension(cc%nvalpr,cc%n) :: Kstatic_proc

passtep=0
iterrk=0

do while ((passtep .eq. 0) .and. (iterrk .le. cc%paramdt(1)))

  iterrk=iterrk+1  
  
  call lin_static_freesurface(sk_proc,phi_proc,Kstatic_proc,cm,cc)
  call diffu(diffu_proc,phi_proc,th_proc,p_proc,s_proc,a_proc,b_proc,dc_proc,cm,chy,cc)
  call transfer_end_points_1d(diffu_proc,p_proc,p0,pe,d0,de,cc)
  call pressrate_nl1d(k1_proc,p_proc,p0,pe,diffu_proc,d0,de,cm,chy,cc,t)
 
 
  do i=1,cc%nvalpr
        l1_proc(i)=frns_age_p(phi_proc(i),th_proc(i),sk_proc(i),tbp_proc(i),p_proc(i),k1_proc(i),&
        				a_proc(i),b_proc(i),dc_proc(i),s_proc(i),cm)
        m1_proc(i)=grns_age_p(phi_proc(i),th_proc(i),p_proc(i),k1_proc(i),&
        					b_proc(i),dc_proc(i),s_proc(i),cm)
  end do

  call MPI_BARRIER(MPI_COMM_WORLD,code)
!-----------------------------------------------------------------------!
dp_proc=cc%rkc(2,1)*dt*k1_proc
dphi_proc=cc%rkc(2,1)*dt*l1_proc
dth_proc=cc%rkc(2,1)*dt*m1_proc
ddt=cc%rka(2)*dt

  call lin_static_freesurface(sk_proc,phi_proc+dphi_proc,Kstatic_proc,cm,cc)
  call diffu(diffu_proc,phi_proc+dphi_proc,th_proc+dth_proc,p_proc+dp_proc,&
  			s_proc,a_proc,b_proc,dc_proc,cm,chy,cc)
  call transfer_end_points_1d(diffu_proc,p_proc+dp_proc,p0,pe,d0,de,cc)
  call pressrate_nl1d(k2_proc,p_proc+dp_proc,p0,pe,diffu_proc,d0,de,cm,chy,cc,t+ddt)
 
 
  do i=1,cc%nvalpr
        l2_proc(i)=frns_age_p(phi_proc(i)+dphi_proc(i),th_proc(i)+dth_proc(i),sk_proc(i),tbp_proc(i),&
					p_proc(i)+dp_proc(i),k2_proc(i),&
        				a_proc(i),b_proc(i),dc_proc(i),s_proc(i),cm)
        m2_proc(i)=grns_age_p(phi_proc(i)+dphi_proc(i),th_proc(i)+dth_proc(i),p_proc(i)+dp_proc(i),k2_proc(i),&
        					b_proc(i),dc_proc(i),s_proc(i),cm)
  end do

  call MPI_BARRIER(MPI_COMM_WORLD,code)
!-----------------------------------------------------------------------!
dp_proc=cc%rkc(3,1)*dt*k1_proc+cc%rkc(3,2)*dt*k2_proc
dphi_proc=cc%rkc(3,1)*dt*l1_proc+cc%rkc(3,2)*dt*l2_proc
dth_proc=cc%rkc(3,1)*dt*m1_proc+cc%rkc(3,2)*dt*m2_proc
ddt=cc%rka(3)*dt

  call lin_static_freesurface(sk_proc,phi_proc+dphi_proc,Kstatic_proc,cm,cc)
  call diffu(diffu_proc,phi_proc+dphi_proc,th_proc+dth_proc,p_proc+dp_proc,&
  			s_proc,a_proc,b_proc,dc_proc,cm,chy,cc)
  call transfer_end_points_1d(diffu_proc,p_proc+dp_proc,p0,pe,d0,de,cc)
  call pressrate_nl1d(k3_proc,p_proc+dp_proc,p0,pe,diffu_proc,d0,de,cm,chy,cc,t+ddt)
 
  do i=1,cc%nvalpr
        l3_proc(i)=frns_age_p(phi_proc(i)+dphi_proc(i),th_proc(i)+dth_proc(i),sk_proc(i),tbp_proc(i),&
					p_proc(i)+dp_proc(i),k3_proc(i),&
        				a_proc(i),b_proc(i),dc_proc(i),s_proc(i),cm)
        m3_proc(i)=grns_age_p(phi_proc(i)+dphi_proc(i),th_proc(i)+dth_proc(i),p_proc(i)+dp_proc(i),k3_proc(i),&
        					b_proc(i),dc_proc(i),s_proc(i),cm)
  end do

  call MPI_BARRIER(MPI_COMM_WORLD,code)
!-----------------------------------------------------------------------!
dp_proc=cc%rkc(4,1)*dt*k1_proc+cc%rkc(4,2)*dt*k2_proc+cc%rkc(4,3)*dt*k3_proc  
dphi_proc=cc%rkc(4,1)*dt*l1_proc+cc%rkc(4,2)*dt*l2_proc+cc%rkc(4,3)*dt*l3_proc
dth_proc=cc%rkc(4,1)*dt*m1_proc+cc%rkc(4,2)*dt*m2_proc+cc%rkc(4,3)*dt*m3_proc
ddt=cc%rka(4)*dt

  call lin_static_freesurface(sk_proc,phi_proc+dphi_proc,Kstatic_proc,cm,cc)
 call diffu(diffu_proc,phi_proc+dphi_proc,th_proc+dth_proc,p_proc+dp_proc,&
  			s_proc,a_proc,b_proc,dc_proc,cm,chy,cc)
  call transfer_end_points_1d(diffu_proc,p_proc+dp_proc,p0,pe,d0,de,cc)
  call pressrate_nl1d(k4_proc,p_proc+dp_proc,p0,pe,diffu_proc,d0,de,cm,chy,cc,t+ddt)
 
  do i=1,cc%nvalpr
        l4_proc(i)=frns_age_p(phi_proc(i)+dphi_proc(i),th_proc(i)+dth_proc(i),sk_proc(i),tbp_proc(i),&
					p_proc(i)+dp_proc(i),k4_proc(i),&
        				a_proc(i),b_proc(i),dc_proc(i),s_proc(i),cm)
        m4_proc(i)=grns_age_p(phi_proc(i)+dphi_proc(i),th_proc(i)+dth_proc(i),p_proc(i)+dp_proc(i),k4_proc(i),&
        					b_proc(i),dc_proc(i),s_proc(i),cm)
  end do

  call MPI_BARRIER(MPI_COMM_WORLD,code)


!-----------------------------------------------------------------------!
dp_proc=cc%rkc(5,1)*dt*k1_proc+cc%rkc(5,2)*dt*k2_proc+cc%rkc(5,3)*dt*k3_proc+cc%rkc(5,4)*dt*k4_proc
dphi_proc=cc%rkc(5,1)*dt*l1_proc+cc%rkc(5,2)*dt*l2_proc+cc%rkc(5,3)*dt*l3_proc+cc%rkc(5,4)*dt*l4_proc
dth_proc=cc%rkc(5,1)*dt*m1_proc+cc%rkc(5,2)*dt*m2_proc+cc%rkc(5,3)*dt*m3_proc+cc%rkc(5,4)*dt*m4_proc
ddt=cc%rka(5)*dt

  call lin_static_freesurface(sk_proc,phi_proc+dphi_proc,Kstatic_proc,cm,cc)
 call diffu(diffu_proc,phi_proc+dphi_proc,th_proc+dth_proc,p_proc+dp_proc,&
  			s_proc,a_proc,b_proc,dc_proc,cm,chy,cc)
  call transfer_end_points_1d(diffu_proc,p_proc+dp_proc,p0,pe,d0,de,cc)
  call pressrate_nl1d(k5_proc,p_proc+dp_proc,p0,pe,diffu_proc,d0,de,cm,chy,cc,t+ddt)
 
  do i=1,cc%nvalpr
        l5_proc(i)=frns_age_p(phi_proc(i)+dphi_proc(i),th_proc(i)+dth_proc(i),sk_proc(i),tbp_proc(i),&
					p_proc(i)+dp_proc(i),k5_proc(i),&
        				a_proc(i),b_proc(i),dc_proc(i),s_proc(i),cm)
        m5_proc(i)=grns_age_p(phi_proc(i)+dphi_proc(i),th_proc(i)+dth_proc(i),p_proc(i)+dp_proc(i),k5_proc(i),&
        					b_proc(i),dc_proc(i),s_proc(i),cm)
  end do

  call MPI_BARRIER(MPI_COMM_WORLD,code)


!-----------------------------------------------------------------------!
dp_proc=cc%rkc(6,1)*dt*k1_proc+cc%rkc(6,2)*dt*k2_proc+cc%rkc(6,3)*dt*k3_proc+cc%rkc(6,4)*dt*k4_proc+cc%rkc(6,5)*dt*k5_proc
dphi_proc=cc%rkc(6,1)*dt*l1_proc+cc%rkc(6,2)*dt*l2_proc+cc%rkc(6,3)*dt*l3_proc+cc%rkc(6,4)*dt*l4_proc+cc%rkc(6,5)*dt*l5_proc
dth_proc=cc%rkc(6,1)*dt*m1_proc+cc%rkc(6,2)*dt*m2_proc+cc%rkc(6,3)*dt*m3_proc+cc%rkc(6,4)*dt*m4_proc+cc%rkc(6,5)*dt*m5_proc
ddt=cc%rka(6)*dt

  call lin_static_freesurface(sk_proc,phi_proc+dphi_proc,Kstatic_proc,cm,cc)
 call diffu(diffu_proc,phi_proc+dphi_proc,th_proc+dth_proc,p_proc+dp_proc,&
  			s_proc,a_proc,b_proc,dc_proc,cm,chy,cc)
  call transfer_end_points_1d(diffu_proc,p_proc+dp_proc,p0,pe,d0,de,cc)
  call pressrate_nl1d(k6_proc,p_proc+dp_proc,p0,pe,diffu_proc,d0,de,cm,chy,cc,t+ddt)
 
  do i=1,cc%nvalpr
        l6_proc(i)=frns_age_p(phi_proc(i)+dphi_proc(i),th_proc(i)+dth_proc(i),sk_proc(i),tbp_proc(i),&
					p_proc(i)+dp_proc(i),k6_proc(i),&
        				a_proc(i),b_proc(i),dc_proc(i),s_proc(i),cm)
        m6_proc(i)=grns_age_p(phi_proc(i)+dphi_proc(i),th_proc(i)+dth_proc(i),p_proc(i)+dp_proc(i),k6_proc(i),&
        					b_proc(i),dc_proc(i),s_proc(i),cm)
  end do

  call MPI_BARRIER(MPI_COMM_WORLD,code)


!-----------------------------------------------------------------------!

  do i=1,cc%nvalpr
	

     e_proc(i)=dt*abs(cc%rkr(1)*l1_proc(i)+cc%rkr(3)*l3_proc(i)+cc%rkr(4)*l4_proc(i)+cc%rkr(5)*l5_proc(i)+cc%rkr(6)*l6_proc(i))
     e_proc(i+cc%nvalpr)=dt*abs(cc%rkr(1)*m1_proc(i)+cc%rkr(3)*m3_proc(i)+cc%rkr(4)*m4_proc(i)&
     										+cc%rkr(5)*m5_proc(i)+cc%rkr(6)*m6_proc(i))
     e_proc(i+2*cc%nvalpr)=dt*abs(cc%rkr(1)*k1_proc(i)+cc%rkr(3)*k3_proc(i)+cc%rkr(4)*k4_proc(i)&
     										+cc%rkr(5)*k5_proc(i)+cc%rkr(6)*k6_proc(i))   										
  end do

  emax_proc=maxval(e_proc)


  call MPI_GATHER(emax_proc,1,MPI_DOUBLE_PRECISION,emax,1,MPI_DOUBLE_PRECISION,0,&
                   MPI_COMM_WORLD,code)

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

p_proc=p_proc+dt0*(cc%rkf(1)*k1_proc+cc%rkf(3)*k3_proc+cc%rkf(4)*k4_proc+cc%rkf(5)*k5_proc+cc%rkf(6)*k6_proc)
phi_proc=phi_proc+dt0*(cc%rkf(1)*l1_proc+cc%rkf(3)*l3_proc+cc%rkf(4)*l4_proc+cc%rkf(5)*l5_proc+cc%rkf(6)*l6_proc)
th_proc=th_proc+dt0*(cc%rkf(1)*m1_proc+cc%rkf(3)*m3_proc+cc%rkf(4)*m4_proc+cc%rkf(5)*m5_proc+cc%rkf(6)*m6_proc)

  call MPI_BARRIER(MPI_COMM_WORLD,code)

end subroutine rk4f_2d_freesurface_age_pnl
