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


program fault_2d_cr_slip

use mod_cst_fault_rns
use netcdf
implicit none
integer :: iter,i,k
integer :: irec,irecp,irecloc,irecevt,recpok,recok,reclocok,isevt
integer :: nelevt_proc,nelevt_proc0
integer, dimension(4) :: iterex
integer, dimension(:), allocatable :: nelevt,vrecpoints
integer, dimension(:), allocatable :: rkr,rks,ad
real*8 :: t,dt,dt0,dtrec,dtrecp,dtrecloc
real*8 :: start
real*8 :: vmax,vmax_proc
real*8 :: t0evt,dtevt,dt0evt,tev,m0_proc,m_proc,du_proc,dtau_proc
real*8 :: x0,y0,xi,yi
real*8, dimension(4) :: tex,vmaxex
real*8, dimension(:), allocatable :: vmax_allproc
real*8, dimension(:), allocatable :: m,du,dtau
real*8, dimension(:), allocatable :: phi,th,u,tau
real*8, dimension(:), allocatable :: a_proc,b_proc,dc_proc,s_proc
real*8, dimension(:), allocatable :: phi_proc,phiex_proc,th_proc,u_proc,tau_proc
real*8, dimension(:), allocatable :: v0_proc,kappa_proc
real*8, dimension(:), allocatable :: sk_proc,tbp_proc
real*8, dimension(:), allocatable :: xa_proc,ya_proc,u0_proc,tau0_proc
real*8, dimension(:), allocatable :: xa,ya
type(pmeca) :: cm
type(pcomput) :: cc
type(prec) :: cr
character(len=10) :: numf


!-------------------------------------------------------------------------------------------------------------------!
!--Parallel Initialization------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------!
 call MPI_INIT(code)
 call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,code)    
 call MPI_COMM_RANK(MPI_COMM_WORLD,rang,code)   

!-------------------------------------------------------------------------------------------------------------------!
!--Debut de l'horloge-----------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------!
 if (rang .eq. 0) then 
    start = MPI_WTIME()
 end if                

!-------------------------------------------------------------------------------------------------------------------!
!--Parameters-------------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------!
 call param_fault_rns(cm,cc,cr)

!-------------------------------------------------------------------------------------------------------------------!
!--Remove previous results and open output files--------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------!
call renewoutfiles(cm,cr)

!-------------------------------------------------------------------------------------------------------------------!
!--Allocation-------------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------!
allocate(a_proc(cc%nvalpr),b_proc(cc%nvalpr),s_proc(cc%nvalpr),dc_proc(cc%nvalpr))
allocate(phi_proc(cc%nvalpr),phiex_proc(cc%nvalpr),th_proc(cc%nvalpr),u_proc(cc%nvalpr),tau_proc(cc%nvalpr))
allocate(sk_proc(cc%nvalpr),tbp_proc(cc%nvalpr))
allocate(v0_proc(2*cc%nvalpr),kappa_proc(2*cc%nvalpr))

if (cr%catal .eq. 1) then
   allocate(xa_proc(cc%nvalpr),ya_proc(cc%nvalpr),u0_proc(cc%nvalpr),tau0_proc(cc%nvalpr))
   if (rang .eq. 0) then
      allocate(xa(nprocs),ya(nprocs),nelevt(nprocs),m(nprocs),du(nprocs),dtau(nprocs))
   end if
end if

allocate(rkr(nprocs),rks(nprocs),ad(nprocs))

 if (rang .eq. 0) then
    allocate(phi(cc%n),th(cc%n),u(cc%n),tau(cc%n))
    allocate(vmax_allproc(nprocs))
 end if

allocate(vrecpoints(cr%nploc))
!-------------------------------------------------------------------------------------------------------------------!
!--Find components to be recorded-----------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------!
if (cr%qloc .eq. 1) then 
   do i=1,cr%nploc
      call x2ci_2d_rg(vrecpoints(i),cc,cr%xr(i))
   end do
end if

!-------------------------------------------------------------------------------------------------------------------!
!--Initialisation-Normalisation-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------!
 call init_rns_2d_rg(a_proc,b_proc,dc_proc,s_proc,phi_proc,th_proc,u_proc,cm,cc,cr)
 t=0.0
 tex(1)=0.0
 tex(2)=0.0
 tex(3)=0.0
 tex(4)=0.0
 irec=0
 irecp=0
 irecevt=0
 dtrec=0.0
 dtrecp=0.0
 recpok=0
 recok=0
 iterex(1)=0
 iterex(2)=0
 iterex(3)=0
 iterex(4)=0
 isevt=0
 dtevt=0.0
 
 call normalize_rns(a_proc,b_proc,dc_proc,s_proc,u_proc,phi_proc,th_proc,cm,cc,cr)
 
    call MPI_GATHER(phi_proc,cc%nvalpr,MPI_DOUBLE_PRECISION,phi,cc%nvalpr,MPI_DOUBLE_PRECISION,0,&
                   MPI_COMM_WORLD,code)
   call MPI_GATHER(th_proc,cc%nvalpr,MPI_DOUBLE_PRECISION,th,cc%nvalpr,MPI_DOUBLE_PRECISION,0,&
                   MPI_COMM_WORLD,code)
                   
    if (rang .eq. 0) then
      vmax=maxval(exp(phi(:)))

       dt=0.01/vmax
         
      vmaxex(1)=vmax
      vmaxex(2)=vmax
      vmaxex(3)=vmax
      vmaxex(4)=vmax
   end if
   call MPI_BCAST(vmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,code)
   call MPI_BCAST(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,code)

   call MPI_BARRIER(MPI_COMM_WORLD,code)
   
   !------------------------------------------------------------------------------------------------------!
   !--Compute modes of elastic kernel------------------------------------------------------------!
      !------------------------------------------------------------------------------------------------------!
 call elastic_kernel_2d_cr(kappa_proc,cc,cm)
 
 if (nprocs .eq. 2) then
 	rks(1)=0
 	rks(2)=1 
  	rkr(1)=0
 	rkr(2)=1 
  	ad(1)=1+cc%nvalpr
 	ad(2)=1 
 elseif (nprocs .eq. 1) then
 	rks(1)=0
 	rkr(1)=0
 	ad(1)=1+cc%nvalpr/2
 else 
 	do k=0,nprocs/2-1
         rks(2*k+1)=2*k
         rks(2*k+2)=2*k+1
         rkr(2*k+1)=nprocs/4+k
         rkr(2*k+2)=nprocs/4+k
         ad(2*k+1)=1
         ad(2*k+2)=1+cc%nvalpr
  	end do
  end if
  

 
 call MPI_BARRIER(MPI_COMM_WORLD,code)               

!-------------------------------------------------------------------------------------------------------------------!
!--Copy parameters/initiation files---------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------!

if (rang .eq. 0) then
   call system('cp parametres_fault_rns.cfg '//trim(cr%path_res))
      if (cc%methinit .eq. 1) then
      call system('cp '//trim(cc%path_init)//'init.nc '//trim(cr%path_res))
   end if
end if



!-------------------------------------------------------------------------------------------------------------------!
!--Iterations-------------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------!
      if (rang .eq. 0) then
          print*,' Starting iterations '
      end if
      
      
      
do iter=1,cc%niter

  if (mod(iter,cc%niterscreen) .eq. 0) then
      if (rang .eq. 0) then
                 print*,' Iteration ',iter,' ,log vmax = ',log10(vmax*cm%vb0),'  (m/s), time t = ',&
         t*cm%dc0/(365.0*24.0*3600.0*cm%vb0),' years'
      end if
   end if





  !------------------------------------------------------------------------------------------!
  !--Conditions for variable record----------------------------------------------------------!
  !------------------------------------------------------------------------------------------!
   vmax_proc=maxval(exp(phi_proc(:)))
   call MPI_GATHER(vmax_proc,1,MPI_DOUBLE_PRECISION,vmax_allproc,1,MPI_DOUBLE_PRECISION,0,&
                   MPI_COMM_WORLD,code)
   if (rang .eq. 0) then
	vmax=maxval(vmax_allproc(:))
   end if
   call MPI_BCAST(vmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,code)
   
   call MPI_BARRIER(MPI_COMM_WORLD,code)
   
  call varecok_2d(cr,iter,iterex,vmax,vmaxex,t,tex,dtrec,dtrecp,dtrecloc,recok,reclocok,recpok)
  !---------------------------------------------------------------------------------------------------------!
  !--Gather Variables for recording-------------------------------------------------------------------------!
  !---------------------------------------------------------------------------------------------------------!
  if ((recok .eq. 1) .or. (recpok .eq. 1)) then

     do i=1,cc%nvalpr
        tau_proc(i)=(cm%mu0+cm%r*a_proc(i)*phi_proc(i)+b_proc(i)*th_proc(i))*s_proc(i)
     end do

     call MPI_GATHER(tau_proc,cc%nvalpr,MPI_DOUBLE_PRECISION,tau,cc%nvalpr,MPI_DOUBLE_PRECISION,0,&
                      MPI_COMM_WORLD,code)
     call MPI_GATHER(u_proc,cc%nvalpr,MPI_DOUBLE_PRECISION,u,cc%nvalpr,MPI_DOUBLE_PRECISION,0,&
                      MPI_COMM_WORLD,code)
     call MPI_GATHER(phi_proc,cc%nvalpr,MPI_DOUBLE_PRECISION,phi,cc%nvalpr,MPI_DOUBLE_PRECISION,0,&
                   MPI_COMM_WORLD,code)
   call MPI_GATHER(th_proc,cc%nvalpr,MPI_DOUBLE_PRECISION,th,cc%nvalpr,MPI_DOUBLE_PRECISION,0,&
                   MPI_COMM_WORLD,code)

  end if
  !---------------------------------------------------------------------------------------------------------!
  !--Record spatial averages+extreme values-----------------------------------------------------------------!
  !---------------------------------------------------------------------------------------------------------!
  if (recok .eq. 1) then
     call recqmoy(recok,cm,cc,cr,phi,th,u,tau,irec,t,dtrec,vmax)
  end if
  !---------------------------------------------------------------------------------------------------------!
  !--Record local values------------------------------------------------------------------------------------!
  !---------------------------------------------------------------------------------------------------------!
  if (reclocok .eq. 1) then
      call recqloc(reclocok,vrecpoints,cm,cc,cr,phi,th,u,tau,irecloc,t,dtrecloc)
   end if
  !-------------------------------------------------------------------------------------------------!
  !--Record profiles--------------------------------------------------------------------------------!
  !-------------------------------------------------------------------------------------------------!
  if (recpok .eq. 1) then
     call recprofils(recpok,cm,cr,cc,phi,th,u,tau,irecp,t,dtrecp)
  end if
  !--------------------------------------------------------------------------------------------------------!
  !--Catalogue Record--------------------------------------------------------------------------------------!
  !--------------------------------------------------------------------------------------------------------!
  if (cr%catal .eq. 1) then   
      call reccatal_2d_rg(cc,cm,cr,phi_proc,th_proc,tau_proc,u_proc,s_proc,a_proc,b_proc,irecevt,t,vmax,isevt,&
      			 t0evt,dt0evt,dtevt,m0_proc,&
                    m_proc,xa_proc,ya_proc,u0_proc,tau0_proc,du_proc,dtau_proc,nelevt_proc,nelevt_proc0,&
                    x0,y0,xa,ya,nelevt,xi,yi,m,du,dtau,tev)
  end if        
  
  call MPI_BARRIER(MPI_COMM_WORLD,code)
 !--------------------------------------------------------------------------------------------------!
 !--Variables update--------------------------------------------------------------------------------!
 !--------------------------------------------------------------------------------------------------!      
 
 phiex_proc(:)=phi_proc(:)
                  
call rk4f_2d_cr_slip(phi_proc,th_proc,tbp_proc,kappa_proc,v0_proc,rkr,rks,ad,&
                  a_proc,b_proc,dc_proc,s_proc,cc,cm,dt,dt0)

  t=t+dt0
  dtrec=dtrec+dt0
  dtrecp=dtrecp+dt0
  dtrecloc=dtrecloc+dt0
  dtevt=dtevt+dt0


  do i=1,cc%nvalpr
     u_proc(i)=u_proc(i)+0.5*dt0*(exp(phi_proc(i))+exp(phiex_proc(i)))
  end do


end do

 if (rang .eq. 0) then
   statusnetcdf=nf90_close(cr%ncqmoyid)
   do i=1,cr%nploc
      statusnetcdf=nf90_close(cr%ncqlocid(i))
   end do
   statusnetcdf=nf90_close(cr%nceqkcatid)
end if


!-------------------------------------------------------------------------------------------------------------------! 
!--End Parallel-----------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------!
 call MPI_FINALIZE(code)

end program fault_2d_cr_slip

include 'param_fault_rns.f90'
include 'init_rns_2d_rg.f90'
include 'renewoutfiles.f90'
include 'normalize_rns.f90'
include 'ci2x_2d_rg.f90'
include 'x2ci_2d_rg.f90'
include 'frns_slip.f90'
include 'grns_slip.f90'
include 'elastic_kernel_2d_cr.f90'
include 'lin_static_stress_cr.f90'
include 'rk4f_2d_cr_slip.f90'
include 'tbp_rns_2d_cr.f90'
include 'varecok_2d.f90'
include 'recqmoy.f90'
include 'recqloc.f90'
include 'recprofils.f90'
include 'reccatal_2d_rg.f90'


