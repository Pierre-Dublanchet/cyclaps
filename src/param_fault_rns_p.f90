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


subroutine param_fault_rns_p(cm,cc,cr,chy)

use mod_cst_fault_rns
use mod_cst_fault_rns_p
use cfgio_mod
implicit none
integer :: i,j
real*8, dimension(2) :: testrec
logical :: repexist
type(pmeca) :: cm
type(pcomput) :: cc
type(prec) :: cr
type(phydro) :: chy
type(cfg_t) :: cfg

 !-----------------------------!
 !-Read parameters from file---!
 !-----------------------------!
 cfg=parse_cfg("parametres_fault_rns.cfg")


call cfg%get("section 1", "simlab", cm%simlab)
call cfg%get("section 1", "slip_mode", cm%slipmode)
call cfg%get("section 1", "young_mod", cm%Y)
call cfg%get("section 1", "poisson_ratio", cm%nu)
call cfg%get("section 1", "rho_r", cm%rho)
call cfg%get("section 1", "slab_thickness", cm%H)
call cfg%get("section 1", "ext_shear_stressing", cm%paramdsext)
call cfg%get("section 1", "ref_fric_coeff", cm%mu0)
call cfg%get("section 1", "dl_coeff", cm%alpha)
call cfg%get("section 1", "vplate", cm%vb)
call cfg%get("section 1", "vstar", cm%vb0)
call cfg%get("section 1", "vsis", cm%vsis)
call cfg%get("section 1", "nasp",cm%n_asp)
call cfg%get("section 1", "a",cm%a_asp)
call cfg%get("section 1", "b",cm%b_asp)
call cfg%get("section 1", "dcasp",cm%dc_asp)
call cfg%get("section 1", "xcasp",cm%xc_asp)
call cfg%get("section 1", "ycasp",cm%yc_asp)
call cfg%get("section 1", "Rasp",cm%R_asp)
call cfg%get("section 1", "sigasp",cm%sig_asp)
call cfg%get("section 1", "viasp",cm%vi_asp)
call cfg%get("section 1", "thiasp",cm%thi_asp)
call cfg%get("section 1", "uiasp",cm%ui_asp)




call cfg%get("section 1", "permea",chy%k)
call cfg%get("section 1", "porosity",chy%phi)
call cfg%get("section 1", "compress",chy%cf)
call cfg%get("section 1", "viscosity",chy%etaf)
call cfg%get("section 1", "rho_f",chy%rhof)
call cfg%get("section 1", "paraminj", chy%paraminj)

call cfg%get("section 1", "nx",cc%nx)
call cfg%get("section 1", "ny",cc%ny)
call cfg%get("section 1", "dx",cc%dx)
call cfg%get("section 1", "dy",cc%dy)
call cfg%get("section 1", "niter",cc%niter)
call cfg%get("section 1", "nit_screen",cc%niterscreen)
call cfg%get("section 1", "meth_init",cc%methinit)
call cfg%get("section 1", "paramdt",cc%paramdt)
call cfg%get("section 1", "pathinit",cc%path_init)

call cfg%get("section 1", "qcat",cr%catal)
call cfg%get("section 1", "qmoy",cr%qmoy)
call cfg%get("section 1", "qprof",cr%profil)
call cfg%get("section 1", "qloc",cr%qloc)
call cfg%get("section 1", "qprofhvc",cr%profilhvc)
call cfg%get("section 1", "vfrec",cr%vfrec)
call cfg%get("section 1", "tfrec",cr%tfrec)
call cfg%get("section 1", "dtfrec",cr%dtfrec)
call cfg%get("section 1", "nitrec",cr%nitrec)
call cfg%get("section 1", "pathres",cr%path_res)
call cfg%get("section 1", "nploc",cr%nploc)
call cfg%get("section 1", "xrloc",cr%xr)
call cfg%get("section 1", "yrloc",cr%yr)
call cfg%get("section 1", "hboxm",cr%hboxm)



call cfg%finalize()



!---------------------------------------------!
!--Reference hydraulic diffusivity------------!
!---------------------------------------------!
chy%D=chy%k/(chy%phi*chy%cf*chy%etaf)
 
!---------------------------------------------!
 !--Runge Kutta Felhberg Coefficients--!
 !---------------------------------------------!
  cc%rka(1)=0.0
  cc%rka(2)=1.0/4.0
  cc%rka(3)=3.0/8.0
  cc%rka(4)=12.0/13.0
  cc%rka(5)=1.0
  cc%rka(6)=1.0/2.0
  
  cc%rkc(2,1)=1.0/4.0
  cc%rkc(3,1)=3.0/32.0
  cc%rkc(3,2)=9.0/32.0
  cc%rkc(4,1)=1932.0/2197.0
  cc%rkc(4,2)=-7200.0/2197.0
  cc%rkc(4,3)=7296.0/2197.0
  cc%rkc(5,1)=439.0/216.0
  cc%rkc(5,2)=-8.0
   cc%rkc(5,3)=3680.0/513.0
   cc%rkc(5,4)=-845.0/4104.0
   cc%rkc(6,1)=-8.0/27.0
   cc%rkc(6,2)=2.0
   cc%rkc(6,3)=-3544.0/2565.0
   cc%rkc(6,4)=1859.0/4104.0
   cc%rkc(6,5)=-11.0/40.0
   
   cc%rkr(1)  =  1.0/360.0
   cc%rkr(3)  = -128.0/4275.0
  cc%rkr(4)  = -2197.0/75240.0
 cc%rkr(5)  =  1.0/50.0
 cc%rkr(6)  =  2.0/55.0

 cc%rkf(1)  =  16.0/135.0
 cc%rkf(3)  =  6656.0/12825.0
 cc%rkf(4)  =  28561.0/56430.0
 cc%rkf(5)  = -9.0/50.0
 cc%rkf(6)  = 2.0/55.0

 !-----------------------------!
 !-Total number of points------!
 !-----------------------------!
 cc%n=cc%nx*cc%ny
 select case (cm%simlab)
      case ('2d_cr_slip_press','2d_cr_aging_press','2d_cr_slip_pnl','2d_cr_aging_pnl','2d_cr_regslip_pnl','2d_cr_regage_pnl')
         cc%n1=2*cc%n
      case default
         cc%n1=cc%n
 end select
 
 !-------------------------------!
 !--slip mode factor------!
 !----------------------------!
 if (cm%slipmode .eq. 2) then
     cm%sm_factor=1-cm%nu
     elseif (cm%slipmode .eq. 3) then
        cm%sm_factor=1
 end if

 if (rang .eq. 0) then

 !-----------------------------!
 !-Print parameters------------!
 !-----------------------------!

 do i=1,2
    print*,'--------------------------------------------------------------------------'
 end do
 print*,'-----------------------CYCLAPS (earthquake CYCLe simulator for multiple AsPerities fault Systems)-----------------'
 print*,'-------------------------------P. Dublanchet et al. (February 2024)-----------'
 do i=1,2
    print*,'--------------------------------------------------------------------------'
 end do

select case (cm%simlab)
   case ('2d_infperiodic_aging_pnl')
      print*,'1d (mode II or III) fault between 2d elastic slabs with replication along strike'
      print*,'classical rate-and-state friction with aging law'
      print*,'numerical resolution of non-linear pore pressure diffusion'
   case ('2d_infperiodic_slip_pnl')
      print*,'1d (mode II or III) fault between 2d elastic slabs with replication along strike'
      print*,'classical rate-and-state friction with slip law'
      print*,'numerical resolution of non-linear pore pressure diffusion'
   case ('2d_infperiodic_aging_press')
      print*,'1d (mode II or III) fault between 2d elastic slabs with replication along strike'
      print*,'classical rate-and-state friction with aging law'
      print*,'imposed pore pressure history'
   case ('2d_infperiodic_slip_press')
      print*,'1d (mode II or III) fault between 2d elastic slabs with replication along strike'
      print*,'classical rate-and-state friction with slip law'
      print*,'imposed pore pressure history'
   case ('2d_freesurface_aging_pnl')
      print*,'1d strike-slip fault in a semi inf elastic half space with free surface'
      print*,'classical rate-and-state friction with aging law'
      print*,'numerical resolution of non-linear pore pressure diffusion'
   case ('2d_freesurface_slip_pnl')
      print*,'1d strike-slip fault in a semi inf elastic half space with free surface'
      print*,'classical rate-and-state friction with slip law'
      print*,'numerical resolution of non-linear pore pressure diffusion'
   case ('2d_freesurface_aging_press')
      print*,'1d strike-slip fault in a semi inf elastic half space with free surface'
      print*,'classical rate-and-state friction with aging law'
      print*,'imposed pore pressure history'
   case ('2d_freesurface_slip_press')
      print*,'1d strike-slip fault in a semi inf elastic half space with free surface'
      print*,'classical rate-and-state friction with slip law'
      print*,'imposed pore pressure history'
   case ('2d_cr_aging_pnl')
      print*,'1d (mode II or III) fault between 2d semi inf elastic half spaces, Cochard & Rice Kernel (no replication)'
      print*,'classical rate-and-state friction with aging law'
      print*,'numerical resolution of non-linear pore pressure diffusion'
   case ('2d_cr_slip_pnl')
      print*,'1d (mode II or III) fault between 2d semi inf elastic half spaces, Cochard & Rice Kernel (no replication)'
      print*,'classical rate-and-state friction with slip law'
      print*,'numerical resolution of non-linear pore pressure diffusion'
   case ('2d_cr_regage_pnl')
      print*,'1d (mode II or III) fault between 2d semi inf elastic half spaces, Cochard & Rice Kernel (no replication)'
      print*,'regularized rate-and-state friction with aging law'
      print*,'numerical resolution of non-linear pore pressure diffusion'
   case ('2d_cr_regslip_pnl')
      print*,'1d (mode II or III) fault between 2d semi inf elastic half spaces, Cochard & Rice Kernel (no replication)'
      print*,'regularized rate-and-state friction with slip law'
      print*,'numerical resolution of non-linear pore pressure diffusion'
   case ('2d_cr_aging_press')
      print*,'1d (mode II or III) fault between 2d semi inf elastic half spaces, Cochard & Rice Kernel (no replication)'
      print*,'classical rate-and-state friction with aging law'
      print*,'imposed pore pressure history'
   case ('2d_cr_slip_press')
      print*,'1d (mode II or III) fault between 2d semi inf elastic half spaces, Cochard & Rice Kernel (no replication)'
      print*,'classical rate-and-state friction with slip law'
      print*,'imposed pore pressure history'
   case ('3d_infperiodic_aging_pnl')
      print*,'2d fault between 3d elastic slabs with replication along strike and depth'
      print*,'classical rate-and-state friction with aging law'
      print*,'numerical resolution of non-linear pore pressure diffusion'
   case ('3d_infperiodic_slip_pnl')
      print*,'2d fault between 3d elastic slabs with replication along strike and depth'
      print*,'classical rate-and-state friction with slip law'
      print*,'numerical resolution of non-linear pore pressure diffusion'
   case ('3d_infperiodic_aging_press')
      print*,'2d fault between 3d elastic slabs with replication along strike and depth'
      print*,'classical rate-and-state friction with aging law'
      print*,'imposed pore pressure history'
   case ('3d_infperiodic_slip_press')
      print*,'2d fault between 3d elastic slabs with replication along strike and depth'
      print*,'classical rate-and-state friction with slip law'
      print*,'imposed pore pressure history'
   case ('3d_infperiodic_aging_pnl_rs')
      print*,'2d fault between 3d elastic slabs with replication along strike and depth'
      print*,'classical rate-and-state friction with aging law'
      print*,'numerical resolution of non-linear pore pressure diffusion with real injection scenario'
end select


print*,' Young modulus (GPa): ',cm%Y/1e9
print*,' Poisson ratio: ',cm%nu
print*,' Density (kg/m3): ',cm%rho
print*,' Slab thickness (m): ',cm%H
print*,'  Fault Length (m): ',cc%nx*cc%dx,' by ',cc%ny*cc%dy
print*,' Constant friction coefficient: ',cm%mu0 
print*,' Dieterich-Linker coefficient: ',cm%alpha
print*,' Tectonic slip rate (m/s): ',cm%vb

print*,' Number of points in x: ',cc%nx,' Number of points in y: ',cc%ny,' Total number of points: ',cc%n
print*,' Number of processors: ',nprocs
print*,' Number of points per processor: ',cc%n/nprocs
print*,' Computational cell size: ',cc%dx,' and ',cc%dy
print*,' Number of iterations: ',cc%niter

if (cc%methinit .eq. 0) then
    print*,' Initial conditions defined in the code '
elseif (cc%methinit .eq. 1) then
    print*,' Initial conditions read from init files '
end if
print*,'Initial conditions read in: ',cc%path_init

print*,' Runge-Kutta Fehlberg algorithm used '
print*,'Maximum number of RK iterations:',cc%paramdt(1)
print*,' Error tolerance: ',cc%paramdt(2)
print*,' Safety factor: ',cc%paramdt(3)
print*,' Max dt (s): ',cc%paramdt(4)
print*,' Min dt (s): ',cc%paramdt(5)
 
 if (cr%catal .eq. 1) then 
    print*,' Catalogue produced '
 end if
 if (cr%qmoy .eq. 1) then
    print*,' Average variables recorded '
 end if
 if (cr%qloc .eq. 1) then
       print*,' Local variables recorded '
 end if
  if (cr%profil .eq. 1) then
       print*,' Variable profils (1d fault)/maps (2d fault) recorded '
 end if
  if (cr%profilhvc .eq. 1) then
       print*,' Variable hor. and vert. profiles (2d fault) recorded '
 end if

 if (cr%vfrec(1) .gt. 0.0) then
    print*,' Record averages when vmax is * or / by: ',cr%vfrec(1)
 end if
 if (cr%vfrec(2) .gt. 0.0) then
    print*,' Record profils when vmax is * or / by: ',cr%vfrec(2)
 end if
 if (cr%vfrec(3) .gt. 0.0) then
    print*,' Record local variables when vmax is * or / by: ',cr%vfrec(3)
 end if
 if (cr%vfrec(4) .gt. 0.0) then
    print*,' Record hor/vert profiles when vmax is * or / by: ',cr%vfrec(4)
 end if

 if (cr%tfrec(1) .gt. 0.0) then
    print*,' Record averages when time is * or / by: ',cr%tfrec(1)
 end if
 if (cr%tfrec(2) .gt. 0.0) then
    print*,' Record profils when time is * or / by: ',cr%tfrec(2)
 end if
  if (cr%tfrec(3) .gt. 0.0) then
    print*,' Record local variables when time is * or / by: ',cr%tfrec(3)
 end if
 if (cr%tfrec(4) .gt. 0.0) then
    print*,' Record hor/vert profiles when time is * or / by: ',cr%tfrec(4)
 end if

 if (cr%dtfrec(1) .gt. 0.0) then
    print*,' Record averages every ',cr%dtfrec(1),'  seconds' 
 end if
 if (cr%dtfrec(2) .gt. 0.0) then
    print*,' Record profils every ',cr%dtfrec(2),' seconds'
 end if
  if (cr%dtfrec(3) .gt. 0.0) then
    print*,' Record local variables every ',cr%dtfrec(3),'  seconds' 
 end if
 if (cr%dtfrec(4) .gt. 0.0) then
    print*,' Record hor/vert profiles every ',cr%dtfrec(4),' seconds'
 end if

 if (cr%nitrec(1) .gt. 0.0) then
    print*,' Record averages every ',cr%nitrec(1),' iteration' 
 end if
 if (cr%nitrec(2) .gt. 0.0) then
    print*,' Record profils every ',cr%nitrec(2),' iteration'
 end if
  if (cr%nitrec(3) .gt. 0.0) then
    print*,' Record local variables every ',cr%nitrec(3),' iteration' 
 end if
 if (cr%nitrec(4) .gt. 0.0) then
    print*,' Record hor/vert profiles every ',cr%nitrec(4),' iteration'
 end if

 print*,' Results stored in: ',cr%path_res

  do i=1,2
    print*,'--------------------------------------------------------------------------'
 end do

 !-----------------------------!
 !-Print error messages--------!
 !-----------------------------!
 testrec(1)=0.0
 testrec(2)=0.0
 do i=1,4
    select case (i)
      case (1)
        do j=1,2
           if (cr%vfrec(j) .gt. 0.0) then
              testrec(j)=testrec(j)+cr%vfrec(j)/cr%vfrec(j)
           end if
        end do
      case (2)
        do j=1,2
           if (cr%tfrec(j) .gt. 0.0) then
              testrec(j)=testrec(j)+cr%tfrec(j)/cr%tfrec(j)
           end if
        end do
      case (3)
        do j=1,2
           if (cr%dtfrec(j) .gt. 0.0) then
              testrec(j)=testrec(j)+cr%dtfrec(j)/cr%dtfrec(j)
           end if
        end do
      case (4)
        do j=1,2
           if (cr%nitrec(j) .gt. 0.0) then
              testrec(j)=testrec(j)+cr%nitrec(j)/cr%nitrec(j)
           end if
        end do
     end select
 end do

 if ((testrec(1) .gt. 1.0) .or. (testrec(2) .gt. 1.0)) then 
    print*,' error: more than two recording methods have been used'
    !stop
 end if

 inquire(file=trim(cc%path_init), exist=repexist)
 if (.not. repexist) then
    print*,'error: ',trim(cc%path_init),' does not exist'
    !stop
 end if

 inquire(file=trim(cr%path_res), exist=repexist)
 if (.not. repexist) then
    print*,'error: ',trim(cr%path_res),' does not exist'
    !stop
 end if

 if (nprocs .gt. cc%n) then
    print*,' error: more processors than number of points '
    stop
 end if


 end if

 !-----------------------------!
 !-Data splitting--------------!
 !-----------------------------! 
 cc%nvalpr=cc%n/nprocs
 


end subroutine param_fault_rns_p

