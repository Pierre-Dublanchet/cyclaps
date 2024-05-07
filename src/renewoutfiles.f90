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


subroutine renewoutfiles(cm,cr)


use mod_cst_fault_rns
use netcdf
implicit none
type(pmeca) :: cm
type(prec) :: cr
integer :: i
character(len=10) :: numf


!-------------------------------------------------------------------------------------------------------------------!
!--Remove previous results------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------!
if (rang .eq. 0) then
   call system('rm -f '//trim(cr%path_res)//'parametres_fault_rns.cfg')
   call system('rm -f '//trim(cr%path_res)//'qmoy.nc')
   call system('rm -f '//trim(cr%path_res)//'qloc*.nc')
   call system('rm -f '//trim(cr%path_res)//'maps*.nc')
   call system('rm -f '//trim(cr%path_res)//'profilsx*.nc')
   call system('rm -f '//trim(cr%path_res)//'profilsy*.nc')
   call system('rm -f '//trim(cr%path_res)//'earthquake_catalogue.nc')
   call system('rm -f '//trim(cr%path_res)//'rupture_times_eqk1.nc')
end if

!-------------------------------------------------------------------------------------------------------------------!
!--Open output files------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------!

if (rang .eq. 0) then


   !--------------------------------------------------------!
   !--Open file for average/extremal variables timeseries---!
   !--------------------------------------------------------!
   statusnetcdf=nf90_create(trim(cr%path_res)//'qmoy.nc', NF90_CLASSIC_MODEL, cr%ncqmoyid)
            
   statusnetcdf=nf90_def_dim(cr%ncqmoyid, "time series",NF90_UNLIMITED,cr%qmoy_dimid)

   statusnetcdf=nf90_def_var(cr%ncqmoyid, "time",NF90_DOUBLE,[cr%qmoy_dimid],cr%qmoy_varid(1))
   statusnetcdf=nf90_def_var(cr%ncqmoyid, "time delai",NF90_DOUBLE,[cr%qmoy_dimid],cr%qmoy_varid(2))
   statusnetcdf=nf90_def_var(cr%ncqmoyid, "mean slip rate",NF90_DOUBLE,[cr%qmoy_dimid],cr%qmoy_varid(3))
   statusnetcdf=nf90_def_var(cr%ncqmoyid, "max slip rate",NF90_DOUBLE,[cr%qmoy_dimid],cr%qmoy_varid(4))
   statusnetcdf=nf90_def_var(cr%ncqmoyid, "mean slip",NF90_DOUBLE,[cr%qmoy_dimid],cr%qmoy_varid(5))
   statusnetcdf=nf90_def_var(cr%ncqmoyid, "mean shear stress",NF90_DOUBLE,[cr%qmoy_dimid],cr%qmoy_varid(6))
   statusnetcdf=nf90_def_var(cr%ncqmoyid, "mean state",NF90_DOUBLE,[cr%qmoy_dimid],cr%qmoy_varid(7))
   if ((cm%simlab .eq. '3d_infperiodic_regage') .or. (cm%simlab .eq. '3d_infperiodic_regslip')) then
   statusnetcdf=nf90_def_var(cr%ncqmoyid, "mean slip rate vw",NF90_DOUBLE,[cr%qmoy_dimid],cr%qmoy_varid(8))
   end if


   statusnetcdf=nf90_enddef(cr%ncqmoyid)

   !--------------------------------------------------------!
   !--Open file for local variables timeseries--------------!
   !--------------------------------------------------------!
   allocate(cr%ncqlocid(cr%nploc),cr%qloc_dimid(cr%nploc),cr%qloc_varid(cr%nploc,8))
   do i=1,cr%nploc
      select case (i)
         case (1:9)
            write(numf,'(I1)') i
         case (10:99)
            write(numf,'(I2)') i
      end select
      statusnetcdf=nf90_create(trim(cr%path_res)//'qloc'//trim(numf)//'.nc', NF90_CLASSIC_MODEL, cr%ncqlocid(i))

      statusnetcdf=nf90_def_dim(cr%ncqlocid(i), "time series",NF90_UNLIMITED,cr%qloc_dimid(i))

      statusnetcdf=nf90_def_var(cr%ncqlocid(i), "time",NF90_DOUBLE,[cr%qloc_dimid(i)],cr%qloc_varid(i,1))
      statusnetcdf=nf90_def_var(cr%ncqlocid(i), "time delai",NF90_DOUBLE,[cr%qloc_dimid(i)],cr%qloc_varid(i,2))
      statusnetcdf=nf90_def_var(cr%ncqlocid(i), "slip",NF90_DOUBLE,[cr%qloc_dimid(i)],cr%qloc_varid(i,3))
      statusnetcdf=nf90_def_var(cr%ncqlocid(i), "state",NF90_DOUBLE,[cr%qloc_dimid(i)],cr%qloc_varid(i,4))
      statusnetcdf=nf90_def_var(cr%ncqlocid(i), "slip rate",NF90_DOUBLE,[cr%qloc_dimid(i)],cr%qloc_varid(i,5))
      statusnetcdf=nf90_def_var(cr%ncqlocid(i), "shear stress",NF90_DOUBLE,[cr%qloc_dimid(i)],cr%qloc_varid(i,6))
      if ((cm%simlab(len_trim(cm%simlab)-2:len_trim(cm%simlab)) .eq. 'ess') &
      	.or. (cm%simlab(len_trim(cm%simlab)-2:len_trim(cm%simlab)) .eq. 'pnl')) then
		statusnetcdf=nf90_def_var(cr%ncqlocid(i), "pore pressure",NF90_DOUBLE,[cr%qloc_dimid(i)],cr%qloc_varid(i,7))
      end if
      if ((cm%simlab .eq. '2d_cr_regage_pnl') .or. (cm%simlab .eq. '2d_cr_regslip_pnl')) then
		statusnetcdf=nf90_def_var(cr%ncqlocid(i), "darcy vel",NF90_DOUBLE,[cr%qloc_dimid(i)],cr%qloc_varid(i,8))
      end if

      statusnetcdf=nf90_enddef(cr%ncqlocid(i))
   end do 

   !--------------------------------------------------------!
   !--Open earthquake catalogue file------------------------!
   !--------------------------------------------------------!
   statusnetcdf=nf90_create(trim(cr%path_res)//'earthquake_catalogue.nc', NF90_CLASSIC_MODEL, cr%nceqkcatid)

   statusnetcdf=nf90_def_dim(cr%nceqkcatid, "time series",NF90_UNLIMITED,cr%eqkcat_dimid)

   statusnetcdf=nf90_def_var(cr%nceqkcatid, "onset time",NF90_DOUBLE,[cr%eqkcat_dimid],cr%eqkcat_varid(1))
   statusnetcdf=nf90_def_var(cr%nceqkcatid, "onset time delai",NF90_DOUBLE,[cr%eqkcat_dimid],cr%eqkcat_varid(2))
   statusnetcdf=nf90_def_var(cr%nceqkcatid, "event duration",NF90_DOUBLE,[cr%eqkcat_dimid],cr%eqkcat_varid(3))
   statusnetcdf=nf90_def_var(cr%nceqkcatid, "x initiation",NF90_DOUBLE,[cr%eqkcat_dimid],cr%eqkcat_varid(4))
   statusnetcdf=nf90_def_var(cr%nceqkcatid, "y initiation",NF90_DOUBLE,[cr%eqkcat_dimid],cr%eqkcat_varid(5))
   statusnetcdf=nf90_def_var(cr%nceqkcatid, "x barycenter",NF90_DOUBLE,[cr%eqkcat_dimid],cr%eqkcat_varid(6))
   statusnetcdf=nf90_def_var(cr%nceqkcatid, "y barycenter",NF90_DOUBLE,[cr%eqkcat_dimid],cr%eqkcat_varid(7))
   statusnetcdf=nf90_def_var(cr%nceqkcatid, "number of elements",NF90_DOUBLE,[cr%eqkcat_dimid],cr%eqkcat_varid(8))
   statusnetcdf=nf90_def_var(cr%nceqkcatid, "coseismic moment",NF90_DOUBLE,[cr%eqkcat_dimid],cr%eqkcat_varid(9))
   statusnetcdf=nf90_def_var(cr%nceqkcatid, "coseismic stress drop",NF90_DOUBLE,[cr%eqkcat_dimid],cr%eqkcat_varid(10))
   statusnetcdf=nf90_def_var(cr%nceqkcatid, "coseismic slip",NF90_DOUBLE,[cr%eqkcat_dimid],cr%eqkcat_varid(11))
   if ((cm%simlab .eq. '3d_infperiodic_regage') .or. (cm%simlab .eq. '3d_infperiodic_regslip')) then
   		statusnetcdf=nf90_def_var(cr%nceqkcatid, "shear stress init",NF90_DOUBLE,[cr%eqkcat_dimid],cr%eqkcat_varid(12))
   		statusnetcdf=nf90_def_var(cr%nceqkcatid, "shear stress final",NF90_DOUBLE,[cr%eqkcat_dimid],cr%eqkcat_varid(13))
   end if

   statusnetcdf=nf90_enddef(cr%nceqkcatid)
end if


end subroutine renewoutfiles