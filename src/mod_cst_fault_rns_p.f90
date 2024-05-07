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


module mod_cst_fault_rns_p

use mpi
implicit none
include 'debug.h'
include 'stat.h'



type phydro
  integer :: rang_inj		   !-rank of injection node	
  integer :: i_inj_proc		   !-index of injection node	
  real*8 :: k                  !-Fault Permeability
  real*8 :: phi                !-Fault Porosity
  real*8 :: cf                 !-Fluid/matrix effective compressibility
  real*8 :: etaf               !--Fluid viscosity
  real*8 :: rhof               !-Fluid density
  !real*8, dimension(7) :: paraminj                 !-Normalized injection parameters
  real*8, dimension(:),allocatable :: paraminj                 !-Normalized injection parameters
  real*8 :: D             !-Hydraulic diffusivity
end type phydro



end module mod_cst_fault_rns_p
