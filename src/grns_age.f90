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


function grns_age(phi,th,dc)

use mod_cst_fault_rns
implicit none
real*8 :: phi,th,dc,grns_age
 
grns_age=exp(-th)-exp(phi)/dc


end function grns_age
