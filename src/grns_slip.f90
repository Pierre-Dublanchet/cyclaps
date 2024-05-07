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


function grns_slip(phi,th,dc)

use mod_cst_fault_rns
implicit none
real*8 :: phi,th,dc,grns_slip
 
 grns_slip=-exp(phi)*(phi+th-log(dc))/dc

end function grns_slip
