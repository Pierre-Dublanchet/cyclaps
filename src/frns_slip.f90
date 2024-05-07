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


function frns_slip(phi,th,sk,tbp,a,b,dc,s,cm)

use mod_cst_fault_rns
implicit none
type(pmeca) :: cm
real*8 :: phi,th,sk,tbp,a,b,dc,s,frns_slip
 
frns_slip=tbp+sk+b*s*exp(phi)*(phi+th-log(dc))/dc
frns_slip=frns_slip/(cm%r*a*s+cm%beta*exp(phi))


end function frns_slip
