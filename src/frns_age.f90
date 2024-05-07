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


function frns_age(phi,th,sk,tbp,a,b,dc,s,cm)

use mod_cst_fault_rns
implicit none
type(pmeca) :: cm
real*8 :: phi,th,sk,tbp,a,b,dc,s,frns_age
 
frns_age=tbp+sk-b*s*(exp(-th)-exp(phi)/dc)
frns_age=frns_age/(cm%r*a*s+cm%beta*exp(phi))

end function frns_age
