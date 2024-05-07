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


function frns_slip_p(phi,th,sk,tbp,p,pp,a,b,dc,s,cm)

use mod_cst_fault_rns
implicit none
type(pmeca) :: cm
real*8 :: phi,th,sk,tbp,p,pp,a,b,dc,s,frns_slip_p
 
!--For rate-and-state--!
frns_slip_p=tbp+sk+pp*(cm%mu0-cm%alpha+cm%r*a*phi+b*th-b*log(dc))
frns_slip_p=frns_slip_p+b*(s-p)*exp(phi)*(phi+th-log(dc))/dc
frns_slip_p=frns_slip_p/(cm%r*a*(s-p)+cm%beta*exp(phi))


end function frns_slip_p
