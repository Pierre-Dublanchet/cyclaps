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


function frns_age_p(phi,th,sk,tbp,p,pp,a,b,dc,s,cm)

use mod_cst_fault_rns
implicit none
type(pmeca) :: cm
real*8 :: phi,th,sk,tbp,p,pp,a,b,dc,s,frns_age_p
 
frns_age_p=tbp+sk+pp*(cm%mu0-cm%alpha+cm%r*a*phi+b*th-b*log(dc))
frns_age_p=frns_age_p-b*(s-p)*(exp(-th)-exp(phi)/dc)
frns_age_p=frns_age_p/(cm%r*a*(s-p)+cm%beta*exp(phi))

end function frns_age_p
