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


function frns_regage(phi,th,sk,tbp,a,b,dc,s,cm)

use mod_cst_fault_rns
implicit none
type(pmeca) :: cm
real*8 :: phi,th,sk,tbp,a,b,dc,s,frns_regage,F,fric

!F=cm%mu0+cm%r*a*phi+b*th-b*log(dc)
!fric=cm%r*a*log(0.5*exp(F/(cm%r*a))+sqrt(1+0.25*exp(2*F/(cm%r*a))))
 
!frns_regage=(tbp+sk)*sqrt(1+4*exp(-2*F/(cm%r*a)))-b*s*(exp(-th)-exp(phi)/dc)
!frns_regage=frns_regage/(cm%r*a*s+cm%beta*exp(phi)*sqrt(1+4*exp(-2*F/(cm%r*a))))


F=cm%r*a*phi+b*th-b*log(dc)
fric=cm%r*a*(log(0.5*exp(F/(cm%r*a))+sqrt(exp(-2*cm%mu0)+0.25*exp(2*F/(cm%r*a))))+cm%mu0)
 
frns_regage=(tbp+sk)*sqrt(1+4*exp(-2*(F+cm%mu0)/(cm%r*a)))-b*s*(exp(-th)-exp(phi)/dc)
frns_regage=frns_regage/(cm%r*a*s+cm%beta*exp(phi)*sqrt(1+4*exp(-2*(F+cm%mu0)/(cm%r*a))))

end function frns_regage
