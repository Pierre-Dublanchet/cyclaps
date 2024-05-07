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


function grns_age_p(phi,th,p,pp,b,dc,s,cm)

use mod_cst_fault_rns
implicit none
type(pmeca) :: cm
real*8 :: phi,th,p,pp,b,dc,s,grns_age_p
 
grns_age_p=cm%alpha*pp/(b*(s-p)) 
grns_age_p=grns_age_p+exp(-th)-exp(phi)/dc

end function grns_age_p
