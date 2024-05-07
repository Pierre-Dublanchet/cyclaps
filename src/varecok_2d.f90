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


subroutine varecok_2d(cr,iter,iterex,vmax,vmaxex,t,tex,dtrec,dtrecp,dtrecloc,recok,reclocok,recpok)

use mod_cst_fault_rns
implicit none
integer :: iter,recok,recpok,reclocok
integer, dimension(3) :: iterex
real*8 :: t,vmax,dtrec,dtrecp,dtrecloc
real*8, dimension(3) :: tex,vmaxex
type(prec) :: cr

if (cr%qmoy .eq. 1) then
     if (iter .eq. 1) then
        recok=1
        vmaxex(1)=vmax
     end if
     if (cr%vfrec(1) .gt. 0.0) then
        !if ((vmax .ge. vmaxex(1)*cr%vfrec(1)) .or. (vmax .le. vmaxex(1)/cr%vfrec(1)) .or. (dtrec .ge. cr%dtfrec(1))) then
        if ((vmax .ge. vmaxex(1)*cr%vfrec(1)) .or. (vmax .le. vmaxex(1)/cr%vfrec(1))) then
           recok=1
           vmaxex(1)=vmax
        end if
     end if
     if (cr%tfrec(1) .gt. 0.0) then
        if (t .ge. tex(1)*cr%tfrec(1)) then
           recok=1
           tex(1)=t
        end if
     end if
     if (cr%dtfrec(1) .gt. 0.0) then
        if (dtrec .ge. cr%dtfrec(1)) then
           recok=1          
        end if
     end if
     if (cr%nitrec(1) .gt. 0.0) then
        if (iter-iterex(1) .ge. cr%nitrec(1)) then
           recok=1
           iterex(1)=iter
        end if
     end if
  end if

  if (cr%profil .eq. 1) then
     if (iter .eq. 1) then
        recpok=1
        vmaxex(2)=vmax
     end if
     if (cr%vfrec(2) .gt. 0.0) then
        !if ((vmax .ge. vmaxex(2)*cr%vfrec(2)) .or. (vmax .le. vmaxex(2)/cr%vfrec(2)) .or. (dtrecp .ge. cr%dtfrec(2))) then
        if ((vmax .ge. vmaxex(2)*cr%vfrec(2)) .or. (vmax .le. vmaxex(2)/cr%vfrec(2))) then
           recpok=1
           vmaxex(2)=vmax
        end if
     end if
     if (cr%tfrec(2) .gt. 0.0) then
        if (t .ge. tex(2)*cr%tfrec(2)) then
           recpok=1
           tex(2)=t
        end if
     end if
     if ((cr%dtfrec(2) .gt. 0.0))then
        if (dtrecp .ge. cr%dtfrec(2)) then
           recpok=1
        end if
     end if
     if (cr%nitrec(2) .gt. 0.0) then
        if (iter-iterex(2) .ge. cr%nitrec(2)) then
           recpok=1
           iterex(2)=iter
        end if
     end if
  end if

  if (cr%qloc .eq. 1) then
     if (iter .eq. 1) then
        reclocok=1
        vmaxex(3)=vmax
     end if
     if (cr%vfrec(3) .gt. 0.0) then
        !if ((vmax .ge. vmaxex(3)*cr%vfrec(3)) .or. (vmax .le. vmaxex(3)/cr%vfrec(3)) .or. (dtrecloc .ge. cr%dtfrec(3))) then
        if ((vmax .ge. vmaxex(3)*cr%vfrec(3)) .or. (vmax .le. vmaxex(3)/cr%vfrec(3))) then
           reclocok=1
           vmaxex(3)=vmax
        end if
     end if
     if (cr%tfrec(3) .gt. 0.0) then
        if (t .ge. tex(3)*cr%tfrec(3)) then
           reclocok=1
           tex(3)=t
        end if
     end if
     if (cr%dtfrec(3) .gt. 0.0) then
        if (dtrecloc .ge. cr%dtfrec(3)) then
           reclocok=1          
        end if
     end if
     if (cr%nitrec(3) .gt. 0.0) then
        if (iter-iterex(3) .ge. cr%nitrec(3)) then
           reclocok=1
           iterex(3)=iter
        end if
     end if
  end if




end subroutine varecok_2d
