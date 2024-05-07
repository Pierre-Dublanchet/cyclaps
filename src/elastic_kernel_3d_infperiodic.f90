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


subroutine elastic_kernel_3d_infperiodic(kappa_proc,cc,cm)

use mod_cst_fault_rns
implicit none
type(pmeca) :: cm
type(pcomput) :: cc
integer :: i,j
real*8 :: kappa,k,l,m,t,alpha
real*8, dimension(cc%ny,cc%nvalpr/cc%ny) :: kappa_proc

alpha=(3-4*cm%nu)

do i=1,cc%nvalpr/cc%ny
       if (i+rang*cc%nvalpr/cc%ny .le. cc%nx/2) then
        k=(i+rang*cc%nvalpr/cc%ny-1)*2*pi/(cc%nx*cc%dx)
       else
        k=(cc%nx+1-i-rang*cc%nvalpr/cc%ny)*2*pi/(cc%nx*cc%dx)
       end if
       do j=1,cc%ny
              if (j .le. cc%ny/2) then
            l=(j-1)*2*pi/(cc%ny*cc%dy)
              else
            l=(cc%ny+1-j)*2*pi/(cc%ny*cc%dy)
      	      end if
         m=sqrt(k**2+l**2)
         t=(1.0-exp(-2*cm%H*m))/(1.0+exp(-2*cm%H*m))

               if ((k .eq. 0.0) .and. (l .eq. 0.0)) then
            kappa=0.0
               else
            kappa=(4*(cm%H**2)*(m**2)+(alpha-1)**2)*(k**2)*(t**3)+m*cm%H*(alpha+1)*(l**2)*(t**2-1)&
                                     -(4*(cm%H**2)*(k**2)*(m**2)+alpha*(alpha+1)*(m**2)+(alpha+1)*(k**2))*t
            kappa=kappa/(m*t*(alpha+1)*(cm%H*m*(t**2-1)-alpha*t))
               end if
         !kappa_proc(j,i)=-0.5*cm%Y*kappa/(2*(1+cm%nu))
         kappa_proc(j,i)=-0.5*kappa
         


       end do
end do

end subroutine elastic_kernel_3d_infperiodic
