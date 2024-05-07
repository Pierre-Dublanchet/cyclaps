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


module mod_cst_fault_rns

use mpi
implicit none
include 'debug.h'
include 'stat.h'
integer :: rang,code,nprocs
integer, dimension(MPI_STATUS_SIZE) :: statut
integer :: statusnetcdf
real*8, parameter :: pi=3.141592653589793
complex*16, parameter :: minus_one=-1.0
complex*16 :: z

type pmeca
  character(len=30) :: simlab
  integer :: slipmode              !-Slip mode for 2D problems 	
  real*8 :: Y                  !-Young's Modulus 
  real*8 :: nu                 !-Poisson's Ratio
  real*8 :: rho               !-Density
  real*8 :: H                  !-Slab Thickness 
  real*8, dimension(:),allocatable :: paramdsext  !--External stressing parameters for variable external stressing (1): homogeneous constant stressing rate (2): amplitude of shear stress perturbation (3): radius of stressed region (4): duration of transient stressing, (5): x (6): y coordinates of the center of the stressed region
  real*8 :: mu0                !-Constant Friction Coefficient
  real*8 :: alpha              !-Dieterich-Linker Coefficient alpha
  real*8 :: vb,vb0                 !-Tectonic plate motion rate and reference speed
  real*8 :: vsis              !-Radiative slip rate (eqk detection)
  real*8 :: mu              !-Shear modulus
  real*8 :: r                  !-Ratio of reference rate-and-state parameters r=a0/b0
  real*8 :: beta             !-Normalized radiation damping
  real*8 :: a0				!-Reference value for a parameter
  real*8 :: b0				!-Reference value for b parameter 
  real*8 :: dc0		     !-Reference value for dc
  real*8 :: sig0	          !-Reference value for normal stress
  real*8 :: sm_factor     !-slip mode factor for 2D applications (sm=1 for mode III, sm=1-nu for mode II)
  integer :: n_asp        !-number of asperities
  real*8, dimension(:), allocatable :: a_asp !-rns direct effect parameter a (1: VS barriers, 2-21: VW asperities)
  real*8, dimension(:), allocatable :: b_asp !-rns state effect parameter b (1: VS barriers, 2-21: VW asperities)
  real*8, dimension(:), allocatable :: dc_asp !-rns critical slip dc (1: VS barriers, 2-21: VW asperities)
  real*8, dimension(:), allocatable :: xc_asp !-x coordinate of asperity center (1: VS barriers, 2-21: VW asperities)
  real*8, dimension(:), allocatable :: yc_asp !-y coordinate of asperity center (1: VS barriers, 2-21: VW asperities)
  real*8, dimension(:), allocatable :: R_asp !-Radius of asperities (1: VS barriers, 2-21: VW asperities)
  real*8, dimension(:), allocatable :: sig_asp !-Normal stress on asperities (1: VS barriers, 2-21: VW asperities)
  real*8, dimension(:), allocatable :: vi_asp !-Initial slip rate on asperities (1: VS barriers, 2-21: VW asperities)
  real*8, dimension(:), allocatable :: thi_asp !-Initial state variable on asperities (1: VS barriers, 2-21: VW asperities)
  real*8, dimension(:), allocatable :: ui_asp !-Initial slip on asperities (1: VS barriers, 2-21: VW asperities)
end type pmeca

type pcomput
  integer :: nx                             !-Number of points in x
  integer :: ny                             !-Number of points in y
  integer :: n                              !-Total number of points
  integer*8 :: n1                           !-int8 version of number of points
  integer :: nvalpr                         !-Number of points per processor
  integer :: niter                          !-Number of Iterations
  integer :: niterscreen                    !-Print Evolution Frequency (print evolution every niterscreen iteration)
  integer :: methinit                       !-Method for initial conditions definition: self defined (0), read from file (1)
  real*8 :: dx,dy                           !-Computational Cell Size
  real*8, dimension(:),allocatable :: paramdt           !-Time Step Size Control Parameters:  paramdt(1)=maximum num of RK iterations, paramdt(2)=error tolerance, paramdt(3)=safety factor, paramdt(4)=max time step, paramdt(5)=min time step
  real*8, dimension(6) :: rka       !-Parameters a of Runge Kutta iterations
  real*8, dimension(6,5) :: rkc       !-Parameters c of Runge Kutta iterations
  real*8, dimension(6) :: rkr       !-Parameters r of Runge Kutta iterations
  real*8, dimension(6) :: rkf       !-Parameters f of Runge Kutta iterations
  character(len=50) :: path_init            !-Path to initial values directory
end type pcomput

type prec
  integer :: catal,qmoy,profil,qloc,profilhvc                          !-Indicator for catalogue (catal), averages (qmoy) and profils (profils) recording
  integer :: nploc                                       !-Number of local variables to monitor
  integer :: ncqmoyid,nceqkcatid,qmoy_dimid,eqkcat_dimid
  integer, dimension(8) :: qmoy_varid
  integer, dimension(15) :: eqkcat_varid
  integer, dimension(:), allocatable :: ncqlocid,qloc_dimid
  integer, dimension(:,:), allocatable :: qloc_varid
  integer, dimension(:), allocatable :: nitrec                        !-Parameters for recording frequency of averages quantities (1) and profils (2)
  real*8, dimension(:), allocatable :: vfrec,tfrec,dtfrec             !-Parameters for recording frequency of averages quantities (1) and profils (2)
  character(len=70) :: path_res                          !-Path to results directory
  real*8, dimension(:), allocatable :: xr,yr                         !-Coordinates of loc where local variables are monitored
  real*8 :: hboxm
end type prec

end module mod_cst_fault_rns
