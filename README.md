---
author:
- Pierre Dublanchet
bibliography:
- biblio.bib
date: May 2024
title: |
  CYCLAPS\
  earthquake CYCLe simulator for Asperities under Poroelastic Stressing
---

# Introduction

CYCLAPS (earthquake CYCLe simulator for Asperities under Poroelastic
Stressing) is an earthquake simulator based on rate-and-state friction
and quasi-dynamic elasticity. It allows to compute slip, slip rate and
shear stress history on 1d or 2d faults, embedded in 2d or 3d elastic
media respectively, and undergoing slow tectonic loading. CYCLAPS is an
asperity model specifically designed to simulate slip on faults with
frictional and (normal and shear) stress heterogeneity. Arbitrary
external normal and shear stress perturbation can be implemented.
Poro-elastic coupling can be handled, with linear or non-lienar fluid
diffusion restricted to the fault (impermeable bulk). CYCLAPS is
parrallelized.

# List of fault models

1.  **fault_2d_infperiodic_aging** : 1d (mode II or III) fault between
    2d elastic slabs of finite thickness. Spectral boundary integral
    approach (replication along strike). Aging law.

2.  **fault_2d_infperiodic_aging_pnl** : same with non linear pore
    pressure diffusion.

3.  **fault_2d_infperiodic_aging_press** : same with imposed pore
    pressure history.

4.  **fault_2d_infperiodic_slip** : 1d (mode II or III) fault between 2d
    elastic slabs of finite thickness. Spectral boundary integral
    approach (replication along strike). Slip law.

5.  **fault_2d_infperiodic_slip_pnl** : same with non linear pore
    pressure diffusion.

6.  **fault_2d_infperiodic_slip_press** : same with imposed pore
    pressure history.

7.  **fault_2d_freesurface_aging** : 1d strike-slip fault in a semi
    infinite elastic half space with free surface, aging law.

8.  **fault_2d_freesurface_aging_pnl** : same with non linear pore
    pressure diffusion.

9.  **fault_2d_freesurface_aging_press** : same with imposed pore
    pressure history.

10. **fault_2d_freesurface_slip** : 1d strike-slip fault in a semi
    infinite elastic half space with free surface, slip law.

11. **fault_2d_freesurface_slip_pnl** : same with non linear pore
    pressure diffusion.

12. **fault_2d_freesurface_slip_press** : same with imposed pore
    pressure history.

13. **fault_2d_cr_aging** : 1d (mode II or III) fault between 2d semi
    infinite elastic half spaces. [@Cochard1997] spectral boundary
    integral approach (no replication along strike). Aging law.

14. **fault_2d_cr_aging_pnl** : same with non linear pore pressure
    diffusion.

15. **fault_2d_cr_aging_press** : same with imposed pore pressure
    history.

16. **fault_2d_cr_slip** : 1d (mode II or III) fault between 2d semi
    infinite elastic half spaces. [@Cochard1997] spectral boundary
    integral approach (no replication along strike). Slip law.

17. **fault_2d_cr_slip_pnl** : same with non linear pore pressure
    diffusion.

18. **fault_2d_cr_slip_press** : same with imposed pore pressure
    history.

19. **fault_2d_cr_regage_pnl** : same with non linear pore pressure
    diffusion and regularized rate-and-state, aging law.

20. **fault_2d_cr_regslip_pnl** : same with non linear pore pressure
    diffusion and regularized rate-and-state, slip law.

21. **fault_3d_infperiodic_aging** : 2d fault between 3d elastic slabs
    of finite thickness. Spectral boundary integral approach
    (replication along depth and strike). Aging law.

22. **fault_3d_infperiodic_aging_pnl** : same with non linear pore
    pressure diffusion.

23. **fault_3d_infperiodic_aging_press** : same with imposed pore
    pressure history.

24. **fault_3d_infperiodic_slip** : 2d fault between 3d elastic slabs of
    finite thickness. Spectral boundary integral approach (replication
    along depth and strike). Slip law.

25. **fault_3d_infperiodic_slip_pnl** : same with non linear pore
    pressure diffusion.

26. **fault_3d_infperiodic_slip_press** : same with imposed pore
    pressure history.

27. **fault_3d_infperiodic_regage** : 2d fault between 3d elastic slabs
    of finite thickness. Spectral boundary integral approach
    (replication along depth and strike). Regularized rate-and-state,
    aging law.

28. **fault_3d_infperiodic_regslip** : same with regularized
    rate-and-state, slip law.

# Model description

## Geometry and constitutive equations

The fault geometry considered is a planar (1d or 2d) fault embedded in
an elastic medium (2d or 3d, finite or infinite), as depicted in Figure
[1](#fig1){reference-type="ref" reference="fig1"}. The fault is the
$z=0$ plane. For 2d configurations, the variables only depend on $x$
coordinate, not on $y$. The fault simulated is a frictional interface,
loaded by a (possibly heterogeneous) lithostatic normal stress
$\sigma(x,y)$ (or $\sigma(x)$ for 2d). Fluid diffusion inside the fault
leads to a pore pressure $p(x,y,t)$, so that the effective normal stress
is $\sigma_e=\sigma-p$. A constant slip rate $v_p$ is imposed either at
a distance $\pm H$ from the fault (fault_xxx_infperiodic_xxx), or within
the plane $z=0$, around the frictional domain (fault_xxx_cr_xxx,
fault_xxx_freesurface_xxx), forcing shear slip $\delta$ on the fault in
the $x$ direction.

In the special case of freesurface configuration
(fault_2d_freesurface_xxx) the model is a vertical strike slip fault
(dip angle $\beta=90^{\circ}$), with slip occuring in the $y$ direction,
$x>0$ is the depth. The freesurface is situated at $x=-L/2$.

Slip is resisted on the fault by rate-and-state friction
[@Dieterich1979; @Marone1998], with a possible normal stress dependence
on the state variable as formulated by [@Linker1992]. For standard
rate-and-state friction, the friction coefficient $f$ writes:

$$f=f_0 + a\ln{\frac{v}{v^*}}+b\ln{\frac{v^*\theta}{d_c}},$$ For
regularized rate-and-state friction (fault_xxx_regage_xxx and
fault_xxx_regslip_xxx), $f$ is given by:
$$f=a \sinh^{-1}\left[\frac{v}{2v^*}\exp{\left(\frac{f_0+b\ln{v^*\theta/d_c}}{a}\right)}\right]$$
where $f_0$, $a$, $b$ and $d_c$ are the rate-and-state parameters, $v^*$
a reference slip rate, $v$ the slip rate, and $\theta$ the state
variable. $a$, $b$ and $d_c$ can be heterogeneous along the fault and
depend on both $x$ and $y$. $f_0$ and $v^*$ are constant.

The distribution of $a$, $b$, $d_c$, $\sigma$, of the initial slip
$\delta$, slip rate $v$ and state variable $\theta$ can either be
specified in the form of circular patches with uniform values
(asperities). Traditionally, velocity weakening (VW $a-b<0$) patches are
distributed on a velocity strengthening background (VS $a-b>0$), as
illustrated in Figure [1](#fig1){reference-type="ref" reference="fig1"}
(right). The other option is to specify these values at each point of
the computational grid, and to provide the matrices as input (see input
section for further details).

The following state evolution laws are used:

-    aging law (fault_xxx_aging_xxx):
    $$\frac{d\theta}{dt} = 1- \frac{v\theta}{d_c}-\alpha\frac{\theta}{b\sigma_e}\frac{d\sigma_e}{dt},$$

-   slip law (fault_xxx_slip_xxx):
    $$\frac{d\theta}{dt} = -\frac{v\theta}{d_c}\ln{\frac{v\theta}{d_c}}-\alpha\frac{\theta}{b\sigma_e}\frac{d\sigma_e}{dt},$$

where $\alpha$ is a constant coefficient, and $\sigma_e$ is the
effective normal stress $\sigma-p$.

The pore pressure evolution within the fault $p$ can either be imposed
(fault_xxx_press). In this case, the pore pressure and pore pressure
rate histories have to be coded in the routines press_1d(2d).f90 and
pressrate_1d(2d).f90 respectively.

The pore pressure can also be solved numerically (fault_xxx_pnl). In
this case, $p$ obeys the following non-linear diffusion equation:

$$\frac{\partial p}{\partial t} = \nabla \left(\frac{k}{\eta_f \phi C^*} \nabla p\right) + s(t) \delta_D(x-x_i,y-y_i),$$
where $k$ is the fault permeability, $\eta_f$ the fluid viscosity,
$\phi$ the fault porosity and $C^*$ the effective (fluid+pore space)
compressibility, $x_i$, $y_i$ the coordinates of a ponctual fluid
source, with time evolution given by $s(t)$. $\delta_D$ is the dirac
delta function. The permeability can be time, space, stress, slip, or
slip-rate dependent. The routine diffu.f90 can be used to define the
resulting hydraulic diffusivity law. The source term $s(t)$ has to be
defined in the routines pressrate_nl1d(2d).f90.

The fault slip evolution is computed assuming a quasi-static balance of
the form:

$$f\sigma_e = \tau_0 + \kappa * \delta -\frac{\mu}{2c_s} v,$$ where
$\tau_0(x,y,t)$ incorporates the initial shear stress (imposed by the
initial slip rate and state variable), and a possible external shear
stressing, arising either from the boundary conditions, or from another
mechanism. An external transient shear stress perturbation can be
imposed, but has to be coded in routines tbp_xxx.f90. $\kappa(x,y)$ is
the stress interaction kernel, accounting for the stress redistribution
associated with slip along the fault. $\kappa(x,y)$ depends on geometry
of the fault and the boundary conditions. The convolution
$\kappa * \delta$ is either computed using a spectral approach
(infperiodic), a spectral approach avoiding the replication of the fault
(cr, following [@Cochard1997]), or in the space domain (freesurface).
$\mu$ is the shear modulus of the elastic medium, and $c_s$ the shear
wave speed of the elastic medium.

<figure id="fig1">
<div class="center">
<img src="./img/aspmodel.png" />
</div>
<figcaption>Fault model geometry (left) and asperity structure of the
fault (right). VW: velocity weakening, VS: velocity
strengthening.</figcaption>
</figure>

## Discretization

The constitutive equations are solved using finite differences. For 3d
configurations (2d fault), the fault plane is discretized in
$n_x \times n_y$ rectangular cells of size $\Delta x \times \Delta y$
ordered columnwise (Figure [2](#fig2){reference-type="ref"
reference="fig2"}). For 2d configurations (1d fault), $n_y=1$.

<figure id="fig2">
<div class="center">
<img src="./img/xymodel.png" style="width:50.0%" />
</div>
<figcaption>Fault discretization for 2d faults (3d configurations). The
colorscale indicates the numbering of the computational cells. 1d faults
(2d configuration) follow the same convention, with <span
class="math inline"><em>n</em><sub><em>y</em></sub> = 1</span>.</figcaption>
</figure>

# Third party source code and librairies

## Third party source code

The source code of CYCLAPS includes external source codes (see licenses
notices in licenses):

::: center
  **Name**            **License**    **URL**                                                                                 **Copyright**
  ------------------- -------------- --------------------------------------------------------------------------------------- --------------------------------------------------------------
  FFTE                see licenses   [www.ffte.jp](www.ffte.jp){.uri}                                                        Copyright (c), 2000-2004, 2008-2014, 2020, Daisuke Takahashi
  Special Functions   see licenses   <https://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.html>   Copyright (c) 1996 Shanjie Zhang and Jianming Jin
  CFGIO               MIT            [ https://github.com/pkgpl/cfgio]( https://github.com/pkgpl/cfgio){.uri}                Copyright (c) 2017 Wansoo Ha
:::

See [@Zhang1996] for details about the computation of special functions.

## Third party librairies

The source code of CYCLAPS uses external librairies (see licenses
notices in licenses):

::: center
  **Name**   **License**    **URL**                                           **Copyright**
  ---------- -------------- ------------------------------------------------- ---------------------------------------------------------------------------------
  MPICH      see licenses   <https://www.mpich.org/>                          Copyright (c) 1998--2024, Argonne National Laboratory
  NetCDF     see licenses   <https://www.unidata.ucar.edu/software/netcdf/>   Copyright (c) 1993-2014 University Corporation for Atmospheric Research/Unidata
:::

# Citation

If you use this software, please cite it as:

-   [@Dublanchet2018] for 2d versions without fluid injection

-   [@Dublanchet2019fluid] for 2d version with fluid injection,

-   [@Dublanchet2019] for 3d version.

# Install, compile, execute

## Requirements

These instructions are valid for MacOS and LINUX.

-   MPICH is needed. Details about the installation of MPICH are
    provided here: <https://www.mpich.org/downloads/> .

-   NETCDF library (netcdf-fortran package for MacOS, or libnetcdff-dev
    package for LINUX) is necessary.

## Install

Download and unzip CYCLAPS.zip. In the folder CYCLAPS, you will find:

-   src : folder containing source files for CYCLAPS

-   init: folder containing netcdf initiation files

-   results: folder where model outputs (netcdf files) are stored

-   makefile

## Compile

Go to the root of CYCLAPS folder. Indicate in the makefile the path to
the netcdf and netcdf-fortran librairies (variables NETCDF and NETCDFF),
or create a symbolic link. Then, to compile fault model fault_xxx type:

       >> make exec=fault_xxx

This will create an executable `fault_xxx`

## Execute

Go to the root of CYCLAPS folder. To execute fault model fault_xxx:

       >> mpiexec -n <np> fault_xxx

where `<np>` is the number of processes to use. The number of
computational cells in $x$ and $y$ directions $n_x$ and $n_y$ have to be
multiples of np.

# Input data

## Parameters file

The simulation parameters have to be written in the file
parametres_fault_rns.cfg. parametres_fault_rns.cfg is a key-value
parameter file. More information is provided in Table
[\[tab1\]](#tab1){reference-type="ref" reference="tab1"}. The frictional
properties ($a$, $b$, $d_c$), the lithostatic normal stress $\sigma$,
the initial slip rate and state variable can be either specified in this
file, in the form of circular asperities with uniform properties, or
specified in init files (see section
[7.2](#sec:init){reference-type="ref" reference="sec:init"}).

:::: center
[]{#tab1 label="tab1"}

::: supertabular
\|p0.2\|p0.2\|p0.15\|p0.35\| simlab & label of fault model used (xxx in
fault_xxx) & n.u. & n.c.\
\
slip_mode & 2 for in-plane (mode II), 3 for anti-plane (mode III) & n.u.
& Only for 1d fault in 2d medium\
young_mod & Young's modulus $E$ & Pa & n.c.\
poisson_ratio & Poisson's ratio $\nu$ & n.u. & n.c.\
rho_r & rock density $\rho$ & kg.m$^{-3}$ & n.c.\
slab_thickness & thickness of elastic slabs in contact $H$ & m & only
used for infperiodic configuration (2d or 3d)\
ext_shear_stressing & vector of shear stressing parameters & (see below)
& only if external shear stress perturbation (used in tbp_xxx.f90
functions)\
ext_shear_stressing(1) & constant shear stressing & Pa.s$^{-1}$ & n.c.\
ext_shear_stressing(2) & amplitude of shear stress perturbation & Pa &
n.c.\
ext_shear_stressing(3) & radius of stressed region & m & n.c.\
ext_shear_stressing(4) & duration of transient stressing & s & n.c.\
ext_shear_stressing(5,6) & $x,y$ coordinates of the center of stressed
region & m & n.c.\
ref_fric_coeff & reference friction coefficient $f_0$ & n.u. & n.c.\
dl_coeff & Linker-Dieterich coefficient $\alpha$ (state dependance on
normal stress) & n.u. & n.c.\
vplate & tectonic plate rate $v_p$ & m.s$^{-1}$ & n.c.\
vstar & reference slip rate $v^*$ & m.s$^{-1}$ & n.c.\
vsis & radiative slip rate (for earthquake detection) $v_{sis}$ &
m.s$^{-1}$ & n.c.\
nasp & number of circular VW asperities & n.u. & Only if meth_init=0.
nasp should be smaller than 20\
a & direct effect rate-and-state parameter $a$ (1: VS region, then 1
value/VW asperity) & n.u. & nasp+1 values should be provided\
b & state effect rate-and-state parameter $b$ (1: VS region, then 1
value/VW asperity) & n.u. & nasp+1 values should be provided\
dcasp & critical slip $d_c$ (1: VS region, then 1 value/VW asperity) & m
& nasp+1 values should be provided\
xcasp & $x$ coordinate of asperity center (1: VS region, then 1 value/VW
asperity) & m & nasp+1 values should be provided\
ycasp & $y$ coordinate of asperity center (1: VS region, then 1 value/VW
asperity) & m & nasp+1 values should be provided\
Rasp & radius of asperity (1: VS region, then 1 value/VW asperity) & m &
nasp+1 values should be provided\
sigasp & normal stress on asperity $\sigma$ (1: VS region, then 1
value/VW asperity) & Pa & nasp+1 values should be provided\
viasp & initial slip rate on asperity (1: VS region, then 1 value/VW
asperity) & m.s$^{-1}$ & nasp+1 values should be provided\
thiasp & initial state variable on asperity (1: VS region, then 1
value/VW asperity) & s & nasp+1 values should be provided\
uiasp & initial slip on asperity (1: VS region, then 1 value/VW
asperity) & m & nasp+1 values should be provided\
\
permea & fault reference permeability $k$ & m$^2$ & n.c.\
porosity & fault reference porosity $\phi$ & n.u. & n.c.\
compress & effective (fluid and pore space) compressibility $c^*$ &
Pa$^{-1}$ & n.c.\
viscosity & fluid viscosity $\eta_f$ & Pa.s & n.c.\
rho_f & fluid density $\rho_f$ & kg.m$^{-3}$ & n.c.\
paraminj & vector of injection parameters & (see below) & all the
parameters are not used, depending on the injection scenario considered
(press_xxx.f90 and pressrate_xxx.f90 functions to define the injection
scenario)\
paraminj(1) & injection duration & s & n.c.\
paraminj(2) & $x$ coordinate of injection borehole & m & n.c.\
paraminj(3) & $y$ coordinate of injection borehole & m & n.c.\
paraminj(4) & Darcy velocity at injection borehole & m.s$^{-1}$ & n.c.\
paraminj(5) & Pore pressure at injection borehole & Pa & n.c.\
paraminj(6) & borehole radius & m & n.c.\
\
nx & number of computational cells $n_x$ in the $x$ direction & n.u. &
n.c.\
ny & number of computational cells $n_y$ in the $y$ direction & n.u. &
use ny$=$`<!-- -->`{=html}1 for 2d problems\
dx & computational cell size $\Delta x$ in the $x$ direction & m & n.c.\
dy & computational cell size $\Delta y$ in the $y$ direction & m & use
dy$=$`<!-- -->`{=html}0.0 for 2d problems\
niter & number of iterations & n.u. & n.c.\
nit_screen & print code progression every nit_screen iteration & n.u. &
n.c.\
meth_init & method for initial conditions (0: use the a, b, dc, sigma,
vi, thi defined in this file, 1: use the values defined in init.nc file)
& n.u. & n.c.\
paramdt & vector of time-step control parameters & (see below) & n.c.\
paramdt(1) & maximum number of Runge-Kutta iterations to adapt the time
step & n.u. & n.c.\
paramdt(2) & absolute error tolerance at each time step (error on
$\ln{v}$, $\ln{\theta}$, normalized pore pressure $p/\sigma_0$) & n.u. &
n.c.\
paramdt(3) & safety factor for time step adaptation & n.u. & should be
between $0$ and $1$\
paramdt(4) & maximum possible time step & s & n.c.\
paramdt(5) & minimum possible time step & s & n.c.\
pathinit & path to initial conditions directory (containing the init.nc
file) & n.u. & n.c.\
\
qcat & earthquake catalog production (1: yes, 0: no) & n.u. & file
earthquake_catalog.nc\
qmoy & average and extremal values history recording (1: yes, 0: no) &
n.u. & file qmoy.nc\
qprof & variable maps recording (1: yes, 0: no) & n.u. & files maps\*.nc
(one file = one time)\
qprofhvc & $x$ and $y$ variable profiles recording (1: yes, 0: no) &
n.u. & files profilsx\*.nc and profilsy\*.nc (one file = one time)\
qloc & local variables history recording (1: yes, 0: no) & n.u. & nploc
files qlloc\*.nc (1 file per point on the fault)\
vfrec & slip rate factor used to write outputs & n.u. & used if $>0$.
vfrec(1): for average extremal values, vfrec(2): for variable maps,
vfrec(3): for local variables, vfrec(4) for $x$ and $y$ profils (only
used in 3d geometry). The output is written when the max slip rate is
multiplied or divided by vfrec(i).\
tfrec & time factor used to write outputs & n.u. & same as vfrec, but
for time instead of slip rate\
dtfrec & time step between two output writings & s & used if $>0$. Write
outputs every dtfrec(i) seconds. The 4 components correspond to the same
outputs as for vfrec and tfrec.\
nitrec & number of iterations between two output writings & n.u. & same
as tfrec, but in terms of number of iterations\
pathres & path to results/output directory & n.u. & n.c.\
nploc & number of locations where variables are written & n.u. & used if
qloc$=$`<!-- -->`{=html}1\
xrloc & $x$ coordinates of locations where variable are written & m &
used if qloc$=$`<!-- -->`{=html}1. nploc values should be provided\
yrloc & $y$ coordinates of locations where variable are written & m &
used if qloc$=$`<!-- -->`{=html}1. nploc values should be provided\
hboxm & size of the boundary zone not included in moment and maximum
slip rate computation & m & only for 3d configuration with regularized
friction law (fault_3d_infperiodic_regxxx)\
:::
::::

## Complex heterogeneous fault structure and initial conditions {#sec:init}

For a more complex fault heterogeneity, or if the number of cricular
asperities is larger than 20, ($a$, $b$, $d_c$, $\sigma$, $v_i$ and
$\theta_i$, $u_i$ and $p_i$) distribution can be specified for each of
the $n_x \times n_y$ cell of the fault. The parameter meth_init has to
be set to 1 in the parametres_fault_rns.cfg file. These values should be
provided in a single netcdf file init.nc localized in the ./INIT/
directory. The init.nc file contains a dataset consisting of the 8
variables: a, b, dc, s, vi, ui, thi, vi and pi. Each variable is a
$n_xn_y \times 1$ vector, where the fault cells are ordered columnwise,
as illustrated in Figure [2](#fig2){reference-type="ref"
reference="fig2"}.

# Outputs

CYCLAPS produces five categories of outputs, along with four ways of
controlling the frequency of output writing. The outputs are written in
netcdf files located in the ./RESULTS/ directory The five categories of
outputs are:

1.  Earthquake catalogue (file earthquake_catalog.nc).

2.  Average and extremal values time histories (file qmoy.nc)

3.  Local variables time series (files qloc\*.nc)

4.  Maps of variables (files maps\*.nc)

5.  Profiles of variables along $x$ and $y$ directions (files
    profilsx\*.nc and profilsy\*.nc)

Details about these different outputs are provided in the next
subsections. The writing of outputs in the files can be controlled in
four different ways, except for the earthquake catalog. The four
possibilities are:

1.  one output is written each time the maximum slip rate is mutliplied
    or divided by a factor provided in vfrec parameter of
    parametres_fault_rns.cfg file (see table
    [\[tab1\]](#tab1){reference-type="ref" reference="tab1"})

2.  one output is written each time the absolute time since the start of
    the simulation is multiplied by a factor provided in tfrec parameter
    of parametres_fault_rns.cfg file (see table
    [\[tab1\]](#tab1){reference-type="ref" reference="tab1"})

3.  outputs are written at constant time steps. The time separating two
    output writings is specified in the dtfrec parameter of
    parametres_fault_rns.cfg file (see table
    [\[tab1\]](#tab1){reference-type="ref" reference="tab1"})

4.  outputs are written based on iterations. The number of iterations
    separating two output writings is specified in the nitrec parameter
    of parametres_fault_rns.cfg file (see table
    [\[tab1\]](#tab1){reference-type="ref" reference="tab1"})

The four components of the vectors vfrec, tfrec, dtfrec and nitrec in
the parametres_fault_rns.cfg file correspond to the 4 last categories of
outputs (not the to the earthquake catalog). Details are provided in
table [\[tab1\]](#tab1){reference-type="ref" reference="tab1"}.

## Earthquake catalog

The file earthquake_catalog.nc contains a dataset consisting of 11
variables listed below:

1.  onset time: onset time of earthquakes (in s)

2.  onset time delai: time delai since the last earthquake on the fault
    (in s)

3.  event duration: event duration (in s)

4.  x initiation: $x$ coordinate of the first point involved in the
    earthquake rupture (in m)

5.  y initiation: $x$ coordinate of the first point involved in the
    earthquake rupture (in m)

6.  x barycenter: $x$ coordinate of the barycenter of the earthquake
    rupture (in m)

7.  y barycenter: $x$ coordinate of the barycenter of the earthquake
    rupture (in m)

8.  number of elements: number of computational cells involved in the
    earthquake rupture (n.u.)

9.  coseismic moment: coseismic moment liberated by the earthquake (N.m)

10. coseismic stress drop: coseismic stress drop associated with the
    earthquake, averaged over the rupture area (Pa)

11. coseismic slip: coseismic slip associated with the earthquake,
    averaged over the rupture area (m)

The models fault_3d_infperiodic_regage and fault_3d_infperiodic_regslip
produce two additionnal variables in the earthquake catalog file:

1.  shear stress init: average shear stress in the rupture area before
    the onset of the earthquake (in Pa)

2.  shear stress final: average shear stress in the rupture area right
    after the earthquake (in Pa)

The difference between shear stress final and shear stress init is the
coseismic stress drop.

Each value contained in one of these variables corresponds to one
earthquake. An earthquake occurs when the maximum slip rate exceeds the
radiative threshold vsis defined in the parametres_fault_rns.cgf file.
The earthquake rupture corresponds to all the fault elements having
experienced a slip rate larger than vsis until the maximum slip rate
decreases below vsis.

## Average and extremal values time histories

The file qmoy.nc contains a dataset consisting of 7 variables capturing
different time series. The 7 variables are:

1.  time: absolute time of each variable writing (in s)

2.  time delai: time delai since the last output writing (in s)

3.  mean slip rate: time series of the spatial average of slip rate (in
    m.s$^{-1}$)

4.  max slip rate: time series of the maximum slip rate on the fault (in
    m.s$^{-1}$)

5.  mean slip: time series of the spatial average of slip on the fault
    (in m)

6.  mean shear stress: time series of the spatial average of shear
    stress on the fault (in Pa)

7.  mean state: time series of the spatial average of the state variable
    on the fault (in s)

The models fault_3d_infperiodic_regage and fault_3d_infperiodic_regslip
produce one additionnal variable in qmoy.nc file:

1.  mean slip rate vw: time series of the spatial average of the slip
    rate within the velocity weakening regions of the fault (in
    m.s$^{-1}$)

Note also that only fault_3d_infperiodic_regage and
fault_3d_infperiodic_regslip consider the parameter hboxm in the
computation of mean slip rate (a zone of width hboxm is excluded). See
table [\[tab1\]](#tab1){reference-type="ref" reference="tab1"} for
details.

## Local variables time series

The time series of slip rate, state variable, shear stress, slip and
eventually pore pressure and Darcy velocity can be written for nploc
locations with coordinates specified in xrloc and yrloc of
parametres_fault_rns.cfg (table [\[tab1\]](#tab1){reference-type="ref"
reference="tab1"}). The times series are written in nploc files named
qloc\*.nc (time series of location i is stored in qloci.nc). The files
qloc\*.nc contain datasets consisting of 6 variables. The variables are:

1.  time: absolute time (in s)

2.  time delai: time delai since the last output writing (in s)

3.  slip: total slip time series (in m)

4.  slip rate: slip rate time series (in m.s$^{-1}$)

5.  state: state variable time series (in s)

6.  shear stress: shear stress time series (in Pa)

For models considering fluid injection (fault_xxx_press or
fault_xxx_pnl), qloc\*.nc contain an additional variable:

1.  pore pressure: pore pressure time series (in Pa).

Finally, models fault_2d_cr_regage_pnl and fault_2d_cr_regslip_pnl write
an additional variable:

1.  darcy vel: Darcy velocity time series (in m.s$^{-1}$).

## Maps of variables

The slip rate, state, shear stress, slip, pore pressure values at each
fault cell at given times can be written in the files maps\*.nc. The
file mapsi.nc contains a dataset consisting of 6 variables,
characterizing the fault at time $t_i$. The vriables are:

1.  time: absolute time $t_i$ (in s)

2.  time delai: time delai since the last output writing (in s)

3.  slip: slip distribution (in m)

4.  slip rate: slip rate distribution (in m.s$^{-1}$)

5.  state: state variable distribution (in s)

6.  shear stress: shear stress distribution (in Pa)

For models considering fluid injection (fault_xxx_press or
fault_xxx_pnl), maps\*.nc contain an additional variable:

1.  pore pressure: pore pressure time series (in Pa).

All the variables, except time and time delai are $n_xn_y \times 1$
vectors containing the values of slip, slip rate, shear stress, state,
and pore pressure at each fault cell. The values are stored in the
variables columnwise (Figure [2](#fig2){reference-type="ref"
reference="fig2"}).

## Profiles of variables along $x$ and $y$ directions

For 3d configuration, writing all the values of slip rate, state, shear
stress, slip and pore pressure at a given time in a maps\*.nc file can
be consuming in time and memory. One can instead write only these values
along two profiles: the first along $x$ centered at $y=0$, the second
along $y$ centered at $x=0$. The results at time $t_i$ are stored in two
files: profilsxi.nc and profilsyi.nc. The files profilsx\*.nc and
profilsy\*.nc contain a dataset consisting in the same variables as in
maps\*.nc (see previous section).

Note that this category of output is only relevant for 3d configurations
(fault_3d_xxx), and is not accounted for in 2d configurations
(fault_2d_xxx).

# License

CYCLAPS is distributed under CeCILL-B License (see licenses).


