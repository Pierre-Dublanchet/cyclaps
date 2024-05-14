# CYCLAPS (earthquake CYCLe simulator for Asperities under Poroelastic Stressing)


## Introduction

CYCLAPS (earthquake CYCLe simulator for Asperities under Poroelastic Stressing) is an earthquake simulator based on rate-and-state friction and quasi-dynamic elasticity. It allows to compute slip, slip rate and shear stress history on 1d or 2d faults, embedded in 2d or 3d elastic media respectively, and undergoing slow tectonic loading. CYCLAPS is an asperity model specifically designed to simulate slip on faults with frictional and (normal and shear) stress heterogeneity. Arbitrary external normal and shear stress perturbation can be implemented. Poro-elastic coupling can be handled, with linear or non-lienar fluid diffusion restricted to the fault (impermeable bulk). CYCLAPS is parrallelized.

## List of fault models

1.***fault\_2d\_infperiodic\_aging*** : 1d (mode II or III) fault between 2d elastic slabs of finite thickness. Spectral boundary integral approach (replication along strike). Aging law.

2.***fault\_2d\_infperiodic\_aging\_pnl*** : same with non linear pore pressure diffusion.

3.***fault\_2d\_infperiodic\_aging\_press*** : same with imposed pore pressure history. 

4.***fault\_2d\_infperiodic\_slip*** : 1d (mode II or III) fault between 2d elastic slabs of finite thickness. Spectral boundary integral approach (replication along strike). Slip law.

5.***fault\_2d\_infperiodic\_slip\_pnl*** : same with non linear pore pressure diffusion.

6.***fault\_2d\_infperiodic\_slip\_press*** : same with imposed pore pressure history.

7.***fault\_2d\_freesurface\_aging*** : 1d strike-slip fault in a semi infinite elastic half space with free surface, aging law.

8.***fault\_2d\_freesurface\_aging\_pnl*** : same with non linear pore pressure diffusion.

9.***fault\_2d\_freesurface\_aging\_press*** : same with imposed pore pressure history.

10.***fault\_2d\_freesurface\_slip*** : 1d strike-slip fault in a semi infinite elastic half space with free surface, slip law.

11.***fault\_2d\_freesurface\_slip\_pnl*** : same with non linear pore pressure diffusion.

12.***fault\_2d\_freesurface\_slip\_press*** : same with imposed pore pressure history.

13.***fault\_2d\_cr\_aging*** : 1d (mode II or III) fault between 2d semi infinite elastic half spaces. \cite{Cochard1997} spectral boundary integral approach (no replication along strike). Aging  law.

14.***fault\_2d\_cr\_aging\_pnl*** : same with non linear pore pressure diffusion.

15.***fault\_2d\_cr\_aging\_press*** : same with imposed pore pressure history.

16.***fault\_2d\_cr\_slip*** : 1d (mode II or III) fault between 2d semi infinite elastic half spaces. \cite{Cochard1997} spectral boundary integral approach (no replication along strike). Slip  law.

17.***fault\_2d\_cr\_slip\_pnl*** : same with non linear pore pressure diffusion.

18.***fault\_2d\_cr\_slip\_press*** : same with imposed pore pressure history.

19.***fault\_2d\_cr\_regage\_pnl*** : same with non linear pore pressure diffusion and regularized rate-and-state, aging law.

20.***fault\_2d\_cr\_regslip\_pnl*** : same with non linear pore pressure diffusion and regularized rate-and-state, slip law.

21.***fault\_3d\_infperiodic\_aging*** : 2d fault between 3d elastic slabs of finite thickness. Spectral boundary integral approach (replication along depth and strike). Aging law.

22.***fault\_3d\_infperiodic\_aging\_pnl*** : same with non linear pore pressure diffusion.

23.***fault\_3d\_infperiodic\_aging\_press*** : same with imposed pore pressure history.

24.***fault\_3d\_infperiodic\_slip*** : 2d fault between 3d elastic slabs of finite thickness. Spectral boundary integral approach (replication along depth and strike). Slip law.

25.***fault\_3d\_infperiodic\_slip\_pnl*** : same with non linear pore pressure diffusion.

26.***fault\_3d\_infperiodic\_slip\_press*** : same with imposed pore pressure history.

27.***fault\_3d\_infperiodic\_regage*** : 2d fault between 3d elastic slabs of finite thickness. Spectral boundary integral approach (replication along depth and strike). Regularized rate-and-state,  aging law.

28.***fault\_3d\_infperiodic\_regslip*** : same with regularized rate-and-state, slip law

## Model description

### Geometry and constitutive equations

The fault geometry considered is a planar (1d or 2d) fault embedded in an elastic medium (2d or 3d, finite or infinite), as depicted in Figure \ref{fig1}. The fault is the $z=0$ plane. For 2d configurations, the variables only depend on $x$ coordinate, not on $y$. The fault simulated is a frictional interface, loaded by a (possibly heterogeneous) lithostatic normal stress $\sigma(x,y)$ (or $\sigma(x)$ for 2d). Fluid diffusion inside the fault leads to a pore pressure $p(x,y,t)$, so that the effective normal stress is $\sigma_e=\sigma-p$.  A constant slip rate $v_p$ is imposed either at a distance $\pm H$ from the fault (fault\_xxx\_infperiodic\_xxx), or within the plane $z=0$, around the frictional domain (fault\_xxx\_cr\_xxx, fault\_xxx\_freesurface\_xxx), forcing shear slip $\delta$ on the fault in the $x$ direction.

In the special case of freesurface configuration (fault\_2d\_freesurface\_xxx) the model is a vertical strike slip fault (dip angle $\beta=90^{\circ}$), with slip occuring in the $y$ direction, $x>0$ is the depth. The freesurface is situated at $x=-L/2$.

Slip is resisted on the fault by rate-and-state friction \cite{Dieterich1979,Marone1998}, with a possible normal stress dependence on the state variable as formulated by \cite{Linker1992}. For standard rate-and-state friction, the friction coefficient $f$ writes:

```math
f=f_0 + a\ln{\frac{v}{v^*}}+b\ln{\frac{v^*\theta}{d_c}},
```

 For regularized rate-and-state friction (fault\_xxx\_regage\_xxx and fault\_xxx\_regslip\_xxx), $f$ is given by:

```math
f=a \sinh^{-1}\left[\frac{v}{2v^*}\exp{\left(\frac{f_0+b\ln{v^*\theta/d_c}}{a}\right)}\right]
```

 where $f_0$, $a$, $b$ and $d_c$ are the rate-and-state parameters, $v^\*$ a reference slip rate, $v$ the slip rate, and $\theta$ the state variable. $a$, $b$ and $d_c$ can be heterogeneous along the fault and depend on both $x$ and $y$. $f_0$ and $v^\*$ are constant.
 
 The distribution of $a$, $b$, $d_c$, $\sigma$, of the initial slip $\delta$, slip rate $v$ and state variable $\theta$ can either be specified in the form of circular patches with uniform values (asperities). Traditionally, velocity weakening (VW $a-b<0$) patches are distributed on a velocity strengthening background (VS $a-b>0$), as illustrated in Figure \ref{fig1} (right). The other option is to specify these values at each point of the computational grid, and to provide the matrices as input (see input section for further details).


The following state evolution laws are used:

-aging law (fault\_xxx\_aging\_xxx):

```math
\frac{d\theta}{dt} = 1- \frac{v\theta}{d_c}-\alpha\frac{\theta}{b\sigma_e}\frac{d\sigma_e}{dt},
```

-slip law (fault\_xxx\_slip\_xxx):
```math
\frac{d\theta}{dt} = -\frac{v\theta}{d_c}\ln{\frac{v\theta}{d_c}}-\alpha\frac{\theta}{b\sigma_e}\frac{d\sigma_e}{dt},
```

where $\alpha$ is a constant coefficient, and $\sigma_e$ is the effective normal stress $\sigma-p$.

The pore pressure evolution within the fault $p$ can either be imposed (fault\_xxx\_press). In this case, the pore pressure and pore pressure rate histories have to be coded in the routines press\_1d(2d).f90 and pressrate\_1d(2d).f90 respectively.

The pore pressure can also be solved numerically (fault\_xxx\_pnl). In this case, $p$ obeys the following non-linear diffusion equation:

```math
\frac{\partial p}{\partial t} = \nabla \left(\frac{k}{\eta_f \phi C^*} \nabla p\right) + s(t) \delta_D(x-x_i,y-y_i),
```


where $k$ is the fault permeability, $\eta_f$ the fluid viscosity, $\phi$ the fault porosity and $C^*$ the effective (fluid+pore space) compressibility, $x_i$, $y_i$ the coordinates of a ponctual fluid source, with time evolution given by $s(t)$. $\delta_D$ is the dirac delta function. The permeability can be time, space, stress, slip, or slip-rate dependent. The routine diffu.f90 can be used to define the resulting hydraulic diffusivity law. The source term $s(t)$ has to be defined in the routines pressrate\_nl1d(2d).f90.

The fault slip evolution is computed assuming a quasi-static balance of the form:

```math
f\sigma_e = \tau_0 + \kappa * \delta -\frac{\mu}{2c_s} v,
```


where $\tau_0(x,y,t)$ incorporates the initial shear stress (imposed by the initial slip rate and state variable), and a possible external shear stressing, arising either from the boundary conditions, or from another mechanism. An external transient shear stress perturbation can be imposed, but has to be coded in routines tbp\_xxx.f90. $\kappa(x,y)$ is the stress interaction kernel, accounting for the stress redistribution associated with slip along the fault. $\kappa(x,y)$ depends on geometry of the fault and the boundary conditions. The convolution $\kappa * \delta$ is either computed using a spectral approach (infperiodic), a spectral approach avoiding the replication of the fault (cr, following \cite{Cochard1997}), or in the space domain (freesurface). $\mu$ is the shear modulus of the elastic medium, and $c_s$ the shear wave speed of the elastic medium.

