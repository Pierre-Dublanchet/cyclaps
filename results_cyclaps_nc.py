#%%
import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy.special import erfc
import configparser

# Constants
YEAR = 365 * 24 * 3600
def safe_float(value):
    if isinstance(value, str):
        value = value.strip()
        value = value.replace('d', 'e').replace('D', 'e')
    return float(value)

def readparametres_cfg(rep, filep):
    """Read parameters from a config file (similar to MATLAB's ini2struct)."""
    config = configparser.ConfigParser()
    config.read(os.path.join(rep, filep))
    par = {}
    for section in config.sections():
        for key, value in config[section].items():
            try:
                par[key] = safe_float(value)
            except ValueError:
                par[key] = value.strip()
    # Handle arrays and special cases as in MATLAB
    par['paraminj'] = np.array([
        par.get('toff', 0),
        par.get('xinj', 0),
        par.get('yinj', 0),
        par.get('qinj', 0)
    ])
    par['toff'], par['xinj'], par['yinj'], par['qinj'] = par['paraminj']
    # Load initial conditions from NetCDF
    with Dataset(os.path.join(rep, 'init.nc'), 'r') as nc:
        par['amap'] = nc.variables['a'][:].reshape(int(par['ny']), int(par['nx']))
        par['bmap'] = nc.variables['b'][:].reshape(int(par['ny']), int(par['nx']))
        par['dcmap'] = nc.variables['dc'][:].reshape(int(par['ny']), int(par['nx']))
        par['smap'] = nc.variables['s'][:].reshape(int(par['ny']), int(par['nx']))
        par['vimap'] = nc.variables['vi'][:].reshape(int(par['ny']), int(par['nx']))
        par['thimap'] = nc.variables['thi'][:].reshape(int(par['ny']), int(par['nx']))
        par['pimap'] = np.zeros_like(par['amap'])  # or load if present
    return par

def extractimeseriesnc(rep):
    """Extract time series data from NetCDF."""
    with Dataset(os.path.join(rep, 'qmoy.nc'), 'r') as nc:
        t = nc.variables['time'][:]
        dt = nc.variables['time delai'][:]
        v = nc.variables['mean slip rate'][:]
        vmax = nc.variables['max slip rate'][:]
        th = nc.variables['mean state'][:]
        tau = nc.variables['mean shear stress'][:]
        u = nc.variables['mean slip'][:]
        vmvw = nc.variables['mean slip rate vw'][:] if 'mean slip rate vw' in nc.variables else np.zeros_like(t)
    return t, dt, v, vmax, vmvw, th, tau, u

def plt_av_var(par, rep, kfig):
    """Plot spatial average of variables."""
    t, dt, v, vmax, vmvw, th, tau, u = extractimeseriesnc(rep)
    #mrate = v * par['nx'] * par['dx'] * par['Y'] / (2 * (1 + par['nu']))
    #mratevw = vmvw * par['nx'] * par['dx'] * par['Y'] / (2 * (1 + par['nu']))

    mrate = (v* par['nx']* par['dx']* par['young_mod']/ (2 * (1 + par['poisson_ratio'])))

    mratevw = (vmvw* par['nx']* par['dx']* par['young_mod']/ (2 * (1 + par['poisson_ratio'])))

    plt.figure(kfig + 1)
    plt.subplot(2, 2, 1)
    plt.semilogy(t / YEAR, v, '-k')
    plt.semilogy(t / YEAR, vmax, '-r')
    plt.xlabel('t (years)')
    plt.ylabel('v (m/s)')

    plt.subplot(2, 2, 2)
    plt.semilogy(t / YEAR, th, '-k')
    plt.xlabel('t (years)')
    plt.ylabel(r'$\theta$')

    plt.subplot(2, 2, 3)
    plt.plot(t / YEAR, u, '-k')
    plt.plot(t / YEAR, par['vplate'] * t, '--r')
    plt.xlabel('t (years)')
    plt.ylabel('u (m)')

    plt.subplot(2, 2, 4)
    plt.plot(t / YEAR, tau * 1e-6, '-k')
    plt.xlabel('t (years)')
    plt.ylabel(r'$\tau$ (MPa)')

    plt.figure(kfig + 2)
    plt.semilogy(t / YEAR, mrate, '-k')
    plt.semilogy(t / YEAR, mratevw, '-r')
    plt.xlabel('t (years)')
    plt.ylabel('Moment rate (N.m/s)')
    plt.legend(['VW+VS', 'VW'])

    plt.figure(kfig + 3)
    plt.semilogy(np.arange(len(t)), dt, 'o-k')
    plt.xlabel('t (years)')
    plt.ylabel('time step dt (s)')

    return kfig + 3

def plt_loc_var(par, rep, kfig):
    """Plot local variables time series."""
    tloc, vloc, thloc, tauloc, uloc, ploc, qdarcyloc = extract_local_timeseries_nc(rep, par)
    for i in range(int(par['nploc'])):
        plt.figure(kfig + 1)
        plt.subplot(3, 2, 1)
        plt.semilogy(tloc[:, i] / YEAR, vloc[:, i], '-')
        plt.xlabel('t (years)')
        plt.ylabel('v (m/s)')

        plt.subplot(3, 2, 2)
        plt.semilogy(tloc[:, i] / YEAR, thloc[:, i], '-')
        plt.xlabel('t (years)')
        plt.ylabel(r'$\theta$')

        plt.subplot(3, 2, 3)
        plt.plot(tloc[:, i] / YEAR, uloc[:, i], '-')
        plt.plot(tloc[:, i] / YEAR, par['vplate'] * tloc[:, i], '--r')
        plt.xlabel('t (years)')
        plt.ylabel('u (m)')

        plt.subplot(3, 2, 4)
        plt.plot(tloc[:, i] / YEAR, tauloc[:, i] * 1e-6, '-')
        plt.xlabel('t (years)')
        plt.ylabel(r'$\tau$ (MPa)')

        plt.subplot(3, 2, 5)
        plt.plot(tloc[:, i] / YEAR, ploc[:, i] * 1e-6, '-')
        plt.xlabel('t (years)')
        plt.ylabel('P (MPa)')

        plt.subplot(3, 2, 6)
        plt.plot(tloc[:, i] / YEAR, qdarcyloc[:, i], '-')
        plt.xlabel('t (years)')
        plt.ylabel(r'$q_{Darcy}$ (m/s)')

    return kfig + 1

def plt_hv_profils(par, rep, kfig):
    """Plot HV profiles of results."""
    x = np.linspace(-0.5 * (par['nx'] - 1) * par['dx'], 0.5 * (par['nx'] - 1) * par['dx'], int(par['nx'])) - 0.5 * par['dx']
    y = np.linspace(-0.5 * (par['ny'] - 1) * par['dy'], 0.5 * (par['ny'] - 1) * par['dy'], int(par['ny'])) - 0.5 * par['dy']
    ns = len([f for f in os.listdir(rep) if f.startswith('profilsx') and f.endswith('.nc')])
    ii = [i for i in range(1, ns + 1) if i % 100 == 0 or (i % 10 == 0 and np.max(extractimeseriesnc(rep)[3]) > 1e-4)]

    for k, i in enumerate(ii):
        with Dataset(os.path.join(rep, f'profilsx{i}.nc'), 'r') as nc:
            vx = nc.variables['slip rate'][:]
            ux = nc.variables['slip'][:]
            thx = nc.variables['state'][:]
            taux = nc.variables['shear stress'][:]
            px = nc.variables['pore pressure'][:] if 'pore pressure' in nc.variables else np.zeros_like(vx)
        with Dataset(os.path.join(rep, f'profilsy{i}.nc'), 'r') as nc:
            vy = nc.variables['slip rate'][:]
            uy = nc.variables['slip'][:]
            thy = nc.variables['state'][:]
            tauy = nc.variables['shear stress'][:]
            py = nc.variables['pore pressure'][:] if 'pore pressure' in nc.variables else np.zeros_like(vy)


        px = np.squeeze(px)
        py = np.squeeze(py)
        taux = np.squeeze(taux)
        tauy = np.squeeze(tauy)
        ux = np.squeeze(ux)
        uy = np.squeeze(uy)
        vx = np.squeeze(vx)
        vy = np.squeeze(vy)
        thx = np.squeeze(thx)
        thy = np.squeeze(thy)
        
        if k % 1 == 0 or k == 0:
            plt.figure(kfig + 1)
            plt.subplot(5, 2, 1)
            plt.plot(x, px * 1e-6, '-k')
            plt.xlabel('x (m)')
            plt.ylabel('P (MPa)')

            plt.subplot(5, 2, 2)
            plt.plot(y, py * 1e-6, '-k')
            plt.xlabel('y (m)')
            plt.ylabel('P (MPa)')

            plt.subplot(5, 2, 3)
            plt.plot(x, taux * 1e-6, '-k')
            plt.xlabel('x (m)')
            plt.ylabel('tau (MPa)')

            plt.subplot(5, 2, 4)
            plt.plot(y, tauy * 1e-6, '-k')
            plt.xlabel('y (m)')
            plt.ylabel('tau (MPa)')

            plt.subplot(5, 2, 5)
            plt.semilogy(x, vx, '-k')
            plt.xlabel('x (m)')
            plt.ylabel('v (m/s)')

            plt.subplot(5, 2, 6)
            plt.semilogy(y, vy, '-k')
            plt.xlabel('y (m)')
            plt.ylabel('v (m/s)')

            plt.subplot(5, 2, 7)
            plt.semilogy(x, thx, '-k')
            plt.xlabel('x (m)')
            plt.ylabel(r'$\theta$ (s)')

            plt.subplot(5, 2, 8)
            plt.semilogy(y, thy, '-k')
            plt.xlabel('y (m)')
            plt.ylabel(r'$\theta$ (s)')

            plt.subplot(5, 2, 9)
            plt.plot(x, ux, '-k')
            plt.xlabel('x (m)')
            plt.ylabel('u (m)')

            plt.subplot(5, 2, 10)
            plt.plot(y, uy, '-k')
            plt.xlabel('y (m)')
            plt.ylabel('u (m)')

    return kfig + 1

def plt_maps(par, rep, kfig):
    """Plot 2D maps of results."""
    x = np.linspace(-0.5 * (par['nx'] - 1) * par['dx'], 0.5 * (par['nx'] - 1) * par['dx'], int(par['nx'])) - 0.5 * par['dx']
    y = np.linspace(-0.5 * (par['ny'] - 1) * par['dy'], 0.5 * (par['ny'] - 1) * par['dy'], int(par['ny'])) - 0.5 * par['dy']
    ns = len([f for f in os.listdir(rep) if f.startswith('maps') and f.endswith('.nc')])
    ii = np.arange(1, ns + 1)

    for i in ii:
        with Dataset(os.path.join(rep, f'maps{i}.nc'), 'r') as nc:
            vx = nc.variables['slip rate'][:]
            ux = nc.variables['slip'][:]
            thx = nc.variables['state'][:]
            taux = nc.variables['shear stress'][:]
            px = nc.variables['pore pressure'][:] if 'pore pressure' in nc.variables else np.zeros_like(vx)
            qx = nc.variables['darcy velocity'][:] if 'darcy velocity' in nc.variables else np.zeros_like(vx)
            t = nc.variables['time'][:]

        mv = vx.reshape(int(par['ny']), int(par['nx']))
        mu = ux.reshape(int(par['ny']), int(par['nx']))
        mth = thx.reshape(int(par['ny']), int(par['nx']))
        mtau = taux.reshape(int(par['ny']), int(par['nx']))
        msigmae = par['smap'] - px.reshape(int(par['ny']), int(par['nx']))
        mqdarcy = qx.reshape(int(par['ny']), int(par['nx']))

        plt.figure(1 + kfig)
        plt.clf()
        plt.subplot(3, 2, 1)
        plt.imshow(msigmae * 1e-6, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower', aspect='auto')
        plt.colorbar(label=r'$\sigma-P$ (MPa)')
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')

        plt.subplot(3, 2, 2)
        plt.imshow(mtau * 1e-6, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower', aspect='auto')
        plt.colorbar(label=r'$\tau$ (MPa)')
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')

        plt.subplot(3, 2, 3)
        plt.imshow(mu, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower', aspect='auto')
        plt.colorbar(label=r'$\delta$ (m)')
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')

        plt.subplot(3, 2, 4)
        plt.imshow(np.log10(mv), extent=[x.min(), x.max(), y.min(), y.max()], origin='lower', aspect='auto')
        plt.colorbar(label='log v (m/s)')
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')

        plt.subplot(3, 2, 5)
        plt.imshow(np.log10(mth), extent=[x.min(), x.max(), y.min(), y.max()], origin='lower', aspect='auto')
        plt.colorbar(label=r'log $\theta$ (m/s)')
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')

        plt.subplot(3, 2, 6)
        plt.imshow(mqdarcy, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower', aspect='auto')
        plt.colorbar(label=r'$v_{Darcy}$ (m/s)')
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')

    return kfig + 1

def plt_catalog(par, rep, kfig):
    """Plot earthquake catalog."""
    t0, dt0, tev, xi, yi, xb, yb, s, m0, dtau, taui, tauf, du = extractcatnc(rep, par)
    mw = (2/3) * np.log10(m0) - 6
    rad = np.sqrt(xb**2 + yb**2)
    reqk = np.sqrt(s) if par['dx'] * par['dy'] > 0 else 0.5 * s / par['dx']
    neqk = len(t0)

    plt.figure(1 + kfig)
    plt.subplot(2, 2, 1)
    plt.plot(t0 / YEAR, mw, 'ok')
    plt.xlabel('time $t$ (years)')
    plt.ylabel('magnitude $M_w$')

    plt.subplot(2, 2, 2)
    plt.plot(t0 / YEAR, -dtau * 1e-6, 'ok')
    plt.xlabel('time $t$ (years)')
    plt.ylabel('stress drop $\Delta \tau$ (MPa)')

    plt.subplot(2, 2, 3)
    plt.plot(t0 / YEAR, du, 'ok')
    plt.xlabel('time $t$ (years)')
    plt.ylabel('coseismic slip $\Delta u$ (m)')

    plt.subplot(2, 2, 4)
    plt.plot(t0 / YEAR, tev, 'ok')
    plt.xlabel('time $t$ (years)')
    plt.ylabel('event duration $t_e$ (s)')

    plt.figure(2 + kfig)
    if par['ny'] > 1:
        for i in range(neqk):
            plt.plot(xi[i], yi[i], '+k')
            plt.gca().add_patch(plt.Rectangle((xb[i] - reqk[i], yb[i] - reqk[i]), 2 * reqk[i], 2 * reqk[i], fill=None, color='k'))
            plt.plot(xb[i], yb[i], '+r')
            plt.gca().add_patch(plt.Rectangle((xb[i] - reqk[i], yb[i] - reqk[i]), 2 * reqk[i], 2 * reqk[i], fill=None, color='r'))
        plt.xlim([np.min(xi), np.max(xi)])
        plt.ylim([np.min(yi), np.max(yi)])
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')
    else:
        for i in range(neqk):
            plt.plot([xb[i] - reqk[i], xb[i] + reqk[i]], [t0[i] / YEAR, t0[i] / YEAR], '-k', linewidth=2)
        plt.xlim([np.min(xi), np.max(xi)])
        plt.xlabel('x (m)')
        plt.ylabel('time t (years)')

    plt.figure(3 + kfig)
    plt.loglog(m0, tev, 'ok')
    plt.xlabel('seismic moment $M_0$ (N.m)')
    plt.ylabel('event duration $t_e$ (s)')

    return kfig + 3

def press_solution(q0, alpha, x, t, toff, beta, phi):
    """Analytical pressure solution."""
    n = len(x)
    p = np.zeros(n)
    if t > 0 and t < toff:
        eta = np.abs(x) / (2 * np.sqrt(alpha * t))
        G = np.sqrt(t) * (np.exp(-eta**2) / np.sqrt(np.pi) - eta * erfc(eta))
        p = q0 * G / (beta * phi * np.sqrt(alpha))
    elif t > toff:
        eta1 = np.abs(x) / (2 * np.sqrt(alpha * t))
        G1 = np.sqrt(t) * (np.exp(-eta1**2) / np.sqrt(np.pi) - eta1 * erfc(eta1))
        eta2 = np.abs(x) / (2 * np.sqrt(alpha * (t - toff)))
        G2 = np.sqrt(t - toff) * (np.exp(-eta2**2) / np.sqrt(np.pi) - eta2 * erfc(eta2))
        p = q0 * (G1 - G2) / (beta * phi * np.sqrt(alpha))
    return p

def extract_local_timeseries_nc(rep, par):
    """Extract local time series data from NetCDF."""
    tloc, vloc, thloc, tauloc, uloc, ploc, qdarcyloc = [], [], [], [], [], [], []
    for i in range(int(par['nploc'])):
        with Dataset(os.path.join(rep, f'qloc{i+1}.nc'), 'r') as nc:
            vt = nc.variables['time'][:]
            vdt = nc.variables['time delai'][:]
            vu = nc.variables['slip'][:]
            vth = nc.variables['state'][:]
            vv = nc.variables['slip rate'][:]
            vtau = nc.variables['shear stress'][:]
            vp = nc.variables['pore pressure'][:] if 'pore pressure' in nc.variables else np.zeros_like(vv)
            vq = nc.variables['darcy vel'][:] if 'darcy vel' in nc.variables else np.zeros_like(vv)
        tloc.append(vt)
        vloc.append(vv)
        thloc.append(vth)
        tauloc.append(vtau)
        uloc.append(vu)
        ploc.append(vp)
        qdarcyloc.append(vq)
    return np.array(tloc).T, np.array(vloc).T, np.array(thloc).T, np.array(tauloc).T, np.array(uloc).T, np.array(ploc).T, np.array(qdarcyloc).T

def extractcatnc(rep, par):
    """Extract earthquake catalog from NetCDF."""
    with Dataset(os.path.join(rep, 'earthquake_catalogue.nc'), 'r') as nc:
        t0 = nc.variables['onset time'][:]
        dt0 = nc.variables['onset time delai'][:]
        tev = nc.variables['event duration'][:]
        xi = nc.variables['x initiation'][:]
        yi = nc.variables['y initiation'][:]
        xb = nc.variables['x barycenter'][:]
        yb = nc.variables['y barycenter'][:]
        s = par['dx'] * par['dy'] * nc.variables['number of elements'][:] if par['dy'] * par['dx'] > 0 else par['dx'] * par['dx'] * nc.variables['number of elements'][:]
        m0 = nc.variables['coseismic moment'][:]
        dtau = nc.variables['coseismic stress drop'][:]
        du = nc.variables['coseismic slip'][:]
        taui = nc.variables['shear stress init'][:] if 'shear stress init' in nc.variables else np.zeros_like(t0)
        tauf = nc.variables['shear stress final'][:] if 'shear stress final' in nc.variables else np.zeros_like(t0)
    return t0, dt0, tev, xi, yi, xb, yb, s, m0, dtau, taui, tauf, du
#%%
# Example usage:
if __name__ == "__main__":
    rep = '../RESULTS/'
    filep = 'parametres_fault_rns.cfg'
    par = readparametres_cfg(rep, filep)
    kfig = 0
    #kfig = plt_av_var(par,rep,kfig)
    #kfig = plt_loc_var(par,rep,kfig)
    kfig = plt_hv_profils(par,rep,kfig)
    #kfig = plt_maps(par, rep, kfig)
    #kfig = plt_catalog(par,rep,kfig)
    plt.show()

#%%
