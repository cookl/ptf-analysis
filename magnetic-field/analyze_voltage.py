# This file analyzes the result of a voltage scan
# You can configure what it does using `volt_config.py`
# The input should be a csv file with the following columns:
#

import traceback

from volt_config import *
from consts import *

from numpy import *
import matplotlib as mpl
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import sys, os
from typing import *


cmap_diverging = plt.get_cmap('RdBu')
cmap_normal    = plt.get_cmap('viridis')


coils_for_dim = {
    'x': (3, 4),
    'y': (5, 6),
    'z': (1, 2),
}
up_dn_sd_for_dim = {
    'x': (4, 7, 10),
    'y': (5, 9, 11),
    'z': (6, 8, 12),
}
sup_sdn_ssd_for_dim = {
    'x': (13, 16, 19),
    'y': (14, 18, 20),
    'z': (15, 17, 21),
}



def heatmap_for_property(data, v_min, v_max, coila, coilb, prop_name, diverging=False, is_abs=False, mg=False, is_fit=False):
    fig, ax = plt.subplots(figsize=(8, 6), dpi=192)

    if diverging:
        amax = max(abs(data.max()), abs(data.min()))
        norm = mpl.colors.Normalize(vmin=-amax, vmax=amax)
        mapp = mpl.cm.ScalarMappable(norm, cmap_diverging)
        cmap = cmap_diverging
    elif mg:  # milligauss, for error plot
        norm = mpl.colors.Normalize(vmin=0, vmax=data.max())
        mapp = mpl.cm.ScalarMappable(norm, cmap_normal)
        cmap = cmap_normal
    else:
        if is_abs:
            norm = mpl.colors.LogNorm(vmin=data.min(), vmax=data.max())
        else:
            norm = mpl.colors.Normalize(vmin=data.min(), vmax=data.max())
        mapp = mpl.cm.ScalarMappable(norm, cmap_normal)
        cmap = cmap_normal

    mapp.set_array([])

    extent = (0, upper_lim + 0.5, 0, upper_lim + 0.5 ) if is_fit else (v_min, v_max, v_min, v_max)

    ax.imshow(data.T, cmap=cmap, norm=norm, extent=extent, origin='lower')
    ax.set_xlabel(f"Voltage for Coil {coila} [V]")
    ax.set_ylabel(f"Voltage for Coil {coilb} [V]")
    ax.set_title(prop_name)
    fig.colorbar(mapp)
    return fig, ax


def bilinear(X, m1, m2, b):
    v0, v1 = X
    return (m1*v0) + (m2*v1) + b


def bilinear_fit(data: ndarray, dimension: str, coila: int, coilb: int):
    up, dn, sd = up_dn_sd_for_dim[dimension]
    # sup, sdn, ssd = sup_sdn_ssd_for_dim[dimension]

    v0 = data[:, :, 0].flatten()
    v1 = data[:, :, 1].flatten()

    b_up = data[:, :, up].flatten()
    # s_up = data[:, :, sup].flatten()

    b_dn = data[:, :, dn].flatten()
    # s_dn = data[:, :, sdn].flatten()

    b_sd = data[:, :, sd].flatten()
    # s_sd = data[:, :, ssd].flatten()

    p_up, c_up = curve_fit(bilinear, (v0, v1), b_up)  # , sigma=s_up)
    p_dn, c_dn = curve_fit(bilinear, (v0, v1), b_dn)  # , sigma=s_dn)
    p_sd, c_sd = curve_fit(bilinear, (v0, v1), b_sd)  # , sigma=s_sd)

    e_up = sqrt(diag(c_up))
    e_dn = sqrt(diag(c_dn))
    e_sd = sqrt(diag(c_sd))

    sys.stdout.write(f" {GREEN}Done.{RESET} Parameters:\n")
    for (d, p, e) in (('Up', p_up, e_up), ('Down', p_dn, e_dn), ('Side', p_sd, e_sd)):
        print("    %s% 4s%s:  (% 5.3f ± %5.3f) ⋅ V_%i  +  (% 5.3f ± %5.3f) ⋅ V_%i  +  (% 5.3f ± %5.3f)  =  B%s" % (
            BOLD, d, RESET, p[0], e[0], coila, p[1], e[1], coilb, p[2], e[2], dimension
        ))

    sys.stdout.write(f"  {BOLD}Down{RESET} full compensation:")
    sys.stdout.write(" V_%i  =  %5.3f  +  %5.3f ⋅ V_%i\n" % (coilb, (-p_dn[2]/p_dn[1]), (-p_dn[0]/p_dn[1]), coila))


    v0s = tile(linspace(0, upper_lim + 0.5, 1024), 1024)
    v1s = repeat(linspace(0, upper_lim + 0.5, 1024), 1024)

    bup = bilinear((v0s, v1s), *p_up)
    bdn = bilinear((v0s, v1s), *p_dn)
    bsd = bilinear((v0s, v1s), *p_sd)

    err = abs(bdn)

    imin = argmin(err)

    print(f"  Minimum of fit error at {GREEN}({v0s[imin]}V, {v1s[imin]}V): {MAGENTA}{err[imin]} mG{RESET}.")

    bup.resize(1024, 1024)
    bdn.resize(1024, 1024)
    bsd.resize(1024, 1024)

    return bup, bdn, bsd, (-p_dn[0]/p_dn[1]), (-p_dn[2]/p_dn[1])



if __name__ == "__main__":
    coila, coilb  = coils_for_dim[dimension]
    up, dn, sd    = up_dn_sd_for_dim[dimension]
    sup, sdn, ssd = sup_sdn_ssd_for_dim[dimension]

    print(f"Running analysis for {BOLD}{dimension.upper()}{RESET} (coils {coila} and {coilb}).")

    sys.stdout.write("Loading data...")

    try:
        rdata = loadtxt(scan_out, delimiter=',', skiprows=2, comments='#')
    except Exception:
        sys.stdout.write(f" {BR_RED}Couldn't read CSV.{RESET} Traceback:\n")
        traceback.print_exc()
        exit(1)

    dim = int(sqrt(rdata.shape[0]))
    data = zeros((dim, dim, 22))

    for i in range(dim):
        for j in range(dim):
            data[i, j, :] = rdata[dim * i + j, 2:]

    del rdata

    sys.stdout.write(f" {GREEN}Loaded data{RESET} ({dim} × {dim}).\n")

    # run fits

    sys.stdout.write(f"Running fits...")

    fit_up, fit_dn, fit_sd, fitm, fitb = bilinear_fit(data, dimension, coila, coilb)

    _heatmap = lambda dat, name, diverging=False, is_abs=False, mg=False, is_fit=False: \
        heatmap_for_property(
            dat, lower_lim - resolution / 2, upper_lim + resolution / 2, coila, coilb, name, diverging, is_abs, mg, is_fit
        )

    bstr = 'B_' + dimension

    PLOTS: List[Tuple[str, Callable]] = [
        ('top_b',    lambda data: _heatmap(data[:,:,up],  f'Top ${bstr}$ [G]', True)),
        ('bottom_b', lambda data: _heatmap(data[:,:,dn],  f'Bottom ${bstr}$ [G]', True)),
        ('side_b',   lambda data: _heatmap(data[:,:,sd],  f'Side ${bstr}$ [G]', True)),
        ('s_top_b',  lambda data: _heatmap(data[:,:,sup], f'σ Top ${bstr}$ [G]', True)),
        ('err_sum',  lambda data:
          _heatmap(abs(data[:,:,dn])+abs(data[:,:,sd])+1e-9, f'|Bottom ${bstr}$| + 0.2|Side ${bstr}$|', False, True)
        ),
        ('top_fit',    lambda data: _heatmap(fit_up.T, f'Fit Top ${bstr}$ [G]', True, is_fit=True)),
        ('bottom_fit', lambda data: _heatmap(fit_dn.T, f'Fit Bottom ${bstr}$ [G]', True, is_fit=True)),
        ('side_fit',   lambda data: _heatmap(fit_sd.T, f'Fit Side ${bstr}$ [G]', True, is_fit=True)),
    ]

    err = (abs(data[:,:,dn])).flatten()
    i_min_err = argmin(err)
    # print(data[:,:,0].shape, "    ", i_min_err)
    min_v0    = data[:,:,0].flatten()[i_min_err]
    min_v1    = data[:,:,1].flatten()[i_min_err]
    print(f"Minimum of measured error at {GREEN} ({min_v0}V, {min_v1}V): {MAGENTA}{1000*err[i_min_err]} mG{RESET}.")

    n_plots = len(PLOTS)

    if do_plots:
        for i, (name, closure) in enumerate(PLOTS):
            sys.stdout.write(f"\r{CLEAR_LINE}Running plot {i+1}/{n_plots}: '{name}'")
            fig, ax = closure(data)
            if name == 'bottom_fit':
                ax.plot(
                    [0, upper_lim + 0.5],
                    [fitb, fitm * (upper_lim + 0.5) + fitb],
                    c='black',
                    ls='--',
                    label='Fit'
                )
                ax.set_xlim(0, upper_lim + 0.5)
                ax.set_ylim(0, upper_lim + 0.5)
            # fig.savefig(os.path.join(plot_dir, f'{name}_{dimension}_{run_addl}.pgf'))
            fig.savefig(os.path.join(plot_dir, f'{name}_{dimension}_{run_addl}.png'))

        sys.stdout.write(f"\r{CLEAR_LINE}Plots {GREEN}done{RESET}.")
