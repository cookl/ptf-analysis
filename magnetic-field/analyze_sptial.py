from numpy import *
import matplotlib as mpl
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import os

from scan_config import *


sc_dat = loadtxt(os.path.join(data_dir + f'out_run{run_n:05}.csv'), delimiter=',', skiprows=1)
idx    = (sc_dat[:,0] != 0) & (sc_dat[:,1] != 0)
sc_dat = sc_dat[idx,:]

viridis = plt.get_cmap('viridis')
rd_bu = plt.get_cmap('RdBu')

sc_xs = unique(sc_dat[:, 0])
sc_ys = unique(sc_dat[:, 1])


# print(sc_dat[:,:2])


def dim_plane(dim):
    off_u = 0 if dim == 'x' else 1 if dim == 'y' else 2
    off_d = 0 if dim == 'x' else 2 if dim == 'y' else 1

    sc_top = zeros((len(sc_xs), len(sc_ys)))
    sc_btm = zeros((len(sc_xs), len(sc_ys)))

    inc = True
    lx = len(sc_xs)
    # ly = len(sc_ys)

    for i in range(lx):
        slice_top = sc_dat[i * lx:(i + 1) * lx, 3 + off_u]
        slice_btm = sc_dat[i * lx:(i + 1) * lx, 9 + off_d]
        if not inc:
            slice_top = reversed(slice_top)
            slice_btm = reversed(slice_btm)
        sc_top[i, :] = array(list(slice_top)) * 1000
        sc_btm[i, :] = array(list(slice_btm)) * 1000
        inc = not inc

    diffx = (sc_xs[1] - sc_xs[0]) / 2
    diffy = (sc_ys[1] - sc_ys[0]) / 2
    extent = (sc_xs.min() - diffx, sc_xs.max() + diffx, sc_ys.min() - diffy, sc_ys.max() + diffy)

    # find field differential
    if dim == 'x':
        mean_low = mean(sc_btm[:, 0])
        mean_hi = mean(sc_btm[:, -1])
        print(f"x differential: {mean_hi - mean_low: 03.7f} ({mean_hi} - {mean_low})")
        print(f"         range: {sc_btm.max() - sc_btm.min(): 03.7f}")
        print(f"          mean: {sc_btm.mean()}")
    elif dim == 'y':
        mean_low = mean(sc_btm[0, :])
        mean_hi = mean(sc_btm[-1, :])
        print(f"y differential: {mean_hi - mean_low:03.7f} ({mean_hi} - {mean_low})")
        print(f"         range: {sc_btm.max() - sc_btm.min(): 03.7f}")
        print(f"          mean: {sc_btm.mean()}")
    else:
        mean_hi = mean(sc_top)
        mean_low = mean(sc_btm)
        print(f"z differential: {mean_hi - mean_low: 03.7f} ({mean_hi} - {mean_low})")
        print(f"         range: {sc_top.max() - sc_btm.min(): 03.7f}")
        print(f"          mean: {sc_btm.mean()}")

    return sc_top, sc_btm, extent


x_top, x_btm, x_ext = dim_plane('x')
y_top, y_btm, y_ext = dim_plane('y')
z_top, z_btm, z_ext = dim_plane('z')


#= Plots =#


fig, ax = plt.subplots(figsize=(8,6),dpi=192)

amax = max(abs(x_btm.min()), abs(x_btm.max()))
norm = mpl.colors.Normalize(vmin=-amax, vmax=amax)
mapp = mpl.cm.ScalarMappable(norm, rd_bu)
cmap = rd_bu

mapp.set_array([])

ax.imshow(x_btm, cmap=cmap, norm=norm, extent=x_ext, origin='lower')
ax.set_xlabel(f"X [m]")
ax.set_ylabel(f"Y [m]")
ax.set_title(f"Bottom $B_x$ (Z = 0.4m) [mG] (Run {run_n})")
fig.colorbar(mapp)
fig.savefig(os.path.join(plot_dir, f'scan_{run_n}_x__bottom.png'), transparent=True)

# ---

fig, ax = plt.subplots(figsize=(8,6),dpi=192)

amax = max(abs(y_btm.min()), abs(y_btm.max()))
norm = mpl.colors.Normalize(vmin=-amax, vmax=amax)
mapp = mpl.cm.ScalarMappable(norm, rd_bu)
cmap = rd_bu

mapp.set_array([])

ax.imshow(y_btm, cmap=cmap, norm=norm, extent=y_ext, origin='lower')
ax.set_xlabel(f"X [m]")
ax.set_ylabel(f"Y [m]")
ax.set_title(f"Bottom $B_y$ (Z = 0.4m) [mG] (Run {run_n})")
fig.colorbar(mapp)
fig.savefig(os.path.join(plot_dir, f'scan_{run_n}_y__bottom.png'), transparent=True)

# ---

fig, ax = plt.subplots(figsize=(8,6),dpi=192)

amax = max(abs(z_btm.min()), abs(z_btm.max()))
norm = mpl.colors.Normalize(vmin=-amax, vmax=amax)
mapp = mpl.cm.ScalarMappable(norm, rd_bu)
cmap = rd_bu

mapp.set_array([])

ax.imshow(z_btm, cmap=cmap, norm=norm, extent=z_ext, origin='lower')
ax.set_xlabel(f"X [m]")
ax.set_ylabel(f"Y [m]")
ax.set_title(f"Bottom $B_z$ (Z = 0.4m) [mG] (Run {run_n})")
fig.colorbar(mapp)
fig.savefig(os.path.join(plot_dir, f'scan_{run_n}_z__bottom.png'), transparent=True)


# ---


fig, ax = plt.subplots(figsize=(8,6),dpi=192)

z_diff = z_top - z_btm

amax = max(abs(z_diff.min()), abs(z_diff.max()))
norm = mpl.colors.Normalize(vmin=-amax, vmax=amax)
mapp = mpl.cm.ScalarMappable(norm, rd_bu)
cmap = rd_bu

mapp.set_array([])

ax.imshow(z_diff, cmap=cmap, norm=norm, extent=z_ext, origin='lower')
ax.set_xlabel(f"X [m]")
ax.set_ylabel(f"Y [m]")
ax.set_title(f"Top $B_z$ - Bottom $B_z$ [mG] (Run {run_n})")
fig.colorbar(mapp)
fig.savefig(os.path.join(plot_dir, f'scan_{run_n}_zd.png'), transparent=True)
