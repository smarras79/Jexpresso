#!/usr/bin/env python3
"""Figures of Leake, Linton & Antiochos (2022, ApJ 934:10) from the output of problems/MHD/fluxEmergenceLeake2022.

RUN is an output directory (output/MHD/fluxEmergenceLeake2022/output-<date>, the one holding iter_N.pvtu);
a restarted run is DIR1,DIR2. LABEL=RUN names a curve; without a label the run's :lk22_sim and :lk22_dipole
(read from the user_inputs.jl copied into RUN) give the paper's "theta=n pi/4" or "No dipole" and its colors.

  ic      [LABEL=RUN ...]                     Fig. 1 (right): T (K) and rho (kg m^-3) of the initial atmosphere
  xsec    RUN --var Baxial|JB --times ...     Figs. 4, 6, 7, 12 (left): the cross-section halfway along the rope
                                              (y_r = 0): log10|B_axial| (G) or |J/J0|/|B/B0|
  height  LABEL=RUN ...                       Fig. 8: flux-rope height on x = y = 0 (axial-field maximum and sign
                                              change of the orthogonal horizontal field)
  energy  LABEL=RUN ...                       Fig. 9: free magnetic and kinetic energy above z = 10 L0, in E0
  lorentz RUN --times 130,160                 Fig. 12 (right): vertical Lorentz force, tension and pressure on x = y = 0
  render  RUN --times ...                     Figs. 1 (left), 3, 5: field lines seeded on x = y = 0, colored by
                                              connectivity, Bz at z = 0 (needs pyvista)

Lengths in Mm (L0 = 0.17 Mm), B in G (B0 = 1300 G), time in t0 = 24.9 s, energy in E0 = 6.607e19 J, as the
paper. Needs numpy and matplotlib (plus pyvista for render). Run any command with -h for its options.
"""
import argparse, os, re, sys

sys.dont_write_bytecode = True          # no tools/__pycache__ in the work tree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_leake2014 import (snapshots, pick_times, read_lattice, ddx, current, runs_arg, savefig,  # noqa: E402
                            trapz, lim, atmosphere, AX)

import numpy as np                      # noqa: E402
import matplotlib                       # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt         # noqa: E402

L0_MM, B0_G, T0_K, RHO0 = 0.17, 1300.0, 5637.0, 2.887e-4
ZTOP = 530.5
# Fig. 8 legend: theta = 7, 0, 1 pi/4 greens (stable), 2, 6 blues (intermediate), 3, 4, 5 red/orange/magenta
SIM_COLOR = {0: "#c9b200", 1: "#5fd35f", 7: "#1b8a2e", 2: "#2346e0", 6: "#0b1a6e",
             3: "#e3261e", 4: "#f28e00", 5: "#d523d0"}


# ---------------------------------------------------------------------------------------------------------
# runs: theta and labels
# ---------------------------------------------------------------------------------------------------------
def run_setup(run, override=None):
    """(n, dipole) of RUN from its copy of user_inputs.jl (:lk22_sim, :lk22_dipole); override = n."""
    n, dip = 4, True
    f = os.path.join(run.split(",")[0], "user_inputs.jl")
    if os.path.isfile(f):
        txt = re.sub(r"#[^\n]*", "", open(f).read())          # the keys, not the comments that mention them
        m = re.search(r":lk22_sim\s*=>\s*([-+0-9.]+)", txt)
        if m:
            n = float(m.group(1))
            n = int(n) if n == int(n) else n
        m = re.search(r":lk22_dipole\s*=>\s*(true|false)", txt)
        dip = (m.group(1) == "true") if m else True
    else:
        print(f"   {run}: no user_inputs.jl copy, assuming :lk22_sim => {n}")
    if override is not None:
        n = int(override) if float(override).is_integer() else override
    return n, dip


def theta(n):
    h = np.sqrt(0.5)
    if float(n).is_integer():
        return ((1.0, 0.0), (h, h), (0.0, 1.0), (-h, h), (-1.0, 0.0), (-h, -h), (0.0, -1.0), (h, -h))[int(n) % 8]
    return np.cos(n * np.pi / 4), np.sin(n * np.pi / 4)


def labelled(items, override=None):
    """[(label, run, n, dipole, color, linestyle)]"""
    out = []
    for i, (lab, run) in enumerate(runs_arg(items)):
        n, dip = run_setup(run, override)
        auto = lab == run
        if not dip:
            lab = "No dipole" if auto else lab
            out.append((lab, run, n, dip, "k", ":"))
        else:
            lab = (rf"$\theta={n}\pi/4$" if auto else lab)
            out.append((lab, run, n, dip, SIM_COLOR.get(n, f"C{i}") if float(n).is_integer() else f"C{i}", "-"))
    return out


# ---------------------------------------------------------------------------------------------------------
# ic: Fig. 1 (right)
# ---------------------------------------------------------------------------------------------------------
def cmd_ic(a):
    z = np.linspace(-30.0, ZTOP, 4001)
    T, p = atmosphere(z, a.cz_gradient, a.mf)
    rho = a.mf * p / T
    fig, ax = plt.subplots(figsize=(a.width, 0.7 * a.width))
    ax.semilogy(z * L0_MM, T * T0_K, color="tab:red")
    ax.set_ylabel("Temperature (K)", color="tab:red")
    ax.tick_params(axis="y", colors="tab:red")
    ax2 = ax.twinx()
    ax2.semilogy(z * L0_MM, rho * RHO0, color="tab:blue")
    ax2.set_ylabel(r"Density (kg/m$^3$)", color="tab:blue")
    ax2.tick_params(axis="y", colors="tab:blue")
    for lab, run, n, dip, c, ls in labelled(a.runs):
        t, f = snapshots(run)[0]
        (x, y, zc), F = read_lattice(f, ("T", "ρ"), box=((-1e-3, 1e-3), (-1e-3, 1e-3), (-1e9, 1e9)))
        ax.semilogy(zc * L0_MM, F["T"][:, 0, 0] * T0_K, "o", ms=2.5, color="tab:red", alpha=0.6)
        ax2.semilogy(zc * L0_MM, F["ρ"][:, 0, 0] * RHO0, "o", ms=2.5, color="tab:blue", alpha=0.6)
        print(f"   {lab}: nodal T, rho on x = y = 0 at t = {t:g} overlaid")
    ax.set_xlabel("Height (Mm)")
    ax.set_xlim(z[0] * L0_MM, z[-1] * L0_MM)
    savefig(fig, os.path.join(a.outdir, f"leake22_fig1_atmosphere.{a.format}"), a.dpi)


# ---------------------------------------------------------------------------------------------------------
# helpers: the rope's cross-section plane, the vertical line x = y = 0
# ---------------------------------------------------------------------------------------------------------
def sample_columns(coords, arrs, px, py):
    """Bilinear (x, y) interpolation of arrs[nz, ny, nx] at the columns (px, py): [nz, npts] each."""
    x, y, z = coords
    ix = np.clip(np.searchsorted(x, px) - 1, 0, len(x) - 2)
    iy = np.clip(np.searchsorted(y, py) - 1, 0, len(y) - 2)
    wx = np.clip((px - x[ix]) / (x[ix + 1] - x[ix]), 0, 1)
    wy = np.clip((py - y[iy]) / (y[iy + 1] - y[iy]), 0, 1)
    out = []
    for A in arrs:
        out.append((1 - wx) * (1 - wy) * A[:, iy, ix] + wx * (1 - wy) * A[:, iy, ix + 1]
                   + (1 - wx) * wy * A[:, iy + 1, ix] + wx * wy * A[:, iy + 1, ix + 1])
    return out


def column(coords, arrs):
    """arrs on the vertical line x = y = 0 (interpolated)."""
    return [c[:, 0] for c in sample_columns(coords, arrs, np.array([0.0]), np.array([0.0]))]


def rope_crossings(z, Bperp, zlim, persist=2):
    """Heights where the horizontal field orthogonal to the rope axis turns from -x_r to +x_r going up (the
    right-handed rope's axis: twist along -x_r below it, +x_r above), inside zlim. The sign must hold on
    `persist` nodes either side, which skips the node-to-node zigzag of an under-resolved field."""
    out = []
    n = len(z)
    for k in range(len(z) - 1):
        a, b = Bperp[k], Bperp[k + 1]
        if not (all(Bperp[j] < 0 for j in range(max(0, k - persist + 1), k + 1)) and
                all(Bperp[j] >= 0 for j in range(k + 1, min(n, k + 1 + persist)))):
            continue
        if a < 0 <= b and a != b:
            zc = z[k] + (z[k + 1] - z[k]) * (-a) / (b - a)
            if zlim[0] <= zc <= zlim[1]:
                out.append(zc)
    return out


def axial_maxima(z, Bax, zlim):
    return [z[k] for k in range(1, len(z) - 1)
            if Bax[k] > 0 and Bax[k] >= Bax[k - 1] and Bax[k] > Bax[k + 1] and zlim[0] <= z[k] <= zlim[1]]


def rope_height(z, Bax, Bperp, zlim, persist=2):
    """(method 1: axial-field maximum nearest to method 2, method 2: highest sign change, lower sign change)."""
    zc = rope_crossings(z, Bperp, zlim, persist)
    if not zc:
        return np.nan, np.nan, np.nan
    z2 = max(zc)
    zm = axial_maxima(z, Bax, zlim)
    z1 = min(zm, key=lambda q: abs(q - z2)) if zm else np.nan
    return z1, z2, (min(zc) if len(zc) > 1 else np.nan)


# ---------------------------------------------------------------------------------------------------------
# xsec: Figs. 4, 6, 7, 12 (left)
# ---------------------------------------------------------------------------------------------------------
def cmd_xsec(a):
    n, dip = run_setup(a.run, a.sim)
    c, s = theta(n)
    xr = np.linspace(-a.xr_max, a.xr_max, a.nsamp) / L0_MM            # in L0
    zlim = lim(a.zlim)
    W = a.xr_max / L0_MM + 60.0          # margins reach the next lattice line on the coarse meshes
    box = ((-W, W), (-W, W), ((zlim[0] / L0_MM - 60.0, zlim[1] / L0_MM + 60.0) if zlim else (-1e9, 1e9)))
    snaps = pick_times(snapshots(a.run), a.times)
    panels = []
    for t, f in snaps:
        coords, F = read_lattice(f, ("Bx", "By", "Bz"), box)
        px, py = c * xr, s * xr
        Bx, By, Bz = sample_columns(coords, [F["Bx"], F["By"], F["Bz"]], px, py)
        Bax = -s * Bx + c * By
        logB = np.log10(np.maximum(np.abs(Bax) * B0_G, 1e-30))
        JB = None
        if a.var == "JB" or a.contours:
            jx, jy, jz = current(coords, F)
            J = np.sqrt(jx**2 + jy**2 + jz**2)
            Bm = np.sqrt(F["Bx"]**2 + F["By"]**2 + F["Bz"]**2)
            del jx, jy, jz
            Jp, Bp = sample_columns(coords, [J, Bm], px, py)
            JB = Jp / np.maximum(Bp, 1e-300)
        panels.append((t, coords[2], logB, JB))
        print(f"   {f}: t = {t:g}")
    var = a.var
    cmap = a.cmap or ("jet" if var == "Baxial" or a.style != "fig4" else "Reds")
    vmin, vmax = lim(a.clim) if a.clim else ((0.0, 2.0) if var == "Baxial" else ((0.0, 0.5) if a.style == "fig4" else (0.0, 1.0)))
    nP = len(panels)
    zr = zlim if zlim else (panels[0][1][0] * L0_MM, panels[0][1][-1] * L0_MM)
    pw = a.width / nP
    fig, axs = plt.subplots(1, nP, figsize=(a.width, pw * (zr[1] - zr[0]) / (2 * a.xr_max) + 1.0), squeeze=False)
    for ax, (t, z, logB, JB) in zip(axs[0], panels):
        Fp = logB if var == "Baxial" else JB
        ax.set_rasterization_zorder(0)
        im = ax.pcolormesh(xr * L0_MM, z * L0_MM, Fp, shading="gouraud", cmap=cmap, vmin=vmin, vmax=vmax, zorder=-1)
        if a.contours:
            ax.contour(xr * L0_MM, z * L0_MM, logB, levels=np.linspace(a.contour_range[0], a.contour_range[1], a.contours),
                       colors="k", linewidths=0.5)
        if a.photosphere:                # the paper's Fig. 4: the photospheric plane seen edge-on
            ax.axhspan(-0.6, 0.0, color="0.55", zorder=0)
        ax.set_xlim(-a.xr_max, a.xr_max)
        ax.set_ylim(*zr)
        ax.set_aspect("equal")
        ax.set_xlabel(r"$X_r$ (Mm)")
        ax.set_title(rf"{a.title or ''}$t = {t:g}\,t_0$", fontsize=10, loc="left")
    axs[0][0].set_ylabel("Z (Mm)")
    cb = fig.colorbar(im, ax=axs[0].tolist(), shrink=0.85, pad=0.02)
    cb.set_label(r"$\log_{10}|B_{\rm axial}|$ (G)" if var == "Baxial" else r"$|J/J_0|/|B/B_0|$")
    tag = a.name or f"leake22_xsec_{var}_sim{n}"
    savefig(fig, os.path.join(a.outdir, f"{tag}.{a.format}"), a.dpi)


# ---------------------------------------------------------------------------------------------------------
# height: Fig. 8
# ---------------------------------------------------------------------------------------------------------
def cmd_height(a):
    os.makedirs(a.outdir, exist_ok=True)
    zlim = (a.zmin / L0_MM, a.zmax / L0_MM)
    m = a.margin
    fig, ax = plt.subplots(figsize=(a.width, 0.75 * a.width))
    for lab, run, n, dip, col, ls in labelled(a.runs, a.sim):
        c, s = theta(n)
        rows = []
        for t, f in pick_times(snapshots(run), a.times):
            coords, F = read_lattice(f, ("Bx", "By"), ((-m, m), (-m, m), (-1e9, 1e9)))
            Bx, By = column(coords, [F["Bx"], F["By"]])
            z = coords[2]
            z1, z2, zlow = rope_height(z, -s * Bx + c * By, c * Bx + s * By, zlim, a.persist)
            rows.append((t, z1 * L0_MM, z2 * L0_MM, zlow * L0_MM))
            print(f"   {lab}: t = {t:8.3f}  axial max {z1 * L0_MM:7.3f} Mm, sign change {z2 * L0_MM:7.3f} Mm"
                  + (f", lower {zlow * L0_MM:7.3f} Mm" if np.isfinite(zlow) else ""))
        r = np.array(rows)
        tag = re.sub(r"[^A-Za-z0-9]+", "", lab) or "run"
        np.savetxt(os.path.join(a.outdir, f"leake22_height_{tag}.csv"), r, fmt="%.6g",
                   header="t z_axialmax_Mm z_signchange_Mm z_lower_signchange_Mm   (on x = y = 0)")
        lo, hi = np.fmin(r[:, 1], r[:, 2]), np.fmax(r[:, 1], r[:, 2])
        ax.fill_between(r[:, 0], lo, hi, color=col, alpha=0.3, lw=0)
        ax.plot(r[:, 0], np.nanmean(r[:, 1:3], axis=1), color=col, ls=ls, lw=1.8, label=lab)
        ax.plot(r[:, 0], r[:, 3], color=col, ls=(0, (2, 2)), lw=1.0)
    ax.set_xlabel(r"Time/t0"); ax.set_ylabel("Height (Mm)")
    ax.legend(fontsize=8)
    savefig(fig, os.path.join(a.outdir, f"leake22_fig8_height.{a.format}"), a.dpi)


# ---------------------------------------------------------------------------------------------------------
# energy: Fig. 9 (free magnetic energy from a Neumann potential field, kinetic energy, above z = 10 L0)
# ---------------------------------------------------------------------------------------------------------
def dual_widths(c):
    h = np.empty_like(c)
    h[1:-1] = 0.5 * (c[2:] - c[:-2])
    h[0], h[-1] = 0.5 * (c[1] - c[0]), 0.5 * (c[-1] - c[-2])
    return h


def potential_energy(x, y, z, Bx, By, Bz, tol, maxit):
    """Energy of the potential field with the normal component of B on the six faces of the box, by a
    vertex-centred finite-volume Laplacian on the node lattice (conjugate gradients, Jacobi preconditioned).
    Returns (E_pot, iterations, relative residual, net boundary flux / total boundary |flux|)."""
    hx, hy, hz = dual_widths(x), dual_widths(y), dual_widths(z)
    Cx = (hz[:, None, None] * hy[None, :, None]) / np.diff(x)[None, None, :]
    Cy = (hz[:, None, None] * hx[None, None, :]) / np.diff(y)[None, :, None]
    Cz = (hy[None, :, None] * hx[None, None, :]) / np.diff(z)[:, None, None]
    diag = np.zeros(Bx.shape)
    diag[:, :, :-1] += Cx; diag[:, :, 1:] += Cx
    diag[:, :-1, :] += Cy; diag[:, 1:, :] += Cy
    diag[:-1, :, :] += Cz; diag[1:, :, :] += Cz

    def A(phi):                          # -Laplacian (SPD on the complement of the constants)
        out = np.zeros_like(phi)
        f = Cx * np.diff(phi, axis=2); out[:, :, :-1] -= f; out[:, :, 1:] += f
        f = Cy * np.diff(phi, axis=1); out[:, :-1, :] -= f; out[:, 1:, :] += f
        f = Cz * np.diff(phi, axis=0); out[:-1, :, :] -= f; out[1:, :, :] += f
        return out

    g = np.zeros(Bx.shape)               # outward normal flux B.n through each boundary node's dual face
    area = np.zeros(Bx.shape)
    yz = hz[:, None] * hy[None, :]; xz = hz[:, None] * hx[None, :]; xy = hy[:, None] * hx[None, :]
    g[:, :, -1] += Bx[:, :, -1] * yz; g[:, :, 0] -= Bx[:, :, 0] * yz
    g[:, -1, :] += By[:, -1, :] * xz; g[:, 0, :] -= By[:, 0, :] * xz
    g[-1, :, :] += Bz[-1, :, :] * xy; g[0, :, :] -= Bz[0, :, :] * xy
    area[:, :, -1] += yz; area[:, :, 0] += yz; area[:, -1, :] += xz; area[:, 0, :] += xz
    area[-1, :, :] += xy; area[0, :, :] += xy
    net, tot = g.sum(), np.abs(g).sum()
    g -= net * area / area.sum()         # compatibility: the discrete div B leaves a small net flux
    phi = np.zeros_like(g)
    r = g.copy(); zv = r / diag; p = zv.copy(); rz = (r * zv).sum(); gn = np.sqrt((g * g).sum())
    it, res = 0, 1.0
    for it in range(1, maxit + 1):
        Ap = A(p)
        alpha = rz / (p * Ap).sum()
        phi += alpha * p
        r -= alpha * Ap
        res = np.sqrt((r * r).sum()) / gn
        if res < tol:
            break
        zv = r / diag
        rz_new = (r * zv).sum()
        p = zv + (rz_new / rz) * p
        rz = rz_new
    Epot = 0.5 * (phi * A(phi)).sum()
    return Epot, it, res, net / max(tot, 1e-300)


def cmd_energy(a):
    os.makedirs(a.outdir, exist_ok=True)
    fig, (axm, axk) = plt.subplots(1, 2, figsize=(a.width, 0.42 * a.width))
    zc = a.zc_mm / L0_MM if a.zc_mm is not None else a.zc
    for lab, run, n, dip, col, ls in labelled(a.runs, a.sim):
        rows = []
        for t, f in pick_times(snapshots(run), a.times):
            (x, y, z), F = read_lattice(f, ("ρ", "u", "v", "w", "Bx", "By", "Bz"), ((-1e9, 1e9), (-1e9, 1e9), (zc - 1e-6, 1e9)))
            sl = (slice(None, None, a.stride),) * 3
            if a.stride > 1:                   # subsample, keeping the last plane of each direction
                idx = [np.unique(np.r_[np.arange(0, len(c), a.stride), len(c) - 1]) for c in (z, y, x)]
                sl = np.ix_(*idx)
                x, y, z = x[idx[2]], y[idx[1]], z[idx[0]]
            G = {k: F[k][sl] for k in F}
            V = dual_widths(z)[:, None, None] * dual_widths(y)[None, :, None] * dual_widths(x)[None, None, :]
            Em = 0.5 * ((G["Bx"]**2 + G["By"]**2 + G["Bz"]**2) * V).sum()
            Ek = 0.5 * (G["ρ"] * (G["u"]**2 + G["v"]**2 + G["w"]**2) * V).sum()
            Ep, it, res, net = potential_energy(x, y, z, G["Bx"], G["By"], G["Bz"], a.tol, a.maxit)
            rows.append((t, Em - Ep, Ek, Em, Ep))
            print(f"   {lab}: t = {t:8.3f}  E_free = {Em - Ep:10.4f}  E_kin = {Ek:9.4f}  (E_mag {Em:.4f}, E_pot {Ep:.4f};"
                  f" CG {it} its, residual {res:.1e}; net boundary flux {net:.1e}; z >= {z[0] * L0_MM:.3f} Mm)")
        r = np.array(rows)
        tag = re.sub(r"[^A-Za-z0-9]+", "", lab) or "run"
        np.savetxt(os.path.join(a.outdir, f"leake22_energy_{tag}.csv"), r, fmt="%.8g",
                   header=f"t E_free E_kin E_mag E_pot   (in E0, z >= {zc * L0_MM:.3f} Mm)")
        axm.plot(r[:, 0], r[:, 1], color=col, ls=ls, lw=1.6, label=lab)
        axk.plot(r[:, 0], r[:, 2], color=col, ls=ls, lw=1.6, label=lab)
    axm.set_xlabel("Time/t0"); axm.set_ylabel("Free Magnetic Energy/E0")
    axk.set_xlabel("Time/t0"); axk.set_ylabel("Kinetic Energy/E0")
    axm.legend(fontsize=8); axk.legend(fontsize=8)
    savefig(fig, os.path.join(a.outdir, f"leake22_fig9_energy.{a.format}"), a.dpi)


# ---------------------------------------------------------------------------------------------------------
# lorentz: Fig. 12 (right)
# ---------------------------------------------------------------------------------------------------------
def cmd_lorentz(a):
    n, dip = run_setup(a.run, a.sim)
    c, s = theta(n)
    snaps = pick_times(snapshots(a.run), a.times)
    zr = lim(a.zlim)
    m = a.margin
    fig, axs = plt.subplots(len(snaps), 1, figsize=(a.width, 0.55 * a.width * len(snaps)), squeeze=False)
    for ax, (t, f) in zip(axs[:, 0], snaps):
        coords, F = read_lattice(f, ("Bx", "By", "Bz"), ((-m, m), (-m, m), (-1e9, 1e9)))
        x, y, z = coords
        Bx, By, Bz = F["Bx"], F["By"], F["Bz"]
        tens = Bx * ddx(Bz, x, 2) + By * ddx(Bz, y, 1) + Bz * ddx(Bz, z, 0)       # (B.grad) B_z
        pres = -ddx(0.5 * (Bx**2 + By**2 + Bz**2), z, 0)                         # -d(B^2/2)/dz
        T, Pz, bx, by = column(coords, [tens, pres, Bx, By])
        zz = z * L0_MM
        ax.plot(zz, T + Pz, "k-", lw=1.2, label="Lorentz force in z")
        ax.plot(zz, T, "r--", lw=1.2, label="Magn. tension in z")
        ax.plot(zz, Pz, "b-.", lw=1.2, label="Magn. pressure in z")
        z1, z2, _ = rope_height(z, -s * bx + c * by, c * bx + s * by, (a.zmin / L0_MM, 1e9))
        if np.isfinite(z1):
            ax.axvline(z1 * L0_MM, color="k", ls=":", lw=0.8)
        ax.axhline(0, color="0.7", lw=0.6)
        if zr:
            ax.set_xlim(*zr)
            k = (zz >= zr[0]) & (zz <= zr[1])
            v = np.abs(np.r_[T[k] + Pz[k], T[k], Pz[k]]).max()
            ax.set_ylim(-1.1 * v, 1.1 * v)
        ax.set_xlabel("Z (Mm)"); ax.set_ylabel(r"LF$_z$")
        ax.set_title(rf"$t = {t:g}\,t_0$ (on $x = y = 0$)", fontsize=10, loc="left")
        ax.legend(fontsize=8, loc="lower right")
        print(f"   {f}: t = {t:g}, axial-field maximum at {z1 * L0_MM:.3f} Mm")
    savefig(fig, os.path.join(a.outdir, f"{a.name or 'leake22_fig12_lorentz'}.{a.format}"), a.dpi)


# ---------------------------------------------------------------------------------------------------------
# render: Figs. 1 (left), 3, 5
# ---------------------------------------------------------------------------------------------------------
def surface(m):
    """PolyData of a grid of line cells (pyvista >= 0.47 wants the algorithm named, older ones reject it)."""
    try:
        return m.extract_surface(algorithm=None)
    except TypeError:
        return m.extract_surface()


def cmd_render(a):
    try:
        import pyvista as pv
    except ImportError:
        raise SystemExit("render needs pyvista: pip install pyvista")
    pv.OFF_SCREEN = True
    n, dip = run_setup(a.run, a.sim)
    c, s = theta(n)
    os.makedirs(a.outdir, exist_ok=True)
    for t, f in pick_times(snapshots(a.run), a.times):
        (x, y, z), F = read_lattice(f, ["Bx", "By", "Bz"])
        X, Y, Z0, Z1 = x[-1], y[-1], z[0], z[-1]
        grid = pv.RectilinearGrid(x * L0_MM, y * L0_MM, z * L0_MM)
        grid.point_data["B"] = np.stack([F["Bx"].ravel(), F["By"].ravel(), F["Bz"].ravel()], axis=1)
        grid.point_data["Bz"] = F["Bz"].ravel() * B0_G
        if a.jb:
            jx, jy, jz = current((x, y, z), F)
            grid.point_data["JB"] = (np.sqrt(jx**2 + jy**2 + jz**2)
                                     / np.maximum(np.sqrt(F["Bx"]**2 + F["By"]**2 + F["Bz"]**2), 1e-300)).ravel()
            del jx, jy, jz
        z0, z1, dz = lim(a.seeds)
        S = np.array([(0.0, 0.0, zz) for zz in np.arange(z0, z1 + 1e-9, dz)])
        sl = grid.streamlines_from_source(pv.PolyData(S), vectors="B", integration_direction="both",
                                          initial_step_length=0.5, max_step_length=2.0, max_steps=a.max_steps,
                                          max_length=a.max_length, terminal_speed=1e-12)
        # connectivity (the paper's colors): a line end on a side boundary below --rope-zmax belongs to the
        # rope, any other end (bottom, top, a side above the surface) to the dipole. Each traced piece starts at
        # its seed, so its last point is a line end; a seed traced one way only keeps the seed as the other end.
        # The tracer stops up to a step short of the boundary (steps scale with the local cells), so a piece
        # that left the domain is extended along its last segment to the face it crosses.
        lo = np.array([x[0], y[0], Z0]) * L0_MM; hi = np.array([X, Y, Z1]) * L0_MM
        tolh = 1e-6 * hi[0]

        def exit_point(P):
            if len(P) < 2:
                return P[-1]
            q, d = P[-1], P[-1] - P[-2]
            ts = [((hi[k] if d[k] > 0 else lo[k]) - q[k]) / d[k] for k in range(3) if d[k] != 0]
            return q + max(0.0, min(ts)) * d if ts else q

        def rope_end(e):
            return (abs(e[0]) >= hi[0] - tolh or abs(e[1]) >= hi[1] - tolh) and e[2] < a.rope_zmax
        seed = np.asarray(sl.cell_data["SeedIds"]).astype(int)
        why = np.asarray(sl.cell_data["ReasonForTermination"]).astype(int)
        ends = {}
        for ci in range(sl.n_cells):
            P = sl.get_cell(ci).points
            if len(P):                        # 1 = OUT_OF_DOMAIN; other ends (out of steps, null) are interior
                ends.setdefault(seed[ci], []).append(why[ci] == 1 and rope_end(exit_point(P)))
        kind_of = {sd: (1 if len(e) == 2 and all(e) else (2 if any(e) else 0)) for sd, e in ends.items()}
        sl.cell_data["conn"] = np.array([kind_of.get(sd, 0) for sd in seed], float)
        box = lim(a.box)
        if box:
            sl = sl.clip_box(box, invert=False)
        if not isinstance(sl, pv.PolyData):
            sl = surface(sl)
        pl = pv.Plotter(off_screen=True, window_size=a.size)
        pl.set_background(a.background)
        colors = {0: a.dipole_color, 1: a.rope_color, 2: a.mixed_color}
        for k, col in colors.items():
            part = sl.extract_cells(np.nonzero(sl.cell_data["conn"] == k)[0])
            if part.n_cells == 0:
                continue
            part = surface(part)
            geo = part.tube(radius=a.tube_width) if a.tube_width > 0 else part
            pl.add_mesh(geo, color=col, smooth_shading=True)
            print(f"   {('dipole', 'rope', 'mixed')[k]}: {part.n_cells} line pieces")
        sf = grid.slice(normal="z", origin=(0.0, 0.0, a.surface_z))
        if box:
            sf = sf.clip_box(box, invert=False)
        pl.add_mesh(sf, scalars="Bz", cmap="gray", clim=(-a.bz_clim, a.bz_clim), opacity=a.surface_opacity,
                    scalar_bar_args=dict(title="BZ(G)", color="k" if a.background == "white" else "w"))
        if a.jb:
            nrm = (-s, c, 0.0)                # the cross-section y_r = 0 has normal y_r
            mid = grid.slice(normal=nrm, origin=(0.0, 0.0, 0.0))
            if box:
                mid = mid.clip_box(box, invert=False)
            pl.add_mesh(mid, scalars="JB", cmap="Reds", clim=(0, a.jb_clim), opacity=a.jb_opacity,
                        scalar_bar_args=dict(title="J/B", color="k" if a.background == "white" else "w"))
        if a.view == "top":
            pl.view_xy()
        else:
            pl.view_vector((-s, c, 0.0) if a.view == "end" else (c, s, 0.0), viewup=(0.0, 0.0, 1.0))
            pl.camera.azimuth = a.azimuth
            pl.camera.elevation = a.elevation
        pl.camera.zoom(a.zoom)
        if a.title:
            pl.add_text(f"t = {t:g} t0", font_size=12, color="k" if a.background == "white" else "w")
        out = os.path.join(a.outdir, f"{a.name or 'leake22_render'}_sim{n}_t{t:g}_{a.view}.png")
        pl.screenshot(out)
        pl.close()
        print(f" -> {out}")


# ---------------------------------------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(s, width=6.0):
        s.add_argument("--outdir", default="figs_leake2022")
        s.add_argument("--dpi", type=int, default=200)
        s.add_argument("--width", type=float, default=width, help="figure width in inches")
        s.add_argument("--format", default="png", help="png, pdf, svg, ...")
        s.add_argument("--sim", type=float, default=None, help="override the runs' :lk22_sim (theta = n pi/4)")

    s = sub.add_parser("ic", help="Fig. 1 (right)")
    s.add_argument("runs", nargs="*", help="LABEL=RUN: overlay the nodal T, rho on x = y = 0 at t = 0")
    s.add_argument("--cz-gradient", type=float, default=1.0, help=":lk22_cz_gradient of the deck")
    s.add_argument("--mf", type=float, default=1.0, help=":lk22_mf of the deck (p = rho T/m_f)")
    common(s)
    s.set_defaults(func=cmd_ic)

    s = sub.add_parser("xsec", help="cross-section halfway along the rope (Figs. 4, 6, 7, 12 left)")
    s.add_argument("run")
    s.add_argument("--var", choices=("Baxial", "JB"), default="Baxial")
    s.add_argument("--times", default=None)
    s.add_argument("--xr-max", type=float, default=40.0, help="half-width of the cut in Mm")
    s.add_argument("--zlim", default="-5,65", help="zmin,zmax in Mm")
    s.add_argument("--nsamp", type=int, default=801, help="samples along x_r")
    s.add_argument("--clim"); s.add_argument("--cmap")
    s.add_argument("--style", choices=("fig6", "fig4"), default="fig6", help="fig4: white-red |J|/|B| in 0-0.5")
    s.add_argument("--contours", type=int, default=0, help="N contours of log10|B_axial| (G) (Fig. 12 left)")
    s.add_argument("--contour-range", type=lambda t: tuple(float(v) for v in t.split(",")), default=(0.0, 2.5))
    s.add_argument("--photosphere", action="store_true", help="thin gray band at z = 0 (Fig. 4)")
    s.add_argument("--title", default=None, help="text before the time, e.g. 'Sim 4: '")
    s.add_argument("--name")
    common(s, 12.0)
    s.set_defaults(func=cmd_xsec)

    s = sub.add_parser("height", help="flux-rope height on x = y = 0 (Fig. 8)")
    s.add_argument("runs", nargs="+")
    s.add_argument("--times", default=None)
    s.add_argument("--zmin", type=float, default=-4.5, help="search above this height (Mm)")
    s.add_argument("--zmax", type=float, default=79.0, help="search below this height (Mm; top damping at 78.9)")
    s.add_argument("--margin", type=float, default=6.0, help="read |x|, |y| <= margin (L0) around the axis")
    s.add_argument("--persist", type=int, default=2, help="nodes either side on which a sign change must hold")
    common(s)
    s.set_defaults(func=cmd_height)

    s = sub.add_parser("energy", help="free magnetic and kinetic energy above z = 10 L0 (Fig. 9)")
    s.add_argument("runs", nargs="+")
    s.add_argument("--times", default=None)
    s.add_argument("--zc", type=float, default=10.0, help="bottom of the coronal volume in L0 (the paper: 10)")
    s.add_argument("--zc-mm", type=float, default=None, help="the same in Mm")
    s.add_argument("--stride", type=int, default=1, help="use every N-th node (faster potential-field solve)")
    s.add_argument("--tol", type=float, default=1e-7); s.add_argument("--maxit", type=int, default=20000)
    common(s, 11.0)
    s.set_defaults(func=cmd_energy)

    s = sub.add_parser("lorentz", help="vertical Lorentz force on x = y = 0 (Fig. 12 right)")
    s.add_argument("run")
    s.add_argument("--times", default="130,160")
    s.add_argument("--zlim", default="2,12", help="Mm")
    s.add_argument("--zmin", type=float, default=-4.5, help="search the rope axis above this height (Mm)")
    s.add_argument("--margin", type=float, default=6.0)
    s.add_argument("--name")
    common(s)
    s.set_defaults(func=cmd_lorentz)

    s = sub.add_parser("render", help="3D field lines colored by connectivity (Figs. 1 left, 3, 5)")
    s.add_argument("run")
    s.add_argument("--times", default=None)
    s.add_argument("--seeds", default="-4,88,1.5", help="z0,z1,dz (Mm) of the seeds on x = y = 0")
    s.add_argument("--view", choices=("end", "side", "top"), default="end",
                   help="end: along the rope axis (Figs. 1, 3 top); side: across it; top: from above (Fig. 3 bottom)")
    s.add_argument("--azimuth", type=float, default=-20.0); s.add_argument("--elevation", type=float, default=12.0)
    s.add_argument("--zoom", type=float, default=1.2)
    s.add_argument("--box", default="-47.2,47.2,-47.2,47.2,-5.1,90.2", help="display box in Mm")
    s.add_argument("--rope-zmax", type=float, default=0.0, help="a line end on a side boundary below this height "
                   "(Mm) belongs to the rope")
    s.add_argument("--dipole-color", default="#5b1f73"); s.add_argument("--rope-color", default="#e8c51a")
    s.add_argument("--mixed-color", default="#1f7a7a")
    s.add_argument("--surface-z", type=float, default=0.0); s.add_argument("--bz-clim", type=float, default=100.0)
    s.add_argument("--surface-opacity", type=float, default=0.8)
    s.add_argument("--jb", action="store_true", help="add the |J|/|B| cross-section y_r = 0 (Fig. 3 top)")
    s.add_argument("--jb-clim", type=float, default=0.5); s.add_argument("--jb-opacity", type=float, default=0.6)
    s.add_argument("--tube-width", type=float, default=0.12, help="field-line tube radius in Mm (0: lines)")
    s.add_argument("--max-steps", type=int, default=20000); s.add_argument("--max-length", type=float, default=2000.0)
    s.add_argument("--size", type=lambda t: tuple(int(v) for v in t.split(",")), default=(1400, 1100))
    s.add_argument("--background", default="black")
    s.add_argument("--title", action="store_true")
    s.add_argument("--sim", type=float, default=None)
    s.add_argument("--name")
    s.add_argument("--outdir", default="figs_leake2022")
    s.set_defaults(func=cmd_render)

    argv = sys.argv[1:]                  # let comma lists start with a minus sign: --zlim -5,65
    for i in range(len(argv) - 1, 0, -1):
        if argv[i - 1].startswith("--") and "=" not in argv[i - 1] and argv[i][:1] == "-" and \
           len(argv[i]) > 1 and (argv[i][1].isdigit() or argv[i][1] == "."):
            argv[i - 1:i + 1] = [argv[i - 1] + "=" + argv[i]]
    a = p.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
