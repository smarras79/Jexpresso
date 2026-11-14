#!/usr/bin/env python3
"""
Plots for rg2d_clone.jl runs.

    python3 plot_rg2d.py --logs deck.log nofilter.log ... [--dumps deck nofilter ...] [--out figs]

--logs : log files written by rg2d_clone.jl (stdout); gives min H(t) and max |u|(t)
--dumps: dump prefixes passed as dump=<prefix> to the clone; <prefix>_prev.dat is the
         last finite state before the blow-up, <prefix>_last.dat the state at the stop
"""
import argparse, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_log(fname):
    rows = []
    with open(fname) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.replace("(", " ").replace(")", " ").replace(",", " ").split()
            # t hmin x y umax Hthere maxHu nthin nclip mass
            rows.append([float(p[0]), float(p[1]), float(p[4]), float(p[6]), int(p[7]), int(p[8])])
    r = np.array(rows)
    return dict(t=r[:, 0], hmin=r[:, 1], umax=r[:, 2], hu=r[:, 3], nthin=r[:, 4], nclip=r[:, 5])


def read_dump(fname):
    d = np.loadtxt(fname)
    with open(fname) as f:
        t = float(f.readline().split("=")[1].split()[0])
    return t, d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5], d[:, 6]


def plot_series(logs, out, dxmin=13.81e3, dt=600.0):
    fig, ax = plt.subplots(2, 1, figsize=(8, 6.5), sharex=True)
    for fname in logs:
        L = read_log(fname)
        lab = os.path.basename(fname).replace(".log", "")
        ax[0].semilogy(L["t"], np.maximum(L["hmin"], 1e-4), label=lab)
        ax[1].semilogy(L["t"], np.maximum(L["umax"], 1e-4), label=lab)
    ax[0].set_ylabel("min H  [m]"); ax[0].grid(alpha=.3); ax[0].legend(fontsize=8)
    ax[0].axhline(0.1, color="k", ls=":", lw=.8); ax[0].text(0.5, 0.12, "H_WIND_FULL (deck) = 0.1 m", fontsize=7)
    ax[0].axhline(1e-3, color="k", ls="--", lw=.8); ax[0].text(0.5, 1.3e-3, "clip / thin film 1 mm", fontsize=7)
    ax[1].set_ylabel("max |u|  [m/s]"); ax[1].grid(alpha=.3)
    ax[1].axhline(dxmin / dt, color="r", ls="--", lw=.8)
    ax[1].text(0.5, dxmin / dt * 1.2, f"dx_min/dt = {dxmin/dt:.0f} m/s (CFL = 1)", color="r", fontsize=7)
    ax[1].set_xlabel("t  [days]")
    fig.suptitle("ReducedGravity clone: minimum thickness and maximum velocity")
    fig.tight_layout(); fig.savefig(out, dpi=130); print("wrote", out)


def plot_maps(dumpfile, out, title, dxmin=13.81e3, dt=600.0):
    t, x, y, H, Hu, Hv, u, v = read_dump(dumpfile)
    spd = np.hypot(u, v)
    nx = len(np.unique(np.round(x, 3))); ny = len(np.unique(np.round(y, 3)))
    order = np.lexsort((x, y))
    X = x[order].reshape(ny, nx) / 1e3; Y = y[order].reshape(ny, nx) / 1e3
    Hm = H[order].reshape(ny, nx); Sm = spd[order].reshape(ny, nx)
    fig, ax = plt.subplots(1, 2, figsize=(12, 5.2))
    c0 = ax[0].pcolormesh(X, Y, np.log10(np.maximum(Hm, 1e-3)), shading="nearest", cmap="viridis")
    fig.colorbar(c0, ax=ax[0], label="log10 H [m]")
    ax[0].set_title(f"H  (t = {t/86400:.1f} d, min = {H.min():.2e} m)")
    c1 = ax[1].pcolormesh(X, Y, np.log10(np.maximum(Sm, 1e-4)), shading="nearest", cmap="magma")
    fig.colorbar(c1, ax=ax[1], label="log10 |u| [m/s]")
    i = int(np.argmax(spd))
    ax[1].plot(x[i] / 1e3, y[i] / 1e3, "c+", ms=14, mew=2)
    ax[1].set_title(f"|u|  (max = {spd.max():.2e} m/s at H = {H[i]:.2e} m)")
    for a in ax:
        a.set_xlabel("x [km]"); a.set_ylabel("y [km]"); a.set_aspect("equal")
    fig.suptitle(title); fig.tight_layout(); fig.savefig(out, dpi=130); print("wrote", out)
    # |u| against H, every node: the thin-layer velocities against the Ekman line and the CFL line
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.loglog(np.maximum(H, 1e-4), np.maximum(spd, 1e-6), ".", ms=3, alpha=.5)
    Hs = np.logspace(-3, 2.3, 100)
    ax.loglog(Hs, 0.36 / Hs, "r--", lw=1, label="u = (tau0/(rho f0)) / H  (Ekman transport 0.36 m^2/s)")
    ax.axhline(dxmin / dt, color="k", ls=":", lw=1, label="dx_min/dt (CFL = 1)")
    ax.set_xlabel("H [m]"); ax.set_ylabel("|u| [m/s]"); ax.grid(alpha=.3); ax.legend(fontsize=7)
    ax.set_title(title + "\n|u| vs H, all nodes")
    out2 = out.replace(".png", "_uH.png")
    fig.tight_layout(); fig.savefig(out2, dpi=130); print("wrote", out2)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--logs", nargs="*", default=[])
    ap.add_argument("--dumps", nargs="*", default=[], help="dump prefixes (dump=<prefix> of the clone)")
    ap.add_argument("--out", default=".")
    ap.add_argument("--dxmin", type=float, default=13.81e3, help="smallest LGL node spacing [m]")
    ap.add_argument("--dt", type=float, default=600.0)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.logs:
        plot_series(a.logs, os.path.join(a.out, "series.png"), a.dxmin, a.dt)
    for pre in a.dumps:
        for suffix in ("_prev", "_last"):
            f = pre + suffix + ".dat"
            if os.path.exists(f):
                plot_maps(f, os.path.join(a.out, os.path.basename(pre) + suffix + "_maps.png"),
                          os.path.basename(pre) + suffix, a.dxmin, a.dt)
