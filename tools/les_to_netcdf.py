#!/usr/bin/env python3
"""
Convert Jexpresso LES output to LESICP NetCDF format.

Auto-detects the input type from the file extension:

  XZ cross-section statistics  (.vtu or .pvd)
  ─────────────────────────────────────────────
    Input : les_xz*.vtu  (single time-averaged file)
            les_xz*.pvd  (time-series, averages the last 30 min)
    Output: tables_{case}_statistics_{model}_{name}.nc

  3D instantaneous snapshot  (.pvtu)
  ────────────────────────────────────
    Input : iter_NNN.pvtu  (last time step)
    Output: tables_{case}_snapshots_{model}_{name}.nc

LESICP Cartesian output grid:
  x : [0, 10240] m   dx = 20 m  (513 pts)
  z : [5, 3000]  m   dz = 10 m  (300 pts)  — dims (z, x) or (y, z, x)

Both flat and ridge terrain produce the same Cartesian (x, z) output grid.
For ridge cases the z-coordinates stored in the VTU are already physical AGL
heights (Jexpresso warping function).  A two-step column-by-column
interpolation (z first, then x) is used, matching the official
interpolation_tables.py.  Output points below the terrain surface are NaN.

Usage:
  python3 tools/les_to_netcdf.py /path/to/input [output_folder]

Edit MODEL / NAME / CASE at the top of this file.
"""

# ── User configuration ────────────────────────────────────────────────────────
MODEL = "jexpresso"
NAME  = "marras"       # e.g. "rauchoecker"
CASE  = "flat_u10"     # e.g. "flat_u00", "ridge100_u10", "ridge1000_u00"

import sys
import os
import xml.etree.ElementTree as ET
import numpy as np
import pyvista as pv
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import map_coordinates
import xarray as xr

# ── LESICP Cartesian output grid (flat terrain only) ─────────────────────────
X_OUT = np.arange(0.0, 10241.0, 20.0)         # 0, 20, …, 10240  →  513 pts
Z_OUT = np.arange(5. + 1e-6, 3000., 10.0)     # matches official interpolation_tables.py

TAVG_WINDOW  = 1800.0   # 30-minute averaging window for PVD statistics [s]
ROUND_DIGITS = 4        # decimal places for coordinate deduplication

# Fields read from PVTU snapshots (in output order)
SNAP_FIELDS = ['u', 'v', 'w', 'θ', 'p', 'ρ']
_NC_RENAME  = {'θ': 't', 'ρ': 'rho'}


def _nc_name(key):
    return _NC_RENAME.get(key, key)


# ════════════════════════════════════════════════════════════════════════════
# Shared helpers
# ════════════════════════════════════════════════════════════════════════════

def read_pvd(pvd_path):
    """Return sorted list of (time, abs_path), deduplicated by filename."""
    tree    = ET.parse(pvd_path)
    pvd_dir = os.path.dirname(os.path.abspath(pvd_path))
    seen, entries = set(), []
    for ds in tree.getroot().iter("DataSet"):
        rel = ds.attrib["file"]
        if rel in seen:
            continue
        seen.add(rel)
        entries.append((float(ds.attrib["timestep"]),
                        os.path.join(pvd_dir, rel)))
    entries.sort(key=lambda e: e[0])
    return entries


def _detect_grid_mode(pts):
    """
    Return 'rectilinear' or 'terrain'.

    For a flat mesh every (x, y) column shares the same z-axis, so the
    number of globally unique z values ≈ N / (Nx * Ny).
    For a terrain-following mesh z differs per column → unique_z >> that.
    """
    xr_ = np.round(pts[:, 0] - pts[:, 0].min(), ROUND_DIGITS)
    yr_ = np.round(pts[:, 1] - pts[:, 1].min(), ROUND_DIGITS)
    zr_ = np.round(pts[:, 2],                   ROUND_DIGITS)
    Nx  = len(np.unique(xr_))
    Ny  = len(np.unique(yr_))
    return ('rectilinear'
            if abs(len(np.unique(zr_)) - len(pts) // (Nx * Ny)) <= 1
            else 'terrain')



# ════════════════════════════════════════════════════════════════════════════
# XZ cross-section statistics  (VTU / PVD input)
# ════════════════════════════════════════════════════════════════════════════

def build_xz_grid_info(mesh):
    """
    Analyse the 2-D xz mesh geometry (y=0 plane).

    Rectilinear (flat terrain)
      mode='rectilinear', x_sem, z_sem, flat_ix, flat_iz

    Terrain-following (ridge)
      mode='terrain', x_sem, z_cols (Nx, Nz), sort_idx (Nx, Nz)
      z_cols[i, k] = physical AGL height at (x_sem[i], sigma_level_k)
    """
    pts     = mesh.points
    x_raw   = pts[:, 0]
    z_raw   = pts[:, 2]
    x_shift = x_raw - x_raw.min()

    x_sem_r  = np.round(x_shift, ROUND_DIGITS)
    z_sem_r  = np.round(z_raw,   ROUND_DIGITS)
    x_sem    = np.unique(x_sem_r)
    z_unique = np.unique(z_sem_r)
    Nx, N    = len(x_sem), len(x_raw)

    if Nx * len(z_unique) == N:
        flat_ix = np.searchsorted(x_sem,    x_sem_r)
        flat_iz = np.searchsorted(z_unique, z_sem_r)
        return dict(mode='rectilinear',
                    x_sem=x_sem, z_sem=z_unique,
                    flat_ix=flat_ix, flat_iz=flat_iz)
    else:
        Nz       = N // Nx
        flat_ix  = np.searchsorted(x_sem, x_sem_r)
        z_cols   = np.empty((Nx, Nz))
        sort_idx = np.empty((Nx, Nz), dtype=int)
        for i in range(Nx):
            mask        = flat_ix == i
            order       = np.argsort(z_raw[mask])
            z_cols[i]   = z_raw[mask][order]
            sort_idx[i] = np.where(mask)[0][order]
        h_terrain_out = np.interp(X_OUT, x_sem, z_cols[:, 0])
        return dict(mode='terrain',
                    x_sem=x_sem,
                    z_cols=z_cols, sort_idx=sort_idx,
                    h_terrain_out=h_terrain_out)


def _interp_xz_cartesian(vals, grid_info):
    """
    Rectilinear path: interpolate onto fixed (X_OUT × Z_OUT) Cartesian grid.
    Returns (Nz_out, Nx_out) float64.
    """
    x_sem, z_sem = grid_info['x_sem'], grid_info['z_sem']
    Nx_out, Nz_out = len(X_OUT), len(Z_OUT)
    field_2d = np.empty((len(x_sem), len(z_sem)))
    field_2d[grid_info['flat_ix'], grid_info['flat_iz']] = vals
    XX, ZZ = np.meshgrid(X_OUT, Z_OUT, indexing='ij')
    tgt    = np.column_stack([XX.ravel(), ZZ.ravel()])
    interp = RegularGridInterpolator((x_sem, z_sem), field_2d,
                                      method='linear',
                                      bounds_error=False, fill_value=None)
    return interp(tgt).reshape(Nx_out, Nz_out).T   # (Nz_out, Nx_out)


def _interp_xz_terrain(vals, grid_info):
    """
    Terrain path: column-by-column z-interp then x-interp onto Cartesian grid.
    Mirrors official interpolation_tables.py (griddata per x-column).
    Returns (Nz_out, Nx_out) float64.  Points below terrain surface are NaN.
    """
    x_sem         = grid_info['x_sem']
    z_cols        = grid_info['z_cols']
    sort_idx      = grid_info['sort_idx']
    h_terrain_out = grid_info['h_terrain_out']
    Nx_sem        = len(x_sem)
    Nz_out        = len(Z_OUT)
    Nx_out        = len(X_OUT)

    # Step 1: z-interpolation per x column → (Nx_sem, Nz_out)
    field_z = np.empty((Nx_sem, Nz_out))
    for i in range(Nx_sem):
        field_z[i] = np.interp(Z_OUT, z_cols[i], vals[sort_idx[i]],
                                left=np.nan, right=np.nan)

    # Step 2: x-interpolation per Z_OUT level → (Nx_out, Nz_out)
    result = np.full((Nx_out, Nz_out), np.nan)
    for j in range(Nz_out):
        col   = field_z[:, j]
        valid = np.isfinite(col)
        if valid.sum() < 2:
            continue
        interp_x      = np.interp(X_OUT, x_sem[valid], col[valid],
                                  left=np.nan, right=np.nan)
        above         = Z_OUT[j] > h_terrain_out
        result[:, j]  = np.where(above, interp_x, np.nan)

    return result.T   # (Nz_out, Nx_out)


def run_xz_statistics(inpath, out_dir):
    """Process VTU (single) or PVD (time-averaged) → statistics NetCDF."""
    content = "statistics"
    out_nc  = os.path.join(out_dir,
                           f"tables_{CASE}_{content}_{MODEL}_{NAME}.nc")
    _print_header(content, out_nc)

    single_vtu = inpath.lower().endswith(".vtu")
    if single_vtu:
        print(f"Mode: single VTU  →  {inpath}")
        selected = [(0.0, inpath)]
    else:
        print(f"Mode: PVD series  →  {inpath}")
        entries  = read_pvd(inpath)
        if not entries:
            sys.exit("ERROR: no entries in PVD")
        t_max    = entries[-1][0]
        t_cutoff = t_max - TAVG_WINDOW
        selected = [(t, p) for t, p in entries if t >= t_cutoff]
        print(f"  All snapshots : {len(entries)},  t ∈ [{entries[0][0]:.1f}, {t_max:.1f}] s")
        print(f"  Avg window    : t ≥ {t_cutoff:.1f} s  →  {len(selected)} snapshots")
        if not selected:
            sys.exit("ERROR: no snapshots in averaging window")

    # Build grid info from first snapshot (coordinates constant across time)
    mesh0     = pv.read(selected[0][1])
    grid_info = build_xz_grid_info(mesh0)
    mode      = grid_info['mode']
    varnames  = list(mesh0.point_data.keys())

    print(f"\n  Grid mode  : {mode}  ({len(grid_info['x_sem'])} x-cols)")
    Nz_acc = len(Z_OUT)
    Nx_acc = len(X_OUT)
    if mode == 'rectilinear':
        print(f"  z_sem ∈ [{grid_info['z_sem'].min():.1f}, {grid_info['z_sem'].max():.1f}] m")
    else:
        print(f"  z_cols ∈ [{grid_info['z_cols'].min():.1f}, {grid_info['z_cols'].max():.1f}] m")
        print(f"  Terrain h_max : {grid_info['z_cols'][:, 0].max():.1f} m")
    print(f"  Target     : x({Nx_acc}) × z({Nz_acc})  →  dims (z, x)")
    print(f"  Variables  : {varnames}\n")

    # Accumulate time average
    accum = {v: np.zeros((Nz_acc, Nx_acc)) for v in varnames}
    n_ok  = 0
    for i, (t, fpath) in enumerate(selected):
        print(f"  [{i+1:3d}/{len(selected)}] t={t:8.1f} s  {os.path.basename(fpath)}")
        try:
            mesh = pv.read(fpath)
        except Exception as e:
            print(f"    WARNING: {e}"); continue
        for v in varnames:
            vals = np.asarray(mesh.point_data[v], dtype=np.float64).ravel()
            if mode == 'rectilinear':
                accum[v] += _interp_xz_cartesian(vals, grid_info)
            else:
                accum[v] += _interp_xz_terrain(vals, grid_info)
        n_ok += 1

    if n_ok == 0:
        sys.exit("ERROR: no snapshots could be read")
    print(f"\nAveraged over {n_ok} snapshots.")
    for v in varnames:
        accum[v] /= n_ok

    # Build dataset
    print(f"\nWriting: {out_nc}")
    # Both flat and ridge use Cartesian (z, x) output — NaN below terrain for ridge
    ds = xr.Dataset(
        {v: xr.DataArray(accum[v].astype(np.float32), dims=["z", "x"],
                         attrs={"long_name": v})
         for v in varnames},
        coords={
            "z": xr.DataArray(Z_OUT.astype(np.float32), dims=["z"],
                               attrs={"units": "m", "long_name": "height AGL"}),
            "x": xr.DataArray(X_OUT.astype(np.float32), dims=["x"],
                               attrs={"units": "m", "long_name": "streamwise distance"}),
        },
        attrs=dict(title="LES xz cross-section statistics",
                   model=MODEL, name=NAME, case=CASE, content=content,
                   grid_mode=mode, n_samples=int(n_ok), source=inpath,
                   Conventions="LESICP"),
    )
    ds.to_netcdf(out_nc, engine="scipy")
    _print_footer(out_nc, varnames)


# ════════════════════════════════════════════════════════════════════════════
# 3-D instantaneous snapshot  (PVTU input)
# ════════════════════════════════════════════════════════════════════════════

def _build_3d_rectilinear(pts, fields_dict):
    """
    De-duplicate GLL nodes and map onto a rectilinear (Nx, Ny, Nz) grid.
    Returns x_ax, y_ax, z_ax (1-D, shifted to 0), arrays_3d dict.
    """
    xr_ = np.round(pts[:, 0] - pts[:, 0].min(), ROUND_DIGITS)
    yr_ = np.round(pts[:, 1] - pts[:, 1].min(), ROUND_DIGITS)
    zr_ = np.round(pts[:, 2],                   ROUND_DIGITS)
    x_ax = np.unique(xr_); y_ax = np.unique(yr_); z_ax = np.unique(zr_)
    Nx, Ny, Nz = len(x_ax), len(y_ax), len(z_ax)
    print(f"  Structured grid : {Nx} × {Ny} × {Nz}  (raw pts: {len(xr_)})")
    ix  = np.searchsorted(x_ax, xr_)
    iy  = np.searchsorted(y_ax, yr_)
    iz  = np.searchsorted(z_ax, zr_)
    lin = (ix * Ny + iy) * Nz + iz
    uniq_lin, inv = np.unique(lin, return_inverse=True)
    cnt = np.bincount(inv, minlength=len(uniq_lin)).astype(np.float64)
    arrays_3d = {}
    for name, vals in fields_dict.items():
        acc  = np.bincount(inv, weights=vals.astype(np.float64),
                           minlength=len(uniq_lin))
        full = np.zeros(Nx * Ny * Nz)
        full[uniq_lin] = acc / cnt
        arrays_3d[name] = full.reshape(Nx, Ny, Nz)
    return x_ax, y_ax, z_ax, arrays_3d


def _build_3d_terrain(pts, fields_dict):
    """
    De-duplicate GLL nodes and map onto a curvilinear structured (Nx, Ny, Nz)
    grid where z_cols[i, j, k] holds physical AGL height for column (i, j),
    level k.  Returns x_ax, y_ax, z_cols (Nx, Ny, Nz), arrays_3d dict.
    """
    xr_ = np.round(pts[:, 0] - pts[:, 0].min(), ROUND_DIGITS)
    yr_ = np.round(pts[:, 1] - pts[:, 1].min(), ROUND_DIGITS)
    zr_ = pts[:, 2]
    x_ax = np.unique(np.round(xr_, ROUND_DIGITS))
    y_ax = np.unique(np.round(yr_, ROUND_DIGITS))
    Nx, Ny = len(x_ax), len(y_ax)

    # Sort by (x, y, z) — groups each (x,y) column, z ascending within it
    sk = np.lexsort((zr_, np.round(yr_, ROUND_DIGITS), np.round(xr_, ROUND_DIGITS)))
    xs = np.round(xr_[sk], ROUND_DIGITS)
    ys = np.round(yr_[sk], ROUND_DIGITS)
    zs = zr_[sk]

    # Detect duplicate (x, y, z_round) triplets and average them
    zr_s   = np.round(zs, ROUND_DIGITS)
    is_new = np.ones(len(xs), dtype=bool)
    is_new[1:] = (xs[1:] != xs[:-1]) | (ys[1:] != ys[:-1]) | (zr_s[1:] != zr_s[:-1])
    gid  = np.cumsum(is_new) - 1
    cnt  = np.bincount(gid).astype(np.float64)
    Nz   = int(gid[-1] + 1) // (Nx * Ny)
    print(f"  Structured grid : {Nx} × {Ny} × {Nz}  (raw pts: {len(xr_)})")

    z_cols    = zs[is_new].reshape(Nx, Ny, Nz)
    arrays_3d = {}
    for name, vals in fields_dict.items():
        acc = np.bincount(gid, weights=vals[sk].astype(np.float64))
        arrays_3d[name] = (acc / cnt).reshape(Nx, Ny, Nz)
    return x_ax, y_ax, z_cols, arrays_3d


def _interp_3d_rectilinear(arrays_3d, x_ax, y_ax, z_ax, present):
    """
    Interpolate rectilinear 3D arrays to the Cartesian output grid.
    Returns dict {nc_key: (Ny, Nz_out, Nx_out) float32} — dims (y, z, x).
    """
    Nxo, Ny, Nzo = len(X_OUT), len(y_ax), len(Z_OUT)
    xi  = np.interp(X_OUT, x_ax, np.arange(len(x_ax)))
    zi  = np.interp(Z_OUT, z_ax, np.arange(len(z_ax)))
    XI2, ZI2 = np.meshgrid(xi, zi, indexing='ij')
    c2  = np.array([XI2.ravel(), ZI2.ravel()])
    result = {}
    for key in present:
        arr3d = arrays_3d[key]
        out   = np.empty((Nxo, Ny, Nzo), dtype=np.float32)
        for iy in range(Ny):
            out[:, iy, :] = map_coordinates(
                arr3d[:, iy, :], c2, order=1, mode='nearest', prefilter=False,
            ).reshape(Nxo, Nzo)
        nc_key         = _nc_name(key)
        result[nc_key] = out.transpose(1, 2, 0)   # (Ny, Nzo, Nxo) — dims (y, z, x)
        print(f"  {nc_key:4s}  ∈ [{out.min():.4g}, {out.max():.4g}]")
    return result


def _interp_3d_terrain(arrays_3d, z_cols, x_ax, y_ax, present):
    """
    Terrain path: two-step column-by-column interpolation to Cartesian grid.
    Mirrors official interpolation_tables.py for ridge cases.
    Returns dict {nc_key: (Ny, Nz_out, Nx_out) float32} — dims (y, z, x).
    Points below the terrain surface are NaN.
    """
    Nx, Ny, Nz = z_cols.shape
    Nz_out     = len(Z_OUT)
    Nx_out     = len(X_OUT)

    # Terrain surface height at X_OUT (ridge is 2-D: same in y, use y=0)
    h_terrain_out = np.interp(X_OUT, x_ax, z_cols[:, 0, 0])

    result = {}
    for key in present:
        f3d = arrays_3d[key]   # (Nx, Ny, Nz)
        out = np.empty((Nx_out, Ny, Nz_out), dtype=np.float32)
        for iy in range(Ny):
            # Step 1: z-interp per x column → (Nx, Nz_out)
            field_z = np.empty((Nx, Nz_out))
            for i in range(Nx):
                field_z[i] = np.interp(Z_OUT, z_cols[i, iy], f3d[i, iy],
                                        left=np.nan, right=np.nan)
            # Step 2: x-interp per z level → (Nx_out, Nz_out)
            for j in range(Nz_out):
                col   = field_z[:, j]
                valid = np.isfinite(col)
                if valid.sum() < 2:
                    out[:, iy, j] = np.nan
                    continue
                interp_x      = np.interp(X_OUT, x_ax[valid], col[valid],
                                          left=np.nan, right=np.nan)
                above         = Z_OUT[j] > h_terrain_out
                out[:, iy, j] = np.where(above, interp_x, np.nan)

        nc_key         = _nc_name(key)
        result[nc_key] = out.transpose(1, 2, 0)   # (Ny, Nz_out, Nx_out) — dims (y, z, x)
        print(f"  {nc_key:4s}  ∈ [{np.nanmin(out):.4g}, {np.nanmax(out):.4g}]")
    return result


def run_3d_snapshot(pvtu_path, out_dir):
    """Process PVTU → 3D snapshot NetCDF."""
    content = "snapshots"
    out_nc  = os.path.join(out_dir,
                           f"tables_{CASE}_{content}_{MODEL}_{NAME}.nc")
    _print_header(content, out_nc)

    print(f"\nReading: {pvtu_path}")
    mesh    = pv.read(pvtu_path)
    print(f"  Points  : {mesh.n_points}")
    print(f"  Fields  : {list(mesh.point_data.keys())}")
    pts     = mesh.points
    present = [k for k in SNAP_FIELDS if k in mesh.point_data]
    missing = [k for k in SNAP_FIELDS if k not in mesh.point_data]
    if missing:
        print(f"  Skipping (not in file): {missing}")
    fields  = {k: np.asarray(mesh.point_data[k], dtype=np.float64) for k in present}
    del mesh

    mode = _detect_grid_mode(pts)
    print(f"\nBuilding structured 3D grid  [{mode}]...")
    if mode == 'rectilinear':
        x_ax, y_ax, z_ax, arrays_3d = _build_3d_rectilinear(pts, fields)
        del pts, fields
        print(f"  x ∈ [{x_ax.min():.1f}, {x_ax.max():.1f}] m  ({len(x_ax)} pts)")
        print(f"  y ∈ [{y_ax.min():.1f}, {y_ax.max():.1f}] m  ({len(y_ax)} pts)")
        print(f"  z ∈ [{z_ax.min():.1f}, {z_ax.max():.1f}] m  ({len(z_ax)} pts)")
        Nxo, Ny, Nzo = len(X_OUT), len(y_ax), len(Z_OUT)
        print(f"  Target : x({Nxo}) × y({Ny}) × z({Nzo})  →  dims (y, z, x)\n")
    else:
        x_ax, y_ax, z_cols, arrays_3d = _build_3d_terrain(pts, fields)
        del pts, fields
        print(f"  x ∈ [{x_ax.min():.1f}, {x_ax.max():.1f}] m  ({len(x_ax)} pts)")
        print(f"  y ∈ [{y_ax.min():.1f}, {y_ax.max():.1f}] m  ({len(y_ax)} pts)")
        print(f"  z_terrain ∈ [{z_cols.min():.1f}, {z_cols.max():.1f}] m")
        print(f"  Terrain h_max : {z_cols[:, 0, 0].max():.1f} m")
        Nxo, Ny, Nzo = len(X_OUT), len(y_ax), len(Z_OUT)
        print(f"  Target : x({Nxo}) × y({Ny}) × z({Nzo})  →  dims (y, z, x)\n")

    print("Interpolating...")
    if mode == 'rectilinear':
        data_dict = _interp_3d_rectilinear(arrays_3d, x_ax, y_ax, z_ax, present)
    else:
        data_dict = _interp_3d_terrain(arrays_3d, z_cols, x_ax, y_ax, present)
    del arrays_3d

    y_shifted = (y_ax - y_ax.min()).astype(np.float32)

    # Both flat and ridge use Cartesian (y, z, x) output — NaN below terrain for ridge
    ds = xr.Dataset(
        {nc_key: xr.DataArray(arr, dims=["y", "z", "x"],
                              attrs={"long_name": nc_key})
         for nc_key, arr in data_dict.items()},
        coords={
            "y": xr.DataArray(y_shifted, dims=["y"],
                               attrs={"units": "m", "long_name": "spanwise distance"}),
            "z": xr.DataArray(Z_OUT.astype(np.float32), dims=["z"],
                               attrs={"units": "m", "long_name": "height AGL"}),
            "x": xr.DataArray(X_OUT.astype(np.float32), dims=["x"],
                               attrs={"units": "m", "long_name": "streamwise distance"}),
        },
        attrs=dict(title="3D instantaneous snapshot",
                   model=MODEL, name=NAME, case=CASE, content=content,
                   grid_mode=mode, source=pvtu_path, Conventions="LESICP"),
    )

    print(f"\nWriting: {out_nc}")
    ds.to_netcdf(out_nc, engine="scipy")
    _print_footer(out_nc, list(data_dict.keys()))


# ════════════════════════════════════════════════════════════════════════════
# Entry point
# ════════════════════════════════════════════════════════════════════════════

def _print_header(content, out_nc):
    print(f"model   : {MODEL}")
    print(f"name    : {NAME}")
    print(f"case    : {CASE}")
    print(f"content : {content}")
    print(f"outname : {os.path.basename(out_nc)}")


def _print_footer(out_nc, varnames):
    print(f"\nDone.")
    print(f"  Output : {out_nc}")
    print(f"  Vars   : {', '.join(varnames)}")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(0)

    inpath = sys.argv[1]
    if not os.path.isfile(inpath):
        sys.exit(f"ERROR: file not found: {inpath}")

    out_dir = sys.argv[2] if len(sys.argv) >= 3 \
              else os.path.dirname(os.path.abspath(inpath))
    os.makedirs(out_dir, exist_ok=True)

    ext = os.path.splitext(inpath)[1].lower()
    if ext == ".pvtu":
        run_3d_snapshot(inpath, out_dir)
    elif ext in (".vtu", ".pvd"):
        run_xz_statistics(inpath, out_dir)
    else:
        sys.exit(f"ERROR: unrecognised extension '{ext}'. "
                 f"Expected .vtu, .pvd, or .pvtu")


if __name__ == "__main__":
    main()
