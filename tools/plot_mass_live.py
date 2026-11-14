#!/usr/bin/env python3
"""
Live plot of a Jexpresso mass diagnostic (mass.csv, written by :lmass_diagnostic)
while the simulation runs.

    python3 tools/plot_mass_live.py <mass.csv | run output dir | case output root> [options]

The argument can be the CSV itself, the run's output directory, or a parent
such as output/ShallowWater/ReducedGravity_entrain: the newest mass.csv below
it is followed. The file is read incrementally (only the new lines on each
refresh), so it stays cheap on multi-year runs.

Options:
  --interval S   refresh period in seconds (default 5)
  --png FILE     no window: rewrite FILE at every refresh (headless / over ssh)
  --once         draw once and exit (with --png: write FILE once)

Top panel: rel_change = mass/mass0 - 1 against model time in days; bottom: the
absolute change mass - mass0 (units of the integral: m^3 for shallow water),
with mass0 in the title. Close the window (or Ctrl-C) to stop. Needs numpy and
matplotlib.
"""
import argparse, glob, os, sys, time

import numpy as np


def find_csv(path):
    if os.path.isfile(path):
        return path
    direct = os.path.join(path, "mass.csv")
    if os.path.isfile(direct):
        return direct
    cands = glob.glob(os.path.join(path, "**", "mass.csv"), recursive=True)
    return max(cands, key=os.path.getmtime) if cands else None


class Follower:
    """Reads the new complete lines of a growing CSV on each call."""

    def __init__(self, path):
        self.path, self.pos, self.rest = path, 0, ""
        self.t, self.mass, self.rel = [], [], []

    def poll(self):
        try:
            with open(self.path, "r") as f:
                f.seek(self.pos)
                chunk = f.read()
                self.pos = f.tell()
        except FileNotFoundError:
            return 0
        text = self.rest + chunk
        lines = text.split("\n")
        self.rest = lines.pop()                      # an incomplete last line waits for the next poll
        n = 0
        for line in lines:
            if not line or line.startswith("#") or line.startswith("step"):
                continue
            parts = line.split(",")
            if len(parts) < 5:
                continue
            try:
                self.t.append(float(parts[2])); self.mass.append(float(parts[3])); self.rel.append(float(parts[4]))
                n += 1
            except ValueError:
                continue
        return n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("path")
    ap.add_argument("--interval", type=float, default=5.0)
    ap.add_argument("--png", default=None)
    ap.add_argument("--once", action="store_true")
    a = ap.parse_args()

    import matplotlib
    if a.png:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    csv = find_csv(a.path)
    while csv is None:                               # the run may not have started writing yet
        print(f"waiting for a mass.csv under {a.path} ...", flush=True)
        time.sleep(a.interval)
        csv = find_csv(a.path)
    print(f"following {csv}", flush=True)

    fol = Follower(csv)
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    l0, = ax0.plot([], [], lw=1.2)
    l1, = ax1.plot([], [], lw=1.2, color="C1")
    ax0.axhline(0.0, color="k", lw=0.6)
    ax0.set_ylabel("relative change\nmass / mass0 - 1")
    ax1.axhline(0.0, color="k", lw=0.6)
    ax1.set_ylabel("absolute change\nmass - mass0")
    ax1.set_xlabel("t [days]")
    for ax in (ax0, ax1):
        ax.grid(alpha=0.3)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3), useOffset=False)
    fig.suptitle(os.path.join(*os.path.abspath(csv).split(os.sep)[-3:]), fontsize=10)
    if not a.png:
        plt.ion()
        plt.show(block=False)

    while True:
        new = fol.poll()
        if new and fol.t:
            t = np.asarray(fol.t)
            l0.set_data(t, np.asarray(fol.rel))
            l1.set_data(t, np.asarray(fol.mass) - fol.mass[0])
            for ax in (ax0, ax1):
                ax.relim(); ax.autoscale_view()
            ax0.set_title(f"t = {t[-1]:.2f} d   rel. change = {fol.rel[-1]:+.3e}   "
                          f"mass0 = {fol.mass[0]:.6e}   ({len(t)} samples)", fontsize=9)
            fig.tight_layout()
            if a.png:
                fig.savefig(a.png, dpi=110)
        if a.once:
            if not a.png:
                plt.show(block=True)
            return
        if a.png:
            time.sleep(a.interval)
        else:
            if not plt.fignum_exists(fig.number):
                return
            fig.canvas.draw_idle()
            plt.pause(a.interval)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(0)
