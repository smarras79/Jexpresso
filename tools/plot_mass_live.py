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
  --area A       domain area, for files written before the writer recorded it
                 (e.g. 4e12 for the 2000 km x 2000 km ReducedGravity basin)

Top panel: the relative change rel_change = (mass - mass0)/mass0 against model
time in days. Bottom: the
change as a domain mean, (mass - mass0)/area -- for shallow water the
equivalent change of the layer thickness in m -- when the area is known (the
"# area =" header line, or --area); otherwise the absolute change mass - mass0.
Close the window (or Ctrl-C) to stop. Needs numpy and matplotlib.
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
        self.area, self.q1 = None, None              # from the "# area =" / "# q1 =" header lines

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
            if line.startswith("# area ="):
                try:
                    self.area = float(line.split("=", 1)[1])
                except ValueError:
                    pass
                continue
            if line.startswith("# q1 ="):
                self.q1 = line.split("=", 1)[1].strip()
                continue
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
    ap.add_argument("--area", type=float, default=None)
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
    ax0.set_ylabel("relative change\n(mass - mass0) / mass0")
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
            area = a.area if a.area else fol.area
            dm = np.asarray(fol.mass) - fol.mass[0]
            if area:
                l1.set_data(t, dm / area)
                if fol.q1 in (None, "H"):
                    ax1.set_ylabel("equivalent layer\nthickness change [m]")
                else:
                    ax1.set_ylabel(f"mean change of {fol.q1}\n(mass - mass0)/area")
            else:
                l1.set_data(t, dm)
            for ax in (ax0, ax1):
                ax.relim(); ax.autoscale_view()
            extra = f"   = {dm[-1] / area:+.2e} m of layer" if area and fol.q1 in (None, "H") else ""
            ax0.set_title(f"t = {t[-1]:.2f} d   rel. change = {fol.rel[-1]:+.3e}   "
                          f"mass - mass0 = {dm[-1]:+.3e}{extra}   ({len(t)} samples)", fontsize=9)
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
