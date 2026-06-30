#!/usr/bin/env python3
"""
Thesis figure: checkpoint/restart in a Xyce co-simulation IS possible.

The co-simulation integrates the [0, 20 ms] horizon as 50 windows. Each window is a
SEPARATE Xyce process invocation that restarts from the previous window's checkpoint
file (.OPTIONS RESTART). The real series inductor Ls_d carries its state (current +
BDF history) across each restart. If restart worked, the assembled interface current
must be (a) continuous across the 49 window boundaries and (b) identical to a single
monolithic Xyce run of the same circuit (no windows, no restart).

Inputs (from the repo working dir, one dir up):
  Circuit_solution.prn          windowed/restarted co-sim:  Index TIME V(P) V(NX) I(VMEAS)
  ref_monolithic.cir.prn        monolithic reference:       Index TIME I(LF) V(P)
  Field_waveform_solution.prn   per-window end (sync) pts:  Index TIME V(FIELD) I(FIELD)
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import os

HERE = os.path.dirname(os.path.abspath(__file__))
# Read the snapshotted .prn files (thesis_plots/data/) so the figure is reproducible
# regardless of later co-sim runs overwriting the repo-root .prn files. Re-snapshot with:
#   ./main && cp Circuit_solution.prn ref_monolithic.cir.prn Field_waveform_solution.prn thesis_plots/data/
ROOT = os.path.join(HERE, "data")
DT_WIN = (1/50)/50          # window length = dt_field = 4e-4 s
T_END  = 0.02
NWIN   = round(T_END/DT_WIN)  # 50

def load(fname, tcol, ycol):
    t, y = [], []
    for ln in open(os.path.join(ROOT, fname)):
        p = ln.split()
        if len(p) <= max(tcol, ycol):
            continue
        try:
            t.append(float(p[tcol])); y.append(float(p[ycol]))
        except ValueError:
            continue
    return np.asarray(t), np.asarray(y)

# windowed/restarted co-sim interface current, monolithic reference, sync points
tc, ic = load("Circuit_solution.prn", 1, 4)      # I(VMEAS)
tr, ir = load("ref_monolithic.cir.prn", 1, 2)    # I(LF)
ts, isy = load("Field_waveform_solution.prn", 1, 3)  # window-end field current

bounds = np.arange(1, NWIN) * DT_WIN             # 49 interior restart boundaries
ms_c = tc*1e3; ms_r = tr*1e3; ms_s = ts*1e3      # ms for readability
err = ic - np.interp(tc, tr, ir)
peak = np.max(np.abs(ir))

# zoom window: around the current zero-crossing (steepest dI/dt -> hardest continuity test)
z0, z1 = 9.2e-3, 11.2e-3

dense_rms = 100*np.sqrt(np.mean(err**2))/peak     # RMS of the plotted (dense) error
dense_max = 100*np.max(np.abs(err))/peak           # max of the plotted (dense) error

plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3,
                     "figure.dpi": 120})
fig = plt.figure(figsize=(9, 7.6), constrained_layout=True)
gs = fig.add_gridspec(3, 1, height_ratios=[3, 2.4, 1.2])

# ---- (a) full horizon: restarted co-sim vs monolithic --------------------------
axA = fig.add_subplot(gs[0])
axA.plot(ms_r, ir, color="0.55", lw=3.4, label="monolithic Xyce (single run, no restart)")
axA.plot(ms_c, ic, color="C0", lw=1.1, label="co-sim (50 windows, restarted each)")
axA.plot(ms_s, isy, "o", color="C3", ms=4.5, label="window end = checkpoint / restart point")
axA.set_ylabel("interface current  $I$  [A]")
axA.legend(loc="lower left", fontsize=9, framealpha=0.95)
axA.axvspan(z0*1e3, z1*1e3, color="C1", alpha=0.13)
axA.text(0.5*(z0+z1)*1e3, axA.get_ylim()[1]*0.86, "zoom\n(b)", ha="center",
         va="top", fontsize=8, color="C1")
axA.set_xlim(0, T_END*1e3)
axA.set_title("(a) full 20 ms horizon", fontsize=10, loc="left")

# ---- (b) zoom on restart boundaries: continuity across the checkpoint -----------
axB = fig.add_subplot(gs[1])
m  = (tc >= z0) & (tc <= z1)
mr = (tr >= z0) & (tr <= z1)
ms_ = (ts >= z0) & (ts <= z1)
first = True
for b in bounds:
    if z0 <= b <= z1:
        axB.axvline(b*1e3, color="0.7", lw=0.9, ls="--", zorder=0,
                    label="restart boundary" if first else None)
        first = False
axB.plot(ms_r[mr], ir[mr], color="0.55", lw=3.4, label="monolithic")
axB.plot(ms_c[m], ic[m], color="C0", lw=1.5, label="co-sim (restarted windows)")
axB.plot(ms_s[ms_], isy[ms_], "o", color="C3", ms=7, zorder=5, label="restart point")
axB.set_ylabel("interface current  $I$  [A]")
axB.set_xlabel("time  [ms]")
axB.legend(loc="upper right", fontsize=8.5, framealpha=0.95, ncol=2)
axB.set_xlim(z0*1e3, z1*1e3)
axB.set_title("(b) zoom at the zero-crossing (steepest $dI/dt$)", fontsize=9.5, loc="left")

# ---- (c) error vs monolithic ----------------------------------------------------
axC = fig.add_subplot(gs[2])
axC.plot(ms_c, 100*err/peak, color="C2", lw=0.9)
axC.axhline(0, color="0.6", lw=0.7)
axC.set_ylabel("error\n[% peak]")
axC.set_xlabel("time  [ms]")
axC.set_xlim(0, T_END*1e3)
axC.yaxis.set_major_formatter(mtick.FormatStrFormatter("%.1f"))
axC.annotate("window-1 cold start", xy=(0.14, 1.31), xytext=(2.2, 1.15),
             fontsize=8, va="center",
             arrowprops=dict(arrowstyle="->", color="0.4", lw=0.9))
axC.set_title(f"(c) co-sim $-$ monolithic interface current:  RMS {dense_rms:.2f} %,  max {dense_max:.2f} % of peak (at the window-1 cold start)",
              fontsize=9.3, loc="left")

fig.savefig(os.path.join(HERE, "restart_demo.png"), bbox_inches="tight", dpi=300)
fig.savefig(os.path.join(HERE, "restart_demo.pdf"), bbox_inches="tight")
print("wrote restart_demo.png / .pdf")
print(f"co-sim points={len(tc)}  ref points={len(tr)}  windows={NWIN}  boundaries={len(bounds)}")
print(f"peak={peak:.2f} A  max|err|={100*np.max(np.abs(err))/peak:.3f}%  "
      f"sync-RMS={100*np.sqrt(np.mean((np.interp(ts,tr,ir)-isy)**2))/peak:.3f}%")
