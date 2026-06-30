#!/usr/bin/env python3
"""
Thesis figure: WHY restart is needed -- the negative control, through the ACTUAL co-sim.

Same matched-secant co-sim (Bemf/Rs_d/Ls_d + ROM Bfield + FEM + WR), integrated in 50
windows, but with checkpoint/restart DISABLED (sim_config break_restart=1): each window runs
cold in window-local time, inductor IC=0, only the source phase carried. Without the state
handoff the interface current resets each window -> sawtooth that does NOT track the true
(monolithic) solution. Overlaid with the WORKING restarted co-sim (which reproduces the
monolithic run), this shows restart is both necessary and correct.

Data (snapshots in ./data/, regenerate as noted):
  ref_monolithic.cir.prn   monolithic reference          Index TIME I(LF) V(P)
  Circuit_solution.prn      co-sim WITH restart           Index TIME V(P) V(NX) I(VMEAS)   (./main, break_restart=0)
  broken_cosim.prn          co-sim WITHOUT restart        Index TIME V(P) V(NX) I(VMEAS)   (./main, break_restart=1, WRmaxSteps>=120)
"""
import numpy as np, matplotlib.pyplot as plt, os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
DT_WIN = (1/50)/50
T_END  = 0.02
NWIN   = round(T_END/DT_WIN)

def load(name, tcol, ycol):
    t, y = [], []
    for ln in open(os.path.join(DATA, name)):
        p = ln.split()
        if len(p) <= max(tcol, ycol): continue
        try: t.append(float(p[tcol])); y.append(float(p[ycol]))
        except ValueError: continue
    return np.asarray(t), np.asarray(y)

tr, ir = load("ref_monolithic.cir.prn", 1, 2)   # monolithic I(LF)
tc, ic = load("Circuit_solution.prn", 1, 4)      # WITH restart   I(VMEAS)
tb, ib = load("broken_cosim.prn", 1, 4)          # WITHOUT restart I(VMEAS)
peak = np.max(np.abs(ir))
bounds = np.arange(1, NWIN) * DT_WIN

plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3, "figure.dpi": 120})
fig = plt.figure(figsize=(9, 6.6), constrained_layout=True)
gs = fig.add_gridspec(2, 1, height_ratios=[3, 2.3])

z0, z1 = 9.2e-3, 11.2e-3   # zoom at the zero-crossing (steepest dI/dt), matching restart_demo

axA = fig.add_subplot(gs[0])
axA.plot(tr*1e3, ir, color="0.55", lw=3.4, label="monolithic Xyce (truth)")
axA.plot(tc*1e3, ic, color="C0", lw=1.2, label="co-sim WITH restart")
axA.plot(tb*1e3, ib, color="C3", lw=1.0, label="co-sim WITHOUT restart (cold-start each window)")
axA.axvspan(z0*1e3, z1*1e3, color="C1", alpha=0.13)
axA.set_ylabel("interface current  $I$  [A]")
axA.legend(loc="lower left", fontsize=9, framealpha=0.95)
axA.set_xlim(0, T_END*1e3); axA.set_title("(a) full 20 ms horizon", fontsize=10, loc="left")

axB = fig.add_subplot(gs[1])
first = True
for b in bounds:
    if z0 <= b <= z1:
        axB.axvline(b*1e3, color="0.7", lw=0.9, ls="--", zorder=0,
                    label="window / restart boundary" if first else None); first = False
mr=(tr>=z0)&(tr<=z1); mc=(tc>=z0)&(tc<=z1); mb=(tb>=z0)&(tb<=z1)
axB.plot(tr[mr]*1e3, ir[mr], color="0.55", lw=3.4, label="monolithic")
axB.plot(tc[mc]*1e3, ic[mc], color="C0", lw=1.6, label="WITH restart")
axB.plot(tb[mb]*1e3, ib[mb], color="C3", lw=1.4, label="WITHOUT restart")
axB.set_ylabel("interface current  $I$  [A]"); axB.set_xlabel("time  [ms]")
axB.legend(loc="lower right", fontsize=8.5, framealpha=0.95, ncol=2)
axB.set_xlim(z0*1e3, z1*1e3)
axB.set_title("(b) zoom", fontsize=10, loc="left")

fig.savefig(os.path.join(HERE, "restart_control.png"), bbox_inches="tight", dpi=300)
fig.savefig(os.path.join(HERE, "restart_control.pdf"), bbox_inches="tight")
err_broken = 100*np.max(np.abs(ib - np.interp(tb, tr, ir)))/peak
err_work   = 100*np.max(np.abs(ic - np.interp(tc, tr, ir)))/peak
print("wrote restart_control.png / .pdf")
print(f"peak={peak:.2f} A   max|err| WITHOUT restart = {err_broken:.1f}% of peak   WITH restart = {err_work:.2f}%")
