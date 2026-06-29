#!/usr/bin/env python3
"""
Thesis figure: WHY restart is necessary -- the negative control.

Same coupled RL circuit, integrated in 50 windows, but WITHOUT restart: each window is
solved independently and cold-started (no checkpoint, inductor IC=0). Without the state
handoff the inductor current is forced back to 0 at every window start -> a sawtooth that
does NOT track the true (monolithic) solution. Overlaid with the WORKING restarted co-sim
(which reproduces the monolithic run), this shows restart is both necessary and correct.

Generates the broken data here (50 standalone Xyce runs); reads the working co-sim and the
monolithic reference from the snapshot in ./data/.
"""
import numpy as np, matplotlib.pyplot as plt, subprocess, os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
DT_WIN = (1/50)/50          # 4e-4 s
T_END  = 0.02
NWIN   = round(T_END/DT_WIN)
PI = "3.14159265358979"
# true physical params (R_series+L_series, R_FEM+L_FEM) -- same circuit as the monolithic ref
RS, LS, RF, LF = "0.006", "1.6e-07", "0.00051", "1.6e-07"

def load(path, tcol, ycol):
    t, y = [], []
    for ln in open(path):
        p = ln.split()
        if len(p) <= max(tcol, ycol): continue
        try: t.append(float(p[tcol])); y.append(float(p[ycol]))
        except ValueError: continue
    return np.asarray(t), np.asarray(y)

def gen_broken():
    """50 independent cold-started windows (no restart), concatenated to absolute time."""
    tb, ib = [], []
    for k in range(NWIN):
        t0 = k * DT_WIN
        net = f"""Broken window {k} (no restart, cold inductors IC=0)
Bsrc s 0 V = {{ sin(2*{PI}*50*(time + {t0:.10e})) }}
Rs s m {RS}
Ls m p {LS} IC=0
Rf p n {RF}
Lf n 0 {LF} IC=0
.tran 1e-7 {DT_WIN:.10e} 0 UIC
.print tran I(Lf)
.end
"""
        open(os.path.join(HERE, "_broken_win.cir"), "w").write(net)
        subprocess.run(["Xyce", os.path.join(HERE, "_broken_win.cir")],
                       capture_output=True, cwd=HERE)
        tk, ik = load(os.path.join(HERE, "_broken_win.cir.prn"), 1, 2)
        m = tk <= DT_WIN * 1.0000001
        tb.append(tk[m] + t0); ib.append(ik[m])
    for ext in (".cir", ".cir.prn"):
        p = os.path.join(HERE, "_broken_win" + ext)
        if os.path.exists(p): os.remove(p)
    return np.concatenate(tb), np.concatenate(ib)

# --- data ----------------------------------------------------------------------
tb, ib = gen_broken()
np.savetxt(os.path.join(DATA, "broken_no_restart.txt"),
           np.column_stack([tb, ib]), header="time  I_Lf_no_restart")
tr, ir = load(os.path.join(DATA, "ref_monolithic.cir.prn"), 1, 2)   # monolithic I(LF)
tc, ic = load(os.path.join(DATA, "Circuit_solution.prn"), 1, 4)     # working co-sim I(VMEAS)
peak = np.max(np.abs(ir))
bounds = np.arange(1, NWIN) * DT_WIN

# --- figure --------------------------------------------------------------------
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3, "figure.dpi": 120})
fig = plt.figure(figsize=(9, 6.6), constrained_layout=True)
gs = fig.add_gridspec(2, 1, height_ratios=[3, 2.3])
fig.suptitle("Why restart is needed: with restart the windows reproduce the monolithic run; without, they don't",
             fontsize=12, fontweight="bold")

z0, z1 = 3.0e-3, 6.0e-3   # zoom near the peak: broken sawtooth vs smooth restarted curve

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
axB.set_title("(b) zoom: without restart the inductor current resets to 0 each window (sawtooth); with restart it is continuous",
              fontsize=9.3, loc="left")

fig.savefig(os.path.join(HERE, "restart_control.png"), bbox_inches="tight", dpi=300)
fig.savefig(os.path.join(HERE, "restart_control.pdf"), bbox_inches="tight")
err_broken = 100*np.max(np.abs(ib - np.interp(tb, tr, ir)))/peak
err_work   = 100*np.max(np.abs(ic - np.interp(tc, tr, ir)))/peak
print("wrote restart_control.png / .pdf")
print(f"peak={peak:.2f} A   max|err| WITHOUT restart = {err_broken:.1f}% of peak   "
      f"WITH restart = {err_work:.2f}%")
