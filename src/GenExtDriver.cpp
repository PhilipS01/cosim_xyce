// GenExtDriver.cpp -- IMPLICIT field-circuit coupling via Xyce's General External
// Device (YGENEXT) + GenCouplingSimulator, replacing the partitioned system("Xyce")
// + PWL exchange. The field is a real Xyce device whose flux Q = Lrom*i is a SHARED
// STATE (Xyce integrates dQ/dt itself -> no PWL/DDT staircase) and whose admittance
// dFdx/dQdx is in the Newton Jacobian (-> implicit, breaks the strong-Ls WR wall).
// See ../doc/AppNote-GenExt.pdf (3.3.1 computeXyceVectors, 4.1.x init + simulateUntil,
// 5.1 RLC example).
//
// Step 1 (this file): per-window deferred FEM correction (multirate).
//   - Lrom = L_FEM  -> the inductive admittance is EXACT and implicit in Q (no
//     derivative anywhere). The deferred correction is therefore purely RESISTIVE,
//     injected on the F side (no d/dt -> a time-waveform offset is fine, no staircase).
//   - Per window: sub-step simulateUntil at coupling resolution + getSolution to
//     capture V_p(t), i(t); run the lumped FEM over that V_p; set next window's
//     resistive correction voff(t) = (R_FEM - Rrom)*I_FEM(t).  FEM called ONCE per
//     window (multirate), not per Newton step.

#include <Xyce_config.h>
#include <N_CIR_GenCouplingSimulator.h>
#include <N_DEV_VectorComputeInterface.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Field ROM as a 2-terminal series R-L device (RLC example minus C).
// Device vars: external nodes A(0)=port p, B(1)=ground; internal branch (2) = i.
//   F[A]   =  i ;  F[B] = -i
//   F[br]  =  Rrom*i - (V_A - V_B) + voff(t)        (resistive + deferred FEM correction)
//   Q[br]  =  Lrom*i                                (Xyce: dQ/dt = Lrom*di/dt)
// voff(t): piecewise-linear resistive correction over the current window (interpolated
// from a stored waveform). It is on the F (resistive) side, so no derivative -> safe.
// -----------------------------------------------------------------------------
class FieldRomVCI : public Xyce::Device::VectorComputeInterface
{
public:
  FieldRomVCI(double Rrom, double Lrom, bool saturating = false, double Isat = 100.0)
    : Rrom_(Rrom), Lrom_(Lrom), saturating_(saturating), Isat_(Isat)
  {
    jacStamp.resize(3);
    jacStamp[A_].resize(1);  jacStamp[A_][0] = BR_;
    jacStamp[B_].resize(1);  jacStamp[B_][0] = BR_;
    jacStamp[BR_].resize(3); jacStamp[BR_][0] = A_; jacStamp[BR_][1] = B_; jacStamp[BR_][2] = BR_;
  }

  // Saturation flux model put in Q (the SHARED STATE Xyce differentiates itself):
  //   lambda(i) = Lrom*Isat*atan(i/Isat)  -> dlambda/di = Lrom/(1+(i/Isat)^2) = differential L.
  // So dQ/dt = L_dyn(i)*di/dt is the EXACT nonlinear inductive voltage, and dQdx = L_dyn(i)
  // (the true differential inductance) lands in the Newton Jacobian -> implicit, no derivative
  // computed in C++, no staircase. Linear limit (Isat->inf or saturating_=false): lambda=Lrom*i.
  double flux(double i)  const { return saturating_ ? Lrom_ * Isat_ * std::atan(i / Isat_) : Lrom_ * i; }
  double Ldyn(double i)  const { const double r = i / Isat_; return saturating_ ? Lrom_ / (1.0 + r*r) : Lrom_; }

  // Set the resistive correction waveform voff(t) for the upcoming window (absolute time).
  void setCorrection(const std::vector<double> & t, const std::vector<double> & v)
  { ct_ = t; cv_ = v; }

  // BLACK-BOX mode: per-window probed surrogate. Rrom=R_diff; flux is a QUADRATIC around the
  // operating point i0:  Q = lambda0 + L_diff*(i-i0) + 0.5*curv*(i-i0)^2,  so the differential
  // inductance dQ/di = L_diff + curv*(i-i0) VARIES within the window (captures saturation
  // curvature). curv=0 -> tangent line. lambda0 chosen by the driver for Q-continuity at the
  // window boundary (offset is voltage-irrelevant; a Q jump would spike dQ/dt -> Xyce abort).
  void setProbed(double Rdiff, double Ldiff, double i0, double lambda0, double curv = 0.0)
  { Rrom_ = Rdiff; Lrom_ = Ldiff; curv_ = curv; i0_ = i0; lambda0_ = lambda0; linFlux_ = true; }

  double interpCorrection(double time) const
  {
    if (ct_.empty()) return 0.0;
    if (time <= ct_.front()) return cv_.front();
    if (time >= ct_.back())  return cv_.back();
    // linear search is fine (small grids); could bisect.
    for (size_t j = 1; j < ct_.size(); ++j) {
      if (time <= ct_[j]) {
        const double a = (time - ct_[j-1]) / (ct_[j] - ct_[j-1]);
        return cv_[j-1] + a * (cv_[j] - cv_[j-1]);
      }
    }
    return cv_.back();
  }

  bool computeXyceVectors(std::vector<double> & sV,
                          double time,
                          std::vector<double> & F,
                          std::vector<double> & Q,
                          std::vector<double> & B,
                          std::vector<std::vector<double> > & dFdx,
                          std::vector<std::vector<double> > & dQdx) override
  {
    // BACKWARD-STEP SAFETY: Xyce does adaptive stepping and may REJECT a step (LTE too
    // large), then call us again at an EARLIER time (AppNote 3.6). This device is a PURE
    // FUNCTION of (sV, time): F,Q,dFdx,dQdx are recomputed from the live solution and a
    // pure interpolation voff(time); NOTHING is cached across calls. So a backward jump
    // needs no rollback -- we just return the correct values for whatever time is asked.
    // (We only count jumps for diagnostics; the count never changes the result.)
    if (time < lastTime_ - 1e-18) ++backwardJumps_;
    lastTime_ = time;

    const int n = static_cast<int>(sV.size()); // 3
    F.assign(n, 0.0); Q.assign(n, 0.0); B.clear();
    dFdx.assign(n, std::vector<double>(n, 0.0));
    dQdx.assign(n, std::vector<double>(n, 0.0));

    const double i = sV[BR_];
    const double voff = interpCorrection(time);

    F[A_]  =  i;
    F[B_]  = -i;
    F[BR_] =  Rrom_ * i - (sV[A_] - sV[B_]) + voff;
    // Flux in Q: black-box -> probed quadratic lambda0+L_diff*(i-i0)+0.5*curv*(i-i0)^2;
    //            else analytic lambda(i).
    const double di = i - i0_;
    Q[BR_] =  linFlux_ ? (lambda0_ + Lrom_ * di + 0.5 * curv_ * di * di) : flux(i);

    dFdx[A_][BR_]  =  1.0;
    dFdx[B_][BR_]  = -1.0;
    dFdx[BR_][A_]  = -1.0;
    dFdx[BR_][B_]  =  1.0;
    dFdx[BR_][BR_] =  Rrom_;
    dQdx[BR_][BR_] =  linFlux_ ? (Lrom_ + curv_ * di) : Ldyn(i);  // differential inductance
    return true;
  }

  long backwardJumps() const { return backwardJumps_; } // diagnostics (rejected-step retries)

  std::vector<std::vector<int> > jacStamp;

private:
  static const int A_ = 0, B_ = 1, BR_ = 2;
  double Rrom_, Lrom_;
  bool saturating_;
  double Isat_;
  bool   linFlux_ = false;       // black-box: use probed (quadratic) flux instead of analytic
  double i0_ = 0.0, lambda0_ = 0.0, curv_ = 0.0;
  std::vector<double> ct_, cv_; // correction waveform (time, value)
  mutable double lastTime_ = -1e300; // for backward-jump detection (diagnostic only)
  mutable long backwardJumps_ = 0;
};

// -----------------------------------------------------------------------------
// Lumped "FEM": given the port-voltage samples V_p(t) over a window, integrate the
// true field  V_p = R_FEM*I + L_FEM*dI/dt  with backward Euler, return I_FEM(t).
// (A real FEM would replace this; the interface is V_p in -> I_FEM, lambda out.)
// -----------------------------------------------------------------------------
// Black-box FEM interface: given V_p(t), return I_FEM(t) and the flux linkage lambda_FEM(t).
// The DRIVER treats this as opaque -- it never reads R_FEM/L_FEM/Isat; those live ONLY inside
// here (as a real distributed FEM would). The Jacobian (R_diff, L_diff) is recovered by probing
// these outputs (probeJacobian), not by knowing the parameters.
static std::vector<double> femEvaluate(const std::vector<double> & t,
                                       const std::vector<double> & Vp,
                                       double R_FEM, double L_FEM, double I_start,
                                       bool saturating, double Isat,
                                       std::vector<double> * lambdaOut)
{
  auto flux  = [&](double I){ return saturating ? L_FEM * Isat * std::atan(I / Isat) : L_FEM * I; };
  auto Ldyn  = [&](double I){ const double r = I / Isat; return saturating ? L_FEM / (1.0 + r*r) : L_FEM; };
  const size_t N = t.size();
  std::vector<double> I(N, 0.0);
  I[0] = I_start;
  for (size_t j = 1; j < N; ++j) {
    const double h = t[j] - t[j-1];
    const double lam_prev = flux(I[j-1]);
    double Ij = (Vp[j] + (Ldyn(I[j-1]) / h) * I[j-1]) / (R_FEM + Ldyn(I[j-1]) / h); // linear seed
    if (saturating) {
      for (int it = 0; it < 50; ++it) {
        const double g  = R_FEM * Ij + (flux(Ij) - lam_prev) / h - Vp[j];
        const double gp = R_FEM + Ldyn(Ij) / h;
        const double dI = g / gp;
        Ij -= dI;
        if (std::fabs(dI) <= 1e-12 + 1e-10 * std::fabs(Ij)) break;
      }
    }
    I[j] = Ij;
  }
  if (lambdaOut) { lambdaOut->resize(N); for (size_t j = 0; j < N; ++j) (*lambdaOut)[j] = flux(I[j]); }
  return I;
}

// Probe the field's port Jacobian from BLACK-BOX outputs only (V_p, I_FEM, lambda_FEM over a
// window). R_diff: LS slope of (V_p - dlambda/dt) vs I (field relation V_p = R*I + dlambda/dt).
// Flux model lambda(I): QUADRATIC fit  lambda ~ a0 + a1*I + a2*I^2  (centered for conditioning)
// -> differential inductance dlambda/dI = a1 + 2*a2*I = Lslope0 + curv*I, where Lslope0=a1 and
// curv=2*a2. The quadratic captures the saturation CURVATURE within a window (the tangent line
// = curv 0). dlambda/dt by FD of the FEM's smooth lambda (only to FORM the surrogate, never
// re-differentiated by Xyce). Linear field: a2~0 -> recovers R_FEM, L_FEM; tangent line.
static void probeJacobian(const std::vector<double> & t,
                          const std::vector<double> & Vp,
                          const std::vector<double> & I,
                          const std::vector<double> & lam,
                          double & Rdiff, double & Lslope0, double & curv)
{
  const size_t N = t.size();
  auto lsSlope = [&](const std::vector<double> & x, const std::vector<double> & y) {
    double sx=0, sy=0, sxx=0, sxy=0;
    for (size_t j = 0; j < N; ++j) { sx+=x[j]; sy+=y[j]; sxx+=x[j]*x[j]; sxy+=x[j]*y[j]; }
    const double d = (double)N*sxx - sx*sx;
    return (std::fabs(d) < 1e-300) ? 0.0 : ((double)N*sxy - sx*sy) / d;
  };
  // R_diff (resistive Jacobian)
  std::vector<double> vres(N);
  for (size_t j = 0; j < N; ++j) {
    const double dlamdt = (j == 0) ? (lam[1]-lam[0])/(t[1]-t[0])
                                   : (lam[j]-lam[j-1])/(t[j]-t[j-1]);
    vres[j] = Vp[j] - dlamdt;
  }
  Rdiff = lsSlope(I, vres);

  // Quadratic LS fit lambda ~ c0 + c1*x + c2*x^2 with x = I - Imean (centered -> conditioned).
  double Imean = 0.0; for (double v : I) Imean += v; Imean /= (double)N;
  double S0=N, S1=0,S2=0,S3=0,S4=0, T0=0,T1=0,T2=0;
  for (size_t j = 0; j < N; ++j) {
    const double x = I[j]-Imean, x2=x*x, y=lam[j];
    S1+=x; S2+=x2; S3+=x2*x; S4+=x2*x2; T0+=y; T1+=x*y; T2+=x2*y;
  }
  // 3x3 normal equations [[S0 S1 S2],[S1 S2 S3],[S2 S3 S4]] [c0 c1 c2]^T = [T0 T1 T2]^T
  double A[3][3] = {{S0,S1,S2},{S1,S2,S3},{S2,S3,S4}};
  double b[3] = {T0,T1,T2};
  for (int col=0; col<3; ++col) {
    int piv=col; for(int r=col+1;r<3;++r) if(std::fabs(A[r][col])>std::fabs(A[piv][col])) piv=r;
    for(int c=0;c<3;++c) std::swap(A[col][c],A[piv][c]); std::swap(b[col],b[piv]);
    const double d=A[col][col]; if(std::fabs(d)<1e-300) continue;
    for(int r=col+1;r<3;++r){ const double f=A[r][col]/d; for(int c=col;c<3;++c) A[r][c]-=f*A[col][c]; b[r]-=f*b[col]; }
  }
  double c[3]={0,0,0};
  for(int ii=2;ii>=0;--ii){ double s=b[ii]; for(int cc=ii+1;cc<3;++cc) s-=A[ii][cc]*c[cc]; c[ii]=(std::fabs(A[ii][ii])<1e-300)?0.0:s/A[ii][ii]; }
  // centered: dlambda/dI = c1 + 2*c2*(I-Imean). curv = 2*c2; Lslope0 = a1 (un-centered) = c1 - 2*c2*Imean.
  curv    = 2.0 * c[2];
  Lslope0 = c[1] - 2.0 * c[2] * Imean;
}

// -----------------------------------------------------------------------------
// Config (sim_config.txt, same key=value format as the partitioned solver) drives BOTH the
// field model and the generated netlist. Env vars override (GENEXT_*).
struct GxCfg {
  double frequency = 50.0, amplitude = 1.0;
  double R_series = 6.0e-3, L_series = 1.6e-7;   // circuit-side coupling (Rs, Ls)
  double R_FEM = 5.1e-4, L_FEM = 1.6e-7;         // "true" field
  double R_ROM = 4.59e-4;                        // ROM resistance (known-param device; Lrom=L_FEM)
  double I_sat = 100.0;
  int    nonlin_model = 0;                       // 1 -> magnetic saturation
  int    N_periods = 1, N_field_steps = 50;      // -> t_final, dt_window
  int    blackbox = 0;                           // 1 -> probe the Jacobian (no known field params)
  // Coupling mode. 0 = deferred single pass (one continuous BDF sim; correction time-lagged one
  //   window). 1 = per-window WR (each window a separate restarted sim; voff iterated to the
  //   fixpoint by re-running the window -- demonstrates convergence + rewind; floors above 0 due
  //   to per-window order-1 NOOP restart). 2 = GLOBAL waveform iteration (pure WR: one continuous
  //   BDF sim over the WHOLE horizon per sweep, black-box FEM once per sweep, voff(t) relaxed to
  //   the fixpoint over the whole axis -- keeps BDF history -> beats mode 1's floor; restart trivial
  //   since each sweep starts at t=0).
  int    iterate = 0;
  int    max_iter = 20;                          // cap on per-window WR iterations
  double tol = 1.0e-4;                           // WR converged when ||dvoff|| < tol*||voff|| (fixpoint)
  double theta = 1.0;                            // under-relaxation: voff += theta*(r - voff)
};

static void loadGxCfg(const char * path, GxCfg & c)
{
  FILE * f = std::fopen(path, "r");
  if (!f) { printf("loadGxCfg: '%s' not found, using defaults.\n", path); return; }
  char line[256];
  while (std::fgets(line, sizeof(line), f)) {
    std::string s(line); const size_t h = s.find('#'); if (h != std::string::npos) s = s.substr(0, h);
    for (char & ch : s) if (ch == '=' || ch == ',' || ch == '\t') ch = ' ';
    std::string key; double val; std::istringstream iss(s);
    if (!(iss >> key)) continue; if (!(iss >> val)) continue;
    if      (key == "frequency")                       c.frequency = val;
    else if (key == "amplitude")                       c.amplitude = val;
    else if (key == "R_series")                        c.R_series = val;
    else if (key == "L_series")                        c.L_series = val;
    else if (key == "R_FEM")                           c.R_FEM = val;
    else if (key == "L_FEM")                           c.L_FEM = val;
    else if (key == "R_ROM")                           c.R_ROM = val;
    else if (key == "I_sat")                           c.I_sat = val;
    else if (key == "nonlin_model")                    c.nonlin_model = (int)val;
    else if (key == "N_periods")                       c.N_periods = (int)val;
    else if (key == "N_field_steps_per_source_period") c.N_field_steps = (int)val;
    else if (key == "wr_genext_blackbox")              c.blackbox = (int)val;
    else if (key == "wr_genext_iterate")               c.iterate = (int)val;
    else if (key == "wr_genext_max_iter")              c.max_iter = (int)val;
    else if (key == "wr_genext_tol")                   c.tol = val;
    else if (key == "wr_genext_theta")                 c.theta = val;
  }
  std::fclose(f);
}

// Generate the YGENEXT netlist from config: EMF source, series Rs+Ls, the field device, .tran.
//
// WINDOW RESTART (resume == true): Xyce's NATIVE checkpoint/restart (.OPTIONS RESTART) does NOT
//   round-trip the YGENEXT internal branch current -- verified by inspecting an unpacked checkpoint
//   (the field current ~17 A is absent), consistent with AppNote 3.6 (GenExt internal vars are not
//   first-class persistent state). And .IC on the internal node is ignored unless the OP is skipped.
//   So we reconstruct the window-start state EXPLICITLY (the driver is the checkpoint store): the
//   only memory in this circuit is the (series) inductor current i, so we seed it on BOTH the real
//   inductor (Ls_d IC=i) and the field branch (.IC V(YGENEXT!FIELD_internalnode_0)=i), and use NOOP
//   to skip the operating point and start transient directly from these ICs (AppNote 3.6 workaround;
//   "all zeros except where specified by .IC"). The algebraic node voltages re-solve at the first
//   step. resume == false -> fresh window from t=0 (UIC, all IC=0).
static void writeGxNetlist(const char * path, const GxCfg & c, double t_final, double dt_print,
                           bool resume = false, double iSeed = 0.0)
{
  FILE * f = std::fopen(path, "w");
  if (!f) { fprintf(stderr, "writeGxNetlist: cannot open %s\n", path); return; }
  std::fprintf(f, "WR GenExt netlist (generated from sim_config by GenExtDriver)\n");
  std::fprintf(f, "Bemf  emf 0  V = { %.10g*sin(2*3.14159265358979*%.10g*time) }\n", c.amplitude, c.frequency);
  std::fprintf(f, "Rs_d  emf a  %.10g\n", c.R_series);
  if (c.L_series > 0.0) std::fprintf(f, "Ls_d  a  b  %.10g IC=%.10g\n", c.L_series, resume ? iSeed : 0.0);
  else                  std::fprintf(f, "Rls_d a  b  1e-9\n");   // Ls=0 -> tiny-R short
  std::fprintf(f, "Vmeas b  p   0\n");
  std::fprintf(f, "YGENEXT field p 0\n");
  std::fprintf(f, ".tran %.10g %.10g 0 %s\n", dt_print, t_final, resume ? "NOOP" : "UIC");
  if (resume) std::fprintf(f, ".IC V(YGENEXT!FIELD_internalnode_0)=%.10g\n", iSeed);
  std::fprintf(f, ".print tran V(p) I(Vmeas)\n.end\n");
  std::fclose(f);
}

// Linear interpolation of a (time,value) waveform; flat outside the range. Empty -> 0.
static double interp1(double time, const std::vector<double> & T, const std::vector<double> & V)
{
  if (T.empty()) return 0.0;
  if (time <= T.front()) return V.front();
  if (time >= T.back())  return V.back();
  for (size_t j = 1; j < T.size(); ++j)
    if (time <= T[j]) { const double a = (time - T[j-1]) / (T[j] - T[j-1]); return V[j-1] + a*(V[j]-V[j-1]); }
  return V.back();
}

// Per-call device configuration for runWindow (the VCI is rebuilt fresh each window/iteration,
// since a restart needs a fresh GenCouplingSimulator lifecycle).
struct DevCfg {
  double Rrom, Lrom; bool sat; double Isat;
  bool   probed = false; double pR = 0, pL = 0, pI0 = 0, pLam0 = 0, pCurv = 0;
  std::vector<double> voffT, voffV;   // correction waveform over the window (absolute time)
};
struct WinCap { std::vector<double> ts, vp, ic; bool ok = false; long backJumps = 0; };

// Run ONE window [t0,t1] as a self-contained Xyce lifecycle and capture V_p(t), i(t).
// restart>=0 -> restart from the checkpoint at that time; ckpt>0 -> write checkpoints (so one
// lands at t1 for the next window/iteration). The window-start sample (t0, Vp_start, I_seed) is
// carried in (getSolution is only valid AFTER a simulateUntil; on a restart run the first solve
// lands at the first sub-step, not t0). This is the rewind primitive: to iterate a window we just
// call runWindow again with restart=t0 and an updated dev.voff -> Xyce re-solves the SAME window.
static WinCap runWindow(const char * netlist, const GxCfg & cfg, double t0, double t1,
                        double dt_print, int n_sub, bool resume,
                        double Vp_start, double I_seed, const DevCfg & dev, double iBranchSeed)
{
  WinCap w;
  // resume windows reconstruct the start state from the carried inductor current (see writeGxNetlist).
  writeGxNetlist(netlist, cfg, t1, dt_print, resume, iBranchSeed);

  Xyce::Circuit::GenCouplingSimulator xyce;
  char prog[] = "xyce";
  std::vector<char> nlbuf(netlist, netlist + std::string(netlist).size() + 1);
  char * xargv[2] = { prog, nlbuf.data() };
  if (xyce.initializeEarly(2, xargv) == Xyce::Circuit::Simulator::ERROR) return w;

  std::vector<std::string> names;
  if (!xyce.getDeviceNames("YGENEXT", names) || names.empty()) { xyce.finalize(); return w; }
  const std::string d = names.front();

  FieldRomVCI vci(dev.Rrom, dev.Lrom, dev.sat, dev.Isat);
  if (dev.probed) vci.setProbed(dev.pR, dev.pL, dev.pI0, dev.pLam0, dev.pCurv);
  vci.setCorrection(dev.voffT, dev.voffV);
  xyce.setNumInternalVars(d, 1);
  xyce.setJacStamp(d, vci.jacStamp);
  xyce.setVectorLoader(d, &vci);
  if (xyce.initializeLate() == Xyce::Circuit::Simulator::ERROR) { xyce.finalize(); return w; }

  w.ts.push_back(t0); w.vp.push_back(Vp_start); w.ic.push_back(I_seed);
  for (int s = 1; s <= n_sub; ++s) {
    const double t_req = t0 + (t1 - t0) * (double)s / (double)n_sub;
    double t_done = t0;
    if (!xyce.simulateUntil(t_req, t_done)) { xyce.finalize(); return w; }
    std::vector<double> sV;
    if (xyce.getSolution(d, sV) && sV.size() >= 3) { w.ts.push_back(t_done); w.vp.push_back(sV[0]); w.ic.push_back(sV[2]); }
    if (t_done < t_req - 1e-15) break;
  }
  w.backJumps = vci.backwardJumps();
  xyce.finalize();
  w.ok = (w.ts.size() >= 2);
  return w;
}

int main(int argc, char ** argv)
{
  const char * cfgpath = (argc > 1) ? argv[1] : "sim_config.txt";
  GxCfg cfg; loadGxCfg(cfgpath, cfg);
  // Env overrides (handy for sweeps/tests).
  if (getenv("GENEXT_LFEM"))     cfg.L_FEM = atof(getenv("GENEXT_LFEM"));
  if (getenv("GENEXT_ISAT"))     cfg.I_sat = atof(getenv("GENEXT_ISAT"));
  if (getenv("GENEXT_SAT"))      cfg.nonlin_model = atoi(getenv("GENEXT_SAT"));
  if (getenv("GENEXT_BLACKBOX")) cfg.blackbox = atoi(getenv("GENEXT_BLACKBOX"));
  if (getenv("GENEXT_LS"))       cfg.L_series = atof(getenv("GENEXT_LS"));
  if (getenv("GENEXT_ITERATE"))  cfg.iterate  = atoi(getenv("GENEXT_ITERATE"));
  if (getenv("GENEXT_MAXITER"))  cfg.max_iter = atoi(getenv("GENEXT_MAXITER"));
  if (getenv("GENEXT_TOL"))      cfg.tol      = atof(getenv("GENEXT_TOL"));
  if (getenv("GENEXT_THETA"))    cfg.theta    = atof(getenv("GENEXT_THETA"));

  const double R_FEM = cfg.R_FEM;
  const double L_FEM = cfg.L_FEM;
  const double Rrom  = cfg.R_ROM;   // known-param device resistance (mismatched on purpose)
  const double Lrom  = L_FEM;       // matched inductance (exact implicit) in known-param mode
  const bool   saturating = (cfg.nonlin_model == 1);
  const double Isat = cfg.I_sat;
  const bool   blackbox = (cfg.blackbox != 0);

  const char * netlist = "wr_genext.cir";

  // Source (for the pre-probe; matches the generated netlist Bemf).
  const double freq = cfg.frequency, amp = cfg.amplitude;

  // Time stepping from config: window = one field step; t_final = N_periods source periods.
  const double dt_window = (1.0 / cfg.frequency) / (double)cfg.N_field_steps;
  const double t_final   = (double)cfg.N_periods / cfg.frequency;
  const int    n_sub     = 20;      // sub-steps per window to capture V_p(t) for the FEM

  printf("Config: f=%.1f amp=%.3g Rs=%.3g Ls=%.3g R_FEM=%.3g L_FEM=%.3g sat=%d Isat=%.1f blackbox=%d\n",
         cfg.frequency, cfg.amplitude, cfg.R_series, cfg.L_series, R_FEM, L_FEM, (int)saturating, Isat, (int)blackbox);

  using Clock = std::chrono::steady_clock;
  auto secs = [](Clock::time_point a, Clock::time_point b){
    return std::chrono::duration<double>(b - a).count(); };

  // =========================================================================================
  // (B) GLOBAL WAVEFORM ITERATION -- pure WR over the WHOLE horizon (the thesis-aligned mode).
  // Each sweep: ONE continuous Xyce run [0,t_final] (full adaptive BDF -> keeps history, no per-
  // window order-1 restart) with the current correction voff(t); black-box FEM ONCE over the whole
  // V_p(t) (multirate preserved); relax voff(t) += theta*(r - voff) over the whole axis; re-run from
  // t=0. Iteration-lagged (sweep), not time-lagged. Restart trivial (every sweep starts at t=0 ->
  // the YGENEXT internal-var restart problem never arises). The implicit ROM (Rrom,Lrom in the
  // device Jacobian) is the Robin preconditioner that gives WR its contraction (no strong-Ls wall).
  // =========================================================================================
  if (cfg.iterate == 2) {
    const auto t0wall = Clock::now();
    const double w_src = 2.0 * M_PI * freq;
    const int    Ntot  = std::max(2, cfg.N_periods * cfg.N_field_steps * n_sub); // capture pts / horizon
    const double dtp   = t_final / (double)Ntot;
    printf("GLOBAL-WR mode: max_sweeps=%d tol=%.1e theta=%.3g  (continuous BDF, %d capture pts)\n",
           cfg.max_iter, cfg.tol, cfg.theta, Ntot);

    // Black-box: one whole-horizon linearization, refreshed each sweep (cold start from a pre-probe).
    double curR = 1e-3, curL = 1e-6, curC = 0.0;
    if (blackbox) {
      std::vector<double> pts(n_sub+1), pvp(n_sub+1);
      for (int j=0;j<=n_sub;++j){ pts[j]=dt_window*(double)j/n_sub; pvp[j]=amp*w_src*pts[j]; }
      std::vector<double> plam;
      const std::vector<double> pI = femEvaluate(pts,pvp,R_FEM,L_FEM,0.0,saturating,Isat,&plam);
      double Rd,Ls0,cv; probeJacobian(pts,pvp,pI,plam,Rd,Ls0,cv);
      curR=Rd; curL=Ls0; curC=cv;
      printf("pre-probe: R_diff=%.3e L_diff(0)=%.3e curv=%.3e\n",curR,curL,curC);
    }

    std::vector<double> voffT, voffV;           // whole-horizon correction (empty -> 0)
    std::vector<double> Ifem, lam, r;
    double maxerr = 0.0, lastSS = 0.0; long femEvals = 0; int s = 0;
    double tXyce = 0.0, tFem = 0.0;
    printf("\n  sweep   max|i-I_FEM|   steady-state    ||dvoff||(rel)\n");
    for (s = 0; s < cfg.max_iter; ++s) {
      DevCfg dev;
      dev.Rrom = Rrom; dev.Lrom = Lrom; dev.sat = saturating; dev.Isat = Isat;
      dev.voffT = voffT; dev.voffV = voffV;
      if (blackbox) { dev.probed=true; dev.pR=curR; dev.pL=curL; dev.pCurv=curC; dev.pI0=0.0; dev.pLam0=0.0; }

      const auto ta = Clock::now();
      WinCap cap = runWindow(netlist, cfg, 0.0, t_final, dtp, Ntot, /*resume=*/false,
                             /*Vp_start=*/0.0, /*I_seed=*/0.0, dev, /*iBranchSeed=*/0.0);
      tXyce += secs(ta, Clock::now());
      if (!cap.ok) { fprintf(stderr, "GLOBAL-WR: Xyce run failed (sweep %d)\n", s); return 1; }

      const auto tb = Clock::now();
      Ifem = femEvaluate(cap.ts, cap.vp, R_FEM, L_FEM, 0.0, saturating, Isat, &lam); ++femEvals;
      tFem += secs(tb, Clock::now());

      // Two error metrics: global max (dominated by the t=0 cold start, same point in every mode)
      // and steady-state (worst point over the LAST source period -- the fair comparison number).
      maxerr = 0.0; double ssErr = 0.0; const double tSS = 0.75 * t_final; // last quarter = "steady"
      for (size_t j=0;j<cap.ts.size();++j) {
        const double e = std::fabs(cap.ic[j]-Ifem[j]);
        maxerr = std::max(maxerr, e);
        if (cap.ts[j] >= tSS) ssErr = std::max(ssErr, e);
      }
      lastSS = ssErr;

      double Rsur = Rrom, Lsur = Lrom;
      if (blackbox) { double Rd,Ls0,cv; probeJacobian(cap.ts,cap.vp,Ifem,lam,Rd,Ls0,cv);
                      Rsur=Rd; Lsur=Ls0; curR=Rd; curL=Ls0; curC=cv; }

      const size_t Nc = cap.ts.size();
      r.assign(Nc, 0.0);
      for (size_t j=0;j<Nc;++j) {
        // Order-2 (central) derivative of I_FEM, to MATCH Xyce's BDF2 in the device. Order-1
        // backward FD here leaves a (Lrom/Rrom)*(FD-BDF) transmission floor; central cancels it.
        double dIdt;
        if (j==0)            dIdt = (Ifem[1]-Ifem[0])/(cap.ts[1]-cap.ts[0]);
        else if (j==Nc-1)    dIdt = (Ifem[j]-Ifem[j-1])/(cap.ts[j]-cap.ts[j-1]);
        else                 dIdt = (Ifem[j+1]-Ifem[j-1])/(cap.ts[j+1]-cap.ts[j-1]);
        r[j] = cap.vp[j] - (Rsur*Ifem[j] + Lsur*dIdt);
      }
      std::vector<double> voffNew(Nc);
      double dvoff = 0.0, vscale = 1e-30;
      for (size_t j=0;j<Nc;++j) {
        const double vcur = interp1(cap.ts[j], dev.voffT, dev.voffV);
        voffNew[j] = vcur + cfg.theta*(r[j]-vcur);
        dvoff = std::max(dvoff, std::fabs(voffNew[j]-vcur));
        vscale = std::max(vscale, std::fabs(voffNew[j]));
      }
      voffT = cap.ts; voffV = voffNew;
      printf("  %4d    %12.4e    %12.4e    %10.3e\n", s, maxerr, ssErr, dvoff/vscale);
      if (dvoff < cfg.tol * vscale) break;
    }
    const double tWall = secs(t0wall, Clock::now());
    printf("\nGLOBAL-WR finished: %ld sweeps, max|i-I_FEM|=%.4e (global), %.4e (steady-state)\n",
           femEvals, maxerr, lastSS);
    printf("TIMING: wall=%.3f s  (Xyce=%.3f s, FEM=%.3f s, rest=%.3f s)\n",
           tWall, tXyce, tFem, tWall - tXyce - tFem);
    return 0;
  }

  // =========================================================================================
  // (A) L1 ITERATIVE WR via per-window CHECKPOINT/RESTART (rewind).
  // Each window [t0,t1] is a self-contained Xyce run (runWindow). Within a window the correction
  // voff is iterated to the transmission fixpoint by RE-RUNNING the window from its start state --
  // the "rewind" needed for iteration.
  //   voff^(k+1) = voff^(k) + theta*(r^(k) - voff^(k)),  r = V_p - [Rsur*I_FEM + Lsur*dI_FEM/dt]
  // recomputed on the SAME window (no extrapolation lag). Converges geometrically to the WR
  // fixpoint (||dvoff||->0); the transmission residual then floors at the discretization level
  // (coarse captured grid + per-window order-1 NOOP restart), which WR iterations cannot lower.
  //
  // CHECKPOINT/RESTART NOTE: Xyce's NATIVE .OPTIONS RESTART does NOT round-trip the YGENEXT
  // internal branch current (verified by inspecting an unpacked checkpoint: the field current is
  // absent; cf. AppNote 3.6). So the rewind/state-carry is DRIVER-MANAGED: we checkpoint by reading
  // the window-end state (getSolution -> Vp, i; the only memory is the series inductor current) and
  // restart by reconstructing it in a fresh run via NOOP + Ls IC + .IC on the field internal node
  // (see writeGxNetlist). max_iter=1 reduces to pure restart-based multirate (no in-window iter).
  // =========================================================================================
  if (cfg.iterate == 1) {
    const auto t0wall = Clock::now();
    const double w_src = 2.0 * M_PI * freq;
    printf("ITERATE mode: max_iter=%d tol=%.1e theta=%.3g  (driver-managed checkpoint/restart, NOOP+.IC)\n",
           cfg.max_iter, cfg.tol, cfg.theta);

    double Vp_start = 0.0, I_field_start = 0.0;          // carried across windows
    double iBranch_start = 0.0;                          // circuit branch current at the boundary (.IC seed)
    double curR = 1e-3, curL = 1e-6, curC = 0.0, curI0 = 0.0, curLam0 = 0.0; // black-box flux model
    if (blackbox) {
      std::vector<double> pts(n_sub+1), pvp(n_sub+1);
      for (int j=0;j<=n_sub;++j){ pts[j]=dt_window*(double)j/n_sub; pvp[j]=amp*w_src*pts[j]; }
      std::vector<double> plam;
      const std::vector<double> pI = femEvaluate(pts,pvp,R_FEM,L_FEM,0.0,saturating,Isat,&plam);
      double Rd,Ls0,cv; probeJacobian(pts,pvp,pI,plam,Rd,Ls0,cv);
      curR=Rd; curL=Ls0; curC=cv;
      printf("pre-probe: R_diff=%.3e L_diff(0)=%.3e curv=%.3e\n",curR,curL,curC);
    }

    std::vector<double> voffT, voffV;   // warm-start correction for the upcoming window (predictor)
    printf("\n  window    t_end      iters   max|i-I_FEM|     i_end      I_FEM_end\n");
    long totBack = 0; int win = 0;
    for (double t0 = 0.0; t0 < t_final - 1e-12; t0 += dt_window, ++win) {
      const double t1 = std::min(t0 + dt_window, t_final);
      const bool resume = (win > 0);                    // window 0 fresh; else reconstruct from t0 state

      DevCfg dev;
      dev.Rrom = Rrom; dev.Lrom = Lrom; dev.sat = saturating; dev.Isat = Isat;
      dev.voffT = voffT; dev.voffV = voffV;             // predictor warm start (abs time)

      WinCap cap; std::vector<double> Ifem, lam, r;
      double maxerr = 0.0, Rsur = Rrom, Lsur = Lrom;
      int k = 0;
      for (k = 0; k < cfg.max_iter; ++k) {
        if (blackbox) { dev.probed=true; dev.pR=curR; dev.pL=curL; dev.pCurv=curC; dev.pI0=curI0; dev.pLam0=curLam0; }
        cap = runWindow(netlist, cfg, t0, t1, dt_window/n_sub, n_sub, resume,
                        Vp_start, I_field_start, dev, iBranch_start);
        if (!cap.ok) { fprintf(stderr, "runWindow failed (window %d iter %d)\n", win, k); return 1; }
        totBack += cap.backJumps;

        // Resample onto a FIXED uniform window grid. Xyce's adaptive stepper lands on slightly
        // different t_done points each iteration (different voff -> different LTE), which jitters
        // the FD residual ~1% and prevents the WR fixpoint from settling. A fixed grid makes the
        // residual operator identical every iteration -> voff converges cleanly.
        { std::vector<double> tg(n_sub+1), vpg(n_sub+1), icg(n_sub+1);
          for (int j=0;j<=n_sub;++j){ const double tt=t0+(t1-t0)*(double)j/n_sub;
            tg[j]=tt; vpg[j]=interp1(tt,cap.ts,cap.vp); icg[j]=interp1(tt,cap.ts,cap.ic); }
          cap.ts=tg; cap.vp=vpg; cap.ic=icg; }

        Ifem = femEvaluate(cap.ts, cap.vp, R_FEM, L_FEM, I_field_start, saturating, Isat, &lam);
        maxerr = 0.0;
        for (size_t j=0;j<cap.ts.size();++j) maxerr = std::max(maxerr, std::fabs(cap.ic[j]-Ifem[j]));

        Rsur = Rrom; Lsur = Lrom;
        if (blackbox) { double Rd,Ls0,cv; probeJacobian(cap.ts,cap.vp,Ifem,lam,Rd,Ls0,cv);
                        Rsur=Rd; Lsur=Ls0 + cv*cap.ic.back(); }

        const size_t Nc = cap.ts.size();
        r.assign(Nc, 0.0);
        for (size_t j=0;j<Nc;++j) {
          const double dIdt = (j==0)?(Ifem[1]-Ifem[0])/(cap.ts[1]-cap.ts[0])
                                    :(Ifem[j]-Ifem[j-1])/(cap.ts[j]-cap.ts[j-1]);
          r[j] = cap.vp[j] - (Rsur*Ifem[j] + Lsur*dIdt);
        }

        // WR update voff += theta*(r - voff), and its change ||dvoff|| = the WR convergence measure.
        // (We converge on the voff FIXPOINT, not on the transmission maxerr: maxerr floors at the
        //  discretization level -- coarse captured grid + per-window order-1 NOOP restart -- which
        //  no number of WR iterations can lower. ||dvoff||->0 is the true "WR has converged" test.)
        std::vector<double> voffNew(Nc);
        double dvoff = 0.0, vscale = 1e-30;
        for (size_t j=0;j<Nc;++j) {
          const double vcur = interp1(cap.ts[j], dev.voffT, dev.voffV);
          voffNew[j] = vcur + cfg.theta*(r[j]-vcur);
          dvoff = std::max(dvoff, std::fabs(voffNew[j]-vcur));
          vscale = std::max(vscale, std::fabs(voffNew[j]));
        }
        dev.voffT = cap.ts; dev.voffV = voffNew;

        if (getenv("GENEXT_DEBUG") && win == (getenv("GENEXT_DBGWIN")?atoi(getenv("GENEXT_DBGWIN")):2))
          fprintf(stderr, "   [w%d it%d] transm.maxerr=%.4e  ||dvoff||=%.3e (rel %.2e)  voff[end]=%.4e\n",
                  win, k, maxerr, dvoff, dvoff/vscale, voffNew.back());
        if (dvoff < cfg.tol * vscale) break;            // WR fixpoint reached (voff stopped moving)
      }

      const int niter = (k < cfg.max_iter) ? k+1 : cfg.max_iter;  // iterations actually run
      printf("  %6.4f  %9.3e   %5d   %12.4e   %9.3f  %9.3f\n",
             t1, cap.ts.back(), niter, maxerr, cap.ic.back(), Ifem.back());

      // Predictor warm start for the NEXT window from the converged residual (end value + slope).
      const size_t Nc = r.size();
      const double r_end = r.back();
      double r_slope = 0.0; { const double h = cap.ts[Nc-1]-cap.ts[Nc-2]; if (h>0) r_slope=(r[Nc-1]-r[Nc-2])/h; }
      const double t2 = std::min(t1 + dt_window, t_final);
      voffT.assign(n_sub+1, 0.0); voffV.assign(n_sub+1, 0.0);
      for (int j=0;j<=n_sub;++j){ const double tt=t1+(t2-t1)*(double)j/n_sub; voffT[j]=tt; voffV[j]=r_end+r_slope*(tt-t1); }

      if (blackbox) { const double iB=cap.ic.back(), dB=iB-curI0;
                      const double Qb=curLam0 + curL*dB + 0.5*curC*dB*dB;
                      double Rd,Ls0,cv; probeJacobian(cap.ts,cap.vp,Ifem,lam,Rd,Ls0,cv);
                      curR=Rd; curL=Ls0+cv*iB; curC=cv; curI0=iB; curLam0=Qb; }

      Vp_start = cap.vp.back(); I_field_start = Ifem.back();
      iBranch_start = cap.ic.back();   // circuit branch current -> .IC seed for next window's restart
    }
    printf("\nGenExt ITERATE finished. (backward-jump retries handled: %ld)\n", totBack);
    printf("TIMING: wall=%.3f s\n", secs(t0wall, Clock::now()));
    return 0;
  }

  // ----------------------------- deferred single-pass mode (default) -----------------------
  // Generate the netlist from config (source, Rs, Ls, the YGENEXT field, .tran).
  const auto t0wall = Clock::now();
  writeGxNetlist(netlist, cfg, t_final, dt_window / n_sub);

  Xyce::Circuit::GenCouplingSimulator xyce;
  char prog[] = "xyce";
  std::vector<char> nlbuf(netlist, netlist + std::string(netlist).size() + 1);
  char * xargv[2] = { prog, nlbuf.data() };

  Xyce::Circuit::Simulator::RunStatus rs = xyce.initializeEarly(2, xargv);
  if (rs == Xyce::Circuit::Simulator::ERROR) { fprintf(stderr, "initializeEarly ERROR\n"); return 1; }
  if (rs == Xyce::Circuit::Simulator::DONE)  { return 0; }

  std::vector<std::string> names;
  if (!xyce.getDeviceNames("YGENEXT", names) || names.empty()) {
    fprintf(stderr, "No YGENEXT devices found in %s\n", netlist);
    return 1;
  }
  const std::string dev = names.front();

  FieldRomVCI vci(Rrom, Lrom, saturating, Isat);
  printf("Attaching field ROM to YGENEXT device: %s (Rrom=%.3e Lrom=%.3e sat=%d Isat=%.1f blackbox=%d)\n",
         dev.c_str(), Rrom, Lrom, (int)saturating, Isat, (int)blackbox);
  xyce.setNumInternalVars(dev, 1);
  xyce.setJacStamp(dev, vci.jacStamp);
  xyce.setVectorLoader(dev, &vci);

  rs = xyce.initializeLate();
  if (rs == Xyce::Circuit::Simulator::ERROR) { fprintf(stderr, "initializeLate ERROR\n"); return 1; }

  // --- multirate driver loop ---
  double t = 0.0;
  double I_field_start = 0.0; // FEM window-start current (carried)
  double Vp_start = 0.0;      // port voltage at window start (carried; 0 at t=0 with UIC)
  // Black-box device flux model, tracked for Q-continuity across windows. Q = curLam0 +
  // curL*(i-curI0) + 0.5*curC*(i-curI0)^2. The offset is voltage-irrelevant; only continuity
  // in the CIRCUIT current matters (else a Q jump -> dQ/dt spike -> Xyce abort).
  double curR = 1.0e-3, curL = 1.0e-6, curC = 0.0, curI0 = 0.0, curLam0 = 0.0;

  // COLD-START PRE-PROBE: before the run, probe the FEM on a synthetic V_p ramp (the source's
  // small-t slope, V_p ~ amp*2*pi*f*t) to get an initial Jacobian -> avoids the ~6-window
  // cold-start transient that a generic guess would cause.
  if (blackbox) {
    const double w = 2.0 * M_PI * freq;
    std::vector<double> pts(n_sub + 1), pvp(n_sub + 1);
    for (int j = 0; j <= n_sub; ++j) { pts[j] = dt_window * (double)j / n_sub; pvp[j] = amp * w * pts[j]; }
    std::vector<double> plam;
    const std::vector<double> pI = femEvaluate(pts, pvp, R_FEM, L_FEM, 0.0, saturating, Isat, &plam);
    double Rd, Ls0, cv; probeJacobian(pts, pvp, pI, plam, Rd, Ls0, cv);
    curR = Rd; curL = Ls0; curC = cv; curI0 = 0.0; curLam0 = 0.0;
    printf("pre-probe: R_diff=%.3e L_diff(0)=%.3e curv=%.3e\n", curR, curL, curC);
    vci.setProbed(curR, curL, curI0, curLam0, curC); // initial probed device model (black-box only)
  }
  // (known-param mode keeps the analytic flux(i); setProbed is black-box only.)

  printf("\n  window     t_end      max|i_circ-I_FEM|   i_end      I_FEM_end\n");
  while (t < t_final - 1e-12) {
    const double t0 = t;
    const double t1 = std::min(t0 + dt_window, t_final);

    // Window-start point is carried (getSolution is only valid AFTER a simulateUntil).
    std::vector<double> ts, vp, ic;
    ts.push_back(t0); vp.push_back(Vp_start); ic.push_back(I_field_start);
    for (int s = 1; s <= n_sub; ++s) {
      const double t_req = t0 + (t1 - t0) * (double)s / (double)n_sub;
      double t_done = t0;
      if (!xyce.simulateUntil(t_req, t_done)) { fprintf(stderr, "simulateUntil failed near %.6e\n", t_req); xyce.finalize(); return 1; }
      std::vector<double> sV;
      if (xyce.getSolution(dev, sV) && sV.size() >= 3) { ts.push_back(t_done); vp.push_back(sV[0]); ic.push_back(sV[2]); }
      t = t_done;
      if (t_done < t_req - 1e-15) break; // netlist final time reached early
    }
    if (ts.size() < 2) break;

    // FEM over the captured V_p(t)  (ONE call per window -- multirate). Black box: V_p -> I,lambda.
    std::vector<double> lam;
    const std::vector<double> Ifem = femEvaluate(ts, vp, R_FEM, L_FEM, I_field_start, saturating, Isat, &lam);

    // Transmission error (circuit i vs field I_FEM) on the captured grid.
    double maxerr = 0.0;
    for (size_t j = 0; j < ts.size(); ++j) maxerr = std::max(maxerr, std::fabs(ic[j] - Ifem[j]));

    // Surrogate coefficients for the correction/device. Known-param: (Rrom, Lrom). Black-box:
    // PROBE (R_diff, L_diff) from the FEM outputs and set the device's per-window linearization.
    double Rsur = Rrom, Lsur = Lrom;
    if (blackbox) {
      double Rdiff, Lslope0, curv;
      probeJacobian(ts, vp, Ifem, lam, Rdiff, Lslope0, curv);
      const double iB = ic.back();                     // circuit current at the window boundary
      const double Ldiff_iB = Lslope0 + curv * iB;     // differential inductance dlambda/dI at iB
      Rsur = Rdiff; Lsur = Ldiff_iB;
      // Keep Q continuous in the circuit current at the boundary: new offset = old model's Q at
      // iB (value-irrelevant for voltage; a jump would spike dQ/dt). Quadratic slope+curvature
      // update. curOff/curL/curC track the old model to evaluate its Q at iB.
      const double dB = iB - curI0;
      const double Qb = curLam0 + curL * dB + 0.5 * curC * dB * dB; // old device Q at boundary
      curR = Rdiff; curL = Ldiff_iB; curC = curv; curI0 = iB; curLam0 = Qb;
      vci.setProbed(curR, curL, curI0, curLam0, curC);
    }

    // Deferred correction with a LINEAR PREDICTOR. Residual of the field vs the surrogate along
    // the FEM trajectory:  r(t) = V_p - [Rsur*I_FEM + Lsur*dI_FEM/dt]  (dI_FEM/dt by FD).
    // Known-param (Lsur=L_FEM): reduces to (R_FEM-Rrom)*I_FEM. Black-box: residual of the probe.
    // Extrapolate r across the NEXT window from this window's exit (end value + end slope).
    // Backward FD here (NOT central): deferred is one window time-lagged and not at the WR
    // fixpoint, so its accuracy relies on the backward-FD/lag alignment -- central diff (BDF2-
    // consistent) actually worsens it. Central diff belongs in the GLOBAL mode, which reaches
    // the true fixpoint where matching Xyce's BDF2 cancels the (Lrom/Rrom)*(FD-BDF) floor.
    const size_t Nc = ts.size();
    std::vector<double> r(Nc);
    for (size_t j = 0; j < Nc; ++j) {
      const double dIdt = (j == 0) ? (Ifem[1]-Ifem[0])/(ts[1]-ts[0])
                                   : (Ifem[j]-Ifem[j-1])/(ts[j]-ts[j-1]);
      r[j] = vp[j] - (Rsur * Ifem[j] + Lsur * dIdt);
    }
    const double t_end = ts.back();
    const double r_end = r.back();
    double r_slope = 0.0;
    { const double h = ts[Nc-1] - ts[Nc-2]; if (h > 0.0) r_slope = (r[Nc-1] - r[Nc-2]) / h; }
    const double t_nextEnd = std::min(t_end + dt_window, t_final);
    std::vector<double> nts(n_sub + 1), ncv(n_sub + 1);
    for (int j = 0; j <= n_sub; ++j) {
      const double tt = t_end + (t_nextEnd - t_end) * (double)j / (double)n_sub;
      nts[j] = tt;
      ncv[j] = r_end + r_slope * (tt - t_end);
    }
    vci.setCorrection(nts, ncv);

    printf("  %6.4f  %9.3e   %12.4e    %9.3f  %9.3f\n",
           t1, t, maxerr, ic.back(), Ifem.back());

    I_field_start = Ifem.back(); // carry field current to next window
    Vp_start = vp.back();        // carry port voltage to next window start
    if (t <= t0 + 1e-18) break;  // no progress
  }

  xyce.finalize();
  printf("\nGenExt driver finished. (device backward-jump retries handled: %ld)\n",
         vci.backwardJumps());
  printf("TIMING: wall=%.3f s\n", secs(t0wall, Clock::now()));
  return 0;
}
