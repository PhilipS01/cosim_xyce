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
int main(int argc, char ** argv)
{
  const char * netlist = (argc > 1) ? argv[1] : "wr_genext.cir";

  // Field + ROM parameters. Lrom = L_FEM -> inductance exact & implicit;
  // Rrom deliberately != R_FEM so the deferred resistive correction has work to do.
  const double R_FEM = 5.1e-4;
  const double L_FEM = getenv("GENEXT_LFEM") ? atof(getenv("GENEXT_LFEM")) : 1.6e-7;
  const double Rrom  = 0.9 * R_FEM; // mismatched ROM resistance
  const double Lrom  = L_FEM;       // matched inductance (exact implicit)
  // Magnetic saturation of the field: lambda(I)=L_FEM*Isat*atan(I/Isat). Enable via env
  // GENEXT_SAT=1; Isat via GENEXT_ISAT (default 100 A; peak I ~150 A -> strong saturation).
  const bool   saturating = (getenv("GENEXT_SAT") && atoi(getenv("GENEXT_SAT")) != 0);
  const double Isat = getenv("GENEXT_ISAT") ? atof(getenv("GENEXT_ISAT")) : 100.0;
  // BLACK-BOX mode: the device does NOT use R_FEM/L_FEM; it PROBES (R_diff,L_diff) from the
  // FEM's (V_p -> I_FEM, lambda_FEM) outputs each window. (R_FEM/L_FEM still live inside
  // femEvaluate as the opaque "FEM", but the device/coupling never reads them.)
  const bool   blackbox = (getenv("GENEXT_BLACKBOX") && atoi(getenv("GENEXT_BLACKBOX")) != 0);

  // Source (for the cold-start pre-probe; must match the netlist Bemf).
  const double freq = 50.0, amp = 1.0;

  // Time stepping (must match the netlist .tran stop).
  const double t_final   = 0.02;
  const double dt_window = 4.0e-4;
  const int    n_sub     = 20;      // sub-steps per window to capture V_p(t) for the FEM

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
  }
  vci.setProbed(curR, curL, curI0, curLam0, curC); // initial device model

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
  return 0;
}
