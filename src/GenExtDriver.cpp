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
    Q[BR_] =  flux(i);                 // nonlinear flux as a shared state

    dFdx[A_][BR_]  =  1.0;
    dFdx[B_][BR_]  = -1.0;
    dFdx[BR_][A_]  = -1.0;
    dFdx[BR_][B_]  =  1.0;
    dFdx[BR_][BR_] =  Rrom_;
    dQdx[BR_][BR_] =  Ldyn(i);         // differential inductance dlambda/di in the Jacobian
    return true;
  }

  long backwardJumps() const { return backwardJumps_; } // diagnostics (rejected-step retries)

  std::vector<std::vector<int> > jacStamp;

private:
  static const int A_ = 0, B_ = 1, BR_ = 2;
  double Rrom_, Lrom_;
  bool saturating_;
  double Isat_;
  std::vector<double> ct_, cv_; // correction waveform (time, value)
  mutable double lastTime_ = -1e300; // for backward-jump detection (diagnostic only)
  mutable long backwardJumps_ = 0;
};

// -----------------------------------------------------------------------------
// Lumped "FEM": given the port-voltage samples V_p(t) over a window, integrate the
// true field  V_p = R_FEM*I + L_FEM*dI/dt  with backward Euler, return I_FEM(t).
// (A real FEM would replace this; the interface is V_p in -> I_FEM, lambda out.)
// -----------------------------------------------------------------------------
static std::vector<double> femEvaluate(const std::vector<double> & t,
                                       const std::vector<double> & Vp,
                                       double R_FEM, double L_FEM, double I_start,
                                       bool saturating = false, double Isat = 100.0)
{
  // True field:  V_p = R_FEM*I + dlambda/dt,  lambda(I) = L_FEM*Isat*atan(I/Isat) (sat) or L_FEM*I.
  // Backward Euler:  R_FEM*I_j + (lambda(I_j) - lambda(I_{j-1}))/h - V_p = 0  (Newton for sat).
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
  return I;
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
  printf("Attaching field ROM to YGENEXT device: %s (Rrom=%.3e Lrom=%.3e sat=%d Isat=%.1f)\n",
         dev.c_str(), Rrom, Lrom, (int)saturating, Isat);
  xyce.setNumInternalVars(dev, 1);
  xyce.setJacStamp(dev, vci.jacStamp);
  xyce.setVectorLoader(dev, &vci);

  rs = xyce.initializeLate();
  if (rs == Xyce::Circuit::Simulator::ERROR) { fprintf(stderr, "initializeLate ERROR\n"); return 1; }

  // --- multirate driver loop ---
  double t = 0.0;
  double I_field_start = 0.0; // FEM window-start current (carried)
  double Vp_start = 0.0;      // port voltage at window start (carried; 0 at t=0 with UIC)
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

    // FEM over the captured V_p(t)  (ONE call per window -- multirate).
    const std::vector<double> Ifem = femEvaluate(ts, vp, R_FEM, L_FEM, I_field_start, saturating, Isat);

    // Transmission error (circuit i vs field I_FEM) on the captured grid.
    double maxerr = 0.0;
    for (size_t j = 0; j < ts.size(); ++j) maxerr = std::max(maxerr, std::fabs(ic[j] - Ifem[j]));

    // Deferred RESISTIVE correction with a LINEAR PREDICTOR. Instead of zero-order
    // holding this window's I_FEM into the next window (the deferred lag), extrapolate
    // the field-current trajectory across the NEXT window from this window's exit
    // (end value + end slope). voff(t) = (R_FEM - Rrom)*I_pred(t) over [t1, t1+dt].
    // (Inductive part is exact via Q=Lrom*i, Lrom=L_FEM -> still derivative-free here.)
    const double t_end = ts.back();
    const double I_end = Ifem.back();
    double slope = 0.0;
    {
      const double h = ts[ts.size()-1] - ts[ts.size()-2];
      if (h > 0.0) slope = (Ifem[Ifem.size()-1] - Ifem[Ifem.size()-2]) / h;
    }
    const double t_nextEnd = std::min(t_end + dt_window, t_final);
    std::vector<double> nts(n_sub + 1), ncv(n_sub + 1);
    for (int j = 0; j <= n_sub; ++j) {
      const double tt = t_end + (t_nextEnd - t_end) * (double)j / (double)n_sub;
      const double Ipred = I_end + slope * (tt - t_end);
      nts[j] = tt;
      ncv[j] = (R_FEM - Rrom) * Ipred;
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
