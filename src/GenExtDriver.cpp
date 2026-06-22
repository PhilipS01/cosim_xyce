// GenExtDriver.cpp -- IMPLICIT field-circuit coupling via Xyce's General External
// Device (YGENEXT) + GenCouplingSimulator, replacing the partitioned system("Xyce")
// + PWL exchange. This skeleton drives ONE transient with the field represented as a
// real Xyce device whose flux Q = Lrom*i is a SHARED STATE (Xyce integrates dQ/dt
// itself -> no PWL/DDT staircase) and whose admittance dFdx/dQdx is in the Newton
// Jacobian (-> implicit, breaks the strong-Ls WR wall). See ../doc/AppNote-GenExt.pdf
// (3.3.1 computeXyceVectors, 4.1.x init + simulateUntil, 5.1 RLC example).
//
// This first step: get a window solving through libxyce with the field as a series
// R-L YGENEXT device (the linear ROM surrogate). The multirate FEM correction
// (offset updated once per simulateUntil window) is stubbed (setCorrection) for the
// next step.

#include <Xyce_config.h>
#include <N_CIR_GenCouplingSimulator.h>
#include <N_DEV_VectorComputeInterface.h>

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Field ROM as a 2-terminal series R-L device (the RLC example minus the C).
// Device variables: external nodes A(0), B(1); one internal branch unknown (2) =
// the interface current i. DAE (Xyce solves F(X) + dQ/dt - B = 0):
//   F[A]   =  i                      (current into A)
//   F[B]   = -i
//   F[br]  =  Rrom*i - (V_A - V_B) + vcorr   (branch KVL, resistive + FEM correction)
//   Q[br]  =  Lrom*i                 (Xyce differentiates: dQ/dt = Lrom*di/dt = inductive V)
// vcorr is the deferred FEM residual (lumped scalar here; becomes time-dependent per
// window once the FEM is wired in).
// -----------------------------------------------------------------------------
class FieldRomVCI : public Xyce::Device::VectorComputeInterface
{
public:
  FieldRomVCI(double Rrom, double Lrom)
    : Rrom_(Rrom), Lrom_(Lrom), vcorr_(0.0)
  {
    // jacStamp: which solution vars each row depends on (sparsity).
    jacStamp.resize(3);
    jacStamp[A_].resize(1);  jacStamp[A_][0] = BR_;                 // F[A] depends on i
    jacStamp[B_].resize(1);  jacStamp[B_][0] = BR_;                 // F[B] depends on i
    jacStamp[BR_].resize(3); jacStamp[BR_][0] = A_;                 // F[br] depends on V_A,
    jacStamp[BR_][1] = B_;                                          //   V_B,
    jacStamp[BR_][2] = BR_;                                         //   i
  }

  // Updated once per simulateUntil window from the FEM solve (multirate). Stub for now.
  void setCorrection(double vcorr) { vcorr_ = vcorr; }

  bool computeXyceVectors(std::vector<double> & sV,
                          double /*time*/,
                          std::vector<double> & F,
                          std::vector<double> & Q,
                          std::vector<double> & B,
                          std::vector<std::vector<double> > & dFdx,
                          std::vector<std::vector<double> > & dQdx) override
  {
    const int n = static_cast<int>(sV.size()); // expect 3 (A, B, branch)
    F.assign(n, 0.0);
    Q.assign(n, 0.0);
    B.clear();
    dFdx.assign(n, std::vector<double>(n, 0.0));
    dQdx.assign(n, std::vector<double>(n, 0.0));

    const double i = sV[BR_];

    F[A_]  =  i;
    F[B_]  = -i;
    F[BR_] =  Rrom_ * i - (sV[A_] - sV[B_]) + vcorr_;

    Q[BR_] =  Lrom_ * i;

    dFdx[A_][BR_]  =  1.0;
    dFdx[B_][BR_]  = -1.0;
    dFdx[BR_][A_]  = -1.0;
    dFdx[BR_][B_]  =  1.0;
    dFdx[BR_][BR_] =  Rrom_;

    dQdx[BR_][BR_] =  Lrom_;

    return true;
  }

  std::vector<std::vector<int> > jacStamp;

private:
  static const int A_  = 0;
  static const int B_  = 1;
  static const int BR_ = 2;
  double Rrom_, Lrom_, vcorr_;
};

// -----------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  const char * netlist = (argc > 1) ? argv[1] : "wr_genext.cir";

  Xyce::Circuit::GenCouplingSimulator xyce;

  // initializeEarly takes a Xyce-style command line (progName + netlist).
  char prog[]   = "xyce";
  std::vector<char> nlbuf(netlist, netlist + std::string(netlist).size() + 1);
  char * xargv[2] = { prog, nlbuf.data() };

  Xyce::Circuit::Simulator::RunStatus rs = xyce.initializeEarly(2, xargv);
  if (rs == Xyce::Circuit::Simulator::ERROR) { fprintf(stderr, "initializeEarly ERROR\n"); return 1; }
  if (rs == Xyce::Circuit::Simulator::DONE)  { return 0; }

  // Find the YGENEXT field device(s) and attach the ROM vector loader.
  std::vector<std::string> names;
  if (!xyce.getDeviceNames("YGENEXT", names) || names.empty()) {
    fprintf(stderr, "No YGENEXT devices found in %s\n", netlist);
    return 1;
  }

  // ROM values (will come from sim_config later). Match the toy field.
  const double Rrom = 4.59e-4;
  const double Lrom = 1.44e-7;
  FieldRomVCI vci(Rrom, Lrom);

  for (const std::string & nm : names) {
    printf("Attaching field ROM to YGENEXT device: %s\n", nm.c_str());
    xyce.setNumInternalVars(nm, 1);          // the branch current i
    xyce.setJacStamp(nm, vci.jacStamp);
    xyce.setVectorLoader(nm, &vci);
  }

  rs = xyce.initializeLate();
  if (rs == Xyce::Circuit::Simulator::ERROR) { fprintf(stderr, "initializeLate ERROR\n"); return 1; }

  // --- multirate driver loop: advance Xyce one coupling window at a time ---
  // For this skeleton: a few fixed windows to the netlist final time, vcorr = 0
  // (no FEM yet). Next step: call the FEM once per window and vci.setCorrection(...).
  const double t_final  = 0.02;   // matches the netlist .tran stop
  const double dt_window = 4.0e-4;
  double t = 0.0;
  while (t < t_final - 1e-15) {
    const double t_req = std::min(t + dt_window, t_final);
    double t_done = t;
    const bool ok = xyce.simulateUntil(t_req, t_done);
    if (!ok) { fprintf(stderr, "simulateUntil failed at t_req=%.6e (reached %.6e)\n", t_req, t_done); break; }
    printf("window -> requested %.6e, completed %.6e\n", t_req, t_done);
    if (t_done <= t + 1e-18) break; // no progress -> netlist done
    t = t_done;
    // (next step) FEM(V_p over [t-dt_window, t]) -> vci.setCorrection(...)
  }

  xyce.finalize();
  printf("GenExt driver finished.\n");
  return 0;
}
