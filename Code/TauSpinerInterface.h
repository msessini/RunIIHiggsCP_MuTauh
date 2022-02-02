#ifndef TauSpinerInterface_h
#define TauSpinerInterface_h

#include <vector>

#ifdef USE_TauSpinner
#include "SimpleParticle.h"
using namespace TauSpinner;
#endif

class TauSpinerInterface {

 public:
  TauSpinerInterface();
  ~TauSpinerInterface();

  enum TauSpinerType {Spin=0,UnSpin,FlipSpin,hplus,hminus};
  enum TauSpinerSignalCharge {tauminus=-1,tauplus=1};
#ifdef USE_TauSpinner
  double Get(char* CPstate, int type, SimpleParticle X, SimpleParticle tau, std::vector<SimpleParticle> tau_daughters,SimpleParticle tau2, std::vector<SimpleParticle> tau_daughters2);
  void SetTauSignalCharge(int tsc){signalcharge=tsc;}

 private:
  void Initialize(char* CPstate);
  static int signalcharge;
  static bool initialized;
#endif
};
#endif
