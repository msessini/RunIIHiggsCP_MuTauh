#ifndef ControlSample_h
#define ControlSample_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "SVFitStorage.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/MultiProngTauSolver.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "SimpleFits/FitSoftware/interface/TauA1NuConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/DiTauConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/GlobalEventFit.h"
#include "Objects.h"
class ControlSample : public Selection {

 public:
  ControlSample(TString Name_, TString id_);
  virtual ~ControlSample();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0,
	     PrimeVtx,
	     nGoodTaus,
	     nGoodMuons,
	     MuonIsolation,
	     PairCharge,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

  int TriggerOkDummy, selVertexDummy, selMuon_IsoDummy, selMuon_AntiIsoDummy, selTauDummy, ChargeSumDummy;
  double MTDummy, MvisDummy, TauFLSigmaDummy;

  int Charge;

  double  cMu_pt,  cMu_eta,  cTau_pt,  cTau_eta;
 private:
  // Selection Variables and Histos

  std::vector<TH1D> TauPT;
  std::vector<TH1D> TauE;
  std::vector<TH1D> MuonPT;
  std::vector<TH1D> TauMuMass;
  std::vector<TH1D> TauHPSDecayMode;
  std::vector<TH1D> MissingM;


};
#endif
