#ifndef ZTauMuTauH_h
#define ZTauMuTauH_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "SVFitStorage.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/MultiProngTauSolver.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "SimpleFits/FitSoftware/interface/TauA1NuConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/DiTauConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/GlobalEventFit.h"
#include "ReferenceScaleFactors.h"
#include "ScaleFactor.h"
#include "Objects.h"
#include "PUReweight.h"
#include "tauTrigSFreader.h"
#include "DataMCCorrections.h"


class ZTauMuTauH : public Selection {

 public:
  ZTauMuTauH(TString Name_, TString id_);
  virtual ~ZTauMuTauH();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0,
	     nGoodPairs,
	     LeptonVeto,
	     MuonIsolation,
	     TauIsolation,
	     PairCharge,
	     PairMass,
	     deltaR,
	     MTM,
	     NCuts};
 
 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();
  ReferenceScaleFactors *RSF;

  int TriggerOkDummy, selVertexDummy, selMuon_IsoDummy, selMuon_AntiIsoDummy, selTauDummy, ChargeSumDummy;
  double MTDummy, MvisDummy, TauFLSigmaDummy;

  int Charge;

  double  cMu_pt,  cMu_eta,  cTau_pt,  cTau_eta;
  PUReweight reweight;//(PUReweight::RUN2ANALYSIS);
  //  tauTrigSFreader tauTrgSF;
  DataMCCorrections DataMC_Corr;
  DataMCCorrections DataMC_CorrLeptonIso;

 private:
  // Selection Variables and Histos

  std::vector<TH1D> TauPT;
  std::vector<TH1D> TauE;
 
  std::vector<TH1D> MuonPT;
  std::vector<TH1D> MuonE;
  std::vector<TH1D> TauHPSDecayMode;
  std::vector<TH1D> TauTauMass;

  std::vector<TH1D> dRTauTau;
  std::vector<TH1D> QCDShape;

  std::vector<TH1D> NQCD;

  std::vector<TH1D> NPrimeVtx;
  std::vector<TH1D> PVSVSignificance;
  std::vector<TH1D> SVChi2;
  std::vector<TH1D> SVQuality;
  std::vector<TH2D> SVQualityVsSignificance;


};
#endif
