#ifndef ZTauHTauH_h
#define ZTauHTauH_h

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
#include "Objects.h"
#include "PUReweight.h"
#include "tauTrigSFreader.h"
#include "DataMCCorrections.h"



class ZTauHTauH : public Selection {

 public:
  ZTauHTauH(TString Name_, TString id_);
  virtual ~ZTauHTauH();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0,
	     nGoodPairs,
	     LeptonVeto,
	     FirstTauIsolation,
	     SecondTauIsolation,
	     nGoodMuons,
	     PairCharge,
	     PairMass,
	     deltaR,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();
  ReferenceScaleFactors *RSF;
  //  tauTrigSFreader tauTrgSF;
  tauTrigSFreader tauTrgSF;
  DataMCCorrections DataMC_Corr;


  int TriggerOkDummy, selVertexDummy, selMuon_IsoDummy, selMuon_AntiIsoDummy, selTauDummy, ChargeSumDummy;
  double MTDummy, MvisDummy, TauFLSigmaDummy;

  int Charge;

  double  cMu_pt,  cMu_eta,  cTau_pt,  cTau_eta;
  PUReweight reweight;//(PUReweight::RUN2ANALYSIS);

 private:
  // Selection Variables and Histos

  std::vector<TH1D> Tau1PT;
  std::vector<TH1D> Tau1E;
  std::vector<TH1D> Tau1HPSDecayMode;
 
  std::vector<TH1D> Tau2PT;
  std::vector<TH1D> Tau2E;
  std::vector<TH1D> Tau2HPSDecayMode;
  std::vector<TH1D> TauTauMass;

  std::vector<TH1D> dRTauTau;
  std::vector<TH1D> QCDShape;

  std::vector<TH1D> NQCD;
  std::vector<TH1D> Tau1Isolation;
  std::vector<TH1D> Tau2Isolation;

  std::vector<TH1D> Tau1Eta;
  std::vector<TH1D> Tau2Eta;


  std::vector<TH1D> MET;
  std::vector<TH1D> TauHMass1;
  std::vector<TH1D> TauHMass2;

  std::vector<TH1D> NPrimeVtx;

};
#endif
