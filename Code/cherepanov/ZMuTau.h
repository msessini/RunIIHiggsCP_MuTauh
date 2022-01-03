#ifndef ZMuTau_h
#define ZMuTau_h

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


class ZMuTau : public Selection {

 public:
  ZMuTau(TString Name_, TString id_);
  virtual ~ZMuTau();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {Trigger=0,
	     Id_and_Kin,
	     NPairsFound, 
	     MuonIsolation,
	     TauIsolation, 
	     DiMuon_Veto,
	     LeptonVeto,
	     PairCharge,
	     MTM,
	     // PairMass,
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
  std::vector<TH1D> TauMass;
  std::vector<TH1D> TauPhi;
  std::vector<TH1D> TauEta;
  std::vector<TH1D> Taudz;

  std::vector<TH1D> MuonPT;
  std::vector<TH1D> MuonE;
  std::vector<TH1D> MuonMass;
  std::vector<TH1D> MuonPhi;
  std::vector<TH1D> MuonEta;
  std::vector<TH1D> Muondz;
  std::vector<TH1D> Muondxy;
  std::vector<TH1D> MuonIsol;

  std::vector<TH1D> againstElectronVLooseMVA6;
  std::vector<TH1D> againstElectronLooseMVA6;
  std::vector<TH1D> againstElectronMediumMVA6;
  std::vector<TH1D> againstElectronTightMVA6;
  std::vector<TH1D> againstElectronVTightMVA6;
  std::vector<TH1D> againstMuonLoose3;
  std::vector<TH1D> againstMuonTight3;
  std::vector<TH1D> byCombinedIsolationDeltaBetaCorrRaw3Hits;

  std::vector<TH1D> DiMuonVeto;
  std::vector<TH1D> ExtraLeptonVeto;
  
  std::vector<TH1D> TauHPSDecayMode;
  std::vector<TH1D> TauTauMass;

  std::vector<TH1D> dRTauTau;
  std::vector<TH1D> QCDShape;

  std::vector<TH1D> NQCD;
  std::vector<TH1D> NWJets;
  std::vector<TH1D> NWJetsRelaxed;
  
  std::vector<TH1D> MET;
  std::vector<TH1D> METphi;
  std::vector<TH1D> PUPPImet;
  std::vector<TH1D> PUPPImetphi;
  std::vector<TH1D> TransverseMass;
  
  std::vector<TH1D> NPrimeVtx;
  std::vector<TH1D> NPU;
  std::vector<TH1D> RHO;
  
  std::vector<TH1D> SVChi2;
  std::vector<TH1D> SVQuality;
  std::vector<TH2D> SVQualityVsSignificance;
  std::vector<TH1D> PVSVSignificance;
  std::vector<TH1D> NbJets;


};
#endif
