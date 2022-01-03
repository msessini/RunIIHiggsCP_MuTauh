#ifndef HTauTau_h
#define HTauTau_h

#include "Selection.h"
#include <vector>
#include "TString.h"
//#include "SVFitStorage.h"
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
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"

class HTauTau : public Selection {

 public:
  HTauTau(TString Name_, TString id_);
  virtual ~HTauTau();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {Trigger=0,
	     Id_and_Kin,
	     NPairsFound,
	     Tau1Isolation,
	     Tau2Isolation,
	     LeptonVeto,
	     PairCharge,
	     PairMass,
	     //MTM,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();
  ReferenceScaleFactors *RSF;
  int TriggerOkDummy, selVertexDummy, selMuon_IsoDummy, selMuon_AntiIsoDummy, selTauDummy, ChargeSumDummy;
  double MTDummy, MvisDummy, TauFLSigmaDummy;

  int Charge;

  PUReweight reweight;//(PUReweight::RUN2ANALYSIS);
  DataMCCorrections DataMC_Corr;
  tauTrigSFreader tauTrgSF;

  ClassicSVfit svfitAlgo1;
  //ClassicSVfit svfitAlgo2;
  //  SVFitStorage svfitstorage;


 private:
  // Selection Variables and Histos

  std::vector<TH1D> Tau1PT;
  std::vector<TH1D> Tau1E;
  std::vector<TH1D> Tau1Mass;
  std::vector<TH1D> Tau1Phi;
  std::vector<TH1D> Tau1Eta;
  std::vector<TH1D> Tau1dz;

  std::vector<TH1D> Tau2PT;
  std::vector<TH1D> Tau2E;
  std::vector<TH1D> Tau2Mass;
  std::vector<TH1D> Tau2Phi;
  std::vector<TH1D> Tau2Eta;
  std::vector<TH1D> Tau2dz;

  std::vector<TH1D> Tau1isolation;
  std::vector<TH1D> Tau2isolation;

  std::vector<TH1D> againstElectronVLooseMVA6_Tau1;
  std::vector<TH1D> againstElectronLooseMVA6_Tau1;
  std::vector<TH1D> againstElectronMediumMVA6_Tau1;
  std::vector<TH1D> againstElectronTightMVA6_Tau1;
  std::vector<TH1D> againstElectronVTightMVA6_Tau1;
  std::vector<TH1D> againstMuonLoose3_Tau1;
  std::vector<TH1D> againstMuonTight3_Tau1;
  std::vector<TH1D> byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau1;

  std::vector<TH1D> againstElectronVLooseMVA6_Tau2;
  std::vector<TH1D> againstElectronLooseMVA6_Tau2;
  std::vector<TH1D> againstElectronMediumMVA6_Tau2;
  std::vector<TH1D> againstElectronTightMVA6_Tau2;
  std::vector<TH1D> againstElectronVTightMVA6_Tau2;
  std::vector<TH1D> againstMuonLoose3_Tau2;
  std::vector<TH1D> againstMuonTight3_Tau2;
  std::vector<TH1D> byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau2;

  std::vector<TH1D> ExtraLeptonVeto;
  std::vector<TH1D> Tau2HPSDecayMode;
  std::vector<TH1D> Tau1HPSDecayMode;

  std::vector<TH1D> TauTauVisMass;
  std::vector<TH1D> TauTauTruthMass;
  std::vector<TH1D> TauTauFullMass;
 
  std::vector<TH1D> dRTauTau;
  std::vector<TH1D> QCDShape;

  std::vector<TH1D> NQCD;

  std::vector<TH1D> MET;
  std::vector<TH1D> METphi;
  std::vector<TH1D> PUPPImet;
  std::vector<TH1D> PUPPImetphi;
  std::vector<TH1D> TransverseMass;
  
  std::vector<TH1D> NPrimeVtx;
  std::vector<TH1D> NPU;
  std::vector<TH1D> RHO;
  
  std::vector<TH1D> NbJets;
  
  std::vector<TH1D> h_SVFitMass;
  std::vector<TH1D> h_SVFitStatus;
  std::vector<TH1D> svfTau1E;
  std::vector<TH1D> svfTau2E;

  std::vector<TH1D> Etasvfit;
  std::vector<TH1D> Phisvfit;
  std::vector<TH1D> Thetasvfit;
  
  std::vector<TH1D> Etavis;
  std::vector<TH1D> Phivis;
  std::vector<TH1D> Thetavis;
  
  std::vector<TH1D> Etatruth;
  std::vector<TH1D> Phitruth;
  std::vector<TH1D> Thetatruth;
  
  std::vector<TH1D> TauTauFullPtRes;
  std::vector<TH1D> TauTauFullEtaRes;
  std::vector<TH1D> TauTauFullPhiRes;

  std::vector<TH1D> TauplusFullPtRes;
  std::vector<TH1D> TauplusFullEtaRes;
  std::vector<TH1D> TauplusFullPhiRes;

  std::vector<TH1D> TauminusFullPtRes;
  std::vector<TH1D> TauminusFullEtaRes;
  std::vector<TH1D> TauminusFullPhiRes;

  std::vector<TH1D> TauTauVisPtRes;
  std::vector<TH1D> TauTauVisEtaRes;
  std::vector<TH1D> TauTauVisPhiRes;

  std::vector<TH1D> TauplusVisPtRes;
  std::vector<TH1D> TauplusVisEtaRes;
  std::vector<TH1D> TauplusVisPhiRes;
    
  std::vector<TH1D> TauminusVisPtRes;
  std::vector<TH1D> TauminusVisEtaRes;
  std::vector<TH1D> TauminusVisPhiRes;

  std::vector<TH1D> TauTauFullPtResPull;
  std::vector<TH1D> TauTauFullEtaResPull;
  std::vector<TH1D> TauTauFullPhiResPull;

  std::vector<TH1D> TauplusFullPtResPull;
  std::vector<TH1D> TauplusFullEtaResPull;
  std::vector<TH1D> TauplusFullPhiResPull;

  std::vector<TH1D> TauminusFullPtResPull;
  std::vector<TH1D> TauminusFullEtaResPull;
  std::vector<TH1D> TauminusFullPhiResPull;

  std::vector<TH1D> TauTauVisPtResPull;
  std::vector<TH1D> TauTauVisEtaResPull;
  std::vector<TH1D> TauTauVisPhiResPull;

  std::vector<TH1D> TauplusVisPtResPull;
  std::vector<TH1D> TauplusVisEtaResPull;
  std::vector<TH1D> TauplusVisPhiResPull;
    
  std::vector<TH1D> TauminusVisPtResPull;
  std::vector<TH1D> TauminusVisEtaResPull;
  std::vector<TH1D> TauminusVisPhiResPull;

  //std::vector<TH1D> DRTruth;
  //std::vector<TH1D> DRFull;
  // std::vector<TH1D> DRFullTruth;
  //std::vector<TH1D> DRVisTruth;
  
  std::vector<TH1D> Pi0EnergyRes;
  std::vector<TH1D> Pi0EnergyResPull;

  /*std::vector<TH1D> CTN;
    std::vector<TH1D> CTT;*/
};
#endif
