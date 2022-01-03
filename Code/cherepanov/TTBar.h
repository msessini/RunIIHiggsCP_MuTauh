#ifndef TTBar_h
#define TTBar_h

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
#include "DataMCCorrections.h"
#include "PUReweight.h"
#include "tauTrigSFreader.h"
  
class TTBar : public Selection {

 public:
  TTBar(TString Name_, TString id_);
  virtual ~TTBar();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0,
	     nGoodTaus,
	     nGoodMuons,
	     nBJets,
	     nJets,
	     LeptonVeto,
	     ET,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

  int TriggerOkDummy, selVertexDummy, selMuon_IsoDummy, selMuon_AntiIsoDummy, selTauDummy, ChargeSumDummy;
  double MTDummy, MvisDummy, TauFLSigmaDummy;

  int Charge;

  double  cMu_pt,  cMu_eta,  cTau_pt,  cTau_eta;
  PUReweight reweight;//(PUReweight::RUN2ANALYSIS);
  ReferenceScaleFactors *RSF;
  tauTrigSFreader tauTrgSF;
  DataMCCorrections DataMC_Corr;

 private:
  // Selection Variables and Histos

  std::vector<TH1D> AllJetsBTag;
  std::vector<TH1D> BDiscr1;
  std::vector<TH1D> BDiscr2;
  std::vector<TH1D> BDiscr3;

  std::vector<TH1D> MTM;


 

 std::vector<TH1D> LeadintTauPt;
 std::vector<TH1D> LeadintMuonPt;
 std::vector<TH1D> LeadintJetPt;
 std::vector<TH1D> NPrimeVtx;

};
#endif
