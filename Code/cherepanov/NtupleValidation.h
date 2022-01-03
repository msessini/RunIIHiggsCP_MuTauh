#ifndef NtupleValidation_h
#define NtupleValidation_h

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
class NtupleValidation : public Selection {

 public:
  NtupleValidation(TString Name_, TString id_);
  virtual ~NtupleValidation();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0,PrimeVtx,NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables and Histos
  std::vector<TH1D> DaughtersPt;
  std::vector<TH1D> MissingTEnergy;
  std::vector<TH1D> NumVertices;
  std::vector<TH1D> particletype;
  std::vector<TH1D> taudecaytype;


  std::vector<TH1D> PVSVSignificance;
  std::vector<TH1D> SVchi2;
  std::vector<TH1D> SVMatchingQuality;
  std::vector<TH1D> SVz;
  std::vector<TH1D> PVz;


  std::vector<TH1D> isPairCandOS;
  std::vector<TH1D> Pair_part1Type;
  std::vector<TH1D> Pair_part2Type;
  std::vector<TH1D> OSPairMass;
  std::vector<TH1D> SSPairMass;

  std::vector<TH1D> s12;
  std::vector<TH1D> s13;
  std::vector<TH1D> s23;

  std::vector<TH1D> s12reco;
  std::vector<TH1D> s13reco;
  std::vector<TH1D> s23reco;
  std::vector<TH1D> TauA1PtResolution;
  std::vector<TH1D> TauA1PtResolutionPVSV;


  std::vector<TH1D> OSPionPtResolution;
  std::vector<TH1D> SSPionPtResolution;


  //------------- Truth Gen ---------
  std::vector<TH1D> TruthTauTauMass;



};
#endif
