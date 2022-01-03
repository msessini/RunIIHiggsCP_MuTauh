#ifndef SimFit_h
#define SimFit_h

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
class SimFit : public Selection {

 public:
  SimFit(TString Name_, TString id_);
  virtual ~SimFit();

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
  std::vector<TH1D> TauMuPtResolution;


  std::vector<TH1D> TauA1EResolution;
  std::vector<TH1D> TauMuEResolution;

  std::vector<TH1D> TauA1EtaResolution;
  std::vector<TH1D> TauMuEtaResolution;

  std::vector<TH1D> TauA1PhiResolution;
  std::vector<TH1D> TauMuPhiResolution;





  std::vector<TH1D> TauA1PtResolutionPVSV;


  std::vector<TH1D> OSPionPtResolution;
  std::vector<TH1D> SSPionPtResolution;


  //------------- Truth Gen ---------
  std::vector<TH1D> TruthTauTauMass;
  std::vector<TH1D> EventFitMass;
  std::vector<TH1D> EventFitZPt;
  std::vector<TH1D> EventFitZE;
  std::vector<TH1D> EventFitZPhi;
  std::vector<TH1D> EventFitZEta;

  std::vector<TH1D> EventFitTauA1E;
  std::vector<TH1D> EventFitTauMuE;


  std::vector<TH1D> EventFitZPtResolution;
  std::vector<TH1D> EventFitZEResolution;
  std::vector<TH1D> EventFitZPhiResolution;
  std::vector<TH1D> EventFitZEtaResolution;

  std::vector<TH1D> SimpleFMuPtResolution;
  std::vector<TH1D> SimpleFMuEResolution;
  std::vector<TH1D> SimpleFMuPhiResolution;
  std::vector<TH1D> SimpleFMuEtaResolution;


  std::vector<TH1D> SimpleFA1PtResolution;
  std::vector<TH1D> SimpleFA1EResolution;
  std::vector<TH1D> SimpleFA1PhiResolution;
  std::vector<TH1D> SimpleFA1EtaResolution;


  std::vector<TH1D> TrackPtResolution;
  std::vector<TH1D> TrackEResolution;
  std::vector<TH1D> TrackPhiResolution;
  std::vector<TH1D> TrackEtaResolution;


  std::vector<TH1D> A1PtResolution;
  std::vector<TH1D> A1EResolution;
  std::vector<TH1D> A1PhiResolution;
  std::vector<TH1D> A1EtaResolution;

  std::vector<TH1D> TrackMatchingQuality;
  std::vector<TH2D> MQVsVisA1Resolution;
  std::vector<TH1D> A1VisiblePtResolutionMuA1;
  std::vector<TH1D> A1VisiblePtResolution;
  std::vector<TH1D> TrackVisiblePtResolution;


};
#endif
