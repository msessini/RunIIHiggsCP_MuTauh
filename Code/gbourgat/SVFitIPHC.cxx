#include "SVFitIPHC.h"
#include "TLorentzVector.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "SVFitObject.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "PDG_Var.h"
#include "SkimConfig.h"



#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TVector3.h"
#include "TMath.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/MultiProngTauSolver.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "SimpleFits/FitSoftware/interface/TauA1NuConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/DiTauConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/GlobalEventFit.h"
#include "Objects.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"



SVFitIPHC::SVFitIPHC(TString Name_, TString id_):
  Selection(Name_,id_),
  DataMC_Corr(true,true,false),
  tauTrgSF("tight")
{

}

SVFitIPHC::~SVFitIPHC(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  SVFitIPHC::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==Trigger)             cut.at(Trigger)=1;
    if(i==Id_and_Kin)          cut.at(Id_and_Kin)=1;
    if(i==NPairsFound)         cut.at(NPairsFound)=1;
    if(i==Tau1Isolation)       cut.at(Tau1Isolation)=1.;
    if(i==Tau2Isolation)       cut.at(Tau2Isolation)=1.;
    if(i==LeptonVeto)          cut.at(LeptonVeto)=0.;
    if(i==PairCharge)          cut.at(PairCharge)=1.;
    if(i==PairMass)            cut.at(PairMass)=100.;
    //if(i==MTM)                 cut.at(MTM)=40;
    
  }
  // Setup cut plots
  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;

    if(i==Trigger){
      title.at(i)="At least 1 good pair with Trig+Matching";
      hlabel="At least 1 good pair with Trig+Matching";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Trigger_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Trigger_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
     else if(i==Id_and_Kin){
      title.at(i)="Id and Kinematic";
      hlabel="Number of Event with good particles";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==NPairsFound){
      title.at(i)="Pairs with good DeltaR";
      hlabel="Pairs with good DeltaR";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NPairsFound_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NPairsFound_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
     else if(i==Tau1Isolation){
      title.at(i)="Tau1 Isolation";
      hlabel="Isolation of Tau1";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Tau1Isolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Tau1Isolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==Tau2Isolation){
      title.at(i)="Tau2 Isolation";
      hlabel="Isolation of Tau2";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Tau2Isolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Tau2Isolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==LeptonVeto){
      title.at(i)="Lepton Veto";
      hlabel="Third Lepton Veto  ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_LeptonVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_LeptonVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }
    else if(i==PairCharge){
      title.at(i)="Pair Charge";
      hlabel="is pair OS";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PairCharge_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PairCharge_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }
    else if(i==PairMass){
      title.at(i)="Pair Visible Mass";
      hlabel="M(tau-tau)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PairMass_",htitle,30,0,150,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PairMass_",htitle,30,0,150,hlabel,"Events"));
      }
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");

  Tau1PT=HConfig.GetTH1D(Name+"_Tau1PT","Transverse momentum of selected #tau1 candidate",20,35,80," P_{T}(#tau1), GeV","Events");
  Tau1E=HConfig.GetTH1D(Name+"_Tau1E","Energy of selected #tau1 candidate",20,35.,140," E(#tau1), GeV","Events");
  Tau1Mass=HConfig.GetTH1D(Name+"_Tau1Mass","Mass of selected #tau1 candidate",20,0,2.," M(#tau1), GeV","Events");
  Tau1Phi=HConfig.GetTH1D(Name+"_Tau1Phi","Phi of selected #tau1 candidate",10,-3.5,3.5," #phi(#tau1)","Events");
  Tau1Eta=HConfig.GetTH1D(Name+"_Tau1Eta","Pseudorapidity tau1",15,-2.7,2.7," #eta(#tau1)","Events");
  Tau1dz=HConfig.GetTH1D(Name+"_Tau1dz","Tau1 dz",20,-0.04,0.04,"Tau1 dz","Events");
  Tau1HPSDecayMode=HConfig.GetTH1D(Name+"_Tau1HPSDecayMode","Decay mode of the selected #tau candidate",11,-0.5,10.5," HPS Mode ","Events");

  Tau2PT=HConfig.GetTH1D(Name+"_Tau2PT","Transverse momentum of selected #tau2 candidate",20,30,55," P_{T}(#tau2), GeV","Events");
  Tau2E=HConfig.GetTH1D(Name+"_Tau2E","Energy of selected #tau2 candidate",20,30.,140," E(#tau2), GeV","Events");
  Tau2Mass=HConfig.GetTH1D(Name+"_Tau2Mass","Mass of selected #tau2 candidate",20,0,2.," M(#tau2), GeV","Events");
  Tau2Phi=HConfig.GetTH1D(Name+"_Tau2Phi","Phi of selected #tau2 candidate",10,-3.5,3.5," #phi(#tau2)","Events");
  Tau2Eta=HConfig.GetTH1D(Name+"_Tau2Eta","Pseudorapidity Tau2",15,-2.7,2.7," #eta(#tau2)","Events");
  Tau2dz=HConfig.GetTH1D(Name+"_Tau2dz","Tau2dz",20,-0.04,0.04,"Tau2 dz","Events");
  Tau2HPSDecayMode=HConfig.GetTH1D(Name+"_Tau2HPSDecayMode","Decay mode of the selected #tau candidate",11,-0.5,10.5," HPS Mode ","Events");

  h_SVFitMass = HConfig.GetTH1D(Name+"_SVFitMass","SVFitMass",100,0.,200.,"m_{SVfit}(#tau_{h},#tau_{h})/GeV");
  h_SVFitStatus = HConfig.GetTH1D(Name+"_SVFitStatus", "SVFitStatus", 5, -0.5, 4.5, "Status of SVFit calculation");
  svfTau1E = HConfig.GetTH1D(Name+"_svfTau1E","svFitTau1E",40,20.,120.,"E_{SVfit}(#tau_{h}1)/GeV");
  svfTau2E = HConfig.GetTH1D(Name+"_svfTau2E","svFitTau2E",40,20.,120.,"E_{SVfit}(#tau_{h}2)/GeV");

  
  TauTauFullPtRes=HConfig.GetTH1D(Name+"_TauTauFullPtRes","Resolution of the full TauTau Pt",100,-50.,50.,"Full TauTau Pt resol, GeV","Events");
  TauTauFullEtaRes=HConfig.GetTH1D(Name+"_TauTauFullEtaRes","Resolution of the full TauTau Eta",100,-5.,5.,"Full TauTau Eta resol, GeV","Events");
  TauTauFullPhiRes=HConfig.GetTH1D(Name+"_TauTauFullPhiRes","Resolution of the full TauTau Phi",100,-3.,3.,"Full TauTau Phi resol, GeV","Events");

  TauplusFullPtRes=HConfig.GetTH1D(Name+"_TauplusFullPtRes","Resolution of the full Tau+ Pt",100,-50.,50.,"Full Tau+ Pt resol, GeV","Events");
  TauplusFullEtaRes=HConfig.GetTH1D(Name+"_TauplusFullEtaRes","Resolution of the full Tau+ Eta",100,-5.,5.,"Full Tau+ Eta resol, GeV","Events");
  TauplusFullPhiRes=HConfig.GetTH1D(Name+"_TauplusFullPhiRes","Resolution of the full Tau+ Phi",100,-3.,3.,"Full Tau+ Phi resol, GeV","Events");

  TauminusFullPtRes=HConfig.GetTH1D(Name+"_TauminusFullPtRes","Resolution of the full Tau- Pt",100,-50.,50.,"Full Tau- Pt resol, GeV","Events");
  TauminusFullEtaRes=HConfig.GetTH1D(Name+"_TauminusFullEtaRes","Resolution of the full Tau- Eta",50,-2.,2.,"Full Tau- Eta resol, GeV","Events");
  TauminusFullPhiRes=HConfig.GetTH1D(Name+"_TauminusFullPhiRes","Resolution of the full Tau- Phi",100,-3.,3.,"Full Tau- Phi resol, GeV","Events");


  TauTauVisPtRes=HConfig.GetTH1D(Name+"_TauTauVisPtRes","Resolution of the visible TauTau Pt",100,-50.,50.,"Visible TauTau Pt resol, GeV","Events");
  TauTauVisEtaRes=HConfig.GetTH1D(Name+"_TauTauVisEtaRes","Resolution of the visible TauTau Eta",100,-5.,5.,"Visible TauTau Eta resol, GeV","Events");
  TauTauVisPhiRes=HConfig.GetTH1D(Name+"_TauTauVisPhiRes","Resolution of the visible TauTau Phi",100,-3.,3.,"Visible TauTau Phi resol, GeV","Events");

  TauplusVisPtRes=HConfig.GetTH1D(Name+"_TauplusVisPtRes","Resolution of the visible Tau+ Pt",100,-50.,50.,"Visible Tau+ Pt resol, GeV","Events");
  TauplusVisEtaRes=HConfig.GetTH1D(Name+"_TauplusVisEtaRes","Resolution of the visible Tau+ Eta",100,-5.,5.,"Visible Tau+ Eta resol, GeV","Events");
  TauplusVisPhiRes=HConfig.GetTH1D(Name+"_TauplusVisPhiRes","Resolution of the visible Tau+ Phi",100,-3.,3.,"Visible Tau+ Phi resol, GeV","Events");
    
  TauminusVisPtRes=HConfig.GetTH1D(Name+"_TauminusVisPtRes","Resolution of the visible Tau- Pt",100,-50.,50.,"Visible Tau- Pt resol, GeV","Events");
  TauminusVisEtaRes=HConfig.GetTH1D(Name+"_TauminusVisEtaRes","Resolution of the visible Tau- Eta",50,-2.,5.,"Visible Tau- Eta resol, GeV","Events");
  TauminusVisPhiRes=HConfig.GetTH1D(Name+"_TauminusVisPhiRes","Resolution of the visible Tau- Phi",100,-3.,3.,"Visible Tau- Phi resol, GeV","Events");



  TauTauFullPtResPull=HConfig.GetTH1D(Name+"_TauTauFullPtResPull","Pull Plot of the full TauTau Pt",50,-10.,10.,"Full TauTau Pt pull plot, GeV","Events");
  TauTauFullEtaResPull=HConfig.GetTH1D(Name+"_TauTauFullEtaResPull","Pull Plot of the full TauTau Eta",100,-5.,5.,"Full TauTau Eta pull plot, GeV","Events");
  TauTauFullPhiResPull=HConfig.GetTH1D(Name+"_TauTauFullPhiResPull","Pull Plot of the full TauTau Phi",100,-3.,3.,"Full TauTau Phi pull plot, GeV","Events");

  TauplusFullPtResPull=HConfig.GetTH1D(Name+"_TauplusFullPtResPull","Pull Plot of the full Tau+ Pt",50,-10.,10.,"Full Tau+ Pt pull plot, GeV","Events");
  TauplusFullEtaResPull=HConfig.GetTH1D(Name+"_TauplusFullEtaResPull","Pull Plot of the full Tau+ Eta",100,-5.,5.,"Full Tau+ Eta pull plot, GeV","Events");
  TauplusFullPhiResPull=HConfig.GetTH1D(Name+"_TauplusFullPhiResPull","Pull Plot of the full Tau+ Phi",100,-3.,3.,"Full Tau+ Phi pull plot, GeV","Events");

  TauminusFullPtResPull=HConfig.GetTH1D(Name+"_TauminusFullPtResPull","Pull Plot of the full Tau- Pt",50,-2.,2.,"Full Tau- Pt pull plot, GeV","Events");
  TauminusFullEtaResPull=HConfig.GetTH1D(Name+"_TauminusFullEtaResPull","Pull Plot of the full Tau- Eta",50,-2.,2.,"Full Tau- Eta pull plot, GeV","Events");
  TauminusFullPhiResPull=HConfig.GetTH1D(Name+"_TauminusFullPhiResPull","Pull Plot of the full Tau- Phi",100,-3.,3.,"Full Tau- Phi pull plot, GeV","Events");


  TauTauVisPtResPull=HConfig.GetTH1D(Name+"_TauTauVisPtResPull","Pull Plot of the visible TauTau Pt",50,-10.,10.,"Visible TauTau Pt pull plot, GeV","Events");
  TauTauVisEtaResPull=HConfig.GetTH1D(Name+"_TauTauVisEtaResPull","Pull Plot of the visible TauTau Eta",100,-5.,5.,"Visible TauTau Eta pull plot, GeV","Events");
  TauTauVisPhiResPull=HConfig.GetTH1D(Name+"_TauTauVisPhiResPull","Pull Plot of the visible TauTau Phi",100,-3.,3.,"Visible TauTau Phi pull plot, GeV","Events");

  TauplusVisPtResPull=HConfig.GetTH1D(Name+"_TauplusVisPtResPull","Pull Plot of the visible Tau+ Pt",50,-10.,10.,"Visible Tau+ Pt pull plot, GeV","Events");
  TauplusVisEtaResPull=HConfig.GetTH1D(Name+"_TauplusVisEtaResPull","Pull Plot of the visible Tau+ Eta",100,-5.,5.,"Visible Tau+ Eta pull plot, GeV","Events");
  TauplusVisPhiResPull=HConfig.GetTH1D(Name+"_TauplusVisPhiResPull","Pull Plot of the visible Tau+ Phi",100,-3.,3.,"Visible Tau+ Phi pull plot, GeV","Events");
    
  TauminusVisPtResPull=HConfig.GetTH1D(Name+"_TauminusVisPtResPull","Pull Plot of the visible Tau- Pt",50,-10.,10.,"Visible Tau- Pt pull plot, GeV","Events");
  TauminusVisEtaResPull=HConfig.GetTH1D(Name+"_TauminusVisEtaResPull","Pull Plot of the visible Tau- Eta",50,-2.,2.,"Visible Tau- Eta pull plot, GeV","Events");
  TauminusVisPhiResPull=HConfig.GetTH1D(Name+"_TauminusVisPhiResPull","Pull Plot of the visible Tau- Phi",100,-3.,3.,"Visible Tau- Phi pull plot, GeV","Events");


  Pi0EnergyRes=HConfig.GetTH1D(Name+"_Pi0EnergyRes","Energy resolution of Pi0",100,-50.,50.,"Energy resolution of Pi0, GeV","Events");
  Pi0EnergyResPull=HConfig.GetTH1D(Name+"_Pi0EnergyResPull","Energy Pull Plot of Pi0",100,-50.,50.,"Energy Pull Plot of Pi0, GeV","Events");

  Selection::ConfigureHistograms();   //   do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);  // do not remove
}

 

void  SVFitIPHC::Store_ExtraDist(){

  //every new histo should be addedd to Extradist1d vector, just push it back;
  Extradist1d.push_back(&Tau1PT);
  Extradist1d.push_back(&Tau1E);
  Extradist1d.push_back(&Tau1Mass);
  Extradist1d.push_back(&Tau1Phi);
  Extradist1d.push_back(&Tau1Eta);
  Extradist1d.push_back(&Tau1dz);
  Extradist1d.push_back(&Tau1HPSDecayMode);
  
  Extradist1d.push_back(&Tau2PT);
  Extradist1d.push_back(&Tau2E);
  Extradist1d.push_back(&Tau2Mass);
  Extradist1d.push_back(&Tau2Phi);
  Extradist1d.push_back(&Tau2Eta);
  Extradist1d.push_back(&Tau2dz);
  Extradist1d.push_back(&Tau2HPSDecayMode);

  Extradist1d.push_back(&h_SVFitMass);
  Extradist1d.push_back(&h_SVFitStatus);
  Extradist1d.push_back(&svfTau1E);
  Extradist1d.push_back(&svfTau2E);


  Extradist1d.push_back(&TauTauFullPtRes);
  Extradist1d.push_back(&TauTauFullEtaRes);
  Extradist1d.push_back(&TauTauFullPhiRes);

  Extradist1d.push_back(&TauplusFullPtRes);
  Extradist1d.push_back(&TauplusFullEtaRes);
  Extradist1d.push_back(&TauplusFullPhiRes);

  Extradist1d.push_back(&TauminusFullPtRes);
  Extradist1d.push_back(&TauminusFullEtaRes);
  Extradist1d.push_back(&TauminusFullPhiRes);

  Extradist1d.push_back(&TauTauVisPtRes);
  Extradist1d.push_back(&TauTauVisEtaRes);
  Extradist1d.push_back(&TauTauVisPhiRes);

  Extradist1d.push_back(&TauplusVisPtRes);
  Extradist1d.push_back(&TauplusVisEtaRes);
  Extradist1d.push_back(&TauplusVisPhiRes);
    
  Extradist1d.push_back(&TauminusVisPtRes);
  Extradist1d.push_back(&TauminusVisEtaRes);
  Extradist1d.push_back(&TauminusVisPhiRes);

  Extradist1d.push_back(&TauTauFullPtResPull);
  Extradist1d.push_back(&TauTauFullEtaResPull);
  Extradist1d.push_back(&TauTauFullPhiResPull);

  Extradist1d.push_back(&TauplusFullPtResPull);
  Extradist1d.push_back(&TauplusFullEtaResPull);
  Extradist1d.push_back(&TauplusFullPhiResPull);

  Extradist1d.push_back(&TauminusFullPtResPull);
  Extradist1d.push_back(&TauminusFullEtaResPull);
  Extradist1d.push_back(&TauminusFullPhiResPull);

  Extradist1d.push_back(&TauTauVisPtResPull);
  Extradist1d.push_back(&TauTauVisEtaResPull);
  Extradist1d.push_back(&TauTauVisPhiResPull);

  Extradist1d.push_back(&TauplusVisPtResPull);
  Extradist1d.push_back(&TauplusVisEtaResPull);
  Extradist1d.push_back(&TauplusVisPhiResPull);
    
  Extradist1d.push_back(&TauminusVisPtResPull);
  Extradist1d.push_back(&TauminusVisEtaResPull);
  Extradist1d.push_back(&TauminusVisPhiResPull);

  Extradist1d.push_back(&Pi0EnergyRes);
  Extradist1d.push_back(&Pi0EnergyResPull);

}

void  SVFitIPHC::doEvent()  { //  Method called on every event
  unsigned int t;                // sample type, you may manage in your further analysis, if needed
  int id(Ntp->GetMCID());  //read event ID of a sample
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;} //  gives a warning if list of samples in Histo.txt and SkimSummary.log do not coincide 
  bool trig=0;
  std::vector<int> TauIndex ;
  std::vector<int> TriggerIndexVector ;
  std::vector<TString>  MatchedTriggerNames;
  value.at(Trigger)=0;
  MatchedTriggerNames.push_back("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v");
  MatchedTriggerNames.push_back("HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v");
  TriggerIndexVector=Ntp->GetVectorTriggers(MatchedTriggerNames);

  for(unsigned int itrig = 0; itrig < TriggerIndexVector.size(); itrig++){
    if(Ntp->TriggerAccept(TriggerIndexVector.at(itrig))){
      trig=1;
    }
  }
  for(unsigned int iDaughter=0;   iDaughter  <  Ntp->NDaughters() ;iDaughter++ ) {
    if(Ntp->CHECK_BIT(Ntp->Daughters_trgMatched(iDaughter),29) || Ntp->CHECK_BIT(Ntp->Daughters_trgMatched(iDaughter),27)) //Number in the list of triggers
      {
	TauIndex.push_back(iDaughter);
     }
  }
  if(trig && TauIndex.size()>1)value.at(Trigger)=1;
  pass.at(Trigger)=(value.at(Trigger)==cut.at(Trigger));
  value.at(Id_and_Kin)=0;
  int goodTau_counter=0;
  std::vector<int> thirdLeptonCounter;
  std::vector<int> goodTausIndex;
  for(unsigned int iDaughter=0;   iDaughter  < TauIndex.size()  ;iDaughter++ ) {
    if(Ntp->tauBaselineSelection(iDaughter,36., 2.1, 0,0)){
      goodTausIndex.push_back(TauIndex.at(iDaughter));
      goodTau_counter++;
      }
  }
  unsigned int Tau1= -1;
  unsigned int Tau2= -1;
  std::vector<int>  Sorted;
  std::vector<int>  PairsIndexTemp;
  std::vector<int>  PairsIndexTau1Temp;
  std::vector<int>  PairsIndexTau2Temp;
  int j=0;
  if(goodTau_counter>1)value.at(Id_and_Kin)=1;
  pass.at(Id_and_Kin)=(value.at(Id_and_Kin) == cut.at(Id_and_Kin));
  value.at(NPairsFound)=0;

  for(int unsigned ipair =0; ipair <goodTausIndex.size(); ipair++)
    {
      for(int unsigned jpair =1; jpair <goodTausIndex.size(); jpair++)
	{
	  if(jpair>ipair)
	    {
	      if((Ntp->Daughters_P4(goodTausIndex.at(ipair)).DeltaR(Ntp->Daughters_P4(goodTausIndex.at(jpair))))>0.5){
		PairsIndexTemp.push_back(j);
		PairsIndexTau1Temp.push_back(goodTausIndex.at(ipair));
		PairsIndexTau2Temp.push_back(goodTausIndex.at(jpair));
		j++;
	      }
	    }
	}
    }
  if(PairsIndexTemp.size()>0)value.at(NPairsFound)=1;
  pass.at(NPairsFound)=(value.at(NPairsFound)==cut.at(NPairsFound));
  if(pass.at(NPairsFound))
    {
      Sorted = Ntp->SortPair(PairsIndexTemp,PairsIndexTau1Temp,PairsIndexTau2Temp);
      Tau1=PairsIndexTau1Temp.at(Sorted.back());
      Tau2=PairsIndexTau2Temp.at(Sorted.back());
      value.at(Tau1Isolation)=0;
      value.at(Tau1Isolation) = (Ntp->isIsolatedTau(Tau1,"Tight"));
      pass.at(Tau1Isolation) = value.at(Tau1Isolation);
      value.at(Tau2Isolation)=0;
      value.at(Tau2Isolation) = (Ntp->isIsolatedTau(Tau2,"Tight"));
      pass.at(Tau2Isolation) = value.at(Tau2Isolation);

      
      value.at(LeptonVeto)=0;
      for(unsigned int iDaughter=0;   iDaughter  <  Ntp->NDaughters() ;iDaughter++ ) {
	if((iDaughter!=Tau1)&&(iDaughter!=Tau2)){
	  if(Ntp->ElectronVeto(iDaughter) || Ntp->MuonVeto(iDaughter))thirdLeptonCounter.push_back(iDaughter);
	}
      }
      value.at(LeptonVeto) = thirdLeptonCounter.size()>0;
      pass.at(LeptonVeto) = (value.at(LeptonVeto)==cut.at(LeptonVeto));

      value.at(PairCharge)=0;
      bool isOS=false;
      isOS=((Ntp->Daughters_charge(Tau1)/abs(Ntp->Daughters_charge(Tau1))) != (Ntp->Daughters_charge(Tau2)/abs(Ntp->Daughters_charge(Tau2))));
      if(isOS)value.at(PairCharge) = 1;
      pass.at(PairCharge) = value.at(PairCharge);
      value.at(PairMass) = 999.;
      value.at(PairMass)=(Ntp->Daughters_P4(Tau1)+Ntp->Daughters_P4(Tau2)).M();
      pass.at(PairMass) = (value.at(PairMass) < cut.at(PairMass));
    }
  double wobs=1;
  double w=1;

  std::vector<unsigned int> exclude_cuts;
  exclude_cuts.push_back(Tau1Isolation);
  exclude_cuts.push_back(Tau2Isolation);
  exclude_cuts.push_back(PairCharge);
  bool IsQCDEvent = false;
  if(passAllBut(exclude_cuts)){
    if(pass.at(PairCharge)){
      if((Ntp->isIsolatedTau(Tau1,"Medium") && !Ntp->isIsolatedTau(Tau2,"Tight") && Ntp->isIsolatedTau(Tau2,"Loose")) || (Ntp->isIsolatedTau(Tau2,"Medium") && !Ntp->isIsolatedTau(Tau1,"Tight") && Ntp->isIsolatedTau(Tau1,"Loose"))){
	if(id == DataMCType::Data){
	  //QCDShape.at(t).Fill(1,w);
	  t=HConfig.GetType(DataMCType::QCD);
	  IsQCDEvent = true;
	}
      }
    }
  }

  if(IsQCDEvent){ pass.at(PairCharge)= true;pass.at(Tau2Isolation)= true;pass.at(Tau1Isolation)=true;}

  bool status=AnalysisCuts(t,w,wobs);  // boolean that say whether your event passed critera defined in pass vector. The whole vector must be true for status = true  

  ///////////////////////////////////////////////////////////
  // Analyse events which passed selection
  if(status) {
    TLorentzVector Tau1P4 = Ntp->TauP4_Corrected(Tau1);
    TLorentzVector Tau2P4 = Ntp->TauP4_Corrected(Tau2);
      //---------  svfit ---------------------
    std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
    classic_svFit::MeasuredTauLepton lep1(1, Tau1P4.Pt(), Tau1P4.Eta(),  Tau1P4.Phi(), Tau1P4.M());
    classic_svFit::MeasuredTauLepton lep2(1, Tau2P4.Pt(), Tau2P4.Eta(),  Tau2P4.Phi(), Tau2P4.M());
  
    measuredTauLeptons.push_back(lep1);
    measuredTauLeptons.push_back(lep2);
    TMatrixD metcov(2,2);
    double metx = Ntp->MET()*cos(Ntp->METphi());
    double mety = Ntp->MET()*sin(Ntp->METphi());
   
    metcov[0][0] = Ntp->PFMETCov00();
    metcov[1][0] = Ntp->PFMETCov01();
    metcov[0][1] = Ntp->PFMETCov10();
    metcov[1][1] = Ntp->PFMETCov11();

    svfitAlgo1.addLogM_fixed(true,5.0);
    svfitAlgo1.setDiTauMassConstraint(-1.0);
    svfitAlgo1.integrate(measuredTauLeptons,metx,mety, metcov );

    if(svfitAlgo1.isValidSolution()){
      double higgsmass  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svfitAlgo1.getHistogramAdapter())->getMass();
      h_SVFitMass.at(t).Fill(higgsmass,w);
    }
    ClassicSVfit svfitAlgo2;
    svfitAlgo2.setHistogramAdapter(new classic_svFit::TauTauHistogramAdapter());
    svfitAlgo2.addLogM_fixed(true, 5.);
    svfitAlgo2.integrate(measuredTauLeptons,metx,mety, metcov );

    classic_svFit::LorentzVector tau1P4 = static_cast<classic_svFit::TauTauHistogramAdapter*>(svfitAlgo2.getHistogramAdapter())->GetFittedTau1LV();
    classic_svFit::LorentzVector tau2P4 = static_cast<classic_svFit::TauTauHistogramAdapter*>(svfitAlgo2.getHistogramAdapter())->GetFittedTau2LV();


    svfTau1E.at(t).Fill(tau1P4.E(),w);
    svfTau2E.at(t).Fill(tau2P4.E(),w);
    

    Tau1PT.at(t).Fill(Tau1P4.Pt(),w);
    Tau1E.at(t).Fill(Tau1P4.E(),w);
    Tau1Mass.at(t).Fill(Tau1P4.M(),w);
    Tau1Phi.at(t).Fill(Tau1P4.Phi(),w);
    Tau1Eta.at(t).Fill(Tau1P4.Eta(),w);
    Tau1dz.at(t).Fill(Ntp->dz(Tau1),w);
    Tau1HPSDecayMode.at(t).Fill(Ntp->decayMode(Tau1),w);

    Tau2PT.at(t).Fill(Tau2P4.Pt(),w);
    Tau2E.at(t).Fill(Tau2P4.E(),w);
    Tau2Mass.at(t).Fill(Tau2P4.M(),w);
    Tau2Phi.at(t).Fill(Tau2P4.Phi(),w);
    Tau2Eta.at(t).Fill(Tau2P4.Eta(),w);
    Tau2dz.at(t).Fill(Ntp->dz(Tau2),w);
    Tau2HPSDecayMode.at(t).Fill(Ntp->decayMode(Tau2),w);

    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> Tauplussvfit;
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> Tauminussvfit;
    TLorentzVector Tauplusvis;
    TLorentzVector Tauminusvis;
    TLorentzVector Tauplustruth;
    TLorentzVector Tauminustruth;

    if(Ntp->Daughters_charge(Tau1)>0)
      {
	Tauplussvfit=tau1P4;
	Tauminussvfit=tau2P4;
	Tauplusvis=Tau1P4;
	Tauminusvis=Tau2P4;
      }
    else
      {
	Tauplussvfit=tau2P4;
	Tauminussvfit=tau1P4;
	Tauplusvis=Tau2P4;
	Tauminusvis=Tau1P4;
      }
    if(id == 10110333 || id == 10110433|| id == 10130533|| id ==10210333|| id == 10210433|| id == 10230533|| id ==10310333 || id ==10330533 || id ==10410433 || id == 10410333|| id == 10430533|| id == 30530533)
      {
	TLorentzVector Tau1Truth; 
	TLorentzVector Tau2Truth;
	TLorentzVector TruthDecayFromTau1;
	TLorentzVector TruthDecayFromTau2; 
	std::vector<TLorentzVector> Pions1;
	std::vector<TLorentzVector> Pions2;
	bool decay=0;

	if(Ntp->CheckDecayID(1,3)){
	  Tau1Truth=Ntp->GetTruthTauLV(1,0);
	  Tau2Truth=Ntp->GetTruthTauLV(3,1);
	  TruthDecayFromTau1=Ntp->GetTruthTauProductLV(1,11,0);
	  TruthDecayFromTau2=Ntp->GetTruthTauProductLV(3,211,1);decay=1;
	}
	if(Ntp->CheckDecayID(1,4)){
	  Tau1Truth=Ntp->GetTruthTauLV(1,0);
	  Tau2Truth=Ntp->GetTruthTauLV(4,1);
	  TruthDecayFromTau1=Ntp->GetTruthTauProductLV(1,11,0);
	  Pions2=Ntp->GetTruthPionsFromRho(1);
	  TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1);decay=1;
	  if((sqrt((Tau1P4.Eta()-Tau1Truth.Eta())*(Tau1P4.Eta()-Tau1Truth.Eta())+(Tau1P4.Phi()-Tau1Truth.Phi())*(Tau1P4.Phi()-Tau1Truth.Phi())))<0.5 && (sqrt((Tau2P4.Eta()-Tau2Truth.Eta())*(Tau2P4.Eta()-Tau2Truth.Eta())+(Tau2P4.Phi()-Tau2Truth.Phi())*(Tau2P4.Phi()-Tau2Truth.Phi())))<0.5)
	    {
	      Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau2).E()-Pions2.at(1).E(),w);
	      Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau2).E()-Pions2.at(1).E())/Pions2.at(1).E(),w);
	    }
	  else{ Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau1).E()-Pions2.at(1).E(),w);Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau1).E()-Pions2.at(1).E())/Pions2.at(1).E(),w);}
	}
	if(Ntp->CheckDecayID(1,5)){
	  Tau1Truth=Ntp->GetTruthTauLV(1,0);
	  Tau2Truth=Ntp->GetTruthTauLV(5,1);
	  TruthDecayFromTau1=Ntp->GetTruthTauProductLV(1,11,0);
	  Pions2=Ntp->GetTruthPionsFromA1(1);
	  TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1)+Pions2.at(2);decay=1;

	}

	if(Ntp->CheckDecayID(2,3)){
	  Tau1Truth=Ntp->GetTruthTauLV(2,0);
	  Tau2Truth=Ntp->GetTruthTauLV(3,1);
	  TruthDecayFromTau1=Ntp->GetTruthTauProductLV(2,13,0);
	  TruthDecayFromTau2=Ntp->GetTruthTauProductLV(3,211,1);decay=1;
	}
	if(Ntp->CheckDecayID(2,4)){
	  Tau1Truth=Ntp->GetTruthTauLV(2,0);
	  Tau2Truth=Ntp->GetTruthTauLV(4,1);
	  TruthDecayFromTau1=Ntp->GetTruthTauProductLV(2,13,0);
	  Pions2=Ntp->GetTruthPionsFromRho(1);
	  TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1);decay=1;
	  if((sqrt((Tau1P4.Eta()-Tau1Truth.Eta())*(Tau1P4.Eta()-Tau1Truth.Eta())+(Tau1P4.Phi()-Tau1Truth.Phi())*(Tau1P4.Phi()-Tau1Truth.Phi())))<0.5 && (sqrt((Tau2P4.Eta()-Tau2Truth.Eta())*(Tau2P4.Eta()-Tau2Truth.Eta())+(Tau2P4.Phi()-Tau2Truth.Phi())*(Tau2P4.Phi()-Tau2Truth.Phi())))<0.5)
	    {
	      Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau2).E()-Pions2.at(1).E(),w);
	      Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau2).E()-Pions2.at(1).E())/Pions2.at(1).E(),w);
	    }
	  else{ Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau1).E()-Pions2.at(1).E(),w);Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau1).E()-Pions2.at(1).E())/Pions2.at(1).E(),w);}

	}
	if(Ntp->CheckDecayID(2,5)){
	  Tau1Truth=Ntp->GetTruthTauLV(2,0);
	  Tau2Truth=Ntp->GetTruthTauLV(5,1);
	  TruthDecayFromTau1=Ntp->GetTruthTauProductLV(2,13,0);
	  Pions2=Ntp->GetTruthPionsFromA1(1);
	  TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1)+Pions2.at(2);decay=1;

	}

	if(Ntp->CheckDecayID(3,3)){
	  Tau1Truth=Ntp->GetTruthTauLV(3,0);
	  Tau2Truth=Ntp->GetTruthTauLV(3,1);
	  TruthDecayFromTau1=Ntp->GetTruthTauProductLV(3,211,0);
	  TruthDecayFromTau2=Ntp->GetTruthTauProductLV(3,211,1);decay=1;

	}
	if(Ntp->CheckDecayID(3,5)){
	  Tau1Truth=Ntp->GetTruthTauLV(3,0);
	  Tau2Truth=Ntp->GetTruthTauLV(5,1);
	  TruthDecayFromTau1=Ntp->GetTruthTauProductLV(3,211,0);
	  Pions2=Ntp->GetTruthPionsFromA1(1);
	  TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1)+Pions2.at(2);decay=1;
	}


	if(Ntp->CheckDecayID(4,4)){
	  Tau1Truth=Ntp->GetTruthTauLV(4,0);
	  Tau2Truth=Ntp->GetTruthTauLV(4,1);
	  Pions1=Ntp->GetTruthPionsFromRho(0);
	  TruthDecayFromTau1=Pions1.at(0)+Pions1.at(1);
	  Pions2=Ntp->GetTruthPionsFromRho(1);
	  TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1);decay=1;
	  if((sqrt((Tau1P4.Eta()-Pions1.at(1).Eta())*(Tau1P4.Eta()-Pions1.at(1).Eta())+(Tau1P4.Phi()-Pions1.at(1).Phi())*(Tau1P4.Phi()-Pions1.at(1).Phi())))<0.5 && (sqrt((Tau2P4.Eta()-Pions2.at(1).Eta())*(Tau2P4.Eta()-Pions2.at(1).Eta())+(Tau2P4.Phi()-Pions2.at(1).Phi())*(Tau2P4.Phi()-Pions2.at(1).Phi())))<0.5)
	    {
	      Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau1).E()-Pions1.at(1).E(),w);
	      Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau2).E()-Pions2.at(1).E(),w);
	      Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau1).E()-Pions1.at(1).E())/Pions1.at(1).E(),w);
	      Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau2).E()-Pions2.at(1).E())/Pions2.at(1).E(),w);
	    }
	  else
	    {
	      Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau1).E()-Pions2.at(1).E(),w);
	      Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau2).E()-Pions1.at(1).E(),w);
	      Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau1).E()-Pions2.at(1).E())/Pions2.at(1).E(),w);
	      Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau2).E()-Pions1.at(1).E())/Pions1.at(1).E(),w);
	    }

	}
	if(Ntp->CheckDecayID(4,3)){
	  Tau1Truth=Ntp->GetTruthTauLV(4,0);
	  Tau2Truth=Ntp->GetTruthTauLV(3,1);
	  Pions1=Ntp->GetTruthPionsFromRho(0);
	  TruthDecayFromTau1=Pions1.at(0)+Pions1.at(1);
	  TruthDecayFromTau2=Ntp->GetTruthTauProductLV(3,211,1);decay=1;
	  if((sqrt((Tau1P4.Eta()-Tau1Truth.Eta())*(Tau1P4.Eta()-Tau1Truth.Eta())+(Tau1P4.Phi()-Tau1Truth.Phi())*(Tau1P4.Phi()-Tau1Truth.Phi())))<0.5 && (sqrt((Tau2P4.Eta()-Tau2Truth.Eta())*(Tau2P4.Eta()-Tau2Truth.Eta())+(Tau2P4.Phi()-Tau2Truth.Phi())*(Tau2P4.Phi()-Tau2Truth.Phi())))<0.5)
	    {
	      Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau1).E()-Pions1.at(1).E(),w);
	      Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau1).E()-Pions1.at(1).E())/Pions1.at(1).E(),w);
	    }
	  else{ Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau2).E()-Pions1.at(1).E(),w);Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau2).E()-Pions1.at(1).E())/Pions1.at(1).E(),w);}

	}
	if(Ntp->CheckDecayID(4,5)){
	  Tau1Truth=Ntp->GetTruthTauLV(4,0);
	  Tau2Truth=Ntp->GetTruthTauLV(5,1);
	  Pions1=Ntp->GetTruthPionsFromRho(0);
	  TruthDecayFromTau1=Pions1.at(0)+Pions1.at(1);
	  Pions2=Ntp->GetTruthPionsFromA1(1);
	  TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1)+Pions2.at(2);decay=1;
	  if((sqrt((Tau1P4.Eta()-Tau1Truth.Eta())*(Tau1P4.Eta()-Tau1Truth.Eta())+(Tau1P4.Phi()-Tau1Truth.Phi())*(Tau1P4.Phi()-Tau1Truth.Phi())))<0.5 && (sqrt((Tau2P4.Eta()-Tau2Truth.Eta())*(Tau2P4.Eta()-Tau2Truth.Eta())+(Tau2P4.Phi()-Tau2Truth.Phi())*(Tau2P4.Phi()-Tau2Truth.Phi())))<0.5)
	    {
	      Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau1).E()-Pions1.at(1).E(),w);
	      Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau1).E()-Pions1.at(1).E())/Pions1.at(1).E(),w);
	    }
	  else{ Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau2).E()-Pions1.at(1).E(),w);Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau2).E()-Pions1.at(1).E())/Pions1.at(1).E(),w);}

	}


	if(Ntp->CheckDecayID(5,5)){
	  Tau1Truth=Ntp->GetTruthTauLV(5,0);
	  Tau2Truth=Ntp->GetTruthTauLV(5,1);
	  Pions1=Ntp->GetTruthPionsFromA1(0);
	  TruthDecayFromTau1=Pions1.at(0)+Pions1.at(1)+Pions1.at(2);
	  Pions2=Ntp->GetTruthPionsFromA1(1);
	  TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1)+Pions2.at(2);decay=1;
	}
	if (decay==1)
	  {
	    TLorentzVector TruthZ    = Tau1Truth+Tau2Truth;
	    //double visiblePtTruth = (TruthDecayFromTau1 + TruthDecayFromTau2).Pt();
	    TLorentzVector VisDecayfromTauplus;
	    TLorentzVector VisDecayfromTauminus;
	  
	    // Condition on DeltaR but DeltaR method doesn't exist for a TLorentzVector of type PxPyPzE4D
	    if((sqrt((Tauplussvfit.Eta()-Tau1Truth.Eta())*(Tauplussvfit.Eta()-Tau1Truth.Eta())+(Tauplussvfit.Phi()-Tau1Truth.Phi())*(Tauplussvfit.Phi()-Tau1Truth.Phi())))<0.5 && (sqrt((Tauminussvfit.Eta()-Tau2Truth.Eta())*(Tauminussvfit.Eta()-Tau2Truth.Eta())+(Tauminussvfit.Phi()-Tau2Truth.Phi())*(Tauminussvfit.Phi()-Tau2Truth.Phi())))<0.5)
	      {
		Tauplustruth=Tau1Truth;
		Tauminustruth=Tau2Truth;
		VisDecayfromTauplus=TruthDecayFromTau1;
		VisDecayfromTauminus=TruthDecayFromTau2;
	      }
	    else
	      {
		Tauplustruth=Tau2Truth;
		Tauminustruth=Tau1Truth;
		VisDecayfromTauplus=TruthDecayFromTau2;
		VisDecayfromTauminus=TruthDecayFromTau1;
	      }
	    

	    TauTauFullPtRes.at(t).Fill((Tauplussvfit+Tauminussvfit).Pt()-(Tauplustruth+Tauminustruth).Pt(),w);   
	    TauTauFullEtaRes.at(t).Fill((Tauplussvfit+Tauminussvfit).Eta()-(Tauplustruth+Tauminustruth).Eta(),w);
	    TauTauFullPhiRes.at(t).Fill(Ntp->DeltaPhi((Tauplussvfit+Tauminussvfit).Phi(),(Tauplustruth+Tauminustruth).Phi()),w);

	    TauTauVisPtRes.at(t).Fill((Tauplusvis+Tauminusvis).Pt()-(TruthDecayFromTau1 + TruthDecayFromTau2).Pt(),w);   
	    TauTauVisEtaRes.at(t).Fill((Tauplusvis+Tauminusvis).Eta()-(TruthDecayFromTau1 + TruthDecayFromTau2).Eta(),w);
	    TauTauVisPhiRes.at(t).Fill(Ntp->DeltaPhi((Tauplusvis+Tauminusvis).Phi(),(TruthDecayFromTau1 + TruthDecayFromTau2).Phi()),w);


	    TauTauFullPtResPull.at(t).Fill(((Tauplussvfit+Tauminussvfit).Pt()-(Tauplustruth+Tauminustruth).Pt())/(Tauplustruth+Tauminustruth).Pt(),w);   
	    TauTauFullEtaResPull.at(t).Fill(((Tauplussvfit+Tauminussvfit).Eta()-(Tauplustruth+Tauminustruth).Eta())/(Tauplustruth+Tauminustruth).Eta(),w);
	    TauTauFullPhiResPull.at(t).Fill(Ntp->DeltaPhi((Tauplussvfit+Tauminussvfit).Phi(),(Tauplustruth+Tauminustruth).Phi())/(Tauplustruth+Tauminustruth).Phi(),w);
	    
	    TauTauVisPtResPull.at(t).Fill(((Tauplusvis+Tauminusvis).Pt()-(TruthDecayFromTau1+TruthDecayFromTau2).Pt())/(TruthDecayFromTau1+TruthDecayFromTau2).Pt(),w);   
	    TauTauVisEtaResPull.at(t).Fill(((Tauplusvis+Tauminusvis).Eta()-(TruthDecayFromTau1+TruthDecayFromTau2).Eta())/(TruthDecayFromTau1+TruthDecayFromTau2).Eta(),w);
	    TauTauVisPhiResPull.at(t).Fill(Ntp->DeltaPhi((Tauplusvis+Tauminusvis).Phi(),(TruthDecayFromTau1+TruthDecayFromTau2).Phi())/(TruthDecayFromTau1+TruthDecayFromTau2).Phi(),w);


	    TauplusFullPtRes.at(t).Fill(Tauplussvfit.Pt()-Tauplustruth.Pt(),w);
	    TauplusFullEtaRes.at(t).Fill(Tauplussvfit.Eta()-Tauplustruth.Eta(),w);
	    TauplusFullPhiRes.at(t).Fill(Ntp->DeltaPhi(Tauplussvfit.Phi(),Tauplustruth.Phi()),w);
	    
	    TauminusFullPtRes.at(t).Fill(Tauminussvfit.Pt()-Tauminustruth.Pt(),w);
	    TauminusFullEtaRes.at(t).Fill(Tauminussvfit.Eta()-Tauminustruth.Eta(),w);
	    TauminusFullPhiRes.at(t).Fill(Ntp->DeltaPhi(Tauminussvfit.Phi(),Tauminustruth.Phi()),w);
	    
	    TauplusFullPtResPull.at(t).Fill((Tauplussvfit.Pt()-Tauplustruth.Pt())/Tauplustruth.Pt(),w);
	    TauplusFullEtaResPull.at(t).Fill((Tauplussvfit.Eta()-Tauplustruth.Eta())/Tauplustruth.Eta(),w);
	    TauplusFullPhiResPull.at(t).Fill(Ntp->DeltaPhi(Tauplussvfit.Phi(),Tauplustruth.Phi())/Tauplustruth.Phi(),w);

	    TauminusFullPtResPull.at(t).Fill((Tauminussvfit.Pt()-Tauminustruth.Pt())/Tauminustruth.Pt(),w);
	    TauminusFullEtaResPull.at(t).Fill((Tauminussvfit.Eta()-Tauminustruth.Eta())/Tauminustruth.Eta(),w);
	    TauminusFullPhiResPull.at(t).Fill(Ntp->DeltaPhi(Tauminussvfit.Phi(),Tauminustruth.Phi())/Tauminustruth.Phi(),w);
	   

	    TauplusVisPtRes.at(t).Fill(Tauplusvis.Pt()-VisDecayfromTauplus.Pt(),w);
	    TauplusVisEtaRes.at(t).Fill(Tauplusvis.Eta()-VisDecayfromTauplus.Eta(),w);
	    TauplusVisPhiRes.at(t).Fill(Ntp->DeltaPhi(Tauplusvis.Phi(),VisDecayfromTauplus.Phi()),w);

	    TauminusVisPtRes.at(t).Fill(Tauminusvis.Pt()-VisDecayfromTauminus.Pt(),w);
	    TauminusVisEtaRes.at(t).Fill(Tauminusvis.Eta()-VisDecayfromTauminus.Eta(),w);
	    TauminusVisPhiRes.at(t).Fill(Ntp->DeltaPhi(Tauminusvis.Phi(),VisDecayfromTauminus.Phi()),w);

	    TauplusVisPtResPull.at(t).Fill((Tauplusvis.Pt()-VisDecayfromTauplus.Pt())/VisDecayfromTauplus.Pt(),w);
	    TauplusVisEtaResPull.at(t).Fill((Tauplusvis.Eta()-VisDecayfromTauplus.Eta())/VisDecayfromTauplus.Eta(),w);
	    TauplusVisPhiResPull.at(t).Fill(Ntp->DeltaPhi(Tauplusvis.Phi(),VisDecayfromTauplus.Phi())/VisDecayfromTauplus.Phi(),w);

	    TauminusVisPtResPull.at(t).Fill((Tauminusvis.Pt()-VisDecayfromTauminus.Pt())/VisDecayfromTauminus.Pt(),w);
	    TauminusVisEtaResPull.at(t).Fill((Tauminusvis.Eta()-VisDecayfromTauminus.Eta())/VisDecayfromTauminus.Eta(),w);
	    TauminusVisPhiResPull.at(t).Fill(Ntp->DeltaPhi(Tauminusvis.Phi(),VisDecayfromTauminus.Phi())/VisDecayfromTauminus.Phi(),w);
	  }
      }
  }
}
//  This is a function if you want to do something after the event loop
void  SVFitIPHC::Finish() {

  Selection::Finish();
}
