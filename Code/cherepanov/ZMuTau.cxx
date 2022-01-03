#include "ZMuTau.h"
#include "TLorentzVector.h"
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




ZMuTau::ZMuTau(TString Name_, TString id_):
  Selection(Name_,id_),
  cMu_pt(20),
  cMu_eta(2.1),
  cTau_pt(20),
  cTau_eta(2.1),
  DataMC_Corr(true,true,false),
  DataMC_CorrLeptonIso(true,false,true)
  // tauTrgSF("vtight")
{
  // TString basedir = "";
  // basedir = (TString)std::getenv("workdir")+"/Code/CommonFiles/";
  // myScaleFactor->init_ScaleFactor(basedir+"LeptonEfficiencies/Muon/Run2016BtoH/Muon_IsoMu24_OR_TkIsoMu24_2016BtoH_eff.root");

  ChargeSumDummy = -999;
  selMuon_IsoDummy = 999.;
}

ZMuTau::~ZMuTau(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  ZMuTau::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==Trigger)             cut.at(Trigger)=1.;
    if(i==Id_and_Kin)          cut.at(Id_and_Kin)=1.;
    if(i==NPairsFound)         cut.at(NPairsFound)=1.;
    if(i==MuonIsolation)       cut.at(MuonIsolation)=0.15;  
    if(i==TauIsolation)        cut.at(TauIsolation)=1.;
    if(i==DiMuon_Veto)         cut.at(DiMuon_Veto)=0.;
    if(i==LeptonVeto)          cut.at(LeptonVeto)=0.;
    if(i==PairCharge)          cut.at(PairCharge)=1.;
    //    if(i==PairMass)            cut.at(PairMass)=115.;
    if(i==MTM)                 cut.at(MTM)=40.;
    //if(i==deltaR)              cut.at(deltaR)=3.3;
  }
  // Setup cut plots
  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    // if(i==PrimeVtx){
    //   title.at(i)="Number of Prime Vertices $(N>$";
    //   title.at(i)+=cut.at(PrimeVtx);
    //   title.at(i)+=")";
    //   htitle=title.at(i);
    //   htitle.ReplaceAll("$","");
    //   htitle.ReplaceAll("\\","#");
    //   hlabel="Number of Prime Vertices";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PrimeVtx_",htitle,51,-0.5,50.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PrimeVtx_",htitle,51,-0.5,50.5,hlabel,"Events"));
    // }
    if(i==Trigger){
      title.at(i)="At least 1 good pair with Trig+Matching";
      hlabel="At least 1 good pair with Trig+Matching";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Trigger_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Trigger_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==Id_and_Kin){
      title.at(i)="Id and Kinematic";
      hlabel="Has good particles";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==NPairsFound){
      title.at(i)="At least 1 pair with good Delta R";
      hlabel="At least 1 pair with good Delta R";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NPairsFound_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NPairsFound_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==MuonIsolation){
      title.at(i)="Muon Isolation";
      hlabel="Muon Isolation";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonIsolation_",htitle,14,0.,0.35,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonIsolation_",htitle,14,0.,0.35,hlabel,"Events"));
    }
    else if(i==TauIsolation){
      title.at(i)="Pairs with good tau isolation";
      hlabel="Pairs with good tau isolation";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauIsolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauIsolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
     else if(i==DiMuon_Veto){
      title.at(i)="Di-Muon Veto";
      hlabel="Di-Muon";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_DiMuonVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_DiMuonVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }
    else if(i==LeptonVeto){
      title.at(i)="Lepton Veto";
      hlabel="Third Lepton";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_LeptonVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_LeptonVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }
    else if(i==PairCharge){
      title.at(i)="Pair Charge";
      hlabel="is pair OS";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PairCharge_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PairCharge_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    /* else if(i==PairMass){
      title.at(i)="Pair Visible Mass";
      hlabel="M(tau-tau)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PairMass_",htitle,15,0,220,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PairMass_",htitle,15,0,220,hlabel,"Events"));
      }*/

     else if(i==MTM){
      title.at(i)="Missing Transverse Mass";
      hlabel="Missing Transverse Mass";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MTM_",htitle,10,0,100,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MTM_",htitle,10,0,100,hlabel,"Events"));
      }
  } 
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");

  TauPT=HConfig.GetTH1D(Name+"_TauPT","Transverse momentum of selected #tau candidate",25,0,80," P_{T}(#tau), GeV","Events");
  TauE=HConfig.GetTH1D(Name+"_TauE","Energy of selected #tau candidate",20,15,150," E(#tau), GeV","Events");
  TauMass=HConfig.GetTH1D(Name+"_TauMass","Mass of selected #tau candidate",10,0,2.," M(#tau), GeV","Events");
  TauPhi=HConfig.GetTH1D(Name+"_TauPhi","Phi of selected #tau candidate",10,-3.14,3.14," #phi(#tau)","Events");
  TauEta=HConfig.GetTH1D(Name+"_TauEta","Pseudorapidity tau",20,-2.7,2.7," #eta(#tau)","Events");
  Taudz=HConfig.GetTH1D(Name+"_Taudz","Taudz",10,-0.12,0.12,"Taudz","Events");

  MuonPT=HConfig.GetTH1D(Name+"_MuonPT","Transverse momentum of selected #mu candidate",25,0,100," P_{T}(#mu), GeV","Events");
  MuonE=HConfig.GetTH1D(Name+"_MuonE","Energy of selected #mu candidate",20,20,150," E(#mu), GeV","Events");
  MuonMass=HConfig.GetTH1D(Name+"_MuonMass","Mass of selected #mu candidate",10,0.084,0.133," M(#mu), GeV","Events");
  MuonPhi=HConfig.GetTH1D(Name+"_MuonPhi","Phi of selected #mu candidate",10,-3.14,3.14," #phi(#mu)","Events");
  MuonEta=HConfig.GetTH1D(Name+"_MuonEta","Pseudorapidity muon",10,-2.5,2.5," #eta(#mu)","Events");
  Muondz=HConfig.GetTH1D(Name+"_Muondz","Muondz",10,-0.15,0.15,"Muondz","Events");
  Muondxy=HConfig.GetTH1D(Name+"_Muondxy","Muondxy",10,-0.04,0.04,"Muondxy","Events");
  MuonIsol=HConfig.GetTH1D(Name+"_MuonIsol","MuonIsol",20,0,0.15,"iso(#mu)","Events");
  

  againstElectronVLooseMVA6=HConfig.GetTH1D(Name+"_againstElectronVLooseMVA6","againstElectronVLooseMVA6",2,-0.5,1.5,"againstElectronVLooseMVA6","Events");
  againstElectronLooseMVA6=HConfig.GetTH1D(Name+"_againstElectronLooseMVA6","againstElectronLooseMVA6",2,-0.5,1.5,"againstElectronLooseMVA6","Events");
  againstElectronMediumMVA6=HConfig.GetTH1D(Name+"_againstElectronMediumMVA6","againstElectronMediumMVA6",2,-0.5,1.5,"againstElectronMediumMVA6","Events");
  againstElectronTightMVA6=HConfig.GetTH1D(Name+"_againstElectronTightMVA6","againstElectronTightMVA6",2,-0.5,1.5,"againstElectronTightMVA6","Events");
  againstElectronVTightMVA6=HConfig.GetTH1D(Name+"_againstElectronVTightMVA6","againstElectronVTightMVA6",2,-0.5,1.5,"againstElectronVTightMVA6","Events");
  againstMuonLoose3=HConfig.GetTH1D(Name+"_againstMuonLoose3","againstMuonLoose3",2,-0.5,1.5,"againstMuonLoose3","Events");
  againstMuonTight3=HConfig.GetTH1D(Name+"_againstMuonTight3","againstMuonTight3",2,-0.5,1.5,"againstMuonTight3","Events");
  byCombinedIsolationDeltaBetaCorrRaw3Hits=HConfig.GetTH1D(Name+"_byCombinedIsolationDeltaBetaCorrRaw3Hits","byCombinedIsolationDeltaBetaCorrRaw3Hits",10,0,20,"byCombinedIsolationDeltaBetaCorrRaw3Hits","Events");
  

  DiMuonVeto=HConfig.GetTH1D(Name+"_DiMuonVeto","DiMuonVeto",2,-0.5,1.5,"DiMuonVeto","Events");
  ExtraLeptonVeto=HConfig.GetTH1D(Name+"_ExtraLeptonVeto","ExtraLeptonVeto",2,-0.5,1.5,"ExtraLeptonVeto","Events");
  
  TauHPSDecayMode=HConfig.GetTH1D(Name+"_TauHPSDecayMode","Decay mode of the selected #tau candidate",11,-0.5,10.5," HPS Mode ","Events");
  QCDShape=HConfig.GetTH1D(Name+"_QCDShape","QCDShape",2,-0.5,1.5,"QCD Shape","");

  TauTauMass=HConfig.GetTH1D(Name+"_TauTauMass","Visible invariant mass of a tau pair",10,18,90," M(#tau#tau), GeV","Events");
  
  SVChi2=HConfig.GetTH1D(Name+"_SVChi2","SV  #chi^{2}",10,0,18,"#chi^{2}","Events");
  SVQuality=HConfig.GetTH1D(Name+"_SVQuality","Track mathicn #DeltaR",10,0,2,"#Sigma#Delta R","Events");
  SVQualityVsSignificance=HConfig.GetTH2D(Name+"_SVQualityVsSignificance","Track mathicn #DeltaR vs significance",25,0,3,31,-0.5,5.5,"","Events");
  NWJets=HConfig.GetTH1D(Name+"_NWJets","NWJets",4,0.5,4.5,"NWJets in ABCD","Events");
  NWJetsRelaxed=HConfig.GetTH1D(Name+"_NWJetsRelaxed","NWJetsRelaxed",2,0.5,2.5,"NWJetsRelaxed in Low and High MT","Events");
  NQCD=HConfig.GetTH1D(Name+"_NQCD","NQCD",4,0.5,4.5,"NQCD in ABCD","Events");
  PVSVSignificance=HConfig.GetTH1D(Name+"_PVSVSignificance"," PV-SV significance for tau decay mode = 10",10,-0.5,5.5," pType","Events");
  dRTauTau=HConfig.GetTH1D(Name+"_dRTauTau","#Delta R",20,0.,4.," #Delta R","Events");

  MET=HConfig.GetTH1D(Name+"_MET","MET",20,0,200,"MET, GeV","Events");
  METphi=HConfig.GetTH1D(Name+"_METphi","METphi",10,-3.14,3.14,"METphi","Events");
  PUPPImet=HConfig.GetTH1D(Name+"_PUPPImet","PUPPImet",10,0,250,"PUPPImet, GeV","Events");
  PUPPImetphi=HConfig.GetTH1D(Name+"_PUPPImetphi","PUPPImetphi",10,-3.7,3.7,"PUPPImetphi","Events");
  TransverseMass=HConfig.GetTH1D(Name+"_TransverseMass","TransverseMass, GeV",20,0,40,"TransverseMass","Events");
  
  NPrimeVtx=HConfig.GetTH1D(Name+"_NPrimeVtx","NPrimeVtx",10,0,60,"N vtx","Events");
  NPU=HConfig.GetTH1D(Name+"_npu","npu",10,0,50,"N pu","Events");
  RHO=HConfig.GetTH1D(Name+"_rho","rho",10,0,35,"rho","Events");
  
  NbJets=HConfig.GetTH1D(Name+"_NbJets","NbJets",17,0,17,"Number of jets","Events");

  Selection::ConfigureHistograms();   //   do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);  // do not remove
}

void  ZMuTau::Store_ExtraDist(){

  //every new histo should be addedd to Extradist1d vector, just push it back;
  Extradist1d.push_back(&TauPT);
  Extradist1d.push_back(&TauE);
  Extradist1d.push_back(&TauMass);
  Extradist1d.push_back(&TauPhi);
  Extradist1d.push_back(&TauEta);
  Extradist1d.push_back(&Taudz);

  Extradist1d.push_back(&MuonPT);
  Extradist1d.push_back(&MuonE);
  Extradist1d.push_back(&MuonMass);
  Extradist1d.push_back(&MuonPhi);
  Extradist1d.push_back(&MuonEta);
  Extradist1d.push_back(&Muondz);
  Extradist1d.push_back(&Muondxy);
  Extradist1d.push_back(&MuonIsol);

  Extradist1d.push_back(&againstElectronVLooseMVA6);
  Extradist1d.push_back(&againstElectronLooseMVA6);
  Extradist1d.push_back(&againstElectronMediumMVA6);
  Extradist1d.push_back(&againstElectronTightMVA6);
  Extradist1d.push_back(&againstElectronVTightMVA6);
  Extradist1d.push_back(&againstMuonLoose3);
  Extradist1d.push_back(&againstMuonTight3);
  Extradist1d.push_back(&byCombinedIsolationDeltaBetaCorrRaw3Hits);
  //Extradist1d.push_back(&byIsolationMVA3oldDMwLTraw);
  //Extradist1d.push_back(&byIsolationMVA3newDMwLTraw);

  Extradist1d.push_back(&DiMuonVeto);
  Extradist1d.push_back(&ExtraLeptonVeto);

  Extradist1d.push_back(&TauHPSDecayMode);

  Extradist1d.push_back(&dRTauTau);
  Extradist1d.push_back(&TauTauMass);
  Extradist1d.push_back(&SVChi2);
  Extradist1d.push_back(&SVQuality);
  Extradist2d.push_back(&SVQualityVsSignificance);
  Extradist1d.push_back(&QCDShape);
  Extradist1d.push_back(&NQCD);
  Extradist1d.push_back(&NWJets);
  Extradist1d.push_back(&NWJetsRelaxed);

  Extradist1d.push_back(&MET);
  Extradist1d.push_back(&METphi);
  Extradist1d.push_back(&PUPPImet);
  Extradist1d.push_back(&PUPPImetphi);
  Extradist1d.push_back(&TransverseMass);

  Extradist1d.push_back(&NPrimeVtx);
  Extradist1d.push_back(&NPU);
  Extradist1d.push_back(&RHO);

  Extradist1d.push_back(&PVSVSignificance);

  Extradist1d.push_back(&NbJets);

}

void  ZMuTau::doEvent(){ //  Method called on every event
  unsigned int t;                // sample type, you may manage in your further analysis, if needed
  int id(Ntp->GetMCID());  //read event ID of a sample
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}  //  gives a warining if list of samples in Histo.txt  and SkimSummary.log do not coincide 
  //  std::cout<<"------------------ New Event -----------------------"<<std::endl;
  Charge = ChargeSumDummy;
  bool trig=0;
  bool match=0;
  value.at(Trigger)=0;
  std::vector<int> MuonIndex;
  std::vector<int> TriggerIndexVector ;
  std::vector<TString>  MatchedTriggerNames;
  MatchedTriggerNames.push_back("HLT_IsoMu24_v");
  TriggerIndexVector=Ntp->GetVectorTriggers(MatchedTriggerNames);
  
  for(unsigned int itrig = 0; itrig < TriggerIndexVector.size(); itrig++){
    if(Ntp->TriggerAccept(TriggerIndexVector.at(itrig))){
      trig=1;
    }
  }
  
  for(unsigned int iDaughter=0;   iDaughter  <  Ntp->NDaughters() ;iDaughter++ ) {
    if(Ntp->CHECK_BIT(Ntp->Daughters_trgMatched(iDaughter),5))
      {
	MuonIndex.push_back(iDaughter);
	match=1;
     }
  }
  if(trig && match)value.at(Trigger)=1;
  pass.at(Trigger)=(value.at(Trigger)==cut.at(Trigger)); 
  bool goodTau=false;
  bool goodMuon=false;
  std::vector<int> goodMuonsIndex;
  std::vector<int> DiMuon;
  std::vector<int> goodTausIndex;
  value.at(Id_and_Kin)=0;
  for(unsigned int iDaughter=0;   iDaughter  <  Ntp->NDaughters() ;iDaughter++ ) {  // loop over all daughters in the event
      if(Ntp->tauBaselineSelection(iDaughter,30., 2.3, 0,1)){
	goodTausIndex.push_back(iDaughter);
	goodTau=true;
      } 
  }
  for(unsigned int index=0;index<MuonIndex.size();index++)
    {
      if(Ntp->isMediumGoodMuon(MuonIndex.at(index))){
	if(Ntp->muonBaselineSelection(MuonIndex.at(index), 25.,2.1,2)){
	  goodMuonsIndex.push_back(MuonIndex.at(index));
	  goodMuon=true;
	}
      }
    }
  if(goodTau && goodMuon)value.at(Id_and_Kin)=1;
  pass.at(Id_and_Kin)=(value.at(Id_and_Kin) == cut.at(Id_and_Kin));
  unsigned int Muon= 999;
  unsigned int Tau= 999;
  std::vector<int>  Sorted;
  std::vector<int>  PairsIndex;
  std::vector<int>  PairsIndexTau;
  std::vector<int>  PairsIndexMuon;
  int i=0;
  value.at(NPairsFound)=0;
  for(unsigned int ipair =0; ipair <goodTausIndex.size(); ipair++)
    {
      for(unsigned int jpair =0; jpair <goodMuonsIndex.size(); jpair++)
	{
	  if((Ntp->Daughters_P4(goodTausIndex.at(ipair)).DeltaR(Ntp->Daughters_P4(goodMuonsIndex.at(jpair))))>0.5){
	    PairsIndex.push_back(i);
	    PairsIndexTau.push_back(goodTausIndex.at(ipair));
	    PairsIndexMuon.push_back(goodMuonsIndex.at(jpair));
	    value.at(NPairsFound)=1;
	    i++;
	  }
	}
    }
  TLorentzVector MuonP4;
  TLorentzVector TauP4;
  pass.at(NPairsFound)=(value.at(NPairsFound)==cut.at(NPairsFound));
  if(pass.at(NPairsFound)) {
    Sorted = Ntp->SortPair(PairsIndex,PairsIndexTau,PairsIndexMuon);
    Muon=PairsIndexMuon.at(Sorted.back());
    Tau=PairsIndexTau.at(Sorted.back());
    MuonP4 = Ntp->Daughters_P4(Muon);
    TauP4 = Ntp->Daughters_P4(Tau);
    value.at(TauIsolation) = Ntp->isIsolatedTau(Tau,"Tight");
    value.at(MuonIsolation) = Ntp->combreliso(Muon);
    value.at(DiMuon_Veto)=(Ntp->DiMuonVeto());
    std::vector<int> thirdLeptonCounter;
    for(unsigned int iDaughter=0;   iDaughter  <  Ntp->NDaughters() ;iDaughter++ ) {  // loop over all daughters in the event
      if((iDaughter!=Tau)&&(iDaughter!=Muon)){
	if(Ntp->ElectronVeto(iDaughter) || Ntp->MuonVeto(iDaughter))thirdLeptonCounter.push_back(iDaughter);
      }
    }
    value.at(LeptonVeto) = thirdLeptonCounter.size()>0;
    value.at(PairCharge)=0;
    bool isOS=false;
    isOS=((Ntp->Daughters_charge(Muon)/abs(Ntp->Daughters_charge(Muon))) != (Ntp->Daughters_charge(Tau)/abs(Ntp->Daughters_charge(Tau))));
    if(isOS)value.at(PairCharge) = 1;

    //value.at(PairMass) = 999.;
    value.at(MTM) = 999.;
    value.at(MTM)=Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi());
    //value.at(PairMass)=(MuonP4+TauP4).M();
    //pass.at(PairMass) = (value.at(PairMass) <= cut.at(PairMass));
  } 
  pass.at(MTM) = (value.at(MTM) <= cut.at(MTM));
  pass.at(MuonIsolation) = (value.at(MuonIsolation) < cut.at(MuonIsolation));
  //pass.at(againstEleMu)=(value.at(againstEleMu)==cut.at(againstEleMu));
  pass.at(TauIsolation) = (value.at(TauIsolation) == cut.at(TauIsolation));
  pass.at(DiMuon_Veto) = ( value.at(DiMuon_Veto) ==  cut.at(DiMuon_Veto));
  pass.at(LeptonVeto) = ( value.at(LeptonVeto) ==  cut.at(LeptonVeto));
  pass.at(PairCharge) = (value.at(PairCharge) == cut.at(PairCharge));

  // Here you can defined different type of weights you want to apply to events.
  double wobs=1;
  double w=1;
  double MuonSF(1);
  double MuonIso(1);


  if(!Ntp->isData() && id!=DataMCType::QCD && id!=DataMCType::W_lnu){
    //    w *= reweight.weight(2016,26,Ntp->PUNumInteractions());
    w *= reweight.PUweightHTT(Ntp->npu());
      //std::cout<<" pu weigh HTT  "<< reweight.PUweightHTT(Ntp->npu())<<std::endl;
    
    if(!Ntp->isData() && pass.at(NPairsFound) && id==33){
      w *= 0.95;  // Tau ID  correction
       MuonSF = DataMC_Corr.get_ScaleFactor(Ntp->Daughters_P4(Muon).Pt(),Ntp->Daughters_P4(Muon).Eta());
       MuonIso =  DataMC_CorrLeptonIso.get_ScaleFactor(Ntp->Daughters_P4(Muon).Pt(),Ntp->Daughters_P4(Muon).Eta());
      //      std::cout<<"MuonSF   "<<MuonSF <<"muon iso   "<<MuonIso <<std::endl;
    }
    w*=MuonSF;
    w*=MuonIso;

    // std::cout<<"weight "<< MuonSF*MuonIso*0.95 <<std::endl;

  TLorentzVector genMomentum(0,0,0,0);
  if( id == 33){
    for(unsigned int imc=0; imc < Ntp->NGenParts(); imc++){
      // if((fabs(Ntp->Genpart_pdg(imc)) ==11 || fabs(Ntp->Genpart_pdg(imc)) ==13)   &&  Ntp->CHECK_BIT(Ntp->Genpart_flags(imc),5)  && Ntp->Genpart_status(imc) ==2){
      //   if(Ntp->Genpart_P4(imc).Pt() > 8){
      //     // std::cout<<" pdgid   "<< Ntp->Genpart_pdg(imc)<< "  px  " << Ntp->Genpart_P4(imc).Px() << " status   "<<Ntp->Genpart_status(imc) <<" index  "  << imc<<std::endl;
      //     // if(Ntp->Genpart_ZMothInd(imc)!=-1) 	std::cout<<"Mother    pdgid   "<< Ntp->Genpart_pdg(Ntp->Genpart_ZMothInd(imc))<< "  px  " << Ntp->Genpart_P4(Ntp->Genpart_ZMothInd(imc)).Px() <<std::endl;
      //   }
      // }
      if(fabs(Ntp->Genpart_pdg(imc)) ==15   &&  Ntp->CHECK_BIT(Ntp->Genpart_flags(imc),0)&& Ntp->Genpart_status(imc) ==2) {
	if(Ntp->Genpart_P4(imc).Pt() > 8){
	  genMomentum+=Ntp->Genpart_P4(imc);
	}
      }
    }
  }
  if( id == 30){
    for(unsigned int imc=0; imc < Ntp->NGenParts(); imc++){
      if((fabs(Ntp->Genpart_pdg(imc)) ==11 || fabs(Ntp->Genpart_pdg(imc)) ==13) && Ntp->Genpart_status(imc) ==1  ){
	if(Ntp->Genpart_P4(imc).Pt() > 8){
	  genMomentum+=Ntp->Genpart_P4(imc);
	}
      }
    }
  }
  float zptw(1);
  if(genMomentum.Pt()!=0 && genMomentum.M() > 75 && genMomentum.M() < 120){
    zptw = DataMC_Corr.ZPTWeight(genMomentum.M(),genMomentum.Pt());
  }
  w*=zptw;

  w*=Ntp->MC_weight(); //generator weight

  //  std::cout<<"zpt "<< zptw<<std::endl;

  }


 // QCD ABCD BG Method
  /*******************
   *        |        *
   *    C   |    D   *  MuAntiIso
   *        |        *       
   * ----------------*------ 
   *        |        *       
   *    A   |    B   *  MuIso   
   *        |        *
   *******************
   *   OS   |    SS
   *
   *
   */


// W+Jets ABCD BG Method
  /*******************
   *        |        *
   *    C   |    D   *  High MT
   *        |        *       
   * ----------------*------ 
   *        |        *       
   *    A   |    B   *  Low MT   
   *        |        *
   *******************
   *   OS   |    SS
   *
   *
   */
  std::vector<unsigned int> exclude_cuts;
  std::vector<unsigned int> relaxed_cuts;
  std::vector<unsigned int> WJetsexclude_cuts;
  exclude_cuts.push_back(PairCharge);
  exclude_cuts.push_back(MuonIsolation);
  relaxed_cuts.push_back(TauIsolation);
  relaxed_cuts.push_back(MuonIsolation);
  relaxed_cuts.push_back(MTM);
  WJetsexclude_cuts.push_back(PairCharge);
  WJetsexclude_cuts.push_back(MTM);
  
  // std::cout<<" before  " << pass.at(TriggerOk) << "    " <<   pass.at(PrimeVtx) << "    " <<  pass.at(nGoodPairs)<< "    " <<   pass.at(FirstTauIsolation) << "    " <<  pass.at(SecondTauIsolation) << "    " <<  pass.at(nGoodMuons) << "    " <<  pass.at(PairCharge) << "  passAllBut  " << passAllBut(exclude_cuts) <<std::endl;
  if(passAllBut(relaxed_cuts)){
    MuonP4 = Ntp->Daughters_P4(Muon);
    TauP4 = Ntp->Daughters_P4(Tau);
    if(Ntp->combreliso(Muon)<0.3 && Ntp->isMediumGoodTau(Tau) && Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi())<40.){
      NWJetsRelaxed.at(t).Fill(1.,w);
    }
    else if(Ntp->combreliso(Muon)<0.3 && Ntp->isMediumGoodTau(Tau) && Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi())>70.){
      NWJetsRelaxed.at(t).Fill(2.,w);
    }
  } 
  if(passAllBut(WJetsexclude_cuts)) {
    MuonP4 = Ntp->Daughters_P4(Muon);
    TauP4 = Ntp->Daughters_P4(Tau);
    if(pass.at(PairCharge)) {
      if(Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi())<40.) {
	NWJets.at(t).Fill(1.,w); //A Low
      }
      if(Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi())>70.){
	NWJets.at(t).Fill(3.,w); //A High
      }
    }
    if(!pass.at(PairCharge)){
      if(Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi())<40.){
	NWJets.at(t).Fill(2.,w); //B Low
      }
      if(Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi())>70.){
	NWJets.at(t).Fill(4.,w); //B High
      }
    }
  }
  if(passAllBut(exclude_cuts)){
    //    for(unsigned int ia=0; ia<pass.size(); ia++){         std::cout<<" ia  "<< ia <<  "   pass  " <<pass.at(ia) << std::endl;    }
    // if(pass.at(FirstTauIsolation) && pass.at(SecondTauIsolation)){
    //   if(pass.at(PairCharge)){
    // 	NQCD.at(t).Fill(1.,w); //A
    //   }  
    //   if(!pass.at(PairCharge)){
    // 	NQCD.at(t).Fill(2.,w); //B
    //   }
    //   if(Ntp->isIsolatedTau(TauIndex_1,"Medium") && Ntp->isIsolatedTau(TauIndex_2,"Loose")){
    // 	if(pass.at(PairCharge)){
    // 	  NQCD.at(t).Fill(3.,w); //С
    // 	}
    // 	if(!pass.at(PairCharge)){
    // 	  NQCD.at(t).Fill(4.,w); //В
    // 	}  
    //   }
    // }
    if(pass.at(PairCharge)){
      if(pass.at(MuonIsolation) ){
	NQCD.at(t).Fill(1.,w); //A
      }
      if(!pass.at(MuonIsolation)){
	NQCD.at(t).Fill(3.,w); //C
      }
    }
    
    if(!pass.at(PairCharge)){
      if(pass.at(MuonIsolation)){
	NQCD.at(t).Fill(2.,w); //B
      }
      if(!pass.at(MuonIsolation)){
	NQCD.at(t).Fill(4.,w); //D
      }
    }
  }

  bool IsQCDEvent = false;
  //  if(passAllBut(exclude_cuts)){
    if(!pass.at(PairCharge)){
     	if(id == DataMCType::Data){
	  QCDShape.at(t).Fill(1,w);
	  t=HConfig.GetType(DataMCType::QCD);
	  IsQCDEvent = true;
	}
    }
    //  }
  if(IsQCDEvent){    pass.at(PairCharge)= true;}
     

  bool status=AnalysisCuts(t,w,wobs);  // boolean that say whether your event passed critera defined in pass vector. The whole vector must be true for status = true
  ///////////////////////////////////////////////////////////
  // Analyse events  which passed selection
  if(status){
    Sorted = Ntp->SortPair(PairsIndex,PairsIndexTau,PairsIndexMuon);
    Muon=PairsIndexMuon.at(Sorted.back());
    Tau=PairsIndexTau.at(Sorted.back());
    MuonP4 = Ntp->Daughters_P4(Muon);
    TauP4 = Ntp->Daughters_P4(Tau);
    double pvx(0);
    pvx =  Ntp->npv();
    // if(id == DataMCType::Data) pvx =  Ntp->npv();
     if(id !=DataMCType::Data && id !=DataMCType::QCD && id !=DataMCType::W_lnu)	  pvx = Ntp->PUNumInteractions();
     if(Ntp->PFTau_hassecondaryVertex(Tau) && Ntp->isPVCovAvailable()){
       if(Ntp->PFTau_secondaryVertex_TracksMatchingQuality(Tau) < 0.01)   PVSVSignificance.at(t).Fill( Ntp->PFTau_FlightLength_significance(Ntp->PVtx(),Ntp->PFTau_TIP_primaryVertex_cov(), Ntp->PFTau_secondaryVertex_pos(Tau), Ntp->PFTau_TIP_secondaryVertex_cov(Tau)),w);
       SVChi2.at(t).Fill(Ntp->PFTau_secondaryVertex_vtxchi2(Tau),w);
       SVQuality.at(t).Fill(Ntp->PFTau_secondaryVertex_TracksMatchingQuality(Tau),w);
SVQualityVsSignificance.at(t).Fill(Ntp->PFTau_secondaryVertex_TracksMatchingQuality(Tau),Ntp->PFTau_FlightLength_significance(Ntp->PVtx(),Ntp->PFTau_TIP_primaryVertex_cov(), Ntp->PFTau_secondaryVertex_pos(Tau), Ntp->PFTau_TIP_secondaryVertex_cov(Tau)));
     }
  NPrimeVtx.at(t).Fill(pvx,w);
  NPU.at(t).Fill(Ntp->npu(),w);
  RHO.at(t).Fill(Ntp->rho(),w);
  
  std::vector<int> thirdLepton;

  TauPT.at(t).Fill(TauP4.Pt(),w);
  TauE.at(t).Fill(TauP4.E(),w);
  TauMass.at(t).Fill(TauP4.M(),w);
  TauPhi.at(t).Fill(TauP4.Phi(),w);
  TauEta.at(t).Fill(TauP4.Eta(),w);
  Taudz.at(t).Fill(Ntp->dz(Tau),w);
  

  MuonPT.at(t).Fill(MuonP4.Pt(),w);
  MuonE.at(t).Fill(MuonP4.E(),w);
  MuonMass.at(t).Fill(MuonP4.M(),w);
  MuonPhi.at(t).Fill(MuonP4.Phi(),w);
  MuonEta.at(t).Fill(MuonP4.Eta(),w);
  Muondz.at(t).Fill(Ntp->dz(Muon),w);
  Muondxy.at(t).Fill(Ntp->dxy(Muon),w);
  MuonIsol.at(t).Fill(Ntp->combreliso(Muon),w);

  
  againstElectronVLooseMVA6.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau),Ntp->Bit_againstElectronVLooseMVA6),w);
  againstElectronLooseMVA6.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau),Ntp->Bit_againstElectronLooseMVA6),w);
  againstElectronMediumMVA6.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau),Ntp->Bit_againstElectronMediumMVA6),w);
  againstElectronTightMVA6.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau),Ntp->Bit_againstElectronTightMVA6),w);
  againstElectronVTightMVA6.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau),Ntp->Bit_againstElectronVTightMVA6),w);
  againstMuonLoose3.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau),Ntp->Bit_againstMuonLoose3),w);
  againstMuonTight3.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau),Ntp->Bit_againstMuonTight3),w);
  byCombinedIsolationDeltaBetaCorrRaw3Hits.at(t).Fill(Ntp->Daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits(Tau),w);
  
  DiMuonVeto.at(t).Fill(Ntp->DiMuonVeto(),w);
  for(unsigned int iDaughter=0;   iDaughter  <  Ntp->NDaughters() ;iDaughter++ ) {
	if((iDaughter!=Tau)&&(iDaughter!=Muon)){
	  if(Ntp->ElectronVeto(iDaughter) || Ntp->MuonVeto(iDaughter))thirdLepton.push_back(iDaughter);
	}
      }
  if(thirdLepton.size()>0)ExtraLeptonVeto.at(t).Fill(1.,w);
  else ExtraLeptonVeto.at(t).Fill(0.,w);
  TauHPSDecayMode.at(t).Fill(Ntp->decayMode(Tau),w);
  
  TauTauMass.at(t).Fill((MuonP4+TauP4).M(),w);
  dRTauTau.at(t).Fill(MuonP4.DeltaR(TauP4),w);

  MET.at(t).Fill(Ntp->MET(),w);
  METphi.at(t).Fill(Ntp->METphi(),w);
  PUPPImet.at(t).Fill(Ntp->PUPPImet(),w);
  PUPPImetphi.at(t).Fill(Ntp->PUPPImetphi(),w);
  TransverseMass.at(t).Fill(Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi()),w);
  
  int jets_counter=0;
  for(int ijet=0; ijet< Ntp->JetsNumber(); ijet++) {
    if((((Ntp->Jet_P4(ijet)).Pt())>20) && (fabs((Ntp->Jet_P4(ijet).Eta())<4.7))) {
      if((((Ntp->Jet_P4(ijet)).DeltaR(Ntp->Daughters_P4(Muon)))>0.5) && (((Ntp->Jet_P4(ijet)).DeltaR(Ntp->Daughters_P4(Tau)))>0.5))jets_counter++; {
	//if(((Ntp->Jet_P4(ijet).Pt())>20) && (fabs((Ntp->Jet_P4(ijet).Eta())<2.4)) && (((Ntp->Jet_P4(ijet))->(Ntp->bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags'))) > 0.800)) {
	if((((Ntp->jets_nHEF(ijet))<0.99) && ((Ntp->jets_nEmEF(ijet)<0.99)) && ((Ntp->NumConst(ijet))>1)) && (((fabs((Ntp->Jet_P4(ijet).Eta())<=2.4)) && ((Ntp->jets_chHEF(ijet))>0) && ((Ntp->jets_chMult(ijet))>0) && ((Ntp->jets_chEmEF(ijet))<0.99)) || (fabs((Ntp->Jet_P4(ijet).Eta())>2.4))) && (fabs((Ntp->Jet_P4(ijet)).Eta()<=2.7)))jets_counter++;
	else if(((Ntp->jets_nHEF(ijet))<0.98) && ((Ntp->jets_nEmEF(ijet))>0.01) && ((Ntp->jets_neMult(ijet))>2) && (fabs((Ntp->Jet_P4(ijet).Eta())>2.7)) && (fabs((Ntp->Jet_P4(ijet).Eta())<=3.0)))jets_counter++;
	else if(((Ntp->jets_nEmEF(ijet))<0.90) && ((Ntp->jets_neMult(ijet))>10) && (fabs((Ntp->Jet_P4(ijet).Eta())>3.0)))jets_counter++;
      }
    }
  }

  NbJets.at(t).Fill(jets_counter,w);
  
  }
}




//  This is a function if you want to do something after the event loop
void  ZMuTau::Finish(){
  if(mode == RECONSTRUCT){
    std::cout<<" Starting Finish!  " <<std::endl;
    
    std::cout<<"A  Data  "<< NQCD.at(0).GetBinContent(1) << std::endl;
    std::cout<<"B  Data  "<< NQCD.at(0).GetBinContent(2) << std::endl;
    std::cout<<"C  Data  "<< NQCD.at(0).GetBinContent(3) << std::endl;
    std::cout<<"D  Data  "<< NQCD.at(0).GetBinContent(4) << std::endl;
    SkimConfig SC;
    SC.ApplySkimEfficiency(types,Npassed, Npassed_noweight);
    
    double W_ALow=0;
    std::vector<double> ALow;
    double QCD_ALow=0;
    
    std::vector<double> QCD_Integral_B;
    double QCD_IntegralMC_B=0;
    double QCD_Integral_B_Data = 0;
    std::vector<double> BLow;
    double MC_BLow=0;
    double BLow_Data = 0;
    double W_BLow=0;
    double QCD_BLow=0;
    
    std::vector<double> QCD_Integral_C;
    double QCD_IntegralMC_C=0;
    double QCD_Integral_C_Data = 0;
    std::vector<double> AHigh;
    double MC_AHigh=0;
    double AHigh_Data = 0;
    
    std::vector<double> QCD_Integral_D;
    double QCD_IntegralMC_D=0;
    double QCD_Integral_D_Data = 0;
    std::vector<double> BHigh;
    double MC_BHigh=0;
    double BHigh_Data = 0;
    
    double W_MC_ALow=0;
    double W_MC_BLow=0;
    double W_MC_AHigh=0;
    double W_MC_BHigh=0;
    double fBLW=0;
    double fALW=0;
    double W_Factor_HighLow=0;
    double QCD_Integral_C_W;
    double QCD_Integral_D_W;
    //Get Yields in ABCD for QCD Scalefactor
    for(unsigned i=0;i<CrossSectionandAcceptance.size();i++){
      QCD_Integral_B.push_back(NQCD.at(i).GetBinContent(2));
      QCD_Integral_C.push_back(NQCD.at(i).GetBinContent(3));
      QCD_Integral_D.push_back(NQCD.at(i).GetBinContent(4));
      BLow.push_back(NWJets.at(i).GetBinContent(2));
      AHigh.push_back(NWJets.at(i).GetBinContent(3));
      BHigh.push_back(NWJets.at(i).GetBinContent(4));
      if(CrossSectionandAcceptance.at(i)>0 ){
	QCD_Integral_B.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
	QCD_Integral_C.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
	QCD_Integral_D.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
	BLow.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
	AHigh.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
	BHigh.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
      }
      if(HConfig.GetID(i) == DataMCType::W_lnu)
	{
	  QCD_Integral_C_W=QCD_Integral_C.at(i);
	  QCD_Integral_D_W=QCD_Integral_D.at(i);
	  W_MC_ALow+=NWJets.at(i).GetBinContent(1);
	  W_MC_BLow+=BLow.at(i);
	  W_MC_AHigh+=AHigh.at(i);
	  W_MC_BHigh+=BHigh.at(i);
	  W_Factor_HighLow+=(NWJetsRelaxed.at(i).GetBinContent(1))/(NWJetsRelaxed.at(i).GetBinContent(2));
	  fBLW+=W_MC_BLow/W_MC_BHigh;
	  fALW+=W_MC_ALow/W_MC_AHigh;
	}
    }
    for(unsigned i=0;i<CrossSectionandAcceptance.size();i++) {
      if(HConfig.GetID(i) == DataMCType::Data){
	QCD_Integral_B_Data += QCD_Integral_B.at(i);
	QCD_Integral_C_Data += QCD_Integral_C.at(i);
	QCD_Integral_D_Data += QCD_Integral_D.at(i);
	BLow_Data += BLow.at(i);
	AHigh_Data += AHigh.at(i);
	BHigh_Data += BHigh.at(i);
      }
      if(CrossSectionandAcceptance.at(i)>0){
	QCD_IntegralMC_B  += QCD_Integral_B.at(i);
	QCD_IntegralMC_C  += QCD_Integral_C.at(i);
	QCD_IntegralMC_D  += QCD_Integral_D.at(i);
	MC_BLow  += BLow.at(i);
	MC_AHigh  += AHigh.at(i);
	MC_BHigh  += BHigh.at(i);
      }
    }
    //W_ALow=(AHigh_Data-MC_AHigh)*fALW;
    W_BLow=W_MC_BLow;//(BHigh_Data-MC_BHigh)*fBLW;
    double OS2SS=(QCD_Integral_C_Data - QCD_IntegralMC_C-QCD_Integral_C_W )/(QCD_Integral_D_Data - QCD_IntegralMC_D-QCD_Integral_D_W);
    double QCD_AHigh=OS2SS*(BHigh_Data-MC_BHigh);
    double W_AHigh=AHigh_Data-MC_AHigh-QCD_AHigh;
    double W_ScaleFactor=W_AHigh/W_MC_AHigh;
    W_ALow=W_Factor_HighLow*W_AHigh;
    QCD_BLow=BLow_Data-MC_BLow-W_BLow*W_ScaleFactor;
    QCD_ALow=QCD_BLow*OS2SS;
    //double QCD_A=(QCD_Integral_B_Data-W_IntegralMC_B+W_MC_B-W_ScaleFactor*W_MC_B)*OS2SS;

    std::cout << "QCD_Integral_B_Data is: " << QCD_Integral_B_Data << std::endl;
    std::cout << "QCD_Integral_C_Data is: " << QCD_Integral_C_Data << std::endl;
    std::cout << "QCD_Integral_D_Data is: " << QCD_Integral_D_Data << std::endl;
    std::cout << "QCD_IntegralMC_B is: " << QCD_IntegralMC_B << std::endl;
    std::cout << "QCD_IntegralMC_C is: " << QCD_IntegralMC_C << std::endl;
    std::cout << "QCD_IntegralMC_D is: " << QCD_IntegralMC_D << std::endl;
    std::cout << "QCD_AHigh: " << QCD_AHigh << std::endl;
    std::cout << "QCD_ALow: " << QCD_ALow << std::endl;
    std::cout << "QCD_BLow: " << QCD_BLow << std::endl;
    std::cout << "OS/SS QCD Sample: " << OS2SS << std::endl;
    std::cout << "BLow_Data: " <<BLow_Data << std::endl;
    std::cout << "AHigh_Data: " << AHigh_Data << std::endl;
    std::cout << "BHigh_Data: " << BHigh_Data << std::endl;
    std::cout << "MC_BLow: " << MC_BLow << std::endl;
    std::cout << "MC_AHigh: " << MC_AHigh << std::endl;
    std::cout << "MC_BHigh: " << MC_BHigh << std::endl;
    std::cout << "W_ALow: " << W_ALow << std::endl;
    std::cout << "W_BLow: " << W_BLow << std::endl;
    std::cout << "W_AHigh: " << W_AHigh << std::endl;
    std::cout << "W_MC_ALow: " <<W_MC_ALow << std::endl;
    std::cout << "W_MC_AHigh: " <<W_MC_AHigh << std::endl;
    std::cout << "W_MC_BLow: " <<W_MC_BLow << std::endl;
    std::cout << "W_MC_BHigh: " <<W_MC_BHigh << std::endl;
    std::cout << "W_ScaleFactor: " << W_ScaleFactor << std::endl;
    
    
    ScaleAllHistOfType(HConfig.GetType(DataMCType::QCD),QCD_ALow/Nminus0.at(0).at(HConfig.GetType(DataMCType::QCD)).Integral());
    ScaleAllHistOfType(HConfig.GetType(DataMCType::W_lnu),W_ALow/Nminus0.at(0).at(HConfig.GetType(DataMCType::W_lnu)).Integral());
  }

  Selection::Finish();
}
