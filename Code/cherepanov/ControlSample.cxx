#include "ControlSample.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "SVFitObject.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
 
ControlSample::ControlSample(TString Name_, TString id_):
  Selection(Name_,id_),
  cMu_pt(20),
  cMu_eta(2.1),
  cTau_pt(20),
  cTau_eta(2.1)
{
  ChargeSumDummy = -999;
  selMuon_IsoDummy = 999.;
}

ControlSample::~ControlSample(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  ControlSample::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)            cut.at(TriggerOk)=1;
    if(i==nGoodTaus)          cut.at(nGoodTaus)=1;
    if(i==nGoodMuons)       cut.at(nGoodMuons)=1;
    if(i==MuonIsolation)       cut.at(MuonIsolation)=0.3;
    if(i==PairCharge)          cut.at(PairCharge)=0;
    if(i==PrimeVtx)             cut.at(PrimeVtx)=1;
  }
  // Setup cut plots
  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    if(i==PrimeVtx){
      title.at(i)="Number of Prime Vertices $(N>$";
      title.at(i)+=cut.at(PrimeVtx);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of Prime Vertices";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PrimeVtx_",htitle,31,-0.5,30.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PrimeVtx_",htitle,31,-0.5,30.5,hlabel,"Events"));
    }
    else if(i==TriggerOk){
      title.at(i)="Trigger ";
      hlabel="Trigger ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==nGoodTaus){
      title.at(i)="Number of tau leptons ";
      hlabel="Number of tau leptons  ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nGoodTaus_",htitle,10,-0.5,9.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nGoodTaus_",htitle,10,-0.5,9.5,hlabel,"Events"));
    }
    else if(i==nGoodMuons){
      title.at(i)="Number of muons ";
      hlabel="Number of muons  ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nGoodMuons_",htitle,10,-0.5,9.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nGoodMuons_",htitle,10,-0.5,9.5,hlabel,"Events"));
    }
    else if(i==PairCharge){
      title.at(i)="Charge of a pair candidate";
      hlabel="Charge of a pair candidate";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PairCharge_",htitle,3,-1.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PairCharge_",htitle,3,-1.5,1.5,hlabel,"Events"));
    }
    else if(i==MuonIsolation){
      title.at(i)="Muon Isolation  ";
      hlabel="Muon Isolation  ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonIsolation_",htitle,20,0.,4.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonIsolation_",htitle,20,0.,4.5,hlabel,"Events"));
    }


  } 
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");


  TauPT=HConfig.GetTH1D(Name+"_TauPT","Transverse momentum of selected #tau candidate",50,-0.05,60.5," P_{T}(#tau), GeV","Events");
  TauE=HConfig.GetTH1D(Name+"_TauE","Energy of selected #tau candidate",50,-0.05,99.5," E(#tau), GeV","Events");
  MuonPT=HConfig.GetTH1D(Name+"_MuonPT","Transverse momentum of selected #mu candidate",50,-0.05,60.5," P_{T}(#mu), GeV","Events");
  TauHPSDecayMode=HConfig.GetTH1D(Name+"_TauHPSDecayMode","Decay mode of the selected #tau candidate",11,-0.5,10.5," HPS Mode ","Events");
  TauMuMass=HConfig.GetTH1D(Name+"_TauMuMass"," Visible mass of the selected pair candidate",50,40,100," M_{pair}, GeV","Events");
  MissingM=HConfig.GetTH1D(Name+"_MissingM","Missing transverse Mass",100,0,100," M_{T}, GeV","Events");

    Selection::ConfigureHistograms();   //   do not remove
    HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);  // do not remove
}

 

void  ControlSample::Store_ExtraDist(){

  //every new histo should be addedd to Extradist1d vector, just push it back;
  Extradist1d.push_back(&TauPT);
  Extradist1d.push_back(&TauE);
  Extradist1d.push_back(&MuonPT);
  Extradist1d.push_back(&TauMuMass);
  Extradist1d.push_back(&TauHPSDecayMode);
  Extradist1d.push_back(&MissingM);


}

void  ControlSample::doEvent(){ //  Method called on every event
  unsigned int t;                // sample type, you may manage in your further analysis, if needed
  int id(Ntp->GetMCID());  //read event ID of a sample
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}  //  gives a warining if list of samples in Histo.txt  and SkimSummary.log do not coincide 
  //  std::cout<<"------------------ New Event -----------------------"<<std::endl;
  Charge = ChargeSumDummy;

  
  //  int ntau(0); int nmu(0); 
  std::vector<int> goodMuonsIndex;
  std::vector<int> goodTauIndex;
  for(unsigned int iDaugther=0;   iDaugther  <  Ntp->NDaughters() ;iDaugther++ ){  // loop over all daughters in the event
    if(Ntp->isTightGoodTau(iDaugther)){
      if(Ntp->tauBaselineSelection(iDaugther)){
	if(Ntp->Daughters_P4(iDaugther).Pt() > cTau_pt){
	  if(fabs(  Ntp->Daughters_P4(iDaugther).Eta()) < cTau_eta  ){
	    goodTauIndex.push_back(iDaugther) ;  }}}}


    if(Ntp->isMediumGoodMuon(iDaugther)){
      if(Ntp->muonBaselineSelection(iDaugther)){
	if(Ntp->Daughters_P4(iDaugther).Pt() > cMu_pt){
	  if(fabs(  Ntp->Daughters_P4(iDaugther).Eta()) < cMu_eta  ){
	    goodMuonsIndex.push_back(iDaugther) ;  }}}}
  }
    
 



  value.at(PrimeVtx)=Ntp->NVtx();
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));
  
  value.at(TriggerOk)=1;
  pass.at(TriggerOk)=true;
  // std::cout<<"----------------   "<< std::endl;
  // std::cout<<" muon size   "<< goodMuonsIndex.size() <<std::endl;
  // for(unsigned int v = 0;v<goodMuonsIndex.size(); v++ ){
  //   std::cout<<"   muon index  "<<goodMuonsIndex.at(v)<< " pt   " <<Ntp->Daughters_P4(goodMuonsIndex.at(v)).Pt() << "  charge   "<<  Ntp->Daughters_charge(goodMuonsIndex.at(v))<< std::endl; 
  // }

  // std::cout<<" tau size   "<< goodTauIndex.size() <<std::endl;
  // for(unsigned int v = 0;v<goodTauIndex.size(); v++ ){
  //   std::cout<<"   tau index  "<<goodTauIndex.at(v)<< " pt   " <<Ntp->Daughters_P4(goodTauIndex.at(v)).Pt() << "  charge   "<<  Ntp->Daughters_charge(goodTauIndex.at(v))<< std::endl; 
  // }


  value.at(nGoodTaus)=goodTauIndex.size();
  pass.at(nGoodTaus) = (value.at(nGoodTaus) >= cut.at(nGoodTaus));


  value.at(nGoodMuons)=goodMuonsIndex.size();
  pass.at(nGoodMuons) =(value.at(nGoodMuons) >= cut.at(nGoodMuons));


 
  

  value.at(PairCharge) = ChargeSumDummy;
  value.at(MuonIsolation)=selMuon_IsoDummy;


 if(goodMuonsIndex.size()!=0){
   value.at(MuonIsolation)=Ntp->combreliso(goodMuonsIndex.at(0));
   if(goodTauIndex.size()!=0){
     Charge=Ntp->Daughters_charge(goodTauIndex.at(0)) + Ntp->Daughters_charge(goodMuonsIndex.at(0));
     value.at(PairCharge) = Charge;
   }
 }
 pass.at(PairCharge) =(value.at(PairCharge)== cut.at(PairCharge));
 pass.at(MuonIsolation) =(value.at(MuonIsolation) < cut.at(MuonIsolation));


  // Here you can defined different type of weights you want to apply to events. At the moment only PU weight is considered if event is not data
  double wobs=1;
  double w;
  if(!Ntp->isData()){w = Ntp->PUReweight();}
  else{w=1;}

  

  bool status=AnalysisCuts(t,w,wobs);  // boolean that say whether your event passed critera defined in pass vector. The whole vector must be true for status = true
  ///////////////////////////////////////////////////////////
  // Analyse events  which passed selection
  if(status){


  int MuonIndex = goodMuonsIndex.at(0);
  int TauIndex = goodTauIndex.at(0);

  TLorentzVector MuonP4 = Ntp->Daughters_P4(MuonIndex);
  TLorentzVector TauP4 = Ntp->Daughters_P4(TauIndex);


  TauPT.at(t).Fill(Ntp->Daughters_P4(TauIndex).Pt(),w);  // Fill transverse momentum
  TauE.at(t).Fill(Ntp->Daughters_P4(TauIndex).E(),w);  // Fill transverse momentum
  MuonPT.at(t).Fill(Ntp->Daughters_P4(MuonIndex).Pt(),w);  // Fill transverse momentum
  TauHPSDecayMode.at(t).Fill(Ntp->decayMode(TauIndex),w);
  TauMuMass.at(t).Fill(  (Ntp->Daughters_P4(MuonIndex) +Ntp->Daughters_P4(TauIndex) ).M() ,w);
  MissingM.at(t).Fill(Ntp->transverseMass(MuonP4.Pt(),MuonP4.Phi(),Ntp->MET(),Ntp->METphi()),w);

  
  }


}




//  This is a function if you want to do something after the event loop
void  ControlSample::Finish(){
  Selection::Finish();
}





