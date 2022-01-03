#include "ZTauMuTauH.h"
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




ZTauMuTauH::ZTauMuTauH(TString Name_, TString id_):
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

ZTauMuTauH::~ZTauMuTauH(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  ZTauMuTauH::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)           cut.at(TriggerOk)=1;
    if(i==nGoodPairs)          cut.at(nGoodPairs)=0;
    if(i==LeptonVeto)          cut.at(LeptonVeto)=1;
    if(i==MuonIsolation)       cut.at(MuonIsolation)=0.15;
    if(i==TauIsolation)        cut.at(TauIsolation)=1;
    if(i==PairCharge)          cut.at(PairCharge)=1;
    if(i==PairMass)            cut.at(PairMass)=95;
    if(i==MTM)                 cut.at(MTM)=45;
    if(i==deltaR)              cut.at(deltaR)=3.3;
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
    if(i==TriggerOk){
      title.at(i)="Trigger ";
      hlabel="Trigger ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==nGoodPairs){
      title.at(i)="NGoodPairs";
      hlabel="Number of pairs passed baseline ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nGoodPairs_",htitle,10,-0.5,9.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nGoodPairs_",htitle,10,-0.5,9.5,hlabel,"Events"));
    }
    else if(i==LeptonVeto){
      title.at(i)="Lepton Veto";
      hlabel="Third Lepton Veto  ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_LeptonVeto_",htitle,4,-0.5,3.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_LeptonVeto_",htitle,4,-0.5,3.5,hlabel,"Events"));
    }
    else if(i==TauIsolation){
      title.at(i)="Tau Isolation";
      hlabel="Isolation of  Tau";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauIsolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0TauIsolation__",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==MuonIsolation){
      title.at(i)="Muon Isolation";
      hlabel="Muon Isolation";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonIsolation_",htitle,20,0.,0.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonIsolation_",htitle,20,0.,0.5,hlabel,"Events"));
    }
    else if(i==PairCharge){
      title.at(i)="Pair Charge";
      hlabel="is pair OS";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PairCharge_",htitle,5,-2.5,2.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PairCharge_",htitle,5,-2.5,2.5,hlabel,"Events"));
    }
    else if(i==PairMass){
      title.at(i)="Pair Visible Mass";
      hlabel="M(tau-tau)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PairMass_",htitle,50,40,120,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PairMass_",htitle,50,40,120,hlabel,"Events"));
    }

    else if(i==MTM){
      title.at(i)="Missing Transverse Mass";
      hlabel="Missing Transverse Mass";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MTM_",htitle,30,0,100,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MTM_",htitle,30,0,100,hlabel,"Events"));
    }

    else if(i==deltaR){
      title.at(i)="Delta R btw taus";
      hlabel="Delta R";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_deltaR_",htitle,20,0,8,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_deltaR_",htitle,20,0,8,hlabel,"Events"));
    }



  } 
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");


  TauPT=HConfig.GetTH1D(Name+"_TauPT","Transverse momentum of selected #tau candidate",20,27,80.5," P_{T}(#tau), GeV","Events");
  TauE=HConfig.GetTH1D(Name+"_TauE","Energy of selected #tau candidate",25,24.5,99.5," E(#tau), GeV","Events");

  MuonPT=HConfig.GetTH1D(Name+"_MuonPT","Transverse momentum of selected #tau candidate",20,27,60.5," P_{T}(#mu), GeV","Events");
  MuonE=HConfig.GetTH1D(Name+"_MuonE","Energy of selected #tau candidate",25,24.5,99.5," E(#mu), GeV","Events");
  TauHPSDecayMode=HConfig.GetTH1D(Name+"_TauHPSDecayMode","Decay mode of the selected #tau candidate",11,-0.5,10.5," HPS Mode ","Events");
  PVSVSignificance=HConfig.GetTH1D(Name+"_PVSVSignificance"," PV-SV significance for tau decay mode = 10",31,-0.5,5.5," PVSV significance ","Events");

  TauTauMass=HConfig.GetTH1D(Name+"_TauTauMass","Visible invariant mass of a tau pair",39,40 ,200," M(#tau#tau), GeV","Events");
  NQCD=HConfig.GetTH1D(Name+"_NQCD","NQCD",6,0.5,6.5,"NQCD in ABCD","Events");

  QCDShape=HConfig.GetTH1D(Name+"_QCDShape","QCDShape",2,0,2,"QCD Shape","");
  dRTauTau=HConfig.GetTH1D(Name+"_dRTauTau","#Delta R",25,0,5," #Delta R","Events");

  SVChi2=HConfig.GetTH1D(Name+"_SVChi2","SV  #chi^{2}",30,0,25,"#chi^{2}","Events");
  SVQuality=HConfig.GetTH1D(Name+"_SVQuality","Track mathicn #DeltaR",25,0,2,"#Sigma#Delta R","Events");

  NPrimeVtx=HConfig.GetTH1D(Name+"_NPrimeVtx","NPrimeVtx",100,0,100,"N vtx","Events");

  SVQualityVsSignificance=HConfig.GetTH2D(Name+"_SVQualityVsSignificance","Track mathicn #DeltaR vs significance",25,0,3,31,-0.5,5.5,"","Events");


  Selection::ConfigureHistograms();   //   do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);  // do not remove
}

 

void  ZTauMuTauH::Store_ExtraDist(){

  //every new histo should be addedd to Extradist1d vector, just push it back;
  Extradist1d.push_back(&TauPT);
  Extradist1d.push_back(&TauE);

  Extradist1d.push_back(&MuonPT);
  Extradist1d.push_back(&MuonE);
  Extradist1d.push_back(&TauHPSDecayMode);

  Extradist1d.push_back(&dRTauTau);
  Extradist1d.push_back(&TauTauMass);
  Extradist1d.push_back(&QCDShape);
  Extradist1d.push_back(&NQCD);

  Extradist1d.push_back(&NPrimeVtx);
  Extradist1d.push_back(&PVSVSignificance);

  Extradist1d.push_back(&SVChi2);
  Extradist1d.push_back(&SVQuality);
  Extradist2d.push_back(&SVQualityVsSignificance);


}

void  ZTauMuTauH::doEvent(){ //  Method called on every event
  unsigned int t;                // sample type, you may manage in your further analysis, if needed
  int id(Ntp->GetMCID());  //read event ID of a sample
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}  //  gives a warining if list of samples in Histo.txt  and SkimSummary.log do not coincide 
  //  std::cout<<"------------------ New Event -----------------------"<<std::endl;
  Charge = ChargeSumDummy;
  bool PassedTrigger(false);
  //int triggerindex;
  std::vector<int> TriggerIndex; 
  std::vector<int> TriggerIndexVector ;
  std::vector<TString>  MatchedTriggerNames;

  MatchedTriggerNames.push_back("HLT_IsoMu24_v");
  MatchedTriggerNames.push_back("HLT_IsoTkMu24_v");
  TriggerIndexVector=Ntp->GetVectorTriggers(MatchedTriggerNames);
  for(unsigned int itrig = 0; itrig < TriggerIndexVector.size(); itrig++){
    if(Ntp->TriggerAccept(TriggerIndexVector.at(itrig))){  
      //    std::cout<<"  Name  "<< Ntp->TriggerName(TriggerIndexVector.at(itrig)) << "   status   "<< Ntp->TriggerAccept(TriggerIndexVector.at(itrig)) <<std::endl;
      PassedTrigger =Ntp->TriggerAccept(TriggerIndexVector.at(itrig)); }
  }


  //  int ntau(0); int nmu(0); 
  std::vector<int> goodMuonsIndex;
  std::vector<int> thirdLeptonCounter; 
  std::vector<int> goodTauIndex;
  for(unsigned int iDaugther=0;   iDaugther  <  Ntp->NDaughters() ;iDaugther++ ){  // loop over all daughters in the event
    if(Ntp->isTightGoodTau(iDaugther)){
      if(Ntp->tauBaselineSelection(iDaugther,20., 2.1, 4,1)){
	if(Ntp->Daughters_P4(iDaugther).Pt() > cTau_pt){
	  if(fabs(  Ntp->Daughters_P4(iDaugther).Eta()) < cTau_eta  ){
	    goodTauIndex.push_back(iDaugther) ;  }}}}


    if(Ntp->isMediumGoodMuon(iDaugther)){
      if(Ntp->muonBaselineSelection(iDaugther, 25,2.1,2)){
	if(Ntp->Daughters_P4(iDaugther).Pt() > cMu_pt){
	  if(fabs(  Ntp->Daughters_P4(iDaugther).Eta()) < cMu_eta  ){
	    goodMuonsIndex.push_back(iDaugther) ;  }}}}

    if(Ntp->ElectronVeto(iDaugther) || Ntp->MuonVeto(iDaugther)){
      thirdLeptonCounter.push_back(iDaugther);
    }
  } 
    
 
  value.at(TriggerOk)=PassedTrigger;
  pass.at(TriggerOk)=PassedTrigger;

  std::vector<int>  PairsIndex;
  std::vector<int>  SSPairsIndex;

  for( int ipair =0; ipair < Ntp->NPairCandidates(); ipair++){
    //    if(Ntp->isOSCand(ipair))
    {
      if(Ntp->getPairType(Ntp->indexDau1(ipair),Ntp->indexDau2(ipair))==0){
	PairsIndex.push_back(ipair);
      }
      else SSPairsIndex.push_back(ipair);
    }
  }



  std::vector<int> Sorted;
  std::vector<int> SortedSS;
  Sorted= Ntp->SortTauHTauHPair(PairsIndex);
  SortedSS= Ntp->SortTauHTauHPair(SSPairsIndex);

  std::vector<int> SortedPair_PassedBaseline;
  std::vector<int> SortedSSPair_PassedBaseline;

  int MuonIndex= -1;
  int TauIndex= -1;
  bool OSpairFound(false);
  bool SSpairFound(false);

  for(unsigned int ipair =0; ipair < Sorted.size(); ipair++){
    //   bool PairPassBase = (Ntp->tauBaselineSelection(Ntp->indexDau1(Sorted.at(ipair)), 25,2.3,4,1) && Ntp->tauBaselineSelection(Ntp->indexDau2(Sorted.at(ipair)), 25,2.3,4,1));
    
    //    std::cout<<"   "<<Ntp->particleType(Ntp->indexDau1(Sorted.at(ipair))) <<"   "<< Ntp->particleType(Ntp->indexDau2(Sorted.at(ipair)) )<<std::endl; 
    if(Ntp->particleType(Ntp->indexDau1(Sorted.at(ipair)))==0 && Ntp->particleType(Ntp->indexDau2(Sorted.at(ipair)) ) == 2    )
      {
	MuonIndex=Ntp->indexDau1(Sorted.at(ipair));
	TauIndex=Ntp->indexDau2(Sorted.at(ipair));
	OSpairFound = true;
      }
    
    if(Ntp->particleType(Ntp->indexDau1(Sorted.at(ipair)))==2 && Ntp->particleType(Ntp->indexDau2(Sorted.at(ipair)) ) == 0    )
      {
	MuonIndex=Ntp->indexDau2(Sorted.at(ipair));
	TauIndex=Ntp->indexDau1(Sorted.at(ipair));
	OSpairFound = true;
      }

    
    if(OSpairFound){
      if(Ntp->tauBaselineSelection(TauIndex, 20,2.3,4,1) && Ntp->muonBaselineSelection(MuonIndex, 20,2.3,2)) {
	//      if(Ntp->decayMode(TauIndex)==10)
	SortedPair_PassedBaseline.push_back(Sorted.at(ipair));
      }
    }
  }
  
  for(unsigned int ipair =0; ipair < SortedSS.size(); ipair++){

    if(Ntp->particleType(Ntp->indexDau1(SortedSS.at(ipair)))==0 && Ntp->particleType(Ntp->indexDau2(SortedSS.at(ipair)) ) == 2    )
      {
	MuonIndex=Ntp->indexDau1(SortedSS.at(ipair));
	TauIndex=Ntp->indexDau2(SortedSS.at(ipair));
	SSpairFound = true;
      }
    if(Ntp->particleType(Ntp->indexDau1(SortedSS.at(ipair)))==2 && Ntp->particleType(Ntp->indexDau2(SortedSS.at(ipair)) ) == 0    )
      {
	MuonIndex=Ntp->indexDau2(SortedSS.at(ipair));
	TauIndex=Ntp->indexDau1(SortedSS.at(ipair));
	SSpairFound = true;
      }
    if(SSpairFound){
      if(Ntp->tauBaselineSelection(TauIndex, 20,2.3,4,1) && Ntp->muonBaselineSelection(MuonIndex, 20,2.3,2)){
	//      if(Ntp->decayMode(TauIndex)==10)
	SortedSSPair_PassedBaseline.push_back(SortedSS.at(ipair));
      }
    }
  }
  int Muon= -1;
  int Tau= -1;
  int Pair=-1;
	
  value.at(nGoodPairs)=SortedPair_PassedBaseline.size();
  pass.at(nGoodPairs) = (value.at(nGoodPairs) > cut.at(nGoodPairs));
  
  value.at(LeptonVeto) = thirdLeptonCounter.size();
  pass.at(LeptonVeto) = ( value.at(LeptonVeto) ==  value.at(LeptonVeto));

   
  value.at(PairCharge) = ChargeSumDummy;
  value.at(MuonIsolation) = 0;
  value.at(TauIsolation) = 0;
  value.at(deltaR) = 999.;
  value.at(PairMass) = 999.;
  value.at(MTM) = 999.;
  if(pass.at(nGoodPairs)){

    
    if(Ntp->particleType(Ntp->indexDau1(SortedPair_PassedBaseline.back()))==0    && Ntp->particleType(Ntp->indexDau2(SortedPair_PassedBaseline.back()))==2){
      Muon=Ntp->indexDau1(SortedPair_PassedBaseline.back());
      Tau=Ntp->indexDau2(SortedPair_PassedBaseline.back());
    }
    if(Ntp->particleType(Ntp->indexDau1(SortedPair_PassedBaseline.back()))==2    && Ntp->particleType(Ntp->indexDau2(SortedPair_PassedBaseline.back()))==0){
      Muon=Ntp->indexDau2(SortedPair_PassedBaseline.back());
      Tau=Ntp->indexDau1(SortedPair_PassedBaseline.back());
    }
    Pair=SortedPair_PassedBaseline.back();
    //    value.at(PairCharge) = Ntp->Daughters_charge(TauIndex_1) + Ntp->Daughters_charge(TauIndex_2);
    if(Ntp->isOSCand(SortedPair_PassedBaseline.back()))    value.at(PairCharge) =1;
    else  value.at(PairCharge) =0;
    value.at(MTM) = sqrt(2*Ntp->Daughters_P4(Muon).Pt()*Ntp->MET()*(1 - cos(Ntp->Daughters_P4(Muon).Phi() - Ntp->METphi())));
    value.at(MuonIsolation) = Ntp->combreliso(Muon);
    value.at(TauIsolation) = Ntp->isIsolatedTau(Tau,"Tight");

    value.at(deltaR) = Ntp->Daughters_P4(Muon).DeltaR(Ntp->Daughters_P4(Tau));
    value.at(PairMass) = (Ntp->Daughters_P4(Muon) +Ntp->TauP4_Corrected(Tau)).M();
   

    // std::cout<<"  isOS  "<< Ntp->isOSCand(SortedPair_PassedBaseline.back()) <<std::endl;
    // std::cout<< "  1st charge  "<< Ntp->Daughters_charge(TauIndex_1) <<" 2nd charge " << Ntp->Daughters_charge(TauIndex_2) <<std::endl;

  }
  
  pass.at(PairCharge) = (value.at(PairCharge) == cut.at(PairCharge));
  pass.at(MuonIsolation) = (value.at(MuonIsolation) <= cut.at(MuonIsolation));
  pass.at(TauIsolation) = (value.at(TauIsolation) == cut.at(TauIsolation));
  pass.at(deltaR) = (value.at(deltaR) <= cut.at(deltaR) && value.at(deltaR)  > 0.5);
  pass.at(PairMass) = true;//(value.at(PairMass) <= cut.at(PairMass));
  pass.at(MTM) = (value.at(MTM) <= cut.at(MTM));

  // Here you can defined different type of weights you want to apply to events.
  double wobs=1;
  double w=1;
  double MuonSF(1);
  double MuonIso(1);


  if(!Ntp->isData() && id!=DataMCType::QCD){
    //    w *= reweight.weight(2016,26,Ntp->PUNumInteractions());
    w *= reweight.PUweightHTT(Ntp->npu());
      //std::cout<<" pu weigh HTT  "<< reweight.PUweightHTT(Ntp->npu())<<std::endl;
    
    if(!Ntp->isData() && pass.at(nGoodPairs)){
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
   *     TauIsolation
   */

  std::vector<unsigned int> exclude_cuts;
  exclude_cuts.push_back(MuonIsolation);
  exclude_cuts.push_back(PairCharge);
  // std::cout<<" before  " << pass.at(TriggerOk) << "    " <<   pass.at(PrimeVtx) << "    " <<  pass.at(nGoodPairs)<< "    " <<   pass.at(FirstTauIsolation) << "    " <<  pass.at(SecondTauIsolation) << "    " <<  pass.at(nGoodMuons) << "    " <<  pass.at(PairCharge) << "  passAllBut  " << passAllBut(exclude_cuts) <<std::endl;

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

    double pvx(0);
    pvx =  Ntp->npv();
    // if(id == DataMCType::Data) pvx =  Ntp->npv();
     if(id !=DataMCType::Data && id !=DataMCType::QCD)	  pvx = Ntp->PUNumInteractions();

     if(Ntp->PFTau_hassecondaryVertex(Tau) && Ntp->isPVCovAvailable()){
       if(Ntp->PFTau_secondaryVertex_TracksMatchingQuality(Tau) < 0.01)   PVSVSignificance.at(t).Fill( Ntp->PFTau_FlightLength_significance(Ntp->PVtx(),Ntp->PFTau_TIP_primaryVertex_cov(), Ntp->PFTau_secondaryVertex_pos(Tau), Ntp->PFTau_TIP_secondaryVertex_cov(Tau)),w);
      
       SVChi2.at(t).Fill(Ntp->PFTau_secondaryVertex_vtxchi2(Tau),w);
       SVQuality.at(t).Fill(Ntp->PFTau_secondaryVertex_TracksMatchingQuality(Tau),w);
       SVQualityVsSignificance.at(t).Fill(Ntp->PFTau_secondaryVertex_TracksMatchingQuality(Tau),Ntp->PFTau_FlightLength_significance(Ntp->PVtx(),Ntp->PFTau_TIP_primaryVertex_cov(), Ntp->PFTau_secondaryVertex_pos(Tau), Ntp->PFTau_TIP_secondaryVertex_cov(Tau)));
     }



  NPrimeVtx.at(t).Fill(pvx,w);
  TLorentzVector MuonP4 = Ntp->Daughters_P4(Muon);
  TLorentzVector TauP4 = Ntp->TauP4_Corrected(Tau);


  TauPT.at(t).Fill(TauP4.Pt(),w);  // Fill transverse momentum
  TauE.at(t).Fill(TauP4.E(),w);  // Fill transverse momentum

  MuonPT.at(t).Fill(MuonP4.Pt(),w);  // Fill transverse momentum
  MuonE.at(t).Fill(MuonP4.E(),w);  // Fill transverse momentum

  TauHPSDecayMode.at(t).Fill(Ntp->decayMode(Tau),w);
 
  TauTauMass.at(t).Fill((MuonP4+TauP4).M(),w);
  dRTauTau.at(t).Fill(MuonP4.DeltaR(TauP4),w);
  }
}




//  This is a function if you want to do something after the event loop
void  ZTauMuTauH::Finish(){
  if(mode == RECONSTRUCT){
    std::cout<<" Starting Finish!  " <<std::endl;
    
    std::cout<<"A  Data  "<< NQCD.at(0).GetBinContent(1) << std::endl;
    std::cout<<"B  Data  "<< NQCD.at(0).GetBinContent(2) << std::endl;
    std::cout<<"C  Data  "<< NQCD.at(0).GetBinContent(3) << std::endl;
    std::cout<<"D  Data  "<< NQCD.at(0).GetBinContent(4) << std::endl;
    SkimConfig SC;
    SC.ApplySkimEfficiency(types,Npassed, Npassed_noweight);

    std::vector<double> QCD_Integral_B;
    double QCD_IntegralMC_B;
    double QCD_Integral_B_Data_minus_MC = 0;
    
    std::vector<double> QCD_Integral_C;
    double QCD_IntegralMC_C;
    double QCD_Integral_C_Data_minus_MC = 0;
    
    std::vector<double> QCD_Integral_D;
    double QCD_IntegralMC_D;
    double QCD_Integral_D_Data_minus_MC = 0;

    //Get Yields in ABCD for QCD Scalefactor                                                                                                                                                                  
    for(unsigned i=0;i<CrossSectionandAcceptance.size();i++){
      QCD_Integral_B.push_back(NQCD.at(i).GetBinContent(2));
      QCD_Integral_C.push_back(NQCD.at(i).GetBinContent(3));
      QCD_Integral_D.push_back(NQCD.at(i).GetBinContent(4));
      if(CrossSectionandAcceptance.at(i)>0){
	QCD_Integral_B.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
	QCD_Integral_C.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
	QCD_Integral_D.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
      }
    }
    for(unsigned i=0;i<CrossSectionandAcceptance.size();i++){
      if(HConfig.GetID(i) == DataMCType::Data){
	QCD_Integral_B_Data_minus_MC  += QCD_Integral_B.at(i);
	QCD_Integral_C_Data_minus_MC += QCD_Integral_C.at(i);
	QCD_Integral_D_Data_minus_MC += QCD_Integral_D.at(i);
      }
      if(CrossSectionandAcceptance.at(i)>0){
	QCD_IntegralMC_B  += QCD_Integral_B.at(i);
	QCD_IntegralMC_C  += QCD_Integral_C.at(i);
	QCD_IntegralMC_D  += QCD_Integral_D.at(i);
      }
    }

    double OS2SS = (QCD_Integral_C_Data_minus_MC  - QCD_IntegralMC_C )/ (QCD_Integral_D_Data_minus_MC - QCD_IntegralMC_D);
    double QCD_ScaleFactor = QCD_Integral_B_Data_minus_MC *OS2SS;


    std::cout << "OS/SS QCD Sample: " << OS2SS << std::endl;
    std::cout << "Scale Factor for QCD Sample: " << QCD_ScaleFactor << std::endl;
    std::cout << "QCD in B region "<<  QCD_Integral_B_Data_minus_MC <<std::endl;
    std::cout << "QCD_Integral_B_Data_minus_MC is: " << QCD_Integral_B_Data_minus_MC << std::endl;
    std::cout << "QCD_Integral_C_Data_minus_MC is: " << QCD_Integral_C_Data_minus_MC << std::endl;
    std::cout << "QCD_Integral_D_Data_minus_MC is: " << QCD_Integral_D_Data_minus_MC << std::endl;
    std::cout << "QCD_IntegralMC_B is: " << QCD_IntegralMC_B << std::endl;
    std::cout << "QCD_IntegralMC_C is: " << QCD_IntegralMC_C << std::endl;
    std::cout << "QCD_IntegralMC_D is: " << QCD_IntegralMC_D << std::endl;
    ScaleAllHistOfType(HConfig.GetType(DataMCType::QCD),QCD_ScaleFactor/Nminus0.at(0).at(HConfig.GetType(DataMCType::QCD)).Integral());
  }

  Selection::Finish();
}






