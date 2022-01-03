#include "ZTauHTauH.h"
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




ZTauHTauH::ZTauHTauH(TString Name_, TString id_):
  Selection(Name_,id_),
  cMu_pt(20),
  cMu_eta(2.1),
  cTau_pt(20),
  cTau_eta(2.1),
  tauTrgSF("vtight"),
  DataMC_Corr(true)
{
  ChargeSumDummy = -999;
  selMuon_IsoDummy = 999.;
}

ZTauHTauH::~ZTauHTauH(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  ZTauHTauH::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)           cut.at(TriggerOk)=1;
    if(i==nGoodPairs)          cut.at(nGoodPairs)=0;
    if(i==LeptonVeto)          cut.at(LeptonVeto)=0;
    if(i==FirstTauIsolation)   cut.at(FirstTauIsolation)=1;
    if(i==SecondTauIsolation)  cut.at(SecondTauIsolation)=1;
    if(i==nGoodMuons)          cut.at(nGoodMuons)=0;
    if(i==PairCharge)          cut.at(PairCharge)=1;
    if(i==PairMass)            cut.at(PairMass)=115;
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
    else if(i==FirstTauIsolation){
      title.at(i)="First Tau Isolation";
      hlabel="Isolation of First Tau";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_FirstTauIsolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_FirstTauIsolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==SecondTauIsolation){
      title.at(i)="Second Tau Isolation";
      hlabel="Isolation of Second Tau";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_SecondTauIsolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_SecondTauIsolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==nGoodMuons){
      title.at(i)="Number of muons";
      hlabel="Number of muons";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nGoodMuons_",htitle,10,-0.5,9.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nGoodMuons_",htitle,10,-0.5,9.5,hlabel,"Events"));
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

    else if(i==deltaR){
      title.at(i)="Delta R btw taus";
      hlabel="Delta R";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_deltaR_",htitle,20,0,8,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_deltaR_",htitle,20,0,8,hlabel,"Events"));
    }



  } 
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");

  Tau1PT=HConfig.GetTH1D(Name+"_Tau1PT","Transverse momentum of selected #tau candidate",40,0,200," P_{T}(#tau), GeV","Events");
  Tau1E=HConfig.GetTH1D(Name+"_Tau1E","Energy of selected #tau candidate",20,24.5,99.5," E(#tau), GeV","Events");
  Tau1HPSDecayMode=HConfig.GetTH1D(Name+"_Tau1HPSDecayMode","Decay mode of the selected #tau candidate",11,-0.5,10.5," HPS Mode ","Events");

  Tau2PT=HConfig.GetTH1D(Name+"_Tau2PT","Transverse momentum of selected #tau candidate",40,0,200," P_{T}(#tau), GeV","Events");
  Tau2E=HConfig.GetTH1D(Name+"_Tau2E","Energy of selected #tau candidate",20,24.5,99.5," E(#tau), GeV","Events");
  Tau2HPSDecayMode=HConfig.GetTH1D(Name+"_Tau2HPSDecayMode","Decay mode of the selected #tau candidate",11,-0.5,10.5," HPS Mode ","Events");


  TauTauMass=HConfig.GetTH1D(Name+"_TauTauMass","Visible invariant mass of a tau pair",30,0 ,300," M(#tau#tau), GeV","Events");
  NQCD=HConfig.GetTH1D(Name+"_NQCD","NQCD",6,0.5,6.5,"NQCD in ABCD","Events");

  QCDShape=HConfig.GetTH1D(Name+"_QCDShape","QCDShape",2,0,2,"QCD Shape","");
  dRTauTau=HConfig.GetTH1D(Name+"_dRTauTau","#Delta R",25,0,5," #Delta R","Events");

  Tau1Isolation=HConfig.GetTH1D(Name+"_Tau1Isolation","First Tau Isoaltion 1- Loose, 2- Medium, 3 Tight, 4-VTight",5,0.5,5.5," Discrimiantor","Events");
  Tau2Isolation=HConfig.GetTH1D(Name+"_Tau2Isolation","First Tau Isoaltion 1- Loose, 2- Medium, 3 Tight, 4-VTight",5,0.5,5.5," Discrimiantor","Events");

  MET=HConfig.GetTH1D(Name+"_MET","Missing transverse energy",30,0,180,"Missing transverse energy","Events");
  TauHMass1=HConfig.GetTH1D(Name+"_TauHMass1","InvariantMass of tauH",30,0.3,1.7," M, GeV","Events");
  TauHMass2=HConfig.GetTH1D(Name+"_TauHMass2","InvariantMass of tauH",30,0.3,1.7," M, GeV","Events");

  Tau1Eta=HConfig.GetTH1D(Name+"_Tau1Eta","Pseudorapidity tau1",15,-3,3," #eta","Events");
  Tau2Eta=HConfig.GetTH1D(Name+"_Tau2Eta","Pseudorapidity tau2",15,-3,3," #eta","Events");







  NPrimeVtx=HConfig.GetTH1D(Name+"_NPrimeVtx","NPrimeVtx",20,0,40,"N vtx","Events");
    Selection::ConfigureHistograms();   //   do not remove
    HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);  // do not remove
}

 

void  ZTauHTauH::Store_ExtraDist(){ 

  //every new histo should be addedd to Extradist1d vector, just push it back;
  Extradist1d.push_back(&Tau1PT);
  Extradist1d.push_back(&Tau1E);
  Extradist1d.push_back(&Tau1HPSDecayMode);

  Extradist1d.push_back(&Tau2PT);
  Extradist1d.push_back(&Tau2E);
  Extradist1d.push_back(&Tau2HPSDecayMode);

  Extradist1d.push_back(&dRTauTau);
  Extradist1d.push_back(&TauTauMass);
  Extradist1d.push_back(&QCDShape);
  Extradist1d.push_back(&NQCD);
  Extradist1d.push_back(&Tau1Isolation);
  Extradist1d.push_back(&Tau2Isolation);

  Extradist1d.push_back(&MET);
  Extradist1d.push_back(&TauHMass1);
  Extradist1d.push_back(&TauHMass2);


  Extradist1d.push_back(&Tau1Eta);
  Extradist1d.push_back(&Tau2Eta);


  Extradist1d.push_back(&NPrimeVtx);

}

void  ZTauHTauH::doEvent(){ //  Method called on every event
  unsigned int t;                // sample type, you may manage in your further analysis, if needed
  int id(Ntp->GetMCID());  //read event ID of a sample
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}  //  gives a warining if list of samples in Histo.txt  and SkimSummary.log do not coincide 
  //  std::cout<<"------------------ New Event -----------------------"<<std::endl;
  Charge = ChargeSumDummy;

  
  //  int ntau(0); int nmu(0); 
  std::vector<int> goodMuonsIndex;
  std::vector<int> thirdLeptonCounter; 
  std::vector<int> goodTauIndex;
  for(unsigned int iDaugther=0;   iDaugther  <  Ntp->NDaughters() ;iDaugther++ ){  // loop over all daughters in the event
    if(Ntp->isTightGoodTau(iDaugther)){
      if(Ntp->tauBaselineSelection(iDaugther,25., 2.1, 4,1)){
	if(Ntp->Daughters_P4(iDaugther).Pt() > cTau_pt){
	  if(fabs(  Ntp->Daughters_P4(iDaugther).Eta()) < cTau_eta  ){
	    goodTauIndex.push_back(iDaugther) ;  }}}}


    if(Ntp->isMediumGoodMuon(iDaugther)){
      if(Ntp->muonBaselineSelection(iDaugther,20,2.1,2)){
	goodMuonsIndex.push_back(iDaugther) ;  }}
    
    if(Ntp->ElectronVeto(iDaugther) || Ntp->MuonVeto(iDaugther)){
      thirdLeptonCounter.push_back(iDaugther);
    }
  }
    
 
  value.at(TriggerOk)=1;
  pass.at(TriggerOk)=true;

  std::vector<int>  PairsIndex;
  std::vector<int>  SSPairsIndex;

  for(unsigned int ipair =0; ipair < Ntp->NPairCandidates(); ipair++){
    //    if(Ntp->isOSCand(ipair))
    {
      if(Ntp->getPairType(Ntp->indexDau1(ipair),Ntp->indexDau2(ipair))==2){
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

  for(unsigned int ipair =0; ipair < Sorted.size(); ipair++){
    //   bool PairPassBase = (Ntp->tauBaselineSelection(Ntp->indexDau1(Sorted.at(ipair)), 25,2.3,4,1) && Ntp->tauBaselineSelection(Ntp->indexDau2(Sorted.at(ipair)), 25,2.3,4,1));
    if(Ntp->tauBaselineSelection(Ntp->indexDau1(Sorted.at(ipair)), 35,2.3,4,1) && Ntp->tauBaselineSelection(Ntp->indexDau2(Sorted.at(ipair)), 35,2.3,4,1)) 
      SortedPair_PassedBaseline.push_back(Sorted.at(ipair));
  }
  for(unsigned int ipair =0; ipair < SortedSS.size(); ipair++){
    if(Ntp->tauBaselineSelection(Ntp->indexDau1(SortedSS.at(ipair)), 35,2.3,4,1) && Ntp->tauBaselineSelection(Ntp->indexDau2(SortedSS.at(ipair)), 35,2.3,4,1)) 
      SortedSSPair_PassedBaseline.push_back(SortedSS.at(ipair));
  }
  


	
  value.at(nGoodPairs)=SortedPair_PassedBaseline.size();
  pass.at(nGoodPairs) = (value.at(nGoodPairs) > cut.at(nGoodPairs));
  
  value.at(LeptonVeto) = thirdLeptonCounter.size();
  pass.at(LeptonVeto) = ( value.at(LeptonVeto) ==  value.at(LeptonVeto));

  value.at(nGoodMuons)=goodMuonsIndex.size();
  pass.at(nGoodMuons) =(value.at(nGoodMuons) == cut.at(nGoodMuons));
  
  int TauIndex_1= -1;
  int TauIndex_2= -1;
  value.at(PairCharge) = ChargeSumDummy;
  value.at(FirstTauIsolation) = 0;
  value.at(SecondTauIsolation) = 0;
  value.at(deltaR) = 999.;
  value.at(PairMass) = 999.;
  if(pass.at(nGoodPairs)){
    TauIndex_1 = Ntp->indexDau1(SortedPair_PassedBaseline.back());
    TauIndex_2 = Ntp->indexDau2(SortedPair_PassedBaseline.back());
    //    value.at(PairCharge) = Ntp->Daughters_charge(TauIndex_1) + Ntp->Daughters_charge(TauIndex_2);
    if(Ntp->isOSCand(SortedPair_PassedBaseline.back()))    value.at(PairCharge) =1;
    else  value.at(PairCharge) =0;
    value.at(FirstTauIsolation) = Ntp->isIsolatedTau(TauIndex_1,"VTight");
    value.at(SecondTauIsolation) = Ntp->isIsolatedTau(TauIndex_2,"VTight");

    value.at(deltaR) = Ntp->Daughters_P4(TauIndex_1).DeltaR(Ntp->Daughters_P4(TauIndex_2));
    value.at(PairMass) = (Ntp->Daughters_P4(TauIndex_1) +Ntp->Daughters_P4(TauIndex_2)).M();
   

    // std::cout<<"  isOS  "<< Ntp->isOSCand(SortedPair_PassedBaseline.back()) <<std::endl;
    // std::cout<< "  1st charge  "<< Ntp->Daughters_charge(TauIndex_1) <<" 2nd charge " << Ntp->Daughters_charge(TauIndex_2) <<std::endl;

  }
  pass.at(PairCharge) = (value.at(PairCharge) == cut.at(PairCharge));
  pass.at(FirstTauIsolation) = (value.at(FirstTauIsolation) == cut.at(FirstTauIsolation));
  pass.at(SecondTauIsolation) = (value.at(SecondTauIsolation) == cut.at(SecondTauIsolation));
  pass.at(deltaR) = (value.at(deltaR) <= cut.at(deltaR));
  pass.at(PairMass) = (value.at(PairMass) <= cut.at(PairMass));

  // Here you can defined different type of weights you want to apply to events.
  double wobs=1;
  double w=1;
  if(!Ntp->isData() && id!=DataMCType::QCD){
    //    w *= reweight.weight(2016,26,Ntp->PUNumInteractions());
    w *= reweight.PUweightHTT(Ntp->npu());
      //std::cout<<" pu weigh HTT  "<< reweight.PUweightHTT(Ntp->npu())<<std::endl;
    if(!Ntp->isData() && pass.at(nGoodPairs) ){
       double w1 = tauTrgSF.getSF(Ntp->Daughters_P4(TauIndex_1).Pt(),  Ntp->decayMode(TauIndex_1)) ;  //from Luca
       double w2 = tauTrgSF.getSF(Ntp->Daughters_P4(TauIndex_2).Pt(),  Ntp->decayMode(TauIndex_2)) ;
       // w*=w1;
       // w*=w2;
    }
    if(!Ntp->isData() && pass.at(nGoodPairs)){
      //   w *= 0.95;  // Tau ID  correction
    }
  }
  //------------------ Z Pt weights --------------------
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
    //-------------------------  mu/e tau fake rate weights 
    double wAgainstMuon1(1);
    double wAgainstElectron1(1);
    double wAgainstMuon2(1);
    double wAgainstElectron2(1);
    if(id == 33){
      if(pass.at(nGoodPairs)){
	TauIndex_1 = Ntp->indexDau1(SortedPair_PassedBaseline.back());
	TauIndex_2 = Ntp->indexDau2(SortedPair_PassedBaseline.back());
	int matchedIndex1(-1);
	int matchedIndex2(-1);
	double dR1(999);
	double dR2(999);
	for(unsigned int imc=0; imc < Ntp->NGenParts(); imc++){
	  if(fabs(Ntp->Genpart_pdg(imc)) ==11 || fabs(Ntp->Genpart_pdg(imc)) ==13)
	    {
	      if(sqrt(pow(Ntp->TauP4_Corrected(Ntp->indexDau1(SortedPair_PassedBaseline.back())).Phi() - Ntp->Genpart_P4(imc).Phi(),2) + 
		      pow(Ntp->TauP4_Corrected(Ntp->indexDau1(SortedPair_PassedBaseline.back())).Eta() - Ntp->Genpart_P4(imc).Eta(),2)) < dR1)
		{
		  dR1 = sqrt(pow(Ntp->TauP4_Corrected(Ntp->indexDau1(SortedPair_PassedBaseline.back())).Phi() - Ntp->Genpart_P4(imc).Phi(),2) + 
			     pow(Ntp->TauP4_Corrected(Ntp->indexDau1(SortedPair_PassedBaseline.back())).Eta() - Ntp->Genpart_P4(imc).Eta(),2));
		  matchedIndex1=imc;
		}
	      if(sqrt(pow(Ntp->TauP4_Corrected(Ntp->indexDau2(SortedPair_PassedBaseline.back())).Phi() - Ntp->Genpart_P4(imc).Phi(),2) + 
		      pow(Ntp->TauP4_Corrected(Ntp->indexDau2(SortedPair_PassedBaseline.back())).Eta() - Ntp->Genpart_P4(imc).Eta(),2)) < dR2)
		{
		  dR2 = sqrt(pow(Ntp->TauP4_Corrected(Ntp->indexDau2(SortedPair_PassedBaseline.back())).Phi() - Ntp->Genpart_P4(imc).Phi(),2) + 
			     pow(Ntp->TauP4_Corrected(Ntp->indexDau2(SortedPair_PassedBaseline.back())).Eta() - Ntp->Genpart_P4(imc).Eta(),2));
		  matchedIndex2=imc;
		}
	    }
	}
	if(dR1 < 0.2  && matchedIndex1!=-1 ){
	  if(fabs(Ntp->Genpart_pdg(matchedIndex1)) ==13 &&  ( Ntp->CHECK_BIT(Ntp->Genpart_flags(matchedIndex1),0) || Ntp->CHECK_BIT(Ntp->Genpart_flags(matchedIndex1),5)   ) )
	    {
	      wAgainstMuon1 = DataMC_Corr.AgainstMuonDataMCCorrection(Ntp->TauP4_Corrected(TauIndex_1),"AgainstMuonMVATight3");
	    }
	  if(fabs(Ntp->Genpart_pdg(matchedIndex1)) ==11 &&  ( Ntp->CHECK_BIT(Ntp->Genpart_flags(matchedIndex1),0) || Ntp->CHECK_BIT(Ntp->Genpart_flags(matchedIndex1),5)   ) )
	    {
	      wAgainstElectron1 = DataMC_Corr.AgainstElectronDataMCCorrection(Ntp->TauP4_Corrected(TauIndex_1),"AgainstElectronMVATight");
	    }

	}
	if(dR2 < 0.2  && matchedIndex2!=-1 ){
	  if(fabs(Ntp->Genpart_pdg(matchedIndex2)) ==13 &&  ( Ntp->CHECK_BIT(Ntp->Genpart_flags(matchedIndex2),0) || Ntp->CHECK_BIT(Ntp->Genpart_flags(matchedIndex2),5)   ) )
	    {
	      wAgainstMuon2 = DataMC_Corr.AgainstMuonDataMCCorrection(Ntp->TauP4_Corrected(TauIndex_2),"AgainstMuonMVATight3");
	    }
	  if(fabs(Ntp->Genpart_pdg(matchedIndex2)) ==11 &&  ( Ntp->CHECK_BIT(Ntp->Genpart_flags(matchedIndex2),0) || Ntp->CHECK_BIT(Ntp->Genpart_flags(matchedIndex2),5)   ) )
	    {
	      wAgainstElectron2 = DataMC_Corr.AgainstElectronDataMCCorrection(Ntp->TauP4_Corrected(TauIndex_2),"AgainstElectronMVATight");
	    }
	}
      }
    }

     w*=wAgainstMuon1;
     w*=wAgainstElectron1;
     w*=wAgainstMuon2;
     w*=wAgainstElectron2;



 // QCD ABCD BG Method
  /*******************
   *        |        *
   *    C   |    D   *  SS
   *        |        *       S
   * ----------------*------ i
   *        |        *       g
   *    A   |    B   *  OS   n
   *        |        *
   *******************
   *  Iso   | AntiIso
   *
   *     TauIsolation
   */

  std::vector<unsigned int> exclude_cuts;
  exclude_cuts.push_back(FirstTauIsolation);
  exclude_cuts.push_back(SecondTauIsolation);
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
      if(pass.at(FirstTauIsolation) && pass.at(SecondTauIsolation)){
	NQCD.at(t).Fill(1.,w); //A
      }
      if(!Ntp->isIsolatedTau(TauIndex_1,"VTight") && Ntp->isIsolatedTau(TauIndex_2,"Loose")){
	NQCD.at(t).Fill(2.,w); //B
      }
    }
    
    if(!pass.at(PairCharge)){
      if(pass.at(FirstTauIsolation) && pass.at(SecondTauIsolation)){
	NQCD.at(t).Fill(3.,w); //С
      }
      if(!Ntp->isIsolatedTau(TauIndex_1,"VTight") && Ntp->isIsolatedTau(TauIndex_2,"Loose")){
	NQCD.at(t).Fill(4.,w); //D
      }
    }
  }

  bool IsQCDEvent = false;
  //  if(passAllBut(exclude_cuts)){
    if(pass.at(PairCharge)){
      if(!Ntp->isIsolatedTau(TauIndex_1,"VTight") && Ntp->isIsolatedTau(TauIndex_2,"Loose")){
	if(id == DataMCType::Data){
	  QCDShape.at(t).Fill(1,w);
	  t=HConfig.GetType(DataMCType::QCD);
	  IsQCDEvent = true;
	} 
      }
    }
    //  }
  if(IsQCDEvent){    pass.at(PairCharge)= true;pass.at(FirstTauIsolation)= true;pass.at(SecondTauIsolation)= true;}
     

  std::vector<unsigned int> exclude_cuts_ForTauIso;
  exclude_cuts_ForTauIso.push_back(FirstTauIsolation);
  exclude_cuts_ForTauIso.push_back(SecondTauIsolation);

 if(passAllBut(exclude_cuts_ForTauIso)){
   if(Ntp->isIsolatedTau(TauIndex_1,"Loose"))Tau1Isolation.at(t).Fill(1.);
   if(Ntp->isIsolatedTau(TauIndex_1,"Medium"))Tau1Isolation.at(t).Fill(2.);
   if(Ntp->isIsolatedTau(TauIndex_1,"Tight"))Tau1Isolation.at(t).Fill(3.);
   if(Ntp->isIsolatedTau(TauIndex_1,"VTight"))Tau1Isolation.at(t).Fill(4.);
   if(Ntp->isIsolatedTau(TauIndex_2,"Loose"))Tau2Isolation.at(t).Fill(1.);
   if(Ntp->isIsolatedTau(TauIndex_2,"Medium"))Tau2Isolation.at(t).Fill(2.);
   if(Ntp->isIsolatedTau(TauIndex_2,"Tight"))Tau2Isolation.at(t).Fill(3.);
   if(Ntp->isIsolatedTau(TauIndex_2,"VTight"))Tau2Isolation.at(t).Fill(4.);
 }


  bool status=AnalysisCuts(t,w,wobs);  // boolean that say whether your event passed critera defined in pass vector. The whole vector must be true for status = true
  ///////////////////////////////////////////////////////////
  // Analyse events  which passed selection
  if(status){

    double pvx(0);
    pvx =  Ntp->npv();
    // if(id == DataMCType::Data) pvx =  Ntp->npv();
     if(id !=DataMCType::Data && id !=DataMCType::QCD)	  pvx = Ntp->PUNumInteractions();

  NPrimeVtx.at(t).Fill(pvx,w);
  TLorentzVector Tau1P4 = Ntp->Daughters_P4(TauIndex_1);
  TLorentzVector Tau2P4 = Ntp->Daughters_P4(TauIndex_2);


  Tau1PT.at(t).Fill(Tau1P4.Pt(),w);  // Fill transverse momentum
  Tau1E.at(t).Fill(Tau1P4.E(),w);  // Fill transverse momentum
  Tau1HPSDecayMode.at(t).Fill(Ntp->decayMode(TauIndex_1),w);

  Tau2PT.at(t).Fill(Tau2P4.Pt(),w);  // Fill transverse momentum
  Tau2E.at(t).Fill(Tau2P4.E(),w);  // Fill transverse momentum
  Tau2HPSDecayMode.at(t).Fill(Ntp->decayMode(TauIndex_2),w);
  TauTauMass.at(t).Fill((Tau1P4+Tau2P4).M(),w);
  dRTauTau.at(t).Fill(Tau1P4.DeltaR(Tau2P4),w);
  Tau1Eta.at(t).Fill(Tau2P4.Eta(),w); 
  Tau2Eta.at(t).Fill(Tau2P4.Eta(),w); 


  MET.at(t).Fill(Ntp->MET(),w);
  TauHMass1.at(t).Fill(Tau1P4.M(),w); 
  TauHMass2.at(t).Fill(Tau2P4.M(),w); 


  }
}




//  This is a function if you want to do something after the event loop
void  ZTauHTauH::Finish(){
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
    double QCD_ScaleFactor = QCD_Integral_B_Data_minus_MC *OS2SS*1.1;


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






