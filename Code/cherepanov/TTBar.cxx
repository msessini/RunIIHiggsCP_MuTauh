#include "TTBar.h"
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




TTBar::TTBar(TString Name_, TString id_):
  Selection(Name_,id_),
  cMu_pt(20),
  cMu_eta(2.1),
  cTau_pt(30),
  cTau_eta(2.1),
  tauTrgSF("vtight"),
  DataMC_Corr(true)
{
  ChargeSumDummy = -999;
  selMuon_IsoDummy = 999.;
}

TTBar::~TTBar(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  TTBar::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)           cut.at(TriggerOk)=1;
    if(i==LeptonVeto)          cut.at(LeptonVeto)=0;
    if(i==nGoodTaus)           cut.at(nGoodTaus)=1;
    if(i==nGoodMuons)          cut.at(nGoodMuons)=1;
    if(i==nBJets)              cut.at(nBJets)=0;
    if(i==nJets)               cut.at(nJets)=3;
    if(i==ET)                  cut.at(ET)=40;
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
    else if(i==nGoodTaus){
      title.at(i)="NGoodTaus";
      hlabel="Number of good taus";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NGoodTaus_",htitle,10,-0.5,9.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NGoodTaus_",htitle,10,-0.5,9.5,hlabel,"Events"));
    }
    else if(i==nGoodMuons){
      title.at(i)="Number of muons";
      hlabel="Number of muons";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nGoodMuons_",htitle,10,-0.5,9.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nGoodMuons_",htitle,10,-0.5,9.5,hlabel,"Events"));
    }
    else if(i==nBJets){
      title.at(i)="Number of b tagged jets";
      hlabel="Number of b tagged jets";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nBJets_",htitle,10,-0.5,9.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nBJets_",htitle,10,-0.5,9.5,hlabel,"Events"));
    }
    else if(i==nJets){
      title.at(i)="Number of jets";
      hlabel="Number of jets";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_nJets_",htitle,10,-0.5,9.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_nJets_",htitle,10,-0.5,9.5,hlabel,"Events"));
    }
    else if(i==ET){
      title.at(i)="Missing energy";
      hlabel="Missing energy";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ET_",htitle,50,0,150,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ET_",htitle,50,0,150,hlabel,"Events"));
    }
    else if(i==LeptonVeto){
      title.at(i)="Lepton Veto";
      hlabel="Third Lepton Veto  ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_LeptonVeto_",htitle,4,-0.5,3.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_LeptonVeto_",htitle,4,-0.5,3.5,hlabel,"Events"));
    }





  } 
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");





   LeadintTauPt=HConfig.GetTH1D(Name+"_LeadintTauPt","Leading #tau candidate pT",30,0,250,"","Events");
   LeadintMuonPt=HConfig.GetTH1D(Name+"_LeadintMuonPt","Leading #mu candidate pT",30,0,250,"","Events");

   LeadintJetPt=HConfig.GetTH1D(Name+"_LeadintJetPt","Leading jet candidate pT",30,0,250,"","Events");


   MTM=HConfig.GetTH1D(Name+"_MTM","MissingTransverseMass",15,0,300,"","Events");

   BDiscr1=HConfig.GetTH1D(Name+"_BDiscr1","bDiscriminator",20,-2,2,"","Events");
   BDiscr2=HConfig.GetTH1D(Name+"_BDiscr2","bCSVscore",20,0,1,"","Events");
   BDiscr3=HConfig.GetTH1D(Name+"_BDiscr3","pfCombinedMVAV2BJetTags",20,-2,2,"","Events");
   AllJetsBTag=HConfig.GetTH1D(Name+"_AllJetsBTag","AllJetsBTag",20,0,1,"","Events");


   // emissed=HConfig.GetTH1D(Name+"_emissed","emissed",50,0,100,"","Events");
   // emisseduncorr=HConfig.GetTH1D(Name+"_emissed","emissed",50,0,100,"","Events");

// mtxuncorr
// mtyuncorr

 

  
   NPrimeVtx=HConfig.GetTH1D(Name+"_NPrimeVtx","NPrimeVtx",10,0,40,"N vtx","Events");
   Selection::ConfigureHistograms();   //   do not remove
   HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);  // do not remove
}

 

void  TTBar::Store_ExtraDist(){

  //every new histo should be addedd to Extradist1d vector, just push it back;





   Extradist1d.push_back(&BDiscr1);
   Extradist1d.push_back(&BDiscr2);
   Extradist1d.push_back(&BDiscr3);
   Extradist1d.push_back(&AllJetsBTag);

   Extradist1d.push_back(&MTM);
   Extradist1d.push_back(&LeadintTauPt);
   Extradist1d.push_back(&LeadintMuonPt);
   Extradist1d.push_back(&LeadintJetPt);
   Extradist1d.push_back(&NPrimeVtx);



}

void  TTBar::doEvent(){ //  Method called on every event
  unsigned int t;                // sample type, you may manage in your further analysis, if needed
  int id(Ntp->GetMCID());  //read event ID of a sample
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}  //  gives a warining if list of samples in Histo.txt  and SkimSummary.log do not coincide 
  //  std::cout<<"------------------ New Event -----------------------"<<std::endl;
  Charge = ChargeSumDummy;

  
  //  int ntau(0); int nmu(0); 
  std::vector<int> goodMuonsIndex;
  std::vector<int> nBJetsIndex; 
  std::vector<int> nJetsIndex; 
  std::vector<int> nAk8JetsIndex; 
  std::vector<int> goodTauIndex;
  std::vector<int> thirdLeptonCounter; 

  for(unsigned int iDaugther=0;   iDaugther  <  Ntp->NDaughters() ;iDaugther++ ){  // loop over all daughters in the event
    if(Ntp->isTightGoodTau(iDaugther)){
      if(Ntp->tauBaselineSelection(iDaugther,30., 2.1, 4,1)){
	if(Ntp->isIsolatedTau(iDaugther,"Tight")){
	  goodTauIndex.push_back(iDaugther) ;  }}}


    if(Ntp->isMediumGoodMuon(iDaugther)){
      if(Ntp->muonBaselineSelection(iDaugther,27,2.4,2)){
	if(Ntp->combreliso(iDaugther)){
	  goodMuonsIndex.push_back(iDaugther) ;  }}}



    if(Ntp->ElectronVeto(iDaugther) || Ntp->MuonVeto(iDaugther)){
      thirdLeptonCounter.push_back(iDaugther);
    }
  }

   
  for(unsigned int ijet =0; ijet < Ntp->NJets(); ijet++){
    if(Ntp->Jet_P4(ijet).Pt()> 30 && fabs(Ntp->Jet_P4(ijet).Eta()) < 2.5){
      if(Ntp->PFjetID(ijet) == 3){
	if(Ntp->bCSVscore(ijet) > 0.875){
	  nBJetsIndex.push_back(ijet);
	}
      }
    }
  }

  for(unsigned int ijet =0; ijet < Ntp->NAK8Jets(); ijet++){
    if(Ntp->AK8Jet_P4(ijet).Pt()> 10 && fabs(Ntp->AK8Jet_P4(ijet).Eta()) < 2.5){
      // if(Ntp->PFjetID(ijet) == 2){
	//   
	nAk8JetsIndex.push_back(ijet);
	//	}
      }
  }




  for(unsigned int ijet =0; ijet < Ntp->NJets(); ijet++){
    if(Ntp->Jet_P4(ijet).Pt()> 30 && fabs(Ntp->Jet_P4(ijet).Eta()) < 2.5){
      if(Ntp->PFjetID(ijet) == 3){
	nJetsIndex.push_back(ijet);
      }
    }
  }

  value.at(TriggerOk)=1;
  pass.at(TriggerOk)=true;

  value.at(LeptonVeto) = thirdLeptonCounter.size();
  pass.at(LeptonVeto) = ( value.at(LeptonVeto) ==  value.at(LeptonVeto));


  value.at(nGoodMuons)=goodMuonsIndex.size();
  pass.at(nGoodMuons) =(value.at(nGoodMuons) == cut.at(nGoodMuons));

  value.at(nGoodTaus)=goodTauIndex.size();
  pass.at(nGoodTaus) =(value.at(nGoodTaus) == cut.at(nGoodTaus));

  value.at(nBJets)=nBJetsIndex.size();
  pass.at(nBJets) =(value.at(nBJets) > cut.at(nBJets));

  value.at(nJets)=nJetsIndex.size();
  pass.at(nJets) =(value.at(nJets) >=  cut.at(nJets));

  value.at(ET)=Ntp->MET();
  pass.at(ET) =(value.at(ET) >=  cut.at(ET));

  value.at(LeptonVeto) = thirdLeptonCounter.size();
  pass.at(LeptonVeto) = ( value.at(LeptonVeto) ==  value.at(LeptonVeto));


  //============================================================================
  // Here you can defined different type of weights you want to apply to events.
  double wobs=1;
  double w=1;

  if(id == 70){
    w*=Ntp->ttbarPtWeight();
  }


  if(!Ntp->isData() && id!=DataMCType::QCD){
    w *= reweight.PUweightHTT(Ntp->npu());
  }
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
  // std::cout<<" zptw  "<<  zptw <<" err  " << DataMC_Corr.ZPTWeightErr(genMomentum.M(),genMomentum.Pt()) <<std::endl;
  //============================================================================



  bool status=AnalysisCuts(t,w,wobs);  // boolean that say whether your event passed critera defined in pass vector. The whole vector must be true for status = true
  ///////////////////////////////////////////////////////////
  // Analyse events  which passed selection
  if(status){





  // if(id == 30 || id == 33){
  //   std::cout<<" first check event:    "<<std::endl;
  //   for(unsigned int imc=0; imc < Ntp->NMCParticles(); imc++){
  //     std::cout<<" pdgid   "<< Ntp->MCParticle_pdgid(imc)<< "  charge "  << Ntp->MCParticle_charge(imc) <<"  pT  " << Ntp->MCParticle_p4(imc).Px() <<std::endl;
  //   }
  // }





//   for(unsigned int iz =0; iz<Ntp->NMCSignalParticles(); iz++){
//     if(Ntp->MCSignalParticle_Tauidx->at(iz).size()!=0){
//       if(Ntp->MCTau_JAK->at(0) == jak){tauIndex=0; DecayOK = true;}
//       else if(  Ntp->MCTau_JAK->at(1) ==jak ){ tauIndex=1; DecayOK = true;}
//       if(DecayOK){
//         tau = TLorentzVector(Ntp->MCTauandProd_p4->at(Ntp->MCSignalParticle_Tauidx->at(iz).at(tauIndex)).at(0).at(1),Ntp->MCTauandProd_p4->at(Ntp->MCSignalParticle_Tauidx->at(iz).at(tauIndex)).at(0).at(2),
//                              Ntp->MCTauandProd_p4->at(Ntp->MCSignalParticle_Tauidx->at(iz).at(tauIndex)).at(0).at(3),Ntp->MCTauandProd_p4->at(Ntp->MCSignalParticle_Tauidx->at(iz).at(tauIndex)).at(0).at(0));

//       }
//     }
//   }
//   return tau;
// }


  



   
    NPrimeVtx.at(t).Fill(Ntp->npv(),w);
 
    MTM.at(t).Fill(sqrt(2*Ntp->Daughters_P4(goodMuonsIndex.at(0)).Pt()*Ntp->MET()*(1 - cos(Ntp->Daughters_P4(goodMuonsIndex.at(0)).Phi() - Ntp->METphi()))),w);
    LeadintTauPt.at(t).Fill(Ntp->Daughters_P4(goodTauIndex.at(0)).Pt(),w);
    LeadintMuonPt.at(t).Fill(Ntp->Daughters_P4(goodMuonsIndex.at(0)).Pt(),w);
    LeadintJetPt.at(t).Fill(Ntp->Jet_P4(nBJetsIndex.at(0)).Pt(),w);

    BDiscr1.at(t).Fill(Ntp->bDiscriminator(nBJetsIndex.at(0)),w);
    // if(nAk8JetsIndex.size()!=0)
    BDiscr2.at(t).Fill(Ntp->bCSVscore(nBJetsIndex.at(0)),w);
    BDiscr3.at(t).Fill(Ntp->pfCombinedMVAV2BJetTags(nBJetsIndex.at(0)),w);


    //  for(unsigned int alljets =0; alljets < nJetsIndex.size(); alljets++)
    AllJetsBTag.at(t).Fill(Ntp->bCSVscore(nJetsIndex.at(0)),w);
    // std::cout<<" selected leading jet :    "<< Ntp->bCSVscore(nJetsIndex.at(0)) << "  pt  " << Ntp->Jet_P4(nJetsIndex.at(0)).Pt() << " PF ID  " <<Ntp->PFjetID(nJetsIndex.at(0)) <<std::endl;
  
  }

}


//  This is a function if you want to do something after the event loop
void  TTBar::Finish(){
  Selection::Finish();
}






