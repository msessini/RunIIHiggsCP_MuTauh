#include "NtupleValidation.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "SVFitObject.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
 
NtupleValidation::NtupleValidation(TString Name_, TString id_):
  Selection(Name_,id_)
{
}

NtupleValidation::~NtupleValidation(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  NtupleValidation::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)    cut.at(TriggerOk)=1;
    if(i==PrimeVtx)     cut.at(PrimeVtx)=1;
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
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PrimeVtx_",htitle,11,-0.5,10.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PrimeVtx_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i==TriggerOk){
      title.at(i)="Trigger ";
      hlabel="Trigger ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }


  } 
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");


  // Setup Extra Histograms
    NumVertices=HConfig.GetTH1D(Name+"_NumVertices","Number of Primary Vertices",36,-0.05,35.5,"N","Events");
    MissingTEnergy=HConfig.GetTH1D(Name+"_MissingTEnergy","Missing Transverse Energy",50,-0.05,99.5," M_{T}, GeV","Events");
    DaughtersPt=HConfig.GetTH1D(Name+"_DaughtersPt","Transverse momentum of reco particles",50,-0.05,99.5," P_{T}, GeV","Events");
    particletype=HConfig.GetTH1D(Name+"_particletype","Particle Type:0-muon,1-electron,2-tau",4,-0.5,3.5," pType","Events");
    taudecaytype=HConfig.GetTH1D(Name+"_taudecaytype","HPS decay type",11,-0.5,10.5," pType","Events");
    PVSVSignificance=HConfig.GetTH1D(Name+"_PVSVSignificance"," PV-SV significance for tau decay mode = 10",51,-0.5,10.5," pType","Events");
    SVchi2=HConfig.GetTH1D(Name+"_SVchi2"," #chi^{2} of SV",50,0,10,"#chi^{2} of SV","Events");
    SVMatchingQuality=HConfig.GetTH1D(Name+"_SVMatchingQuality"," Matching qualtiy of tracks used to build PV and general tracks",50,0,0.01," Sum #Delta R","Events");

    isPairCandOS=HConfig.GetTH1D(Name+"_isPairCandOS","is pair OS ",2,-0.5,1.5,"OS/SS","Events");
    Pair_part1Type=HConfig.GetTH1D(Name+"_Pair_part1Type","Particle Type:0-muon,1-electron,2-tau",4,-0.5,3.5," pType","Events");
    Pair_part2Type=HConfig.GetTH1D(Name+"_Pair_part2Type","Particle Type:0-muon,1-electron,2-tau",4,-0.5,3.5," pType","Events");

    OSPairMass=HConfig.GetTH1D(Name+"_OSPairMass"," Mass of two OS candidates",50,50,100," M_{pair}, GeV","Events");
    SSPairMass=HConfig.GetTH1D(Name+"_SSPairMass"," Mass of two SS candidates",50,50,100," M_{pair}, GeV","Events");
    TruthTauTauMass=HConfig.GetTH1D(Name+"_TruthTauTauMass"," Truth TauTau Mass",50,50,100," M_{pair}, GeV","Events");
    s12=HConfig.GetTH1D(Name+"_s12"," Truth dalitz s12",50,0,1.5," s_{12}, GeV","Events");
    s13=HConfig.GetTH1D(Name+"_s13","Truth dalitz s13 ",50,0,1.5,"s_{13} , GeV","Events");
    s23=HConfig.GetTH1D(Name+"_s23"," Truth dalitz s23",50,0,1.5," s_{23}, GeV","Events");

    s12reco=HConfig.GetTH1D(Name+"_s12reco","reco dalitz s12",50,0,1.5," s_{12}, GeV","Events");
    s13reco=HConfig.GetTH1D(Name+"_s13reco","reco dalitz s13",50,0,1.5," s_{13}, GeV","Events");
    s23reco=HConfig.GetTH1D(Name+"_s23reco","reco dalitz s23",50,0,1.5," s_{23}, GeV","Events");
    TauA1PtResolution=HConfig.GetTH1D(Name+"_TauA1PtResolution","  pT resolution #tau_{a1}",50,-100,100," #Delta pT, GeV","Events");
    TauA1PtResolutionPVSV=HConfig.GetTH1D(Name+"_TauA1PtResolutionPVSV","  pT resolution #tau_{a1}",50,-100,100," #Delta pT, GeV","Events");

    OSPionPtResolution=HConfig.GetTH1D(Name+"_OSPionPtResolution","OS Pion Pt resolution",50,-5,5," #DeltaP_{T}, GeV","Events");
    SSPionPtResolution=HConfig.GetTH1D(Name+"_SSPionPtResolution","SS Pion Pt resolution",50,-5,5," #DeltaP_{T}, GeV","Events");
    Selection::ConfigureHistograms();   //   do not remove
    HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);  // do not remove
}

 

void  NtupleValidation::Store_ExtraDist(){

  //every new histo should be addedd to Extradist1d vector, just push it back;
   Extradist1d.push_back(&NumVertices);
   Extradist1d.push_back(&MissingTEnergy);
   Extradist1d.push_back(&DaughtersPt);
   Extradist1d.push_back(&particletype);
   Extradist1d.push_back(&taudecaytype);
   Extradist1d.push_back(&PVSVSignificance);
   Extradist1d.push_back(&SVchi2);
   Extradist1d.push_back(&SVMatchingQuality);
   Extradist1d.push_back(&isPairCandOS);
   Extradist1d.push_back(&Pair_part1Type);
   Extradist1d.push_back(&Pair_part2Type);
   Extradist1d.push_back(&OSPairMass);
   Extradist1d.push_back(&SSPairMass);
   Extradist1d.push_back(&TruthTauTauMass);
   Extradist1d.push_back(&TauA1PtResolution);
   Extradist1d.push_back(&TauA1PtResolutionPVSV);

   Extradist1d.push_back(&s12);
   Extradist1d.push_back(&s13);
   Extradist1d.push_back(&s23);


   Extradist1d.push_back(&s12reco);
   Extradist1d.push_back(&s13reco);
   Extradist1d.push_back(&s23reco);
   Extradist1d.push_back(&OSPionPtResolution);
   Extradist1d.push_back(&SSPionPtResolution);


}

void  NtupleValidation::doEvent(){ //  Method called on every event
  unsigned int t;                // sample type, you may manage in your further analysis, if needed
  int id(Ntp->GetMCID());  //read event ID of a sample
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}  //  gives a warining if list of samples in Histo.txt  and SkimSummary.log do not coincide 
  //std::cout<<"------------------ New Event -----------------------"<<std::endl;
  bool PassedTrigger(false);
  int triggerindex;

  std::vector<int>TriggerIndexVector ;
  std::vector<TString>  MatchedTriggerNames;

  MatchedTriggerNames.push_back("LooseIsoPFTau20");
  //  TriggerIndexVector=Ntp->GetVectorCrossTriggers("HLT_IsoMu","LooseIsoPFTau20_v");
  TriggerIndexVector=Ntp->GetVectorTriggers(MatchedTriggerNames);
  
  
  // int triggerMask(0);
  // triggerMask |= (1<< Ntp->getBitOfGivenTrigger("HLT_IsoMu22_eta2p1_v"));
  // PassedTrigger = (( Ntp->triggerbit() & triggerMask ) == triggerMask);
  
  for(int itrig = 0; itrig < TriggerIndexVector.size(); itrig++){
    if(Ntp->TriggerAccept(TriggerIndexVector.at(itrig))){  
      //  std::cout<<"  Name  "<< Ntp->TriggerName(TriggerIndexVector.at(itrig)) << "   status   "<< Ntp->TriggerAccept(TriggerIndexVector.at(itrig)) <<std::endl;
      PassedTrigger =Ntp->TriggerAccept(TriggerIndexVector.at(itrig)); }
  }
  



  // Apply Selection
  // two vectors value and pass  are used to apply selection, store value.at(<A cut variable defined in .h as enumerator>) a quantitity you want to use in your selection.
  // Define in the pass.at(<A cut variable defined in .h as enumerator>) a passing condidtion, true/false
  value.at(PrimeVtx)=Ntp->NVtx();
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));
  
  value.at(TriggerOk)=PassedTrigger;
  pass.at(TriggerOk)=PassedTrigger;
  

  // Here you can defined different type of weights you want to apply to events. At the moment only PU weight is considered if event is not data
  double wobs=1;
  double w;
  if(!Ntp->isData()){w = Ntp->PUReweight();}
  else{w=1;}



  bool status=AnalysisCuts(t,w,wobs);  // boolean that say whether your event passed critera defined in pass vector. The whole vector must be true for status = true
  ///////////////////////////////////////////////////////////
  // Analyse events  which passed selection
  if(status){
    
    //    std::cout<<"  CheckDecayID   "  << Ntp->CheckDecayID(2,5) <<std::endl;//Ntp->printMCDecayChainOfEvent(true,true,true,true);
    MissingTEnergy.at(t).Fill(Ntp->MET(),w);
    NumVertices.at(t).Fill(Ntp->NVtx(),w);


    if(id == 10230533){
      if(Ntp->CheckDecayID(2,5)){
	TLorentzVector TruthTauMu = Ntp->GetTruthTauLV(2);
	TLorentzVector TruthTauA1 = Ntp->GetTruthTauLV(5);
	
      TruthTauTauMass.at(t).Fill((TruthTauMu+TruthTauA1).M(),w);
      TLorentzVector mc_ospion =Ntp->GetTruthPionsFromA1().at(0);
      TLorentzVector mc_ss1pion=Ntp->GetTruthPionsFromA1().at(1);
      TLorentzVector mc_ss2pion=Ntp->GetTruthPionsFromA1().at(2);
      s12.at(t).Fill( (mc_ospion + mc_ss1pion).M(),w);
      s13.at(t).Fill( (mc_ospion + mc_ss2pion).M(),w);
      s23.at(t).Fill( (mc_ss1pion + mc_ss2pion).M(),w);
      }
    }

    bool isMuon(false);
    bool isTau(false);
      int MuonCandidate;
      int TauCandidate;
    for(unsigned int iDaugther=0;   iDaugther  <  Ntp->NDaughters() ;iDaugther++ ){  // loop over all daughters in the event

      DaughtersPt.at(t).Fill(Ntp->Daughters_P4(iDaugther).Pt(),w);  // Fill transverse momentum

      particletype.at(t).Fill(Ntp->particleType(iDaugther),w);

      if(Ntp->particleType(iDaugther) == 0){
	if(Ntp->Muon_TrackParticleHasMomentum(iDaugther)){	 MuonCandidate= iDaugther; isMuon=true;}
      }


    
      if(Ntp->PFTau_TIP_hasA1Momentum(iDaugther) && Ntp->isPVCovAvailable() && Ntp->PFTau_hassecondaryVertex(iDaugther) && Ntp->particleType(iDaugther) == 2) {TauCandidate = iDaugther; isTau=true;}

	taudecaytype.at(t).Fill(Ntp->decayMode(iDaugther),w);
    
	if(Ntp->PFTau_hassecondaryVertex(iDaugther) && Ntp->isPVCovAvailable()){
	  
	  PVSVSignificance.at(t).Fill( Ntp->PFTau_FlightLength_significance(Ntp->PVtx(),Ntp->PFTau_TIP_primaryVertex_cov(), Ntp->PFTau_secondaryVertex_pos(iDaugther), Ntp->PFTau_TIP_secondaryVertex_cov(iDaugther)),w);
	  SVchi2.at(t).Fill(Ntp->PFTau_secondaryVertex_vtxchi2(iDaugther),w);
	  SVMatchingQuality.at(t).Fill(Ntp->PFTau_secondaryVertex_TracksMatchingQuality(iDaugther),w);
	}
    }

    if(isTau&& isMuon){
    if(Ntp->decayMode(TauCandidate)==10 && Ntp->PFtauHasPions(TauCandidate) && Ntp->PFtauHasThreePions(TauCandidate) && id == 10230533 ){
      TLorentzVector SSPion1,SSPion2, OSPion;
      unsigned int ssindex2(0);
      unsigned int OSPionIndex,SSPion1Index,SSPion2Index;
      


      if(Ntp->PFTau_PionsCharge(TauCandidate,0)*Ntp->PFTau_PionsCharge(TauCandidate,1) == 1 && Ntp->PFTau_PionsCharge(TauCandidate,0)*Ntp->PFTau_PionsCharge(TauCandidate, 2) == -1){
	SSPion1 = Ntp->PFTau_PionsP4(TauCandidate,0);
	SSPion2 = Ntp->PFTau_PionsP4(TauCandidate,1);
	OSPion = Ntp->PFTau_PionsP4(TauCandidate,2);
	OSPionIndex  = 2;
	SSPion1Index = 0;
	SSPion2Index = 1;
      }
      if(Ntp->PFTau_PionsCharge(TauCandidate,0)*Ntp->PFTau_PionsCharge(TauCandidate,1) == -1 && Ntp->PFTau_PionsCharge(TauCandidate,0)*Ntp->PFTau_PionsCharge(TauCandidate,2) == 1){
	SSPion1 = Ntp->PFTau_PionsP4(TauCandidate,0);
	SSPion2 = Ntp->PFTau_PionsP4(TauCandidate,2);
	OSPion = Ntp->PFTau_PionsP4(TauCandidate,1);
	OSPionIndex  = 1;
	SSPion1Index = 0;
	SSPion2Index = 2;
	
      }
      if(Ntp->PFTau_PionsCharge(TauCandidate,0)*Ntp->PFTau_PionsCharge(TauCandidate,1) == -1 && Ntp->PFTau_PionsCharge(TauCandidate,0)*Ntp->PFTau_PionsCharge(TauCandidate,2) == -1){
	SSPion1 = Ntp->PFTau_PionsP4(TauCandidate,1);
	SSPion2 = Ntp->PFTau_PionsP4(TauCandidate,2);
	OSPion = Ntp->PFTau_PionsP4(TauCandidate,0);
	OSPionIndex  = 0;
	SSPion1Index = 1;
	SSPion2Index = 2;
      }
      s12reco.at(t).Fill(  (OSPion+ SSPion1).M() ,w);
      s13reco.at(t).Fill(  (OSPion+ SSPion2).M() ,w);
      s23reco.at(t).Fill(  (SSPion1+ SSPion2).M(),w);
      
     
      if(Ntp->CheckDecayID(2,5) && id == 10230533){
	TLorentzVector mc_ospion =Ntp->GetTruthPionsFromA1().at(0);
	TLorentzVector mc_ss1pion=Ntp->GetTruthPionsFromA1().at(1);
	TLorentzVector mc_ss2pion=Ntp->GetTruthPionsFromA1().at(2);
	OSPionPtResolution.at(t).Fill(mc_ospion.Pt() - OSPion.Pt(),w);
	SSPionPtResolution.at(t).Fill( (mc_ss1pion + mc_ss2pion).Pt() - ( SSPion1+ SSPion2).Pt(),w);
	//std::cout<<"  pt Diff    "<< Ntp->Daughters_P4(TauCandidate).Pt() - (OSPion+ SSPion1+ SSPion2).Pt() <<std::endl;
	
	TLorentzVector TruthTauMu = Ntp->GetTruthTauLV(2);
	TLorentzVector TruthTauA1 = Ntp->GetTruthTauLV(5);

	 
	  double phiz = 0.01;
	  
	  LorentzVectorParticle LVPTauA1 =Ntp->PFTau_a1_lvp(TauCandidate) ;
	 
	  TMatrixT<double> METpar(2,1); METpar(0,0) = Ntp->MET()*cos(Ntp->METphi()); METpar(1,0) = Ntp->MET()*sin(Ntp->METphi());
	  TMatrixTSym<double> METCov; METCov.ResizeTo(2,2);

	  METCov[0][0] = Ntp->PFMETCov00();
	  METCov[1][0] = Ntp->PFMETCov01();
	  METCov[0][1] = Ntp->PFMETCov10();
	  METCov[1][1] = Ntp->PFMETCov11();
	  

      //      METCov = MET.significanceMatrix<TMatrixTSym<double> >();                                                                                                                                      
	  PTObject METObj(METpar, METCov);


	  //	  GlobalEventFit EF(Ntp->Muon_TrackParticle(MuonCandidate), LVPTauA1,MET2 , Ntp->PFTau_TIP_primaryVertex_pos(TauCandidate),Ntp->PFTau_TIP_primaryVertex_cov(TauCandidate) );  
	  //	  GlobalEventFit EF(Ntp->Muon_TrackParticle(MuonCandidate), LVPTauA1, phiz, Ntp->PVtx(),Ntp->PFTau_TIP_primaryVertex_cov());
	  GlobalEventFit EF(Ntp->Muon_TrackParticle(MuonCandidate), LVPTauA1,METObj , Ntp->PVtx(),Ntp->PFTau_TIP_primaryVertex_cov());
	  LVPTauA1.LVCov().Print();
	  std::cout<<" vertex  cov "<< std::endl;
	  LVPTauA1.VertexCov().Print();
	
	  GEFObject GEF = EF.Fit();
	  
	  TLorentzVector EventFitTauA1 =GEF.getTauH().LV();
	  TLorentzVector EventFitTauMu =GEF.getTauMu().LV();
	  double sign=  Ntp->PFTau_FlightLength_significance(Ntp->PVtx(),Ntp->PFTau_TIP_primaryVertex_cov(), Ntp->PFTau_secondaryVertex_pos(TauCandidate), Ntp->PFTau_TIP_secondaryVertex_cov(TauCandidate));


	  if(GEF.isValid() && sign >3 )  TauA1PtResolutionPVSV.at(t).Fill(EventFitTauA1.Pt() - TruthTauA1.Pt(),1);
	  if(GEF.isValid()  ){TauA1PtResolution.at(t).Fill(EventFitTauA1.Pt() - TruthTauA1.Pt(),1);
	  
	  // std::cout<<"   Covariance(3,3) "<< LVPTauA1.Covariance(3,3) <<std::endl;
	  // std::cout<<"   Covariance(4,4) "<< LVPTauA1.Covariance(4,4) <<std::endl;
	  // std::cout<<"   Covariance(5,5) "<< LVPTauA1.Covariance(5,5) <<std::endl;
	    std::cout<<"   GEF Status  "<< GEF.isValid() <<EventFitTauA1.Pt()  << " = " <<TruthTauA1.Pt() <<std::endl;
      
      }
	
      }
	
    }
    
    }
  
  }

    
    


    for(unsigned int ipair=0; ipair < Ntp->NPairCandidates(); ipair++){
      isPairCandOS.at(t).Fill(Ntp->isOSCand(ipair),w);
      Pair_part1Type.at(t).Fill(Ntp->particleType(Ntp->indexDau1(ipair)),w);
      Pair_part2Type.at(t).Fill(Ntp->particleType(Ntp->indexDau2(ipair)),w);
      if(Ntp->isOSCand(ipair)==1)      OSPairMass.at(t).Fill((  (Ntp->Daughters_P4(Ntp->indexDau1(ipair)) +  Ntp->Daughters_P4(Ntp->indexDau2(ipair) ) ).M()  ),w);
      if(Ntp->isOSCand(ipair)==0)      SSPairMass.at(t).Fill((  (Ntp->Daughters_P4(Ntp->indexDau1(ipair)) +  Ntp->Daughters_P4( Ntp->indexDau2(ipair) ) ).M()  ),w);



    }


}





//  This is a function if you want to do something after the event loop
void  NtupleValidation::Finish(){
  Selection::Finish();
}





