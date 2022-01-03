#include "SimFit.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "SVFitObject.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
 
SimFit::SimFit(TString Name_, TString id_):
  Selection(Name_,id_)
{
		std::cout<<"****************************************************************************"<<std::endl;
		std::cout<<"*     event *         lumi *           run *        met *        event ID  *"<<std::endl;
		std::cout<<"****************************************************************************"<<std::endl;

}

SimFit::~SimFit(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }

  Logger(Logger::Info) << "complete." << std::endl;
}

void  SimFit::Configure(){
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
    PVSVSignificance=HConfig.GetTH1D(Name+"_PVSVSignificance"," PV-SV significance for tau decay mode = 10",21,0,4.5," pType","Events");
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
    TauMuPtResolution=HConfig.GetTH1D(Name+"_TauMuPtResolution","  pT resolution #tau_{a1}",50,-100,100," #Delta pT, GeV","Events");

    TauA1EResolution=HConfig.GetTH1D(Name+"_TauA1EResolution","  E resolution #tau_{a1}",50,-100,100," #Delta E, GeV","Events");
    TauMuEResolution=HConfig.GetTH1D(Name+"_TauMuEResolution","  E resolution #tau_{a1}",50,-100,100," #Delta E, GeV","Events");

    TauA1EtaResolution=HConfig.GetTH1D(Name+"_TauA1EtaResolution","  #eta resolution #tau_{a1}",50,-3,3," #Delta #eta, GeV","Events");
    TauMuEtaResolution=HConfig.GetTH1D(Name+"_TauMuEtaResolution","  #eta resolution #tau_{#mu}",50,-3,3," #Delta #eta, GeV","Events");

    TauA1PhiResolution=HConfig.GetTH1D(Name+"_TauA1PhiResolution","  #phi resolution #tau_{a1}",50,-3,3," #Delta #phi, GeV","Events");
    TauMuPhiResolution=HConfig.GetTH1D(Name+"_TauMuPhiResolution","  #phi resolution #tau_{#mu}",50,-3,3," #Delta #phi, GeV","Events");



    TrackPtResolution=HConfig.GetTH1D(Name+"_TrackPtResolution","  pT resolution #tau_{h}",50,-100,100," #Delta pT, GeV","Events");
    TrackEResolution=HConfig.GetTH1D(Name+"_TrackEResolution","  E resolution #tau_{h}",50,-100,100," #Delta E, GeV","Events");
    TrackPhiResolution=HConfig.GetTH1D(Name+"_TrackPhiResolution","  #phi resolution #tau_{h}",50,-3,3," #Delta #phi, GeV","Events");
    TrackEtaResolution=HConfig.GetTH1D(Name+"_TrackEtaResolution","  #eta resolution #tau_{h}",50,-3,3," #Delta #eta, GeV","Events");


    A1PtResolution=HConfig.GetTH1D(Name+"_A1PtResolution","  pT resolution #tau_{a1}",50,-100,100," #Delta pT, GeV","Events");
    A1EResolution=HConfig.GetTH1D(Name+"_A1EResolution","  E resolution #tau_{a1}",50,-100,100," #Delta E, GeV","Events");
    A1PhiResolution=HConfig.GetTH1D(Name+"_A1PhiResolution","  #phi resolution #tau_{a1}",50,-3,3," #Delta #phi, GeV","Events");
    A1EtaResolution=HConfig.GetTH1D(Name+"_A1EtaResolution","  #eta resolution #tau_{a1}",50,-3,3," #Delta #eta, GeV","Events");

    EventFitTauA1E=HConfig.GetTH1D(Name+"_EventFitTauA1E","  E_{#tau#rightarrow a1}",50,20,80," E, GeV","Events");
    EventFitTauMuE=HConfig.GetTH1D(Name+"_EventFitTauMuE","  E_{#tau#rightarrow #mu}",50,20,80," E, GeV","Events");


    TauA1PtResolutionPVSV=HConfig.GetTH1D(Name+"_TauA1PtResolutionPVSV","  pT resolution #tau_{a1}",50,-100,100," #Delta pT, GeV","Events");

    OSPionPtResolution=HConfig.GetTH1D(Name+"_OSPionPtResolution","OS Pion Pt resolution",50,-5,5," #DeltaP_{T}, GeV","Events");
    SSPionPtResolution=HConfig.GetTH1D(Name+"_SSPionPtResolution","SS Pion Pt resolution",50,-5,5," #DeltaP_{T}, GeV","Events");

    EventFitMass=HConfig.GetTH1D(Name+"_EventFitMass","",100,80,180," ","Events");
    EventFitZPt=HConfig.GetTH1D(Name+"_EventFitZPt","",50,0,100," ","Events");
    EventFitZE=HConfig.GetTH1D(Name+"_EventFitZE","",50,0,500," ","Events");
    EventFitZPhi=HConfig.GetTH1D(Name+"_EventFitZPhi","",50,-3.14,3.14," ","Events");
    EventFitZEta=HConfig.GetTH1D(Name+"_EventFitZEta","",50,-3,3," ","Events");
    EventFitZPtResolution=HConfig.GetTH1D(Name+"_EventFitZPtResolution","",50,-50,50," ","Events");
    EventFitZEResolution=HConfig.GetTH1D(Name+"_EventFitZEResolution","",50,-5,5," ","Events");
    EventFitZPhiResolution=HConfig.GetTH1D(Name+"_EventFitZPhiResolution","",50,-3,3," ","Events");
    EventFitZEtaResolution=HConfig.GetTH1D(Name+"_EventFitZEtaResolution","",50,-3,3," ","Events");


    A1VisiblePtResolution=HConfig.GetTH1D(Name+"_A1VisiblePtResolution","  pT resolution a1",50,-100,100," #Delta pT, GeV","Events");
    TrackVisiblePtResolution=HConfig.GetTH1D(Name+"_TrackVisiblePtResolution","  pT resolution track",50,-100,100," #Delta pT, GeV","Events");
    A1VisiblePtResolutionMuA1=HConfig.GetTH1D(Name+"_A1VisiblePtResolutionMuA1","  pT resolution track",50,-100,100," #Delta pT, GeV","Events");
    TrackMatchingQuality=HConfig.GetTH1D(Name+"_TrackMatchingQuality","  Track  Matching Quality",50,0,0.1," dR sum tracks","Events");
    MQVsVisA1Resolution=HConfig.GetTH2D(Name+"_MQVsVisA1Resolution","  Track  Matching Quality",50,0,1,50,-40,40," dR sum tracks","Events");
    Selection::ConfigureHistograms();   //   do not remove
    HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);  // do not remove
}

 

void  SimFit::Store_ExtraDist(){

  //every new histo should be addedd to Extradist1d vector, just push it back;
   Extradist1d.push_back(&TrackMatchingQuality);
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
   Extradist1d.push_back(&TauMuPtResolution);

   Extradist1d.push_back(&EventFitTauA1E);
   Extradist1d.push_back(&EventFitTauMuE);



   Extradist1d.push_back(&TauA1EResolution);
   Extradist1d.push_back(&TauMuEResolution);


   Extradist1d.push_back(&TauA1EtaResolution);
   Extradist1d.push_back(&TauMuEtaResolution);

   Extradist1d.push_back(&TauA1PhiResolution);
   Extradist1d.push_back(&TauMuPhiResolution);



   Extradist1d.push_back(&TauA1PtResolutionPVSV);

   Extradist1d.push_back(&s12);
   Extradist1d.push_back(&s13);
   Extradist1d.push_back(&s23);


   Extradist1d.push_back(&s12reco);
   Extradist1d.push_back(&s13reco);
   Extradist1d.push_back(&s23reco);
   Extradist1d.push_back(&OSPionPtResolution);
   Extradist1d.push_back(&SSPionPtResolution);
   Extradist1d.push_back(&EventFitMass);
   Extradist1d.push_back(&EventFitZPt);
   Extradist1d.push_back(&EventFitZE);
   Extradist1d.push_back(&EventFitZPhi);
   Extradist1d.push_back(&EventFitZEta);


   Extradist1d.push_back(&EventFitZPtResolution);
   Extradist1d.push_back(&EventFitZEResolution);
   Extradist1d.push_back(&EventFitZPhiResolution);
   Extradist1d.push_back(&EventFitZEtaResolution);




   Extradist1d.push_back(&TrackPtResolution);
   Extradist1d.push_back(&TrackEResolution);
   Extradist1d.push_back(&TrackPhiResolution);
   Extradist1d.push_back(&TrackEtaResolution);


   Extradist1d.push_back(&A1PtResolution);
   Extradist1d.push_back(&A1EResolution);
   Extradist1d.push_back(&A1PhiResolution);
   Extradist1d.push_back(&A1EtaResolution);


   Extradist1d.push_back(&A1VisiblePtResolution);
   Extradist1d.push_back(&TrackVisiblePtResolution);
   Extradist1d.push_back(&A1VisiblePtResolutionMuA1);
   Extradist2d.push_back(&MQVsVisA1Resolution);

}

void  SimFit::doEvent(){ //  Method called on every event
  unsigned int t;                // sample type, you may manage in your further analysis, if needed
  int id(Ntp->GetMCID());  //read event ID of a sample
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}  //  gives a warining if list of samples in Histo.txt  and SkimSummary.log do not coincide 

  //  std::cout<<" id  "<<id<<std::endl;
  value.at(PrimeVtx)=Ntp->NVtx();
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));
  
  bool PassedTrigger(false);
  int triggerindex;
  std::vector<int> TriggerIndex; 
  std::vector<int>TriggerIndexVector ;
  std::vector<TString>  MatchedTriggerNames;

  MatchedTriggerNames.push_back("HLT_IsoTkMu24_v");
  //  TriggerIndexVector=Ntp->GetVectorCrossTriggers("HLT_IsoMu","LooseIsoPFTau20_v");
  TriggerIndexVector=Ntp->GetVectorTriggers(MatchedTriggerNames);
   

  //int triggerMask(0);
  //triggerMask |= (1<< Ntp->getBitOfGivenTrigger("HLT_IsoMu22_eta2p1_v"));
  // PassedTrigger = (( Ntp->triggerbit() & triggerMask ) == triggerMask);
  
  for(int itrig = 0; itrig < TriggerIndexVector.size(); itrig++){
    if(Ntp->TriggerAccept(TriggerIndexVector.at(itrig))){  
      //  std::cout<<"  Name  "<< Ntp->TriggerName(TriggerIndexVector.at(itrig)) << "   status   "<< Ntp->TriggerAccept(TriggerIndexVector.at(itrig)) <<std::endl;
      PassedTrigger =Ntp->TriggerAccept(TriggerIndexVector.at(itrig)); }
  }

  value.at(TriggerOk)=true;//
  pass.at(TriggerOk)=true;//PassedTrigger;
  

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
    bool isTauTrack(false);
    bool isTau(false);
    int MuonCandidate;
    int TauTrackCandidate;
    int TauCandidate;
    for(unsigned int iDaugther=0;   iDaugther  <  Ntp->NDaughters() ;iDaugther++ ){  // loop over all daughters in the event

      DaughtersPt.at(t).Fill(Ntp->Daughters_P4(iDaugther).Pt(),w);  // Fill transverse momentum

      particletype.at(t).Fill(Ntp->particleType(iDaugther),w);

      if(Ntp->particleType(iDaugther) == 0){
	if(Ntp->Muon_TrackParticleHasMomentum(iDaugther)){	 MuonCandidate= iDaugther; isMuon=true;}
      }
    
      //   std::cout<<" isTauTrack   "<< Ntp->PFTau_TrackParticleHasMomentum(iDaugther) <<std::endl;
      if(Ntp->particleType(iDaugther)==2){

	//	std::cout<<" decMode  "<<Ntp->decayMode(iDaugther)<< "  isA1   "<< Ntp->PFTau_TIP_hasA1Momentum(iDaugther)<<"  is Track "<< Ntp->PFTau_TrackParticleHasMomentum(iDaugther) <<std::endl;

	if(Ntp->PFTau_TrackParticleHasMomentum(iDaugther)){	 TauTrackCandidate= iDaugther; isTauTrack=true;}

      }


    
      if(Ntp->PFTau_TIP_hasA1Momentum(iDaugther) && Ntp->isPVCovAvailable() && Ntp->PFTau_hassecondaryVertex(iDaugther) && Ntp->particleType(iDaugther) == 2) {TauCandidate = iDaugther; isTau=true;}

    

      //      std::cout<<"isTau   is Track "<< isTau<<"  " << isTauTrack<<std::endl;
      //  Ntp->deb(0); std::cout<<" id  "<< id <<std::endl;
      //  Ntp->deb(1157);
	taudecaytype.at(t).Fill(Ntp->decayMode(iDaugther),w);
	//	Ntp->deb(101);
	if(Ntp->PFTau_hassecondaryVertex(iDaugther) && Ntp->isPVCovAvailable()){
	  //  Ntp->deb(102);
	  if(Ntp->PFTau_secondaryVertex_TracksMatchingQuality(iDaugther) < 0.01)	  PVSVSignificance.at(t).Fill( Ntp->PFTau_FlightLength_significance(Ntp->PVtx(),Ntp->PFTau_TIP_primaryVertex_cov(), Ntp->PFTau_secondaryVertex_pos(iDaugther), Ntp->PFTau_TIP_secondaryVertex_cov(iDaugther)),w);
	  //  Ntp->deb(103);
	  SVchi2.at(t).Fill(Ntp->PFTau_secondaryVertex_vtxchi2(iDaugther),w);
	  // Ntp->deb(104);
	  SVMatchingQuality.at(t).Fill(Ntp->PFTau_secondaryVertex_TracksMatchingQuality(iDaugther),w);
	}
 	
	//	Ntp->deb(1.5);
    }
    
    if(isTau && isMuon){
    if(Ntp->decayMode(TauCandidate)==10  &&  Ntp->PFtauHasPions(TauCandidate) && Ntp->PFtauHasThreePions(TauCandidate) && id == 10230533  ){
      //	Ntp->deb(2);
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


	std::cout<<"reco pions   os,ss1,ss2  "<<std::endl;
	OSPion.Print();
	SSPion1.Print();
	SSPion2.Print();


	std::cout<<"mc pions   os,ss1,ss2  "<<std::endl;
	mc_ospion.Print();
	mc_ss1pion.Print();
	mc_ss2pion.Print();


	OSPionPtResolution.at(t).Fill(mc_ospion.Pt() - OSPion.Pt(),w);
	SSPionPtResolution.at(t).Fill( (mc_ss1pion + mc_ss2pion).Pt() - ( SSPion1+ SSPion2).Pt(),w);
	//std::cout<<"  pt Diff    "<< Ntp->Daughters_P4(TauCandidate).Pt() - (OSPion+ SSPion1+ SSPion2).Pt() <<std::endl;
	
	TLorentzVector TruthTauMu = Ntp->GetTruthTauLV(2);
	TLorentzVector TruthTauA1 = Ntp->GetTruthTauLV(5);
	TrackMatchingQuality.at(t).Fill(Ntp->PFTau_secondaryVertex_TracksMatchingQuality(TauCandidate),1);
	  double phiz = 0.01;
	  TLorentzVector TruthMuon =	Ntp-> GetTruthTauProductLV(2,13);
	  //	  std::cout<<"reco muon "<<Ntp->Daughters_P4(MuonCandidate).Pt()<< "truth    " <<TruthMuon.Pt() <<std::endl;
	  LorentzVectorParticle LVPTauA1 =Ntp->PFTau_a1_lvp(TauCandidate) ;
	  
	  if(Ntp->PFTau_secondaryVertex_TracksMatchingQuality(TauCandidate) < 0.01) A1VisiblePtResolutionMuA1.at(t).Fill(LVPTauA1.LV().Pt() - (mc_ospion+mc_ss1pion+mc_ss2pion).Pt(),1);       
	  MQVsVisA1Resolution.at(t).Fill(Ntp->PFTau_secondaryVertex_TracksMatchingQuality(TauCandidate),LVPTauA1.LV().Pt() - (mc_ospion+mc_ss1pion+mc_ss2pion).Pt());
	
	  TMatrixT<double> METpar(2,1); METpar(0,0) = Ntp->MET()*cos(Ntp->METphi()); METpar(1,0) = Ntp->MET()*sin(Ntp->METphi());
	  TMatrixTSym<double> METCov; METCov.ResizeTo(2,2);

	  METCov[0][0] = Ntp->PFMETCov00();
	  METCov[1][0] = Ntp->PFMETCov01();
	  METCov[0][1] = Ntp->PFMETCov10();
	  METCov[1][1] = Ntp->PFMETCov11();
	  
	  //	Ntp->deb(3);
      //      METCov = MET.significanceMatrix<TMatrixTSym<double> >();                                                                                                                                      
	  PTObject METObj(METpar, METCov);


	  //	  GlobalEventFit EF(Ntp->Muon_TrackParticle(MuonCandidate), LVPTauA1,MET2 , Ntp->PFTau_TIP_primaryVertex_pos(TauCandidate),Ntp->PFTau_TIP_primaryVertex_cov(TauCandidate) );  
	  //	  GlobalEventFit EF(Ntp->Muon_TrackParticle(MuonCandidate), LVPTauA1, phiz, Ntp->PVtx(),Ntp->PFTau_TIP_primaryVertex_cov());
	  GlobalEventFit EF(Ntp->Muon_TrackParticle(MuonCandidate), LVPTauA1, METObj , Ntp->PVtx(),Ntp->PFTau_TIP_primaryVertex_cov());

	  //  LVPTauA1.LVCov().Print();
	  //	  std::cout<<" vertex  cov "<< std::endl;
	  //	  LVPTauA1.VertexCov().Print();
	
	  GEFObject GEF = EF.Fit();
	  EF.SetCorrectPt(false);

	  
	  for(unsigned int iSigParticle =0; iSigParticle < Ntp->NMCSignalParticles() ; iSigParticle++ ){
	    // std::cout<<"pdgId   "<< Ntp-> MCSignalParticle_pdgid(iSigParticle)<< " Px   "<< Ntp->MCSignalParticle_p4(iSigParticle).Px() <<  " =  px2   " << (TruthTauMu+TruthTauA1).Px()  << std::endl;
	    // std::cout<< Ntp->MCSignalParticle_p4(iSigParticle).M() <<  " =  " << (TruthTauMu+TruthTauA1).M()  << std::endl;
	  }



	  double sign=  Ntp->PFTau_FlightLength_significance(Ntp->PVtx(),Ntp->PFTau_TIP_primaryVertex_cov(), Ntp->PFTau_secondaryVertex_pos(TauCandidate), Ntp->PFTau_TIP_secondaryVertex_cov(TauCandidate));

	 
	  if(GEF.isValid()){

	    std::cout<<"*     "<<Ntp->EventNumber()<< " *      " << Ntp->lumi() <<  " *     " << Ntp->RunNumber() << " *     " <<  Ntp->MET()<<" *      "<<id<<"  *"<<std::endl;
	    std::cout<<"_________________________________________________________________________________"<<std::endl;

	    //  std::cout<<" EventNumber:   "<<Ntp->EventNumber()<< " lumi   " << Ntp->lumi() <<  "    run Number  " << Ntp->RunNumber() << "        event id   " <<  id  <<std::endl;
	    //  if (Ntp->EventNumber()==  811625 || Ntp->EventNumber()==    991805 || Ntp->EventNumber() ==  2017040   || Ntp->EventNumber()==  3359188  || Ntp->EventNumber()==  3448786 || Ntp->EventNumber()==   6608279 || Ntp->EventNumber()==   4161423 || Ntp->EventNumber()==   6597584  || Ntp->EventNumber()==   882233  || Ntp->EventNumber()==  2007199  || Ntp->EventNumber()== 66274705   || Ntp->EventNumber()== 6770939   || Ntp->EventNumber()== 6862074   || Ntp->EventNumber()== 6871581  || Ntp->EventNumber()==  9282844 || Ntp->EventNumber()==   3508828  || Ntp->EventNumber()==  3509759  || Ntp->EventNumber()==  3510715  || Ntp->EventNumber()==  4084437  || Ntp->EventNumber()==  6760343)
{

	      
	      // std::cout<<" EventNumber:   "<<Ntp->EventNumber()<< "        event id   " <<  id <<std::endl;

 
	      std::cout<<" Primary vertex:  "<<std::endl;Ntp->PVtx().Print();   std::cout<<std::endl;
	      std::cout<<" Primary Vertex Covariance:  "<<std::endl;Ntp->PFTau_TIP_primaryVertex_cov().Print();std::cout<<std::endl; std::cout<<std::endl; 
	      std::cout<<"--------"<<std::endl;
	      std::cout<<" SecondaryVertex of 3 pions:  "<<std::endl;LVPTauA1.Vertex().Print();std::cout<<std::endl; 
	      std::cout<<" Secondary Vertex  Covariance:  "<<std::endl;LVPTauA1.VertexCov().Print(); std::cout<<std::endl; std::cout<<std::endl; 
	      std::cout<<"--------"<<std::endl;

	      std::cout<<" Muon Lorentz Vector  :  "<<std::endl;Ntp->Daughters_P4(MuonCandidate).Print();  std::cout<<std::endl; 
 	      std::cout<<" Muon Parameters Helix :  "<<std::endl;  Ntp->Muon_TrackParticle(MuonCandidate).getParMatrix().Print();std::cout<<std::endl; 
	      std::cout<<" Muon Helix Paramters Covariance:  "<<std::endl;  Ntp->Muon_TrackParticle(MuonCandidate).getCovMatrix().Print();std::cout<<std::endl; std::cout<<std::endl; 
	      std::cout<<"--------"<<std::endl;
	      std::cout<<" A1 LorentzVecotr:  "<<std::endl;LVPTauA1.LV().Print();std::cout<<std::endl; 
	      std::cout<<" A1 LorentzVector  Covariance:  "<<std::endl;LVPTauA1.LVCov().Print();std::cout<<std::endl; std::cout<<std::endl; 
	      std::cout<<"_________________________________________________________________________________"<<std::endl;
	     }

	    TLorentzVector EventFitTauA1 =GEF.getTauH().LV();
	    TLorentzVector EventFitTauMu =GEF.getTauMu().LV();
	    TLorentzVector EventFitZ = EventFitTauMu + EventFitTauA1;
	    TLorentzVector TruthZ = TruthTauMu+TruthTauA1;
	    //	    std::cout<<" SF   "<< "  taua1 px "<<EventFitTauA1.Px()<< " + taumu  "<<EventFitTauMu.Px() << "  =  "<< GEF.getResonance().LV().Px()<<std::endl;
	    EventFitMass.at(t).Fill(EventFitZ.M(),1);// std::cout<<"  EF Mass =   "<< EventFitZ.M() << "   Res  " << GEF.getResonance().LV().M()<< std::endl;
	    EventFitZPt.at(t).Fill(EventFitZ.Pt(),1);
	    EventFitZE.at(t).Fill(EventFitZ.E(),1);
	    EventFitZPhi.at(t).Fill(EventFitZ.Phi(),1);
	    EventFitZEta.at(t).Fill(EventFitZ.Eta(),1);

	    EventFitZPtResolution.at(t).Fill(EventFitZ.Pt()-TruthZ.Pt(),1);
	    EventFitZEResolution.at(t).Fill((EventFitZ.E()-TruthZ.E())/TruthZ.E(),1);
	    EventFitZPhiResolution.at(t).Fill(EventFitZ.Phi()-TruthZ.Phi(),1);
	    EventFitZEtaResolution.at(t).Fill(EventFitZ.Eta()-TruthZ.Eta(),1);
	    TauA1PtResolution.at(t).Fill(EventFitTauA1.Pt() - TruthTauA1.Pt(),1);
	    TauMuPtResolution.at(t).Fill(EventFitTauMu.Pt() - TruthTauMu.Pt(),1);

	    TauA1EResolution.at(t).Fill(EventFitTauA1.E() - TruthTauA1.E(),1);
	    TauMuEResolution.at(t).Fill(EventFitTauMu.E() - TruthTauMu.E(),1);

	    TauA1EtaResolution.at(t).Fill(EventFitTauA1.Eta() - TruthTauA1.Eta(),1);
	    TauMuEtaResolution.at(t).Fill(EventFitTauA1.Eta() - TruthTauA1.Eta(),1);

	    TauA1PhiResolution.at(t).Fill(EventFitTauA1.Phi() - TruthTauA1.Phi(),1);
	    TauMuPhiResolution.at(t).Fill(EventFitTauA1.Phi() - TruthTauA1.Phi(),1);

	    EventFitTauA1E.at(t).Fill(EventFitTauA1.E(),1);
	    EventFitTauMuE.at(t).Fill(EventFitTauMu.E(),1);



	    if(GEF.isValid() && sign >3 )  TauA1PtResolutionPVSV.at(t).Fill(EventFitTauA1.Pt() - TruthTauA1.Pt(),1);
	   
	  
	    // std::cout<<"   Covariance(3,3) "<< LVPTauA1.Covariance(3,3) <<std::endl;
	    // std::cout<<"   Covariance(4,4) "<< LVPTauA1.Covariance(4,4) <<std::endl;
	    // std::cout<<"   Covariance(5,5) "<< LVPTauA1.Covariance(5,5) <<std::endl;
	    //	    std::cout<<"   GEF Status  "<< GEF.isValid() <<EventFitTauA1.Pt()  << " = " <<TruthTauA1.Pt() <<std::endl;



      
      }
	
      }

    



	
    }
    
    }

    if(isTau && isTauTrack){
      //      std::cout<<"  passed ? "<< id<< std::endl;
      //	Ntp->deb(4);
      if(Ntp->CheckDecayID(3,5) && id == 10330533 || id == 10430533 && Ntp->PFtauHasPions(TauCandidate) && Ntp->PFtauHasThreePions(TauCandidate)){

	TLorentzVector TruthTauTrack(0,0,0,0);
	TLorentzVector TruthTauA1 = Ntp->GetTruthTauLV(5);

	TLorentzVector  TruthPi(0,0,0,0);
	TLorentzVector  TruthPiRho(0,0,0,0);
	if(id ==10330533){ TruthPi= Ntp->GetTruthTauProductLV(3,211);TruthTauTrack = Ntp->GetTruthTauLV(3); }
	if(id ==10430533){ TruthPiRho= Ntp->GetTruthTauProductLV(4,211);TruthTauTrack = Ntp->GetTruthTauLV(4); }



	  if(Ntp->isPFTauTrackAvailable(TauTrackCandidate)){
	    //	    std::cout<<"  id  "<< id<< "   truth pt   "<< TruthTauTrack.Pt() << "   reco pt "<< Ntp->PFTauLeadTrackLV(TauTrackCandidate).Pt() <<  "  dR  " << Ntp->PFTauTrack_deltaR(TauTrackCandidate) <<std::endl;
	    //	    std::cout<<" truth pi  "<< TruthPi.Pt() << "  truth piro "<< TruthPiRho.Pt()<< std::endl;
	  }
  

	LorentzVectorParticle LVPTauA1 =Ntp->PFTau_a1_lvp(TauCandidate) ;
	 
	TMatrixT<double> METpar(2,1); METpar(0,0) = Ntp->MET()*cos(Ntp->METphi()); METpar(1,0) = Ntp->MET()*sin(Ntp->METphi());
	TMatrixTSym<double> METCov; METCov.ResizeTo(2,2);

	METCov[0][0] = Ntp->PFMETCov00();
	METCov[1][0] = Ntp->PFMETCov01();
	METCov[0][1] = Ntp->PFMETCov10();
	METCov[1][1] = Ntp->PFMETCov11();
	  
	TLorentzVector TruthA1LV=	Ntp->GetTruthPionsFromA1().at(0) + Ntp->GetTruthPionsFromA1().at(1) + Ntp->GetTruthPionsFromA1().at(2);
	if(Ntp->PFTau_secondaryVertex_TracksMatchingQuality(TauCandidate) < 0.01)	A1VisiblePtResolution.at(t).Fill(LVPTauA1.LV().Pt() - TruthA1LV.Pt(),1);       
	TrackVisiblePtResolution.at(t).Fill( Ntp->PFTauLeadTrackLV(TauTrackCandidate).Pt() - TruthTauTrack.Pt(),1);                                                            
	PTObject METObj(METpar, METCov);


	GlobalEventFit EF(Ntp->PFTau_TrackParticle(TauTrackCandidate), LVPTauA1,METObj , Ntp->PVtx(),Ntp->PFTau_TIP_primaryVertex_cov());
	//	Ntp->PFTau_TrackParticle(TauTrackCandidate).getParMatrix().Print();
	GEFObject GEF = EF.Fit();
	TLorentzVector EventFitTauA1 =GEF.getTauH().LV();
	TLorentzVector EventFitTauPi =GEF.getTauMu().LV();
	
	//	std::cout<<"   GEF Status   Track "<< GEF.isValid() << "  " <<EventFitTauA1.Pt()  << " = " <<TruthTauA1.Pt() <<std::endl;
	if(GEF.isValid()  && Ntp->PFTau_secondaryVertex_TracksMatchingQuality(TauCandidate) < 0.01 && TruthTauTrack.Px()!=0 ){ 

	  //	  std::cout<<" Mass  "<< (EventFitTauPi+EventFitTauA1 ).M() << "   id  "<< id <<std::endl;
	  TrackPtResolution.at(t).Fill(EventFitTauPi.Pt() - TruthTauTrack.Pt(),1);
	  TrackEResolution.at(t).Fill(EventFitTauPi.E() - TruthTauTrack.E(),1);
	  TrackPhiResolution.at(t).Fill(EventFitTauPi.Phi() - TruthTauTrack.Phi(),1);
	  TrackEtaResolution.at(t).Fill(EventFitTauPi.Eta() - TruthTauTrack.Eta(),1);
	  
	  
	  A1PtResolution.at(t).Fill(EventFitTauA1.Pt() - TruthTauA1.Pt(),1);
	  A1EResolution.at(t).Fill(EventFitTauA1.E() - TruthTauA1.E(),1);
	  A1PhiResolution.at(t).Fill(EventFitTauA1.Phi() - TruthTauA1.Phi(),1);
	  A1EtaResolution.at(t).Fill(EventFitTauA1.Eta() - TruthTauA1.Eta(),1);
	  
	}
      }
    }
    

    
  }
  
    
  
  

  

}





//  This is a function if you want to do something after the event loop
void  SimFit::Finish(){
  Selection::Finish();
}





