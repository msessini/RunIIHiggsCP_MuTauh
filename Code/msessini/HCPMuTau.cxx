#include "HCPMuTau.h"
#include "TLorentzVector.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "SVFitObject.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
//#include "SVfitProvider.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "PDG_Var.h"
#include "SkimConfig.h"
#include "TauSpinerInterface.h"


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
#include "TauAnalysis/ClassicSVfit/interface/FastMTT.h"
#include "TauPolSoftware/TauDecaysInterface/interface/fonction_a1.h"
#include "TauPolSoftware/TauDecaysInterface/interface/SCalculator.h"
#include <algorithm>

HCPMuTau::HCPTauTau(TString Name_, TString id_, char* Channel_, char* CPstate_):
  Selection(Name_,id_)
  //DataMC_Corr(true,true,false),
  //tauTrgSF("tight")
{
  Channel = Channel_;
  CPstate = CPstate_; 
  ChargeSumDummy = -999;
  selMuon_IsoDummy = 999.;
  WorkSpaceFF2016=TFile::Open(((std::string)std::getenv("workdir")+"Code/fake_factors_tt_dRcorr/fakefactors_ws_tt_lite_2016_dR_corr.root").c_str(), "READ");
  wFF2016= (RooWorkspace*)gDirectory->Get("w");
  WorkSpaceFF2016->Close();
  WorkSpaceFF2017=TFile::Open(((std::string)std::getenv("workdir")+"Code/fake_factors_tt_dRcorr/fakefactors_ws_tt_lite_2017_dR_corr.root").c_str(), "READ");
  wFF2017= (RooWorkspace*)gDirectory->Get("w");
  WorkSpaceFF2017->Close();
  WorkSpaceFF2018=TFile::Open(((std::string)std::getenv("workdir")+"Code/fake_factors_tt_dRcorr/fakefactors_ws_tt_lite_2018_dR_corr.root").c_str(), "READ");
  wFF2018= (RooWorkspace*)gDirectory->Get("w");
  WorkSpaceFF2018->Close();
  BDT=new BDTClassification();
  BDT->PreAnalysis();
}

HCPMuTau::~HCPTauTau(){
  for(unsigned int j=0; j<Npassed.size(); j++){
    Logger(Logger::Info) << "Selection Summary before: "
			 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
			 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
  delete wFF2016;
  delete wFF2017;
  delete wFF2018;
}

void  HCPMuTau::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==Trigger)             cut.at(Trigger)=1;
    if(i==Id_and_Kin)            cut.at(Id_and_Kin)=1;
    //if(i==NPairsFound)         cut.at(NPairsFound)=1;
    //if(i==GoodIndex)           cut.at(GoodIndex)=1.;
    //if(i==ZTTMC)                 cut.at(ZTTMC)=1.;
    //if(i==METFilters)            cut.at(METFilters)=1.;
    //if(i==genmatch)              cut.at(genmatch)=1;
    if(i==TausIsolation)         cut.at(TausIsolation)=1;
    if(i==AgainstEleMu)          cut.at(AgainstEleMu)=1;
    //if(i==Tau2Isolation)       cut.at(Tau2Isolation)=1.;
    if(i==LeptonVeto)            cut.at(LeptonVeto)=0;
    if(i==PairCharge)            cut.at(PairCharge)=1.;
    if(i==PairMass)              cut.at(PairMass)=40.;
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
      title.at(i)="Trigger Matching";
      hlabel="At least 1 good pair with Trig+Matching";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Trigger_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Trigger_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    // else if(i==Id_and_Kin){
    //   title.at(i)="Id and Kinematic";
    //   hlabel="Number of Event with good particles";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    // else if(i==NPairsFound){
    //   title.at(i)="Pairs with good DeltaR";
    //   hlabel="Pairs with good DeltaR";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NPairsFound_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NPairsFound_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    //if(i==GoodIndex){
    //title.at(i)="Valid Index";
    //hlabel="Valid Index";
    //Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_GoodIndex_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_GoodIndex_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    // if(i==ZTTMC){
    //   title.at(i)="ZTT MC";
    //   hlabel="ZTT MC";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ZTTMC_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ZTTMC_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    // else if(i==METFilters){
    //   title.at(i)="MET Filters";
    //   hlabel="MET Filters";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_METFilters_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_METFilters_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    else if(i==Id_and_Kin){
      title.at(i)="Id and Kinematic";
      hlabel="Number of Event with good particles";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    // else if(i==genmatch){
    //   title.at(i)="genmatch";
    //   hlabel="genmatch";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_genmatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_genmatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    else if(i==TausIsolation){
      title.at(i)="Taus Isolation";
      hlabel="Isolation of Taus";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TausIsolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TausIsolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==AgainstEleMu){
      title.at(i)="Against Electrons and Muons";
      hlabel="Against Electrons and Muons";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_AgainstEleMu_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_AgainstEleMu_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    // else if(i==Tau2Isolation){
    //   title.at(i)="Tau2 Isolation";
    //   hlabel="Isolation of Tau2";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Tau2Isolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Tau2Isolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    else if(i==LeptonVeto){
      title.at(i)="Third Lepton Veto";
      hlabel="Third Lepton Veto";
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

  polarimetricAcopAngle=HConfig.GetTH1D(Name+"_polarimetricAcopAngle"," ",5,0.,2*TMath::Pi()," "," ");
  polarimetricGEFAcopAngle=HConfig.GetTH1D(Name+"_polarimetricGEFAcopAngle"," ",5,0.,2*TMath::Pi()," "," ");
  decayplaneAcopAngle=HConfig.GetTH1D(Name+"_decayplaneAcopAngle"," ",5,0.,2*TMath::Pi()," "," ");
  impactparameterAcopAngle=HConfig.GetTH1D(Name+"_impactparameterAcopAngle"," ",5,0.,2*TMath::Pi()," "," ");
  DPIPAcopAngle=HConfig.GetTH1D(Name+"_DPIPAcopAngle"," ",5,0.,2*TMath::Pi()," "," ");
  PVIPAcopAngle=HConfig.GetTH1D(Name+"_PVIPAcopAngle"," ",5,0.,2*TMath::Pi()," "," ");

  polarimetricAcopAngleTruth=HConfig.GetTH1D(Name+"_polarimetricAcopAngleTruth"," ",5,0.,2*TMath::Pi()," "," ");
  decayplaneAcopAngleTruth=HConfig.GetTH1D(Name+"_decayplaneAcopAngleTruth"," ",5,0.,2*TMath::Pi()," "," ");
  impactparameterAcopAngleTruth=HConfig.GetTH1D(Name+"_impactparameterAcopAngleTruth"," ",5,0.,2*TMath::Pi()," "," ");
  DPIPAcopAngleTruth=HConfig.GetTH1D(Name+"_DPIPAcopAngleTruth"," ",5,0.,2*TMath::Pi()," "," ");
  PVIPAcopAngleTruth=HConfig.GetTH1D(Name+"_PVIPAcopAngleTruth"," ",5,0.,2*TMath::Pi()," "," ");

  // mtt
  DeltaPhitauMuMTT=HConfig.GetTH1D(Name+"_DeltaPhitauMuMTT","DeltaPhitauMuMTT",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatauMuMTT=HConfig.GetTH1D(Name+"_DeltaEtatauMuMTT","DeltaEtatauMuMTT",50,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtauMuMTT=HConfig.GetTH1D(Name+"_DeltaPtauMuMTT","DeltaPtauMuMTT",50,-1,1,"#Delta P","Events");
  DeltaEtauMuMTT=HConfig.GetTH1D(Name+"_DeltaEtauMuMTT","DeltaEtauMuMTT",50,-1,1,"#Delta E","Events");
  DeltaPhitauHadMTT=HConfig.GetTH1D(Name+"_DeltaPhitauHadMTT","DeltaPhitauHadMTT",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatauHadMTT=HConfig.GetTH1D(Name+"_DeltaEtatauHadMTT","DeltaEtatauHadMTT",50,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtauHadMTT=HConfig.GetTH1D(Name+"_DeltaPtauHadMTT","DeltaPtauHadMTT",50,-1,1,"#Delta P","Events");
  DeltaEtauHadMTT=HConfig.GetTH1D(Name+"_DeltaEtauHadMTT","DeltaEtauHadMTT",50,-1,1,"#Delta E","Events");
  //svfit
  DeltaPhitauMuSVFit=HConfig.GetTH1D(Name+"_DeltaPhitauMuSVFit","DeltaPhitauMuSVFit",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatauMuSVFit=HConfig.GetTH1D(Name+"_DeltaEtatauMuSVFit","DeltaEtatauMuSVFit",50,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtauMuSVFit=HConfig.GetTH1D(Name+"_DeltaPtauMuSVFit","DeltaPtauMuSVFit",50,-1,1,"#Delta P","Events");
  DeltaEtauMuSVFit=HConfig.GetTH1D(Name+"_DeltaEtauMuSVFit","DeltaEtauMuSVFit",50,-1,1,"#Delta E","Events");
  DeltaPhitauHadSVFit=HConfig.GetTH1D(Name+"_DeltaPhitauHadSVFit","DeltaPhitauHadSVFit",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatauHadSVFit=HConfig.GetTH1D(Name+"_DeltaEtatauHadSVFit","DeltaEtatauHadSVFit",100,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtauHadSVFit=HConfig.GetTH1D(Name+"_DeltaPtauHadSVFit","DeltaPtauHadSVFit",50,-1,1,"#Delta P","Events");
  DeltaEtauHadSVFit=HConfig.GetTH1D(Name+"_DeltaEtauHadSVFit","DeltaEtauHadSVFit",50,-1,1,"#Delta E","Events");
  //mixed
  DeltaPhitauMuMixed=HConfig.GetTH1D(Name+"_DeltaPhitauMuMixed","DeltaPhitauMuMixed",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatauMuMixed=HConfig.GetTH1D(Name+"_DeltaEtatauMuMixed","DeltaEtatauMuMixed",50,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtauMuMixed=HConfig.GetTH1D(Name+"_DeltaPtauMuMixed","DeltaPtauMuMixed",50,-1,1,"#Delta P","Events");
  DeltaEtauMuMixed=HConfig.GetTH1D(Name+"_DeltaEtauMuMixed","DeltaEtauMuMixed",50,-1,1,"#Delta E","Events");
  DeltaPhitauHadMixed=HConfig.GetTH1D(Name+"_DeltaPhitauHadMixed","DeltaPhitauHadMixed",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatauHadMixed=HConfig.GetTH1D(Name+"_DeltaEtatauHadMixed","DeltaEtatauHadMixed",50,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtauHadMixed=HConfig.GetTH1D(Name+"_DeltaPtauHadMixed","DeltaPtauHadMixed",50,-1,1,"#Delta P","Events");
  DeltaEtauHadMixed=HConfig.GetTH1D(Name+"_DeltaEtauHadMixed","DeltaEtauHadMixed",50,-1,1,"#Delta E","Events");
  //GEF
  DeltaPhitauMuGEF=HConfig.GetTH1D(Name+"_DeltaPhitauMuGEF","DeltaPhitauMuGEF",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatauMuGEF=HConfig.GetTH1D(Name+"_DeltaEtatauMuGEF","DeltaEtatauMuGEF",50,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtauMuGEF=HConfig.GetTH1D(Name+"_DeltaPtauMuGEF","DeltaPtauMuGEF",50,-1,1,"#Delta P","Events");
  DeltaEtauMuGEF=HConfig.GetTH1D(Name+"_DeltaEtauMuGEF","DeltaEtauMuGEF",50,-1,1,"#Delta E","Events");
  DeltaPhitauHGEF=HConfig.GetTH1D(Name+"_DeltaPhitauHGEF","DeltaPhitauHGEF",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatauHGEF=HConfig.GetTH1D(Name+"_DeltaEtatauHGEF","DeltaEtatauHGEF",50,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtauHGEF=HConfig.GetTH1D(Name+"_DeltaPtauHGEF","DeltaPtauHGEF",50,-1,1,"#Delta P","Events");
  DeltaEtauHGEF=HConfig.GetTH1D(Name+"_DeltaEtauHGEF","DeltaEtauHGEF",50,-1,1,"#Delta E","Events");
  //
  /*DeltaPhitauMuGEF=HConfig.GetTH1D(Name+"_DeltaPhitauMuGEF","DeltaPhitauMuGEF",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatauMuGEF=HConfig.GetTH1D(Name+"_DeltaEtatauMuGEF","DeltaEtatauMuGEF",50,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtauMuGEF=HConfig.GetTH1D(Name+"_DeltaPtauMuGEF","DeltaPtauMuGEF",50,-1,1,"#Delta P","Events");
  DeltaEtauMuGEF=HConfig.GetTH1D(Name+"_DeltaEtauMuGEF","DeltaEtauMuGEF",50,-1,1,"#Delta E","Events");
  DeltaPhitauHGEF=HConfig.GetTH1D(Name+"_DeltaPhitauHGEF","DeltaPhitauHGEF",50,-0.1,0.1,"#Delta#phi","Events");
  DeltaEtatauHGEF=HConfig.GetTH1D(Name+"_DeltaEtatauHGEF","DeltaEtatauHGEF",50,-0.1,0.1,"#Delta#eta","Events");
  DeltaPtauHGEF=HConfig.GetTH1D(Name+"_DeltaPtauHGEF","DeltaPtauHGEF",50,-1,1,"#Delta P","Events");
  DeltaEtauHGEF=HConfig.GetTH1D(Name+"_DeltaEtauHGEF","DeltaEtauHGEF",50,-1,1,"#Delta E","Events");*/

  //REF
  RefX=HConfig.GetTH1D(Name+"_RefX","RefX",50,-1,1,"ref","Events");
  RefY=HConfig.GetTH1D(Name+"_RefY","RefY",50,-1,1,"ref","Events");
  RefZ=HConfig.GetTH1D(Name+"_RefZ","RefZ",50,-1,1,"ref","Events");

  SVfitMTTdR1=HConfig.GetTH1D(Name+"_SVfitMTTdR1","SVfitMTTdR1",50,0,0.06,"#Delta (R) #tau 1","a.u");
  SVfitMTTdR2=HConfig.GetTH1D(Name+"_SVfitMTTdR2","SVfitMTTdR2",50,0,0.06,"#Delta (R) #tau 2","a.u");

  PcorrEtaSVfitMTT1=HConfig.GetTH2D(Name+"_PcorrEtaSVfitMTT1","PcorrEtaSVfitMTT1",50,-1,1,50,-0.2,0.2,"P SVfit","#Eta MTT");
  PcorrPhiSVfitMTT1=HConfig.GetTH2D(Name+"_PcorrPhiSVfitMTT1","PcorrPhiSVfitMTT1",50,-1,1,50,-0.2,0.2,"P SVfit","#Phi MTT");

  PcorrEtaSVfitMTT2=HConfig.GetTH2D(Name+"_PcorrEtaSVfitMTT2","PcorrEtaSVfitMTT2",50,-1,1,50,-0.2,0.2,"P SVfit","#Eta MTT");
  PcorrPhiSVfitMTT2=HConfig.GetTH2D(Name+"_PcorrPhiSVfitMTT2","PcorrPhiSVfitMTT2",50,-1,1,50,-0.2,0.2,"P SVfit","#Phi MTT");

  dRandPcorrEta1=HConfig.GetTH3D(Name+"_dRandPcorrEta1","dRandPcorrEta1",50,-1,1,50,-0.2,0.2,50,0,0.06,"P SVfit","#Eta MTT","#Delta R");
  dRandPcorrPhi1=HConfig.GetTH3D(Name+"_dRandPcorrPhi1","dRandPcorrPhi1",50,-1,1,50,-0.2,0.2,50,0,0.06,"P SVfit","#Phi MTT","#Delta R");

  dRandPcorrEta2=HConfig.GetTH3D(Name+"_dRandPcorrEta2","dRandPcorrEta2",50,-1,1,50,-0.2,0.2,50,0,0.06,"P SVfit","#Eta MTT","#Delta R");
  dRandPcorrPhi2=HConfig.GetTH3D(Name+"_dRandPcorrPhi2","dRandPcorrPhi2",50,-1,1,50,-0.2,0.2,50,0,0.06,"P SVfit","#Phi MTT","#Delta R");

  PullAcopPV=HConfig.GetTH1D(Name+"_PullAcopPV","PullAcopPV",50,-1,1,"pull AcopPV","u.a");
  dR1vsAcopPV=HConfig.GetTH2D(Name+"_dR1vsAcopPV","dR1vsAcopPV",50,0,0.06,50,-1,1,"dR","acop");
  dR2vsAcopPV=HConfig.GetTH2D(Name+"_dR2vsAcopPV","dR2vsAcopPV",50,0,0.06,50,-1,1,"dR","acop");
  P1vsAcopPV=HConfig.GetTH2D(Name+"_P1vsAcopPV","P1vsAcopPV",50,-1,1,50,-1,1,"P svfit","acop");
  P2vsAcopPV=HConfig.GetTH2D(Name+"_P2vsAcopPV","P2vsAcopPV",50,-1,1,50,-1,1,"P svfit","acop");
  Phi1vsAcopPV=HConfig.GetTH2D(Name+"_Phi1vsAcopPV","Phi1vsAcopPV",50,-0.2,0.2,50,-1,1,"phi mtt","acop");
  Eta1vsAcopPV=HConfig.GetTH2D(Name+"_Eta1vsAcopPV","Eta1vsAcopPV",50,-0.2,0.2,50,-1,1,"eta mtt","acop");
  Phi2vsAcopPV=HConfig.GetTH2D(Name+"_Phi2vsAcopPV","Phi2vsAcopPV",50,-0.2,0.2,50,-1,1,"phi mtt","acop");
  Eta2vsAcopPV=HConfig.GetTH2D(Name+"_Eta2vsAcopPV","Eta2vsAcopPV",50,-0.2,0.2,50,-1,1,"eta mtt","acop");

  Fraction1=HConfig.GetTH1D(Name+"_Fraction1","Fraction1",2,0,2," "," ");
  Fraction2=HConfig.GetTH1D(Name+"_Fraction2","Fraction2",2,0,2," "," ");

  Selection::ConfigureHistograms();   //   do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);  // do not remove

}

void  HCPMuTau::Store_ExtraDist(){

  Extradist1d.push_back(&polarimetricAcopAngle);
  Extradist1d.push_back(&polarimetricGEFAcopAngle);
  Extradist1d.push_back(&decayplaneAcopAngle);
  Extradist1d.push_back(&impactparameterAcopAngle);
  Extradist1d.push_back(&DPIPAcopAngle);
  Extradist1d.push_back(&PVIPAcopAngle);

  Extradist1d.push_back(&polarimetricAcopAngleTruth);
  Extradist1d.push_back(&decayplaneAcopAngleTruth);
  Extradist1d.push_back(&impactparameterAcopAngleTruth);
  Extradist1d.push_back(&DPIPAcopAngleTruth);
  Extradist1d.push_back(&PVIPAcopAngleTruth);

  //mtt
  Extradist1d.push_back(&DeltaPhitauMuMTT);
  Extradist1d.push_back(&DeltaEtatauMuMTT);
  Extradist1d.push_back(&DeltaPtauMuMTT);
  Extradist1d.push_back(&DeltaEtauMuMTT);
  Extradist1d.push_back(&DeltaPhitauHadMTT);
  Extradist1d.push_back(&DeltaEtatauHadMTT);
  Extradist1d.push_back(&DeltaPtauHadMTT);
  Extradist1d.push_back(&DeltaEtauHadMTT);
  //svfit
  Extradist1d.push_back(&DeltaPhitauMuSVFit);
  Extradist1d.push_back(&DeltaEtatauMuSVFit);
  Extradist1d.push_back(&DeltaPtauMuSVFit);
  Extradist1d.push_back(&DeltaEtauMuSVFit);
  Extradist1d.push_back(&DeltaPhitauHadSVFit);
  Extradist1d.push_back(&DeltaEtatauHadSVFit);
  Extradist1d.push_back(&DeltaPtauHadSVFit);
  Extradist1d.push_back(&DeltaEtauHadSVFit);
  //mixed
  Extradist1d.push_back(&DeltaPhitauMuMixed);
  Extradist1d.push_back(&DeltaEtatauMuMixed);
  Extradist1d.push_back(&DeltaPtauMuMixed);
  Extradist1d.push_back(&DeltaEtauMuMixed);
  Extradist1d.push_back(&DeltaPhitauHadMixed);
  Extradist1d.push_back(&DeltaEtatauHadMixed);
  Extradist1d.push_back(&DeltaPtauHadMixed);
  Extradist1d.push_back(&DeltaEtauHadMixed);
  //GEF
  Extradist1d.push_back(&DeltaPhitauMuGEF);
  Extradist1d.push_back(&DeltaEtatauMuGEF);
  Extradist1d.push_back(&DeltaPtauMuGEF);
  Extradist1d.push_back(&DeltaEtauMuGEF);
  Extradist1d.push_back(&DeltaPhitauHGEF);
  Extradist1d.push_back(&DeltaEtatauHGEF);
  Extradist1d.push_back(&DeltaPtauHGEF);
  Extradist1d.push_back(&DeltaEtauHGEF);
  //
  /*Extradist1d.push_back(&DeltaPhitauPiGEF);
  Extradist1d.push_back(&DeltaEtatauPiGEF);
  Extradist1d.push_back(&DeltaPtauPiGEF);
  Extradist1d.push_back(&DeltaEtauPiGEF);
  Extradist1d.push_back(&DeltaPhitauHGEF);
  Extradist1d.push_back(&DeltaEtatauHGEF);
  Extradist1d.push_back(&DeltaPtauHGEF);
  Extradist1d.push_back(&DeltaEtauHGEF);*/
  //ref
  Extradist1d.push_back(&RefX);
  Extradist1d.push_back(&RefY);
  Extradist1d.push_back(&RefZ);

  Extradist1d.push_back(&SVfitMTTdR1);
  Extradist1d.push_back(&SVfitMTTdR2);

  Extradist2d.push_back(&PcorrEtaSVfitMTT1);
  Extradist2d.push_back(&PcorrPhiSVfitMTT1);

  Extradist2d.push_back(&PcorrEtaSVfitMTT2);
  Extradist2d.push_back(&PcorrPhiSVfitMTT2);

  Extradist3d.push_back(&dRandPcorrEta1);
  Extradist3d.push_back(&dRandPcorrPhi1);

  Extradist3d.push_back(&dRandPcorrEta2);
  Extradist3d.push_back(&dRandPcorrPhi2);
 
  Extradist1d.push_back(&PullAcopPV);
  Extradist2d.push_back(&dR1vsAcopPV);
  Extradist2d.push_back(&dR2vsAcopPV);
  Extradist2d.push_back(&P1vsAcopPV);
  Extradist2d.push_back(&P2vsAcopPV);
  Extradist2d.push_back(&Phi1vsAcopPV);
  Extradist2d.push_back(&Eta1vsAcopPV);
  Extradist2d.push_back(&Phi2vsAcopPV);
  Extradist2d.push_back(&Eta2vsAcopPV);

  Extradist1d.push_back(&Fraction1);
  Extradist1d.push_back(&Fraction2);
}

void  HCPMuTau::doEvent()  { //  Method called on every event

  unsigned int t;                // sample type, you may manage in your further analysis, if needed
  int id(Ntp->GetMCID());  //read event ID of a sample

  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}  //  gives a warning if list of samples in Histo.txt  and SkimSummary.log do not coincide

  bool isEmbed=(id==36);
  bool trig=false;
  std::vector<int> TriggerIndexVector ;
  std::vector<TString>  MatchedTriggerNames;
  value.at(Trigger)=0;

  int TauMu= -1;
  int TauHad= -1;
  bool GenMatchSelection=false;
  
  TauMu = Ntp->MuIndex();
  TauHad = Ntp->TauIndex();
  int j=0;
  int GenMatch1=6, GenMatch2=6;
  
  //if(!Ntp->isData()|| isEmbed){GenMatch1=Ntp->gen_match_1(0);}
  GenMatch1=Ntp->GetGenMatch(TauHad);
  //else{GenMatch1=6;}
  if(id==30 || id==33 ||id==203)GenMatchSelection=(!(GenMatch1==6) && !(GenMatch1==5)); //ZL
  else if(id==36)GenMatchSelection=(GenMatch1==5);
  else if(id==23 && id==20)GenMatchSelection=false; //WTaunu
  else if(id==201|| id==202)GenMatchSelection=!((GenMatch1==5)||(GenMatch1==6));
  else if(id==71||id==72||id==73||id==74||id==47||id==48||id==49||id==50||id==51||id==52||id==53||id==54||id==55||id==56||id==57||id==58)GenMatchSelection=(!(GenMatch1==6) &&  !(GenMatch1==5)); // VVT
  else if(id==70 ||id==701 ||id==702 ||id==703 )GenMatchSelection=(!(GenMatch1==6) &&  !(GenMatch1==5));
  else GenMatchSelection=true;
  // Here you can defined different type of weights you want to apply to events.
  double wobs=1;
  double w=1;
  double wIDSF1=1.,wIDSF2=1.;
  double wIDSFAntiE1=1.,wIDSFAntiE2=1.;
  double wIDSFAntiMu1=1.,wIDSFAntiMu2=1.;
  double wTrgSF1=1.,wTrgSF2=1.;

  if(!Ntp->isData() && id!=DataMCType::QCD && !isEmbed) w*=Ntp->PUReweight();
  string TES="Nom";
  if((!Ntp->isData() && id!=DataMCType::QCD) || isEmbed) {

    wTrgSF1=Ntp->TriggerSF(TauHad,GenMatch1,TES,"Nom");
    //wTrgSF2=Ntp->TriggerSF(TauHad,GenMatch2,TES,"Nom");

    w*=wTrgSF1;
    w*=wTrgSF2;

    wIDSF1=Ntp->IDSF(TauHad,GenMatch1,TES);
    //wIDSF2=Ntp->IDSF(TauHad,GenMatch2,TES);
    if(!isEmbed){
      wIDSFAntiE1=Ntp->IDSF(TauHad,GenMatch1,TES,"ele");
      //wIDSFAntiE2=Ntp->IDSF(TauHad,GenMatch2,TES,"ele");
      wIDSFAntiMu1=Ntp->IDSF(TauHad,GenMatch1,TES,"mu");
      //wIDSFAntiMu2=Ntp->IDSF(TauHad,GenMatch2,TES,"mu");
    }
    w*=wIDSF1*wIDSF2*wIDSFAntiE1*wIDSFAntiE2*wIDSFAntiMu1*wIDSFAntiMu2;
  }
  TLorentzVector genMomentum(0,0,0,0);
  if(id==33 || id == 10110333 || id == 10110433|| id == 10130533|| id ==10210333|| id == 10210433|| id == 10230533|| id ==10310333 || id ==10330533 || id ==10410433 || id == 10410333|| id == 10430533|| id == 30530533) {
    for(unsigned int imc=0; imc < Ntp->NGenParts(); imc++){
      if(fabs(Ntp->Genpart_pdg(imc)) ==15   &&  Ntp->CHECK_BIT(Ntp->Genpart_flags(imc),0)&& Ntp->Genpart_status(imc) ==2) {
	if(Ntp->Genpart_P4(imc).Pt() > 8){
	  genMomentum+=Ntp->Genpart_P4(imc);
	}
      }
    }
  }
  if( id == 30) {
    for(unsigned int imc=0; imc < Ntp->NGenParts(); imc++){
      if((fabs(Ntp->Genpart_pdg(imc)) ==11 || fabs(Ntp->Genpart_pdg(imc)) ==13) && Ntp->Genpart_status(imc) ==1  ){
	if(Ntp->Genpart_P4(imc).Pt() > 8){
	  genMomentum+=Ntp->Genpart_P4(imc);
	}
      }
    }
  }
  float zptw(1.);
  if(!Ntp->isData() && genMomentum.Pt()!=0 && genMomentum.M() > 75 && genMomentum.M() < 120) {
    zptw = Ntp->ZPtReweight(genMomentum);
  }
  w*=zptw;

  double top_wt = 1.0;
  if(id==70 ||id==701||id==702||id==703){
    for (unsigned i = 0; i < Ntp->NGenParts(); i++){
      unsigned pdgid = abs(Ntp->Genpart_pdg(i));
      if(pdgid == 6 && Ntp->Genpart_flags(i)==8 && Ntp->Genpart_flags(i)==13){ // Flag From Hard Process and isLastCopy
	double pt = Ntp->Genpart_P4(i).Pt();
	pt = std::min(pt, 472.);
	double a = 0.088, b = -0.00087, c = 9.2e-07;
	top_wt *= std::exp(a + b * pt + c * pt*pt);
      }
    }
    top_wt = std::sqrt(top_wt);
    w*=top_wt;
  }
  int imc1=-1,imc2=-1;
  if(isEmbed){
    for(unsigned int imc=0; imc < Ntp->NGenParts(); imc++){
      if(fabs(Ntp->Genpart_pdg(imc))==15 && imc1==-1)imc1=imc;
      else if(fabs(Ntp->Genpart_pdg(imc))==15 && imc1!=-1)imc2=imc;
    }
    w*=Ntp->EmbeddingSelectionSF(imc1,imc2);
  }
  float PUPPImetCorr_px=Ntp->PUPPImet()*cos(Ntp->PUPPImetphi());
  float PUPPImetCorr_py=Ntp->PUPPImet()*sin(Ntp->PUPPImetphi());


  if((!Ntp->isData() && id!=DataMCType::QCD) || isEmbed){

    if (isEmbed && Ntp->MC_weight()>10000.)w*=Ntp->MC_weight()*0.000000001; //problem with 2016
    else w*=Ntp->MC_weight();

  }  //generator weight because negative weights for this samples
  if(Ntp->year()==2016)
    {
      if(id==11)w*=0.2455;
      else if(id==12)w*=0.2727;
      else if(id==45)w*=0.2546;
      else if(id==461)w*=0.2596;
      else if(id==460)w*=0.2425;
    }
  else if(Ntp->year()==2017)
    {
      if(id==11)w*=0.2447;
      else if(id==12)w*=0.2697;
      else if(id==45)w*=0.2514;
      else if(id==461)w*=0.2567;
      else if(id==460)w*=0.2394;
    }
  else if(Ntp->year()==2018)
    {
      if(id==11)w*=0.2446;
      else if(id==12)w*=0.2695;
      else if(id==45)w*=0.2513;
      else if(id==461)w*=0.2563;
      else if(id==460)w*=0.2397;
    }
  w*=Ntp->stitch_weight(isDY1050);


  if(!Ntp->isData() && !isEmbed && (Ntp->year()==2016||Ntp->year()==2017))w*=Ntp->prefiringweight();
  //std::vector<unsigned int> exclude_cuts;
  //exclude_cuts.push_back(TausIsolation);
  classic_svFit::LorentzVector tau1P4;
  classic_svFit::LorentzVector tau2P4;
  classic_svFit::LorentzVector tau1P4mtt;
  classic_svFit::LorentzVector tau2P4mtt;

  TLorentzVector TauMuVis;
  TLorentzVector TauHadVis;
  //if(passAllBut(exclude_cuts)){
  TauMuVis = Ntp->Daughters_P4(TauMu);
  TauHadVis = Ntp->Daughters_P4(TauHad);
  //}
  unsigned long long evt_ = Ntp->EventNumber();
  std::pair<float, int> max_pair;
  std::vector<float> scores = {};
  double PUPPIMET=sqrt(PUPPImetCorr_px*PUPPImetCorr_px+PUPPImetCorr_py*PUPPImetCorr_py);
  double PUPPIMETPhi=(TVector3(PUPPImetCorr_px,PUPPImetCorr_py,0)).Phi();
  std::vector<TLorentzVector> Pions1;
  std::vector<TLorentzVector> Pions2;
  std::vector<double> Pions1Charge;
  std::vector<double> Pions2Charge;

  TVector3 tauPrimaryVertex , tauNoBSPrimaryVertex, tauBSPrimaryVertex, tauNoBSZNominalPrimaryVertex, tauBSZNominalPrimaryVertex, TauminusSecondaryVertex , TauplusSecondaryVertex,tauWithTracksBSZNominalPrimaryVertex,tauWithTracksBSPrimaryVertex;

  TVector3 TauminusDirection , TauplusDirection, TauplusDirectionNoBS, TauplusDirectionBS, TauplusDirectionNoBSZNominal, TauplusDirectionBSZNominal,TauminusDirectionNoBS, TauminusDirectionBS, TauminusDirectionNoBSZNominal, TauminusDirectionBSZNominal;

  double thetaGJ_Tauminus , thetaGJ_Tauplus;
  TLorentzVector a1LV_Tauminus , a1LV_Tauplus, a1LVRefit_Tauminus , a1LVRefit_Tauplus;
  TLorentzVector TauminusPairConstraint, TauplusPairConstraintNoBS, TauplusPairConstraintBS, TauplusPairConstraintNoBSZNominal, TauplusPairConstraintBSZNominal, TauminusPairConstraintNoBS, TauminusPairConstraintBS, TauminusPairConstraintNoBSZNominal, TauminusPairConstraintBSZNominal, TauminusPairConstraintWithTracksBSZNominal,TauminusPairConstraintWithTracksBS, TauplusPairConstraintWithTracksBSZNominal,TauplusPairConstraintWithTracksBS;
  TLorentzVector TauplusSmall, TauplusLarge, TauplusPairConstraint,TauplusPairConstraintMVA,TauplusPairConstraintNoBSNewMVA,TauplusPairConstraintBSNewMVA;
  bool isPlusReal=true, isMinusReal=true;
  std::vector<TLorentzVector> solutions, solutionsNoBS, solutionsBS, solutionsNoBSZNominal, solutionsBSZNominal, solutionsWithTracksBSZNominal, solutionsWithTracksBS;
  TLorentzVector HPairConstraint;

  std::vector<size_t> hashes;
  size_t hash = 0;

  bool isRefitNoBS=true;
  bool isRefitBS=true;
  bool isRefitNoBSZNominal=true;
  bool isRefitBSZNominal=true;
  double Spin_WT=Ntp->TauSpinerGet(TauSpinerInterface::Spin,CPstate);
  double FlipSpin_WT=Ntp->TauSpinerGet(TauSpinerInterface::FlipSpin,CPstate);
  double Wspin;
  double Wflipspin;
  if(id == 11 || id == 12|| id == 45 || id ==461||id == 460)Wspin=w*Spin_WT;
  else Wspin=w;
  if(id == 11 || id == 12|| id == 45 ||id ==461 ||id ==460)Wflipspin=w*FlipSpin_WT;
  else Wflipspin=w;
  TLorentzVector testTruth(0,0,0,0);

  double DEtatauMuMTT=-99, DEtatauHadMTT=-99, DPhitauMuMTT=-99, DPhitauHadMTT=-99;
  double DEtauMuMTT=-99, DEtauHadMTT=-99, DPtauMuMTT=-99, DPtauHadMTT=-99;
  double DEtatauMuSVFit=-99, DEtatauHadSVFit=-99, DPhitauMuSVFit=-99, DPhitauHadSVFit=-99;
  double DEtauMuSVFit=-99, DEtauHadSVFit=-99, DPtauMuSVFit=-99, DPtauHadSVFit=-99;
  double DEtatauMuMixed=-99, DEtatauHadMixed=-99, DPhitauMuMixed=-99, DPhitauHadMixed=-99;
  double DEtauMuMixed=-99, DEtauHadMixed=-99, DPtauMuMixed=-99, DPtauHadMixed=-99;
  double DEtatauMuGEF=-99, DEtatauHGEF=-99, DPhitauMuGEF=-99, DPhitauHGEF=-99;
  double DEtauMuGEF=-99, DEtauHGEF=-99, DPtauMuGEF=-99, DPtauHGEF=-99;

  vector<TLorentzVector> HadPionsTruth;
  vector<double> HadPionsChargeTruth;

  vector<TLorentzVector> HadRefitPions;
  vector<double> HadRefitPionsCharge;

  TVector3 Muon_ref, Muon_refTruth;
  TVector3 Pion_ref, Pion_refTruth;

  SCalculator ScalcPVTruth;
  SCalculator ScalcIPTruth;
  SCalculator ScalcDPTruth;
  SCalculator ScalcDPIPTruth;
  SCalculator ScalcPVIPTruth;
  SCalculator ScalcPVRefitWithTracksBS;
  SCalculator ScalcPV;
  SCalculator ScalcDP;
  SCalculator ScalcIP;
  SCalculator ScalcDPIP;
  SCalculator ScalcPVIP;
  double Acop_PVIPGEF = -99, Acop_PVIP = -99, Acop_DPIP = -99;
  double Acop_PVIPTruth = -99, Acop_DPIPTruth = -99;

  TLorentzVector zeroLV(0,0,0,0);
  std::vector<TLorentzVector> VectZeroLV;
  VectZeroLV.push_back(zeroLV);
  VectZeroLV.push_back(zeroLV);
  VectZeroLV.push_back(zeroLV);

  TLorentzVector TauMuSVFit, TauHadSVFit;
  TLorentzVector TauMuMTT, TauHadMTT;
  TLorentzVector TauMuMixed, TauHadMixed;
  TLorentzVector TauPlusMixed, TauMinusMixed;
  TLorentzVector TauHGEF, TauMuGEF;
  TLorentzVector TauMuTruth, TauHadTruth;
  TLorentzVector TauMinusTruth, TauPlusTruth;
  TLorentzVector MuonP4Truth;
  unsigned int Taumutruth=0;
  unsigned int Tauhadtruth=0;
  unsigned int Tauminustruth=0;
  unsigned int Tauplustruth=0;
  string CHANNEL = string(Channel);
  //
  bool a1 = false;
  bool rho = false;
  bool pion = false;
  bool a1mu = false;
  bool rhomu = false;
  bool pionmu = false;
  bool charge = false, channel = false;
  bool minushad = false, plushad = false;
  //
  bool mutruth = false;
  bool a1truth = false;
  bool rhotruth = false;
  bool piontruth =false;
  bool a1mutruth = false;
  bool rhomutruth = false;
  bool pionmutruth = false;
  bool chargetruth = false, channeltruth = false;
  bool minushadtruth = false, plushadtruth = false;
  //
  bool recoPionsAreOk = false, genPionsAreOk = false;
  bool recoTausAreOk = false, GEFTausAreOk = false, mixedTausAreOk = false, genTausAreOk = false;
  //GENERATED LEVEL
  if(Ntp->MCTau_JAK(0) == 2){
    TauMuTruth = Ntp->MCTau_p4(0);
    TauHadTruth = Ntp->MCTau_p4(1);
    Taumutruth = 0;
    Tauhadtruth = 1;
    mutruth = true;
  }
  else if(Ntp->MCTau_JAK(1) == 2){
    TauMuTruth = Ntp->MCTau_p4(1);
    TauHadTruth = Ntp->MCTau_p4(0);
    Taumutruth = 1;
    Tauhadtruth = 0;
    mutruth = true;
  }

  if(Ntp->MCTau_charge(Tauhadtruth)<0){
    minushadtruth = true;
    Tauminustruth = Tauhadtruth;
    Tauplustruth = Taumutruth;
    TauMinusTruth = TauHadTruth;
    TauPlusTruth = TauMuTruth;
  }
  if(Ntp->MCTau_charge(Tauhadtruth)>0){
    plushadtruth = true;
    Tauplustruth = Tauhadtruth;
    Tauminustruth = Taumutruth;
    TauPlusTruth = TauHadTruth;
    TauMinusTruth = TauMuTruth;
  }
  //3 prongs decay
  if(Ntp->MCTau_JAK(Tauhadtruth) == 5) a1truth=true;
  //1prong + pi0 decay
  if(Ntp->MCTau_JAK(Tauhadtruth) == 4) rhotruth=true;
  //1 prong decay
  if(Ntp->MCTau_JAK(Tauhadtruth) == 3) piontruth=true;

  if(a1truth && mutruth && CHANNEL == "A1MU") a1mutruth=true;
  if(rhotruth && mutruth && CHANNEL == "RHOMU") rhomutruth=true;
  if(piontruth && mutruth && CHANNEL == "PIONMU") pionmutruth=true;
  if(a1mutruth || rhomutruth || pionmutruth) channeltruth = true;
     if(channeltruth){
       if(a1mutruth){
	 if(minushadtruth){
           HadPionsTruth = Ntp->GetTruthPionsFromA1(Tauminustruth);
	   HadPionsChargeTruth.push_back(1);
	   HadPionsChargeTruth.push_back(-1);
	   HadPionsChargeTruth.push_back(-1);

           Muon_refTruth = Ntp->MCTauandProd_Vertex(Tauplustruth,2) - Ntp->MCTauandProd_Vertex(Tauplustruth,0);
           MuonP4Truth = Ntp->GetTruthTauProductLV(2,13,Tauplustruth);
	 }
	 if(plushadtruth){
	   HadPionsTruth = Ntp->GetTruthPionsFromA1(Tauminustruth);
	   HadPionsChargeTruth.push_back(-1);
	   HadPionsChargeTruth.push_back(1);
	   HadPionsChargeTruth.push_back(1);

           Muon_refTruth = Ntp->MCTauandProd_Vertex(Tauminustruth,2) - Ntp->MCTauandProd_Vertex(Tauminustruth,0);
           MuonP4Truth = Ntp->GetTruthTauProductLV(2,13,Tauminustruth);
	 }
       }
       if(rhomutruth){
	 if(minushadtruth){
	   HadPionsTruth = Ntp->GetTruthPionsFromRho(Tauminustruth);
	   HadPionsChargeTruth.push_back(-1);
	   HadPionsChargeTruth.push_back(0);

           Muon_refTruth = Ntp->MCTauandProd_Vertex(Tauplustruth,2) - Ntp->MCTauandProd_Vertex(Tauplustruth,0);
           MuonP4Truth = Ntp->GetTruthTauProductLV(2,13,Tauplustruth);
	 }
	 if(plushadtruth){
	   HadPionsTruth = Ntp->GetTruthPionsFromRho(Tauplustruth);
	   HadPionsChargeTruth.push_back(1);
	   HadPionsChargeTruth.push_back(0);

           Muon_refTruth = Ntp->MCTauandProd_Vertex(Tauminustruth,2) - Ntp->MCTauandProd_Vertex(Tauminustruth,0);
           MuonP4Truth = Ntp->GetTruthTauProductLV(2,13,Tauminustruth);
	 }
       }
       if(pionmutruth){
	 if(minushadtruth){
	   HadPionsChargeTruth.push_back(-1);
           HadPionsTruth.push_back(Ntp->GetTruthTauProductLV(3,211,Tauminustruth));
           Pion_refTruth = Ntp->MCTauandProd_Vertex(Tauminustruth,2) - Ntp->MCTauandProd_Vertex(Tauminustruth,0);

           Muon_refTruth = Ntp->MCTauandProd_Vertex(Tauplustruth,2) - Ntp->MCTauandProd_Vertex(Tauplustruth,0);
           MuonP4Truth = Ntp->GetTruthTauProductLV(2,13,Tauplustruth);
	 }
	 if(plushadtruth){
	   HadPionsChargeTruth.push_back(1);
           HadPionsTruth.push_back(Ntp->GetTruthTauProductLV(3,211,Tauplustruth));
           Pion_refTruth = Ntp->MCTauandProd_Vertex(Tauplustruth,2) - Ntp->MCTauandProd_Vertex(Tauplustruth,0);

           Muon_refTruth = Ntp->MCTauandProd_Vertex(Tauminustruth,2) - Ntp->MCTauandProd_Vertex(Tauminustruth,0);
           MuonP4Truth = Ntp->GetTruthTauProductLV(2,13,Tauminustruth);
	 }
       }
     }
     //RECONSTRUCTED LEVEL
     if(Ntp->Daughters_charge(TauMu)>0 && Ntp->Daughters_charge(TauHad)<0 && TauMuVis.Pt() > 20 && TauHadVis.Pt() > 20)
       {
	 minushad = true;
	 charge = true;
       }
     else if(Ntp->Daughters_charge(TauHad)>0 && Ntp->Daughters_charge(TauMu)<0 && TauMuVis.Pt() > 20 && TauHadVis.Pt() > 20)
       {
         plushad = true;
	 charge = true;
       }
     if(Ntp->MVADM2017(TauHad) == 10 && Ntp->PFTau_hassecondaryVertex(TauHad) && Ntp->PFtauHasThreePions(TauHad)) a1=true;
     if(Ntp->MVADM2017(TauHad) == 1) rho=true;
     if(Ntp->MVADM2017(TauHad) == 0) pion=true;

     if(a1 && CHANNEL == "A1MU") a1mu=true;
     if(rho && CHANNEL == "RHOMU") rhomu=true;
     if(pion && CHANNEL == "PIONMU") pionmu=true;

     if((a1mu || rhomu || pionmu) && charge) channel = true;
     if(channel){

       TMatrixD metcov(2,2);
       double metx = PUPPIMET*cos(PUPPIMETPhi);
       double mety = PUPPIMET*sin(PUPPIMETPhi);

       metcov[0][0] = Ntp->PUPPIMETCov00();
       metcov[1][0] = Ntp->PUPPIMETCov01();
       metcov[0][1] = Ntp->PUPPIMETCov10();
       metcov[1][1] = Ntp->PUPPIMETCov11();

       TMatrixTSym<double> covMET(2);
       covMET[0][0] = Ntp->PUPPIMETCov00();
       covMET[0][1] = Ntp->PUPPIMETCov01();
       covMET[1][0] = Ntp->PUPPIMETCov10();
       covMET[1][1] = Ntp->PUPPIMETCov11();

       TMatrixT<double> MET(2,1);
       MET[0][0] = PUPPIMET*cos(PUPPIMETPhi);
       MET[1][0] = PUPPIMET*sin(PUPPIMETPhi);

       PTObject METInput(MET,covMET);

       FastMTT FastMTTAlgo;
       ClassicSVfit svfitAlgo1(0);
       double higgsmass;
       if(rhomu || pionmu){
       // // //---------  svfit ---------------------
       std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
       classic_svFit::MeasuredTauLepton lep1(classic_svFit::MeasuredTauLepton::kTauToMuDecay, TauMuVis.Pt(), TauMuVis.Eta(),  TauMuVis.Phi(), TauMuVis.M(),0);
       classic_svFit::MeasuredTauLepton lep2(classic_svFit::MeasuredTauLepton::kTauToHadDecay, TauHadVis.Pt(), TauHadVis.Eta(),  TauHadVis.Phi(), TauHadVis.M(),Ntp->MVADM2017(TauHad));

       measuredTauLeptons.push_back(lep1);
       measuredTauLeptons.push_back(lep2);

       svfitAlgo1.setHistogramAdapter(new classic_svFit::TauTauHistogramAdapter());
       svfitAlgo1.addLogM_fixed(true,5.);
       svfitAlgo1.setDiTauMassConstraint(125.35);
       FastMTTAlgo.run(measuredTauLeptons, metx, mety, metcov);
       svfitAlgo1.integrate(measuredTauLeptons,metx,mety, metcov );
       if(svfitAlgo1.isValidSolution()) {

	 higgsmass  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svfitAlgo1.getHistogramAdapter())->getMass();
	 //h_SVFitMass.at(t).Fill(higgsmass,w);
	 tau1P4mtt = FastMTTAlgo.getTau1P4();
	 tau2P4mtt = FastMTTAlgo.getTau2P4();

	 tau1P4 = static_cast<classic_svFit::TauTauHistogramAdapter*>(svfitAlgo1.getHistogramAdapter())->GetFittedTau1LV();
	 tau2P4 = static_cast<classic_svFit::TauTauHistogramAdapter*>(svfitAlgo1.getHistogramAdapter())->GetFittedTau2LV();

	 TauMuSVFit.SetPxPyPzE(tau1P4.x(),tau1P4.y(),tau1P4.z(),tau1P4.t());
	 TauHadSVFit.SetPxPyPzE(tau2P4.x(),tau2P4.y(),tau2P4.z(),tau2P4.t());
	 TauMuMTT.SetPxPyPzE(tau1P4mtt.x(),tau1P4mtt.y(),tau1P4mtt.z(),tau1P4mtt.t()); 
	 TauHadMTT.SetPxPyPzE(tau2P4mtt.x(),tau2P4mtt.y(),tau2P4mtt.z(),tau2P4mtt.t());

	 //MIXED TAUS
	 if(TauMuSVFit!=zeroLV && TauHadSVFit!=zeroLV && TauMuMTT!=zeroLV && TauHadMTT!=zeroLV && TauMuSVFit!=TauHadSVFit && TauMuMTT!=TauHadMTT){
	   //if(TauMinusSVFit.DeltaR(TauMinusMTT) < 0.03 && TauPlusSVFit.DeltaR(TauPlusMTT) < 0.03){
	   TauMuMixed.SetVect(TauMuSVFit.P()*TauMuMTT.Vect().Unit());
	   TauMuMixed.SetE(TauMuSVFit.E());
	   TauHadMixed.SetVect(TauHadSVFit.P()*TauHadMTT.Vect().Unit());
	   TauHadMixed.SetE(TauHadSVFit.E());
	   //}	
	   /*else{
	     TauPlusMixed = TauPlusMTT;
	     TauMinusMixed = TauMinusMTT;
	     }*/
	 }
       }
       }
       float m_sv_;
       float met_;
       float p_tt_;
       float m_vis_;
       float pt_1_;
       float pt_tt_;
       float pt_vis_;
       if(svfitAlgo1.isValidSolution()) m_sv_=higgsmass;
       else m_sv_=-1;
       m_vis_=(TauMuVis+TauHadVis).M();
       met_=PUPPIMET;

       pt_1_=TauMuVis.Pt();
       //pt_2_=Tau2P4.Pt();
       pt_tt_=(TauMuVis+TauHadVis+TLorentzVector(PUPPImetCorr_px,PUPPImetCorr_py,0,0)).Pt();
       pt_vis_=(TauMuVis+TauHadVis).Pt();
       //cout<<jdeta_<<" "<<jpt_1_<<"  "<<m_vis_<<"  "<<met_<<"  "<<mjj_<<"  "<<n_jets_<<endl;

       //BDT->Execute(jdeta_, jpt_1_, m_vis_, met_, mjj_, n_jets_, pt_1_, pt_tt_, pt_vis_, m_sv_, evt_,scores,max_pair);
  
       std::vector<double> PVRefitNoBS_X_temp, PVRefitNoBS_Y_temp, PVRefitNoBS_Z_temp;
       for(unsigned int i=0;i<Ntp->NPVRefitNoBS();i++)
	 {
	   PVRefitNoBS_X_temp.push_back(Ntp->RefitPVNoBS_x(i));
	   PVRefitNoBS_Y_temp.push_back(Ntp->RefitPVNoBS_y(i));
	   PVRefitNoBS_Z_temp.push_back(Ntp->RefitPVNoBS_z(i));
	 }

       std::vector<double> PVRefitBS_X_temp, PVRefitBS_Y_temp, PVRefitBS_Z_temp;
       for(unsigned int i=0;i<Ntp->NPVRefitBS();i++)
	 {
	   PVRefitBS_X_temp.push_back(Ntp->RefitPVBS_x(i));
	   PVRefitBS_Y_temp.push_back(Ntp->RefitPVBS_y(i));
	   PVRefitBS_Z_temp.push_back(Ntp->RefitPVBS_z(i));
	 }


       vector<size_t> VertexHashNoBS1_temp, VertexHashNoBS2_temp;
       for(unsigned int i=0;i<Ntp->NVertexHashNoBS();i++)
	 {
	   VertexHashNoBS1_temp.push_back(Ntp->VertexHashNoBS1(i));
	   VertexHashNoBS2_temp.push_back(Ntp->VertexHashNoBS2(i));

	 }

       vector<size_t> VertexHashBS1_temp, VertexHashBS2_temp;
       for(unsigned int i=0;i<Ntp->NVertexHashBS();i++)
	 {
	   VertexHashBS1_temp.push_back(Ntp->VertexHashBS1(i));
	   VertexHashBS2_temp.push_back(Ntp->VertexHashBS2(i));

	 }

       boost::hash_combine(hash, Ntp->LeptonHash(TauHad));
       boost::hash_combine(hash, Ntp->LeptonHash(TauMu));
       hashes.push_back(hash);
       hash = 0;
       boost::hash_combine(hash, Ntp->LeptonHash(TauMu));
       boost::hash_combine(hash, Ntp->LeptonHash(TauHad));
       hashes.push_back(hash);

       Muon_ref = Ntp->Daughters_pcaRefitPV(TauMu);

       if(a1mu){
	 //GEF
	 if(Ntp->Muon_TrackParticleHasMomentum(TauMu) == true && Ntp->isPVCovAvailable() == true){
	   TMatrixTSym<double> PV_cov = Ntp->PFTau_TIP_primaryVertex_cov();
	   LorentzVectorParticle A1 = Ntp->PFTau_a1_lvp(TauHad);
	   TrackParticle MUON = Ntp->Muon_TrackParticle(TauMu);
	   GlobalEventFit GEF(MUON,A1,METInput,tauBSPrimaryVertex,PV_cov);
	   GEF.setMassConstraint(125.3);
	   GEF.SetCorrectPt(false);
	   GEFObject FitTaus = GEF.Fit();
	   if(FitTaus.Fitconverged()){
	     TauHGEF = FitTaus.getTauH().LV();
	     TauMuGEF = FitTaus.getTauMu().LV();
	   }
	 }
	 //
	 HadRefitPions.push_back(Ntp->PFTauRefit_PionsP4(TauHad,0));
	 HadRefitPions.push_back(Ntp->PFTauRefit_PionsP4(TauHad,1));
	 HadRefitPions.push_back(Ntp->PFTauRefit_PionsP4(TauHad,2));
	 HadRefitPionsCharge.push_back(Ntp->PFTauRefit_PionsCharge(TauHad, 0));
	 HadRefitPionsCharge.push_back(Ntp->PFTauRefit_PionsCharge(TauHad, 1));
	 HadRefitPionsCharge.push_back(Ntp->PFTauRefit_PionsCharge(TauHad, 2));
       }
       if(rhomu){
	 HadRefitPions.push_back(Ntp->ChargedDaughters_P4(TauHad));
	 HadRefitPions.push_back(TauHadVis - Ntp->ChargedDaughters_P4(TauHad));
         if(minushad){
	   HadRefitPionsCharge.push_back(-1);
	   HadRefitPionsCharge.push_back(0);
         }
         if(plushad){
           HadRefitPionsCharge.push_back(1);
           HadRefitPionsCharge.push_back(0);
         }
       }
       if(pionmu){
	 HadRefitPions.push_back(Ntp->ChargedDaughters_P4(TauHad));
	 Pion_ref = Ntp->Daughters_pcaRefitPV(TauHad);
         if(minushad){
           HadRefitPionsCharge.push_back(-1);
         }
         if(minushad){
           HadRefitPionsCharge.push_back(1);
         }
       }
     }

if(HadPionsTruth!=VectZeroLV) genPionsAreOk = true;
if((TauMuTruth != TauHadTruth) && (TauHadTruth != zeroLV) && (TauMuTruth != zeroLV)) genTausAreOk = true;
   
if (std::isnan(Wspin)!=true && GenMatchSelection)
  {
    if(channeltruth)
      {
	//GEN LEVEL PLOTS
	if(genPionsAreOk)
	  {
	    if(a1mutruth){
	     if(minushadtruth){
	      if(genTausAreOk){
 	          Acop_PVIPTruth = ScalcPVIPTruth.AcopAngle_PVIP("a1","muon", TauMinusTruth, -1., HadPionsTruth, HadPionsChargeTruth, TauPlusTruth, MuonP4Truth, Muon_refTruth);
              }
	      Acop_DPIPTruth = ScalcDPIPTruth.AcopAngle_DPIP("a1","muon", HadPionsTruth, MuonP4Truth, Muon_refTruth);
	     }
             if(plushadtruth){
              if(genTausAreOk){
                  //Acop_PVIPTruth = ScalcPVIPTruth.AcopAngle_PVIP("a1","muon", TauPlusTruth, 1., HadPionsTruth, HadPionsChargeTruth, TauMinusTruth, MuonP4Truth, Muon_refTruth);
              }
              Acop_DPIPTruth = ScalcDPIPTruth.AcopAngle_DPIP("a1","muon", HadPionsTruth, MuonP4Truth, Muon_refTruth);
             }
	    }
	    if(rhomutruth){
	      if(genTausAreOk){
                Acop_PVIPTruth = ScalcPVIPTruth.AcopAngle_PVIP("rho","muon", TauHadTruth, 1., HadPionsTruth, HadPionsChargeTruth, TauMuTruth, MuonP4Truth, Muon_refTruth);
	      }
              Acop_DPIPTruth = ScalcDPIPTruth.AcopAngle_DPIP("rho","muon", HadPionsTruth, MuonP4Truth, Muon_refTruth);
	    }
	    if(pionmutruth){
	      if(genTausAreOk){
                Acop_DPIPTruth = ScalcDPIPTruth.AcopAngle_DPIP("pion","muon", HadPionsTruth, MuonP4Truth, Muon_refTruth);
	      }
              //Acop_IPTruth = ScalcIPTruth.AcopAngle_IP("pion","muon", HadPionsTruth, MuonP4Truth, Muon_refTruth);
	    }
	    DPIPAcopAngleTruth.at(t).Fill(Acop_DPIPTruth,Wspin);
	    PVIPAcopAngleTruth.at(t).Fill(Acop_PVIPTruth,Wspin);
	  } //gen pions
      } //gen channel

    if((HadRefitPions!=VectZeroLV) && (HadRefitPions!=VectZeroLV)) recoPionsAreOk = true;
    if((TauHadMixed!=TauMuMixed) && (TauHadMixed!=zeroLV) && (TauMuMixed!=zeroLV)) mixedTausAreOk = true;
    if((TauHGEF!=TauMuGEF) && (TauHGEF!=zeroLV) && (TauMuGEF!=zeroLV)) GEFTausAreOk = true;

    //RECO LEVEL PLOTS
    if(channel)
      {
	if(recoPionsAreOk)
	  {
	    if(a1mu){
	      if(GEFTausAreOk){
		/*if(minushad){
                  Acop_PVIP = ScalcPVIP.AcopAngle_PVIP("a1","muon", TauHGEF, -1., HadRefitPions, HadRefitPionsCharge, TauMuGEF, TauMuVis, Muon_ref);
                }
                if(plushad){
                  Acop_PVIP = ScalcPVIP.AcopAngle_PVIP("a1","muon", TauHGEF, 1., HadRefitPions, HadRefitPionsCharge, TauMuGEF, TauMuVis, Muon_ref);
                }*/
		TLorentzVector ZMF = TauHGEF + TauMuVis;
                TVector3 boost = -ZMF.BoostVector();

                vector<TLorentzVector> tauandprod;
		tauandprod.push_back(TauHGEF);
		for(unsigned int i = 0; i<HadRefitPions.size(); i++){
		  tauandprod.push_back(HadRefitPions.at(i));
		}

		SCalculator S1("a1");
                if(minushad){
		  S1.Configure(tauandprod,ZMF,-1);
                }
                if(plushad){
		  S1.Configure(tauandprod,ZMF,+1);
		}
                //VECTEUR POLA
                TVector3 h = S1.pv();
                TauHGEF.Boost(boost);
                h*=1./h.Mag();
                TVector3 n = TauHGEF.Vect().Unit();
                TVector3 k = (h.Cross(n)).Unit();
                //IMPACT PARAMETER
                TLorentzVector Muon_tlv = TauMuVis;
                //TVector3 Muon_ref = Ntp->Muon_TrackRef(MuIndex) - tauBSPrimaryVertex;
                TVector3 Muon_dir = Muon_tlv.Vect();
                double proj = Muon_ref*Muon_dir/Muon_dir.Mag2();
                TVector3 Muon_IP = Muon_ref-Muon_dir*proj;
                TLorentzVector eta(Muon_IP,0);
                eta.Boost(boost);
                TVector3 etaVec = eta.Vect().Unit();
                Muon_tlv.Boost(boost);
                TVector3 MuVec = Muon_tlv.Vect().Unit();
                TVector3 etaVecTrans = (etaVec - MuVec*(MuVec*etaVec)).Unit();
                Acop_PVIP = TMath::ATan2((k.Cross(etaVecTrans)).Mag(),k*etaVecTrans);
                double sign_PV = (k.Cross(etaVecTrans))*n;
                if (sign_PV>0) Acop_PVIP = 2.0*TMath::Pi() - Acop_PVIP;
                  Acop_PVIP = Acop_PVIP - 0.5*TMath::Pi();
                if (Acop_PVIP<0) Acop_PVIP = Acop_PVIP + 2*TMath::Pi();
		
	      }
	      /*if(mixedTausAreOk){
                if(minushad){
                  Acop_PVIP = ScalcPVIP.AcopAngle_PVIP("a1","muon", TauHadMixed, -1., HadRefitPions, HadRefitPionsCharge, TauMuMixed, TauMuVis, Muon_ref);
                }
                if(plushad){
                  Acop_PVIP = ScalcPVIP.AcopAngle_PVIP("a1","muon", TauHadMixed, 1., HadRefitPions, HadRefitPionsCharge, TauMuMixed, TauMuVis, Muon_ref);
                }
	      }*/
	      Acop_DPIP = ScalcDPIP.AcopAngle_DPIP("a1","muon", HadRefitPions, TauMuVis, Muon_ref);
	    }
	    if(rhomu){
	      if(mixedTausAreOk){
                if(minushad){
                  Acop_PVIP = ScalcPVIP.AcopAngle_PVIP("rho","muon", TauHadMixed, -1., HadRefitPions, HadRefitPionsCharge, TauMuMixed, TauMuVis, Muon_ref);
                }
                if(plushad){
                  Acop_PVIP = ScalcPVIP.AcopAngle_PVIP("rho","muon", TauHadMixed, 1., HadRefitPions, HadRefitPionsCharge, TauMuMixed, TauMuVis, Muon_ref);
                }
	      }
              Acop_DPIP = ScalcDPIP.AcopAngle_DPIP("rho","muon", HadRefitPions, TauMuVis, Muon_ref);
	    }
	    if(pionmu){
	      if(mixedTausAreOk){
                if(minushad){
                  Acop_PVIP = ScalcPVIP.AcopAngle_PVIP("pion","muon", TauHadMixed, -1., HadRefitPions, HadRefitPionsCharge, TauMuMixed, TauMuVis, Muon_ref);
                }
                if(plushad){
                  Acop_PVIP = ScalcPVIP.AcopAngle_PVIP("pion","muon", TauHadMixed, 1., HadRefitPions, HadRefitPionsCharge, TauMuMixed, TauMuVis, Muon_ref);
                }
	      }
	    }
	    //polarimetricGEFAcopAngle.at(t).Fill(Acop_PVIPGEF,Wspin);
	    DPIPAcopAngle.at(t).Fill(Acop_DPIP,Wspin);
	    PVIPAcopAngle.at(t).Fill(Acop_PVIP,Wspin);
	  } //pions
      } //channel

    //PULLS
    RefX.at(t).Fill((Muon_ref.X()-Muon_refTruth.X())/Muon_refTruth.X());
    RefY.at(t).Fill((Muon_ref.X()-Muon_refTruth.Y())/Muon_refTruth.Y());
    RefZ.at(t).Fill((Muon_ref.X()-Muon_refTruth.Z())/Muon_refTruth.Z());

    /*if(mixedTausAreOk){
      SVfitMTTdR1.at(t).Fill(TauMuHad.DeltaR(TauMuSVFit),w);
      SVfitMTTdR2.at(t).Fill(TauHadMTT.DeltaR(TauHadSVFit),w);

      double P1svfit = (TauMuSVFit.P()-TauMuTruth.P())/TauMuTruth.P();
      double P2svfit = (TauHadSVFit.P()-TauHadTruth.P())/TauHadTruth.P();
      double Eta1MTT = (TauMuMTT.Eta()-TauMuTruth.Eta())/TauMuTruth.Eta();
      double Phi1MTT = (TauMuMTT.Phi()-TauMuTruth.Phi())/TauMuTruth.Phi();
      double Eta2MTT = (TauHadMTT.Eta()-TauHadTruth.Eta())/TauHadTruth.Eta();
      double Phi2MTT = (TauHadMTT.Phi()-TauHadTruth.Phi())/TauHadTruth.Phi();
      double DAcopPV = (Acop_PVIP-Acop_PVIPTruth)/Acop_PVIPTruth;
	
      PcorrEtaSVfitMTT1.at(t).Fill(P1svfit,Eta1MTT,w);
      PcorrPhiSVfitMTT1.at(t).Fill(P1svfit,Phi1MTT,w);

      PcorrEtaSVfitMTT2.at(t).Fill(P2svfit,Eta2MTT,w);
      PcorrPhiSVfitMTT2.at(t).Fill(P2svfit,Phi2MTT,w);

      dRandPcorrEta1.at(t).Fill(P1svfit,Eta1MTT,Tau1MTT.DeltaR(Tau1SVFit),w);
      dRandPcorrPhi1.at(t).Fill(P1svfit,Phi1MTT,Tau1MTT.DeltaR(Tau1SVFit),w);

      dRandPcorrEta2.at(t).Fill(P2svfit,Eta2MTT,Tau2MTT.DeltaR(Tau2SVFit),w);
      dRandPcorrPhi2.at(t).Fill(P2svfit,Phi2MTT,Tau2MTT.DeltaR(Tau2SVFit),w);

      PullAcopPV.at(t).Fill((Acop_PV-Acop_PVTruth)/Acop_PVTruth,w);
      dR1vsAcopPV.at(t).Fill(Tau1MTT.DeltaR(Tau1SVFit),DAcopPV,w);
      dR2vsAcopPV.at(t).Fill(Tau2MTT.DeltaR(Tau2SVFit),DAcopPV,w);
      P1vsAcopPV.at(t).Fill(P1svfit,DAcopPV,w);
      P2vsAcopPV.at(t).Fill(P2svfit,DAcopPV,w);
      Phi1vsAcopPV.at(t).Fill(Phi1MTT,DAcopPV,w);
      Eta1vsAcopPV.at(t).Fill(Eta1MTT,DAcopPV,w);
      Phi2vsAcopPV.at(t).Fill(Phi2MTT,DAcopPV,w);
      Eta2vsAcopPV.at(t).Fill(Eta2MTT,DAcopPV,w);

    }*/

    //FastMTT
    if(TauMuMTT!=zeroLV && TauHadMTT!=zeroLV){
      DPhitauHadMTT = (TauHadMTT.Phi()-TauHadTruth.Phi())/TauHadTruth.Phi();
      DEtatauHadMTT = (TauHadMTT.Eta()-TauHadTruth.Eta())/TauHadTruth.Eta();
      DPtauHadMTT = (TauHadMTT.P()-TauHadTruth.P())/TauHadTruth.P();
      DEtauHadMTT = (TauHadMTT.E()-TauHadTruth.E())/TauHadTruth.E();
      DPhitauMuMTT = (TauMuMTT.Phi()-TauMuTruth.Phi())/TauMuTruth.Phi();
      DEtatauMuMTT = (TauMuMTT.Eta()-TauMuTruth.Eta())/TauMuTruth.Eta();
      DPtauMuMTT = (TauMuMTT.P()-TauMuTruth.P())/TauMuTruth.P();
      DEtauMuMTT = (TauMuMTT.E()-TauMuTruth.E())/TauMuTruth.E();
      DeltaPhitauHadMTT.at(t).Fill(DPhitauHadMTT,w);
      DeltaEtatauHadMTT.at(t).Fill(DEtatauHadMTT,w);
      DeltaPtauHadMTT.at(t).Fill(DPtauHadMTT,w);
      DeltaEtauHadMTT.at(t).Fill(DEtauHadMTT,w);
      DeltaPhitauMuMTT.at(t).Fill(DPhitauMuMTT,w);
      DeltaEtatauMuMTT.at(t).Fill(DEtatauMuMTT,w);
      DeltaPtauMuMTT.at(t).Fill(DPtauMuMTT,w);
      DeltaEtauMuMTT.at(t).Fill(DEtauMuMTT,w);
    }
    //SVFit
    if(TauMuSVFit!=zeroLV && TauHadSVFit!=zeroLV){
      DPhitauHadSVFit = (TauHadSVFit.Phi()-TauHadTruth.Phi())/TauHadTruth.Phi();
      DEtatauHadSVFit = (TauHadSVFit.Eta()-TauHadTruth.Eta())/TauHadTruth.Eta();
      DPtauHadSVFit = (TauHadSVFit.P()-TauHadTruth.P())/TauHadTruth.P();
      DEtauHadSVFit = (TauHadSVFit.E()-TauHadTruth.E())/TauHadTruth.E();
      DPhitauMuSVFit = (TauMuSVFit.Phi()-TauMuTruth.Phi())/TauMuTruth.Phi();
      DEtatauMuSVFit = (TauMuSVFit.Eta()-TauMuTruth.Eta())/TauMuTruth.Eta();
      DPtauMuSVFit = (TauMuSVFit.P()-TauMuTruth.P())/TauMuTruth.P();
      DEtauMuSVFit = (TauMuSVFit.E()-TauMuTruth.E())/TauMuTruth.E();
      DeltaPhitauHadSVFit.at(t).Fill(DPhitauHadSVFit,w);
      DeltaEtatauHadSVFit.at(t).Fill(DEtatauHadSVFit,w);
      DeltaPtauHadSVFit.at(t).Fill(DPtauHadSVFit,w);
      DeltaEtauHadSVFit.at(t).Fill(DEtauHadSVFit,w);
      DeltaPhitauMuSVFit.at(t).Fill(DPhitauMuSVFit,w);
      DeltaEtatauMuSVFit.at(t).Fill(DEtatauMuSVFit,w);
      DeltaPtauMuSVFit.at(t).Fill(DPtauMuSVFit,w);
      DeltaEtauMuSVFit.at(t).Fill(DEtauMuSVFit,w);
    }
    //Mixed
    if(TauMuMixed!=zeroLV && TauHadMixed!=zeroLV){
      DPhitauHadMixed = (TauHadMixed.Phi()-TauHadTruth.Phi())/TauHadTruth.Phi();
      DEtatauHadMixed = (TauHadMixed.Eta()-TauHadTruth.Eta())/TauHadTruth.Eta();
      DPtauHadMixed = (TauHadMixed.P()-TauHadTruth.P())/TauHadTruth.P();
      DEtauHadMixed = (TauHadMixed.E()-TauHadTruth.E())/TauHadTruth.E();
      DPhitauMuMixed = (TauMuMixed.Phi()-TauMuTruth.Phi())/TauMuTruth.Phi();
      DEtatauMuMixed = (TauMuMixed.Eta()-TauMuTruth.Eta())/TauMuTruth.Eta();
      DPtauMuMixed = (TauMuMixed.P()-TauMuTruth.P())/TauMuTruth.P();
      DEtauMuMixed = (TauMuMixed.E()-TauMuTruth.E())/TauMuTruth.E();
      DeltaPhitauHadMixed.at(t).Fill(DPhitauHadMixed,w);
      DeltaEtatauHadMixed.at(t).Fill(DEtatauHadMixed,w);
      DeltaPtauHadMixed.at(t).Fill(DPtauHadMixed,w);
      DeltaEtauHadMixed.at(t).Fill(DEtauHadMixed,w);
      DeltaPhitauMuMixed.at(t).Fill(DPhitauMuMixed,w);
      DeltaEtatauMuMixed.at(t).Fill(DEtatauMuMixed,w);
      DeltaPtauMuMixed.at(t).Fill(DPtauMuMixed,w);
      DeltaEtauMuMixed.at(t).Fill(DEtauMuMixed,w);
    }
    //GEF
    /*if(TauplusPairConstraintBS!=zeroLV && TauminusPairConstraintBS!=zeroLV){
      DPhitauHadGEF = (TauplusPairConstraintBS.Phi()-TauPlusTruth.Phi())/TauPlusTruth.Phi();
      DEtatauHadGEF = (TauplusPairConstraintBS.Eta()-TauPlusTruth.Eta())/TauPlusTruth.Eta();
      DPtauHadGEF = (TauplusPairConstraintBS.P()-TauPlusTruth.P())/TauPlusTruth.P();
      DEtauHadGEF = (TauplusPairConstraintBS.E()-TauPlusTruth.E())/TauPlusTruth.E();
      DPhitauMuGEF = (TauminusPairConstraintBS.Phi()-TauMinusTruth.Phi())/TauMinusTruth.Phi();
      DEtatauMuGEF = (TauminusPairConstraintBS.Eta()-TauMinusTruth.Eta())/TauMinusTruth.Eta();
      DPtauMuGEF = (TauminusPairConstraintBS.P()-TauMinusTruth.P())/TauMinusTruth.P();
      DEtauMuGEF = (TauminusPairConstraintBS.E()-TauMinusTruth.E())/TauMinusTruth.E();
      DeltaPhitauHadGEF.at(t).Fill(DPhitauHadGEF,w);
      DeltaEtatauHadGEF.at(t).Fill(DEtatauHadGEF,w);
      DeltaPtauHadGEF.at(t).Fill(DPtauHadGEF,w);
      DeltaEtauHadGEF.at(t).Fill(DEtauHadGEF,w);
      DeltaPhitauMuGEF.at(t).Fill(DPhitauMuGEF,w);
      DeltaEtatauMuGEF.at(t).Fill(DEtatauMuGEF,w);
      DeltaPtauMuGEF.at(t).Fill(DPtauMuGEF,w);
      DeltaEtauMuGEF.at(t).Fill(DEtauMuGEF,w);
    }*/
    if(TauMuGEF !=zeroLV && TauHGEF!=zeroLV){
      DPhitauMuGEF = (TauMuGEF.Phi()-TauMuTruth.Phi())/TauMuTruth.Phi();
      DEtatauMuGEF = (TauMuGEF.Eta()-TauMuTruth.Eta())/TauMuTruth.Eta();
      DPtauMuGEF = (TauMuGEF.P()-TauMuTruth.P())/TauMuTruth.P();
      DEtauMuGEF = (TauMuGEF.E()-TauMuTruth.E())/TauMuTruth.E();
      DPhitauHGEF = (TauHGEF.Phi()-TauHadTruth.Phi())/TauHadTruth.Phi();
      DEtatauHGEF = (TauHGEF.Eta()-TauHadTruth.Eta())/TauHadTruth.Eta();
      DPtauHGEF = (TauHGEF.P()-TauHadTruth.P())/TauHadTruth.P();
      DEtauHGEF = (TauHGEF.E()-TauHadTruth.E())/TauHadTruth.E();
      DeltaPhitauMuGEF.at(t).Fill(DPhitauMuGEF,w);
      DeltaEtatauMuGEF.at(t).Fill(DEtatauMuGEF,w);
      DeltaPtauMuGEF.at(t).Fill(DPtauMuGEF,w);
      DeltaEtauMuGEF.at(t).Fill(DEtauMuGEF,w);
      DeltaPhitauHGEF.at(t).Fill(DPhitauHGEF,w);
      DeltaEtatauHGEF.at(t).Fill(DEtatauHGEF,w);
      DeltaPtauHGEF.at(t).Fill(DPtauHGEF,w);
      DeltaEtauHGEF.at(t).Fill(DEtauHGEF,w);
    }
  } //isnan
} //do event

//  This is a function if you want to do something after the event loop
void HCPMuTau::Finish() {

  if(mode == RECONSTRUCT) {
    SkimConfig SC;
    SC.ApplySkimEfficiency(types,Npassed, Npassed_noweight);

    double norm=1.;

    for(unsigned i=0;i<CrossSectionandAcceptance.size();i++){
      if(CrossSectionandAcceptance.at(i)>0 || HConfig.GetID(i)==36 || HConfig.GetID(i)==20 || HConfig.GetID(i)==23 || HConfig.GetID(i)==30 || HConfig.GetID(i)==33){
        if(CrossSectionandAcceptance.at(i)>0)norm= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
        else norm=1.;

        cout<<"Soustraction du QCD: "<<endl;
        std::cout << "end " << std::endl;
      }
    }
  }
  Selection::Finish(Channel,CPstate);
}
