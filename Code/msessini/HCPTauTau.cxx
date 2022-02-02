#include "HCPTauTau.h"
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

HCPTauTau::HCPTauTau(TString Name_, TString id_, char* Channel_, char* CPstate_):
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

HCPTauTau::~HCPTauTau(){
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

void  HCPTauTau::Configure(){
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
  DeltaPhitau1MTT=HConfig.GetTH1D(Name+"_DeltaPhitau1MTT","DeltaPhitau1MTT",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatau1MTT=HConfig.GetTH1D(Name+"_DeltaEtatau1MTT","DeltaEtatau1MTT",50,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtau1MTT=HConfig.GetTH1D(Name+"_DeltaPtau1MTT","DeltaPtau1MTT",50,-1,1,"#Delta P","Events");
  DeltaEtau1MTT=HConfig.GetTH1D(Name+"_DeltaEtau1MTT","DeltaEtau1MTT",50,-1,1,"#Delta E","Events");
  DeltaPhitau2MTT=HConfig.GetTH1D(Name+"_DeltaPhitau2MTT","DeltaPhitau2MTT",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatau2MTT=HConfig.GetTH1D(Name+"_DeltaEtatau2MTT","DeltaEtatau2MTT",50,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtau2MTT=HConfig.GetTH1D(Name+"_DeltaPtau2MTT","DeltaPtau2MTT",50,-1,1,"#Delta P","Events");
  DeltaEtau2MTT=HConfig.GetTH1D(Name+"_DeltaEtau2MTT","DeltaEtau2MTT",50,-1,1,"#Delta E","Events");
  //svfit
  DeltaPhitau1SVFit=HConfig.GetTH1D(Name+"_DeltaPhitau1SVFit","DeltaPhitau1SVFit",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatau1SVFit=HConfig.GetTH1D(Name+"_DeltaEtatau1SVFit","DeltaEtatau1SVFit",50,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtau1SVFit=HConfig.GetTH1D(Name+"_DeltaPtau1SVFit","DeltaPtau1SVFit",50,-1,1,"#Delta P","Events");
  DeltaEtau1SVFit=HConfig.GetTH1D(Name+"_DeltaEtau1SVFit","DeltaEtau1SVFit",50,-1,1,"#Delta E","Events");
  DeltaPhitau2SVFit=HConfig.GetTH1D(Name+"_DeltaPhitau2SVFit","DeltaPhitau2SVFit",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatau2SVFit=HConfig.GetTH1D(Name+"_DeltaEtatau2SVFit","DeltaEtatau2SVFit",100,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtau2SVFit=HConfig.GetTH1D(Name+"_DeltaPtau2SVFit","DeltaPtau2SVFit",50,-1,1,"#Delta P","Events");
  DeltaEtau2SVFit=HConfig.GetTH1D(Name+"_DeltaEtau2SVFit","DeltaEtau2SVFit",50,-1,1,"#Delta E","Events");
  //mixed
  DeltaPhitau1Mixed=HConfig.GetTH1D(Name+"_DeltaPhitau1Mixed","DeltaPhitau1Mixed",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatau1Mixed=HConfig.GetTH1D(Name+"_DeltaEtatau1Mixed","DeltaEtatau1Mixed",50,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtau1Mixed=HConfig.GetTH1D(Name+"_DeltaPtau1Mixed","DeltaPtau1Mixed",50,-1,1,"#Delta P","Events");
  DeltaEtau1Mixed=HConfig.GetTH1D(Name+"_DeltaEtau1Mixed","DeltaEtau1Mixed",50,-1,1,"#Delta E","Events");
  DeltaPhitau2Mixed=HConfig.GetTH1D(Name+"_DeltaPhitau2Mixed","DeltaPhitau2Mixed",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatau2Mixed=HConfig.GetTH1D(Name+"_DeltaEtatau2Mixed","DeltaEtatau2Mixed",50,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtau2Mixed=HConfig.GetTH1D(Name+"_DeltaPtau2Mixed","DeltaPtau2Mixed",50,-1,1,"#Delta P","Events");
  DeltaEtau2Mixed=HConfig.GetTH1D(Name+"_DeltaEtau2Mixed","DeltaEtau2Mixed",50,-1,1,"#Delta E","Events");
  //GEF
  DeltaPhitau1GEF=HConfig.GetTH1D(Name+"_DeltaPhitau1GEF","DeltaPhitau1GEF",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatau1GEF=HConfig.GetTH1D(Name+"_DeltaEtatau1GEF","DeltaEtatau1GEF",50,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtau1GEF=HConfig.GetTH1D(Name+"_DeltaPtau1GEF","DeltaPtau1GEF",50,-1,1,"#Delta P","Events");
  DeltaEtau1GEF=HConfig.GetTH1D(Name+"_DeltaEtau1GEF","DeltaEtau1GEF",50,-1,1,"#Delta E","Events");
  DeltaPhitau2GEF=HConfig.GetTH1D(Name+"_DeltaPhitau2GEF","DeltaPhitau2GEF",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatau2GEF=HConfig.GetTH1D(Name+"_DeltaEtatau2GEF","DeltaEtatau2GEF",50,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtau2GEF=HConfig.GetTH1D(Name+"_DeltaPtau2GEF","DeltaPtau2GEF",50,-1,1,"#Delta P","Events");
  DeltaEtau2GEF=HConfig.GetTH1D(Name+"_DeltaEtau2GEF","DeltaEtau2GEF",50,-1,1,"#Delta E","Events");
  //
  DeltaPhitauPiGEF=HConfig.GetTH1D(Name+"_DeltaPhitauPiGEF","DeltaPhitauPiGEF",50,-0.2,0.2,"#Delta#phi","Events");
  DeltaEtatauPiGEF=HConfig.GetTH1D(Name+"_DeltaEtatauPiGEF","DeltaEtatauPiGEF",50,-0.2,0.2,"#Delta#eta","Events");
  DeltaPtauPiGEF=HConfig.GetTH1D(Name+"_DeltaPtauPiGEF","DeltaPtauPiGEF",50,-1,1,"#Delta P","Events");
  DeltaEtauPiGEF=HConfig.GetTH1D(Name+"_DeltaEtauPiGEF","DeltaEtauPiGEF",50,-1,1,"#Delta E","Events");
  DeltaPhitauHGEF=HConfig.GetTH1D(Name+"_DeltaPhitauHGEF","DeltaPhitauHGEF",50,-0.1,0.1,"#Delta#phi","Events");
  DeltaEtatauHGEF=HConfig.GetTH1D(Name+"_DeltaEtatauHGEF","DeltaEtatauHGEF",50,-0.1,0.1,"#Delta#eta","Events");
  DeltaPtauHGEF=HConfig.GetTH1D(Name+"_DeltaPtauHGEF","DeltaPtauHGEF",50,-1,1,"#Delta P","Events");
  DeltaEtauHGEF=HConfig.GetTH1D(Name+"_DeltaEtauHGEF","DeltaEtauHGEF",50,-1,1,"#Delta E","Events");

  //REF
  Ref=HConfig.GetTH1D(Name+"_Ref","Ref",50,0,2,"ref","Events");
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

void  HCPTauTau::Store_ExtraDist(){

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
  Extradist1d.push_back(&DeltaPhitau1MTT);
  Extradist1d.push_back(&DeltaEtatau1MTT);
  Extradist1d.push_back(&DeltaPtau1MTT);
  Extradist1d.push_back(&DeltaEtau1MTT);
  Extradist1d.push_back(&DeltaPhitau2MTT);
  Extradist1d.push_back(&DeltaEtatau2MTT);
  Extradist1d.push_back(&DeltaPtau2MTT);
  Extradist1d.push_back(&DeltaEtau2MTT);
  //svfit
  Extradist1d.push_back(&DeltaPhitau1SVFit);
  Extradist1d.push_back(&DeltaEtatau1SVFit);
  Extradist1d.push_back(&DeltaPtau1SVFit);
  Extradist1d.push_back(&DeltaEtau1SVFit);
  Extradist1d.push_back(&DeltaPhitau2SVFit);
  Extradist1d.push_back(&DeltaEtatau2SVFit);
  Extradist1d.push_back(&DeltaPtau2SVFit);
  Extradist1d.push_back(&DeltaEtau2SVFit);
  //mixed
  Extradist1d.push_back(&DeltaPhitau1Mixed);
  Extradist1d.push_back(&DeltaEtatau1Mixed);
  Extradist1d.push_back(&DeltaPtau1Mixed);
  Extradist1d.push_back(&DeltaEtau1Mixed);
  Extradist1d.push_back(&DeltaPhitau2Mixed);
  Extradist1d.push_back(&DeltaEtatau2Mixed);
  Extradist1d.push_back(&DeltaPtau2Mixed);
  Extradist1d.push_back(&DeltaEtau2Mixed);
  //GEF
  Extradist1d.push_back(&DeltaPhitau1GEF);
  Extradist1d.push_back(&DeltaEtatau1GEF);
  Extradist1d.push_back(&DeltaPtau1GEF);
  Extradist1d.push_back(&DeltaEtau1GEF);
  Extradist1d.push_back(&DeltaPhitau2GEF);
  Extradist1d.push_back(&DeltaEtatau2GEF);
  Extradist1d.push_back(&DeltaPtau2GEF);
  Extradist1d.push_back(&DeltaEtau2GEF);
  //
  Extradist1d.push_back(&DeltaPhitauPiGEF);
  Extradist1d.push_back(&DeltaEtatauPiGEF);
  Extradist1d.push_back(&DeltaPtauPiGEF);
  Extradist1d.push_back(&DeltaEtauPiGEF);
  Extradist1d.push_back(&DeltaPhitauHGEF);
  Extradist1d.push_back(&DeltaEtatauHGEF);
  Extradist1d.push_back(&DeltaPtauHGEF);
  Extradist1d.push_back(&DeltaEtauHGEF);
  //ref
  Extradist1d.push_back(&Ref);
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

void  HCPTauTau::doEvent()  { //  Method called on every event

  unsigned int t;                // sample type, you may manage in your further analysis, if needed
  int id(Ntp->GetMCID());  //read event ID of a sample

  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}  //  gives a warning if list of samples in Histo.txt  and SkimSummary.log do not coincide

  bool isEmbed=(id==36);
  bool trig=false;
  std::vector<int> TauIndex ;
  std::vector<int> TriggerIndexVector ;
  std::vector<TString>  MatchedTriggerNames;
  value.at(Trigger)=0;

  int Tau1= -1;
  int Tau2= -1;
  bool GenMatchSelection=false;

  Tau1=Ntp->tau1IndexVect(0);
  Tau2=Ntp->tau2IndexVect(0);
  int j=0;
  int GenMatch1=6, GenMatch2=6;
  
  if(!Ntp->isData()|| isEmbed){GenMatch1=Ntp->gen_match_1(0);GenMatch2=Ntp->gen_match_2(0);}
  else{GenMatch1=6;GenMatch2=6;}
  if(id==30 || id==33 ||id==203)GenMatchSelection=(!(GenMatch1==6 || GenMatch2==6) && !(GenMatch1==5 && GenMatch2==5)); //ZL
  else if(id==36)GenMatchSelection=(GenMatch1==5 && GenMatch2==5);
  else if(id==23 && id==20)GenMatchSelection=false; //WTaunu
  else if(id==201|| id==202)GenMatchSelection=!((GenMatch1==5 && GenMatch2==6)||(GenMatch2==5 && GenMatch1==6));
  else if(id==71||id==72||id==73||id==74||id==47||id==48||id==49||id==50||id==51||id==52||id==53||id==54||id==55||id==56||id==57||id==58)GenMatchSelection=(!(GenMatch1==6 || GenMatch2==6) &&  !(GenMatch1==5 && GenMatch2==5)); // VVT
  else if(id==70 ||id==701 ||id==702 ||id==703 )GenMatchSelection=(!(GenMatch1==6 || GenMatch2==6) &&  !(GenMatch1==5 && GenMatch2==5));
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

    wTrgSF1=Ntp->TriggerSF(Tau1,GenMatch1,TES,"Nom");
    wTrgSF2=Ntp->TriggerSF(Tau2,GenMatch2,TES,"Nom");

    w*=wTrgSF1;
    w*=wTrgSF2;

    wIDSF1=Ntp->IDSF(Tau1,GenMatch1,TES);
    wIDSF2=Ntp->IDSF(Tau2,GenMatch2,TES);
    if(!isEmbed){
      wIDSFAntiE1=Ntp->IDSF(Tau1,GenMatch1,TES,"ele");
      wIDSFAntiE2=Ntp->IDSF(Tau2,GenMatch2,TES,"ele");
      wIDSFAntiMu1=Ntp->IDSF(Tau1,GenMatch1,TES,"mu");
      wIDSFAntiMu2=Ntp->IDSF(Tau2,GenMatch2,TES,"mu");
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

  TLorentzVector Tau1P4;
  TLorentzVector Tau2P4;
  //if(passAllBut(exclude_cuts)){
  Tau1P4 = Ntp->Daughters_P4(Tau1);
  Tau2P4 = Ntp->Daughters_P4(Tau2);
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

  double DEtatau1MTT=-99, DEtatau2MTT=-99, DPhitau1MTT=-99, DPhitau2MTT=-99;
  double DEtau1MTT=-99, DEtau2MTT=-99, DPtau1MTT=-99, DPtau2MTT=-99;
  double DEtatau1SVFit=-99, DEtatau2SVFit=-99, DPhitau1SVFit=-99, DPhitau2SVFit=-99;
  double DEtau1SVFit=-99, DEtau2SVFit=-99, DPtau1SVFit=-99, DPtau2SVFit=-99;
  double DEtatau1Mixed=-99, DEtatau2Mixed=-99, DPhitau1Mixed=-99, DPhitau2Mixed=-99;
  double DEtau1Mixed=-99, DEtau2Mixed=-99, DPtau1Mixed=-99, DPtau2Mixed=-99;
  double DEtatau1GEF=-99, DEtatau2GEF=-99, DPhitau1GEF=-99, DPhitau2GEF=-99;
  double DEtau1GEF=-99, DEtau2GEF=-99, DPtau1GEF=-99, DPtau2GEF=-99;
  double DEtatauPiGEF=-99, DEtatauHGEF=-99, DPhitauPiGEF=-99, DPhitauHGEF=-99;
  double DEtauPiGEF=-99, DEtauHGEF=-99, DPtauPiGEF=-99, DPtauHGEF=-99;

  vector<TLorentzVector> HadPionsTruth_minus;
  vector<TLorentzVector> HadPionsTruth_plus;
  vector<double> HadPionsChargeTruth_minus;
  vector<double> HadPionsChargeTruth_plus;
  vector<TLorentzVector> tauandprodTruthminus;
  vector<TLorentzVector> tauandprodTruthplus;

  vector<TLorentzVector> HadRefitPions_minus;
  vector<double> HadRefitPionsCharge_minus;
  vector<TLorentzVector> HadRefitPions_plus;
  vector<double> HadRefitPionsCharge_plus;

  TVector3 PionMinus_ref, PionPlus_ref, PionMinus_refTruth, PionPlus_refTruth;

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
  double Acop_PVGEF = -99, Acop_PV = -99, Acop_DP = -99, Acop_DPIP = -99, Acop_PVIP = -99, Acop_IP = -99;;
  double Acop_PVTruth = -99, Acop_DPTruth = -99, Acop_DPIPTruth = -99, Acop_PVIPTruth = -99, Acop_IPTruth = -99;

  TLorentzVector zeroLV(0,0,0,0);
  std::vector<TLorentzVector> VectZeroLV;
  VectZeroLV.push_back(zeroLV);
  VectZeroLV.push_back(zeroLV);
  VectZeroLV.push_back(zeroLV);

  TLorentzVector TauPlusSVFit, TauMinusSVFit;
  TLorentzVector TauPlusMTT, TauMinusMTT;
  TLorentzVector TauPlusMixed, TauMinusMixed;
  TLorentzVector Tau2SVFit, Tau1SVFit;
  TLorentzVector Tau2MTT, Tau1MTT;
  TLorentzVector Tau2Mixed, Tau1Mixed;
  TLorentzVector TauHGEF, TauPiGEF;
  TLorentzVector Tauplusvis;
  TLorentzVector Tauminusvis;
  TLorentzVector TauPlusTruth, TauMinusTruth;
  TLorentzVector Tau1Truth, Tau2Truth;
  TLorentzVector TauPiTruth, TauHTruth;
  unsigned int Tauplustruth=0;
  unsigned int Tauminustruth=0;
  unsigned int Tauplus=0;
  unsigned int Tauminus=0;
  string CHANNEL = string(Channel);
  cout<<CHANNEL<<" "<<CPstate<<endl;
  //
  bool a1minus = false, a1plus = false, a1a1 = false;
  bool rhominus = false, rhoplus = false, rhorho =false;
  bool pionminus = false, pionplus = false, pionpion =false;
  bool a1rho = false, a1pion = false, rhopion = false;
  bool minusplus = false, plusminus = false;
  bool charge = false, channel = false;
  //
  bool a1minustruth = false, a1plustruth = false, a1a1truth = false;
  bool rhominustruth = false, rhoplustruth = false, rhorhotruth =false;
  bool pionminustruth = false, pionplustruth = false, pionpiontruth =false;
  bool a1rhotruth = false, a1piontruth = false, rhopiontruth = false;
  bool minusplustruth = false, plusminustruth = false;
  bool chargetruth = false, channeltruth = false;
  //
  bool recoPionsAreOk = false, genPionsAreOk = false;
  bool recoTausAreOk = false, GEFTausAreOk = false, mixedTausAreOk = false, genTausAreOk = false;

  //GENERATED LEVEL
  if(Ntp->MCTau_charge(0)<0 && Ntp->MCTau_charge(1)>0){
    Tauminustruth = 0;
    Tauplustruth = 1;
    chargetruth = true;
  }
  else if(Ntp->MCTau_charge(1)<0 && Ntp->MCTau_charge(0)>0){
    Tauminustruth = 1;
    Tauplustruth = 0;
    chargetruth = true;
  }
  TauMinusTruth = Ntp->MCTau_p4(Tauminustruth);
  TauPlusTruth = Ntp->MCTau_p4(Tauplustruth);

  //3 prongs decay
  if(Ntp->MCTau_JAK(Tauminustruth) == 5) a1minustruth=true;
  if(Ntp->MCTau_JAK(Tauplustruth) == 5) a1plustruth=true;
  //1prong + pi0 decay
  if(Ntp->MCTau_JAK(Tauminustruth) == 4) rhominustruth=true;
  if(Ntp->MCTau_JAK(Tauplustruth) == 4) rhoplustruth=true;
  //1 prong decay
  if(Ntp->MCTau_JAK(Tauminustruth) == 3) pionminustruth=true;
  if(Ntp->MCTau_JAK(Tauplustruth) == 3) pionplustruth=true;

  if(a1minustruth && a1plustruth && CHANNEL == "A1A1") a1a1truth=true;
  if(rhominustruth && rhoplustruth && CHANNEL == "RHORHO") rhorhotruth=true;
  if(pionminustruth && pionplustruth && CHANNEL == "PIONPION") pionpiontruth=true;
  if(((a1minustruth && rhoplustruth) || (a1plustruth && rhominustruth)) && CHANNEL == "A1RHO"){
    a1rhotruth=true;
    if(a1minustruth && rhoplustruth) minusplustruth=true;
    if(a1plustruth && rhominustruth) plusminustruth=true;
  }
  if(((a1minustruth && pionplustruth) || (a1plustruth && pionminustruth)) && CHANNEL == "A1PION"){
    a1piontruth=true;
    if(a1minustruth && pionplustruth) minusplustruth=true;
    if(a1plustruth && pionminustruth) plusminustruth=true;
  }
  if(((rhominustruth && pionplustruth) || (rhoplustruth && pionminustruth)) && CHANNEL == "RHOPION"){
    rhopiontruth=true;
    if(rhominustruth && pionplustruth) minusplustruth=true;
    if(rhoplustruth && pionminustruth) plusminustruth=true;
  }
  if((a1a1truth || rhorhotruth || pionpiontruth || a1rhotruth || a1piontruth || rhopiontruth) && chargetruth) channeltruth = true;
  if(channeltruth){
    if(a1a1truth){
      HadPionsTruth_minus = Ntp->GetTruthPionsFromA1(Tauminustruth);
      HadPionsChargeTruth_minus.push_back(1);
      HadPionsChargeTruth_minus.push_back(-1);
      HadPionsChargeTruth_minus.push_back(-1);

      HadPionsTruth_plus = Ntp->GetTruthPionsFromA1(Tauplustruth);
      HadPionsChargeTruth_plus.push_back(-1);
      HadPionsChargeTruth_plus.push_back(1);
      HadPionsChargeTruth_plus.push_back(1);
    }
    if(rhorhotruth){
      HadPionsTruth_minus = Ntp->GetTruthPionsFromRho(Tauminustruth);
      HadPionsChargeTruth_minus.push_back(-1);
      HadPionsChargeTruth_minus.push_back(0);

      HadPionsTruth_plus = Ntp->GetTruthPionsFromRho(Tauplustruth);
      HadPionsChargeTruth_plus.push_back(1);
      HadPionsChargeTruth_plus.push_back(0);
    }
    if(pionpiontruth){
    cout<<"yes"<<endl;
      HadPionsTruth_minus.push_back(Ntp->GetTruthTauProductLV(3,211,Tauminustruth));
      HadPionsChargeTruth_minus.push_back(-1);
      PionMinus_refTruth = Ntp->MCTauandProd_Vertex(Tauminustruth,2) - Ntp->MCTauandProd_Vertex(Tauminustruth,0);

      HadPionsTruth_plus.push_back(Ntp->GetTruthTauProductLV(3,211,Tauplustruth));
      HadPionsChargeTruth_plus.push_back(1);
      PionPlus_refTruth = Ntp->MCTauandProd_Vertex(Tauplustruth,2) - Ntp->MCTauandProd_Vertex(Tauplustruth,0);
    }
    if(a1a1truth || rhorhotruth || pionpiontruth){
      Tau1Truth = TauMinusTruth; Tau2Truth = TauPlusTruth;
    }
    if(plusminustruth){
      if(a1rhotruth){
        HadPionsTruth_plus = Ntp->GetTruthPionsFromA1(Tauplustruth);
        HadPionsChargeTruth_plus.push_back(-1);
        HadPionsChargeTruth_plus.push_back(1);
        HadPionsChargeTruth_plus.push_back(1);

        HadPionsTruth_minus = Ntp->GetTruthPionsFromRho(Tauminustruth);
        HadPionsChargeTruth_minus.push_back(-1);
        HadPionsChargeTruth_minus.push_back(0);
      }
      if(a1piontruth){
        HadPionsTruth_plus = Ntp->GetTruthPionsFromA1(Tauplustruth);
        HadPionsChargeTruth_plus.push_back(-1);
        HadPionsChargeTruth_plus.push_back(1);
        HadPionsChargeTruth_plus.push_back(1);

        HadPionsTruth_minus.push_back(Ntp->GetTruthTauProductLV(3,211,Tauminustruth));
        HadPionsChargeTruth_minus.push_back(-1);
        PionMinus_refTruth = Ntp->MCTauandProd_Vertex(Tauminustruth,2) - Ntp->MCTauandProd_Vertex(Tauminustruth,0);
      }
      if(rhopiontruth){
        HadPionsTruth_plus = Ntp->GetTruthPionsFromRho(Tauplustruth);
        HadPionsChargeTruth_plus.push_back(1);
        HadPionsChargeTruth_plus.push_back(0);

        HadPionsTruth_minus.push_back(Ntp->GetTruthTauProductLV(3,211,Tauminustruth));
        HadPionsChargeTruth_minus.push_back(-1);
        PionPlus_refTruth = Ntp->MCTauandProd_Vertex(Tauminustruth,2) - Ntp->MCTauandProd_Vertex(Tauminustruth,0);
      }
      Tau1Truth = TauPlusTruth; Tau2Truth = TauMinusTruth;
      TauHTruth = TauPlusTruth; TauPiTruth = TauMinusTruth;
    }
    if(minusplustruth){
      if(a1rhotruth){
        HadPionsTruth_minus = Ntp->GetTruthPionsFromA1(Tauminustruth);
        HadPionsChargeTruth_minus.push_back(1);
        HadPionsChargeTruth_minus.push_back(-1);
        HadPionsChargeTruth_minus.push_back(-1);

        HadPionsTruth_plus = Ntp->GetTruthPionsFromRho(Tauplustruth);
        HadPionsChargeTruth_plus.push_back(1);
        HadPionsChargeTruth_plus.push_back(0);
      }
      if(a1piontruth){
        HadPionsTruth_minus = Ntp->GetTruthPionsFromA1(Tauminustruth);
        HadPionsChargeTruth_minus.push_back(1);
        HadPionsChargeTruth_minus.push_back(-1);
        HadPionsChargeTruth_minus.push_back(-1);

        HadPionsTruth_plus.push_back(Ntp->GetTruthTauProductLV(3,211,Tauplustruth));
        HadPionsChargeTruth_plus.push_back(1);
        PionPlus_refTruth = Ntp->MCTauandProd_Vertex(Tauplustruth,2) - Ntp->MCTauandProd_Vertex(Tauplustruth,0);
      }
      if(rhopiontruth){
        HadPionsTruth_minus = Ntp->GetTruthPionsFromRho(Tauminustruth);
        HadPionsChargeTruth_minus.push_back(-1);
        HadPionsChargeTruth_minus.push_back(0);

        HadPionsTruth_plus.push_back(Ntp->GetTruthTauProductLV(3,211,Tauplustruth));
        HadPionsChargeTruth_plus.push_back(1);
        PionPlus_refTruth = Ntp->MCTauandProd_Vertex(Tauplustruth,2) - Ntp->MCTauandProd_Vertex(Tauplustruth,0);
      }
      Tau1Truth = TauMinusTruth; Tau2Truth = TauPlusTruth;
      TauHTruth = TauMinusTruth; TauPiTruth = TauPlusTruth;
    }
  }
  //RECONSTRUCTED LEVEL
  if(Ntp->Daughters_charge(Tau1)>0 && Ntp->Daughters_charge(Tau2)<0 && Tau1P4.Pt() > 40 && Tau2P4.Pt() > 40)
    {
      Tauplusvis=Tau1P4;
      Tauminusvis=Tau2P4;
      Tauplus=Tau1;
      Tauminus=Tau2;
      charge = true;
    }
  else if(Ntp->Daughters_charge(Tau2)>0 && Ntp->Daughters_charge(Tau1)<0 && Tau1P4.Pt() > 40 && Tau2P4.Pt() > 40)
    {
      Tauplusvis=Tau2P4;
      Tauminusvis=Tau1P4;
      Tauplus=Tau2;
      Tauminus=Tau1;
      charge = true;
    }

  if(Ntp->MVADM2017(Tauminus) == 10 && Ntp->PFTau_hassecondaryVertex(Tauminus) && Ntp->PFtauHasThreePions(Tauminus)) a1minus=true;
  if(Ntp->MVADM2017(Tauplus) == 10 && Ntp->PFTau_hassecondaryVertex(Tauplus) && Ntp->PFtauHasThreePions(Tauplus)) a1plus=true;

  if(Ntp->MVADM2017(Tauminus) == 1) rhominus=true;
  if(Ntp->MVADM2017(Tauplus) == 1) rhoplus=true;

  if(Ntp->MVADM2017(Tauminus) == 0) pionminus=true;
  if(Ntp->MVADM2017(Tauplus) == 0) pionplus=true;

  if(a1minus && a1plus && CHANNEL == "A1A1") a1a1=true;
  if(rhominus && rhoplus && CHANNEL == "RHORHO") rhorho=true;
  if(pionminus && pionplus && CHANNEL == "PIONPION") pionpion=true;
  if(((a1minus && rhoplus) || (a1plus && rhominus)) && CHANNEL == "A1RHO"){
    a1rho=true;
    if(a1minus && rhoplus) minusplus=true;
    if(a1plus && rhominus) plusminus=true;
  }
  if(((a1minus && pionplus) || (a1plus && pionminus)) && CHANNEL == "A1PION"){
    a1pion=true;
    if(a1minus && pionplus) minusplus=true;
    if(a1plus && pionminus) plusminus=true;
  }
  if(((rhominus && pionplus) || (rhoplus && pionminus)) && CHANNEL == "RHOPION"){
    rhopion=true;
    if(rhominus && pionplus) minusplus=true;
    if(rhoplus && pionminus) plusminus=true;
  }
  if((a1a1 || rhorho || pionpion || a1rho || a1pion || rhopion) && charge) channel = true;
  ClassicSVfit svfitAlgo1(0);
  double higgsmass;
  if(channel){
    FastMTT FastMTTAlgo;
    // // //---------  svfit ---------------------
    std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
    classic_svFit::MeasuredTauLepton lep1(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Tau1P4.Pt(), Tau1P4.Eta(),  Tau1P4.Phi(), Tau1P4.M(),Ntp->MVADM2017(Tau1));
    classic_svFit::MeasuredTauLepton lep2(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Tau2P4.Pt(), Tau2P4.Eta(),  Tau2P4.Phi(), Tau2P4.M(),Ntp->MVADM2017(Tau2));

    measuredTauLeptons.push_back(lep1);
    measuredTauLeptons.push_back(lep2);
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

      if(Ntp->Daughters_charge(Tau1)>0)
	{
	  TauPlusSVFit.SetPxPyPzE(tau1P4.x(),tau1P4.y(),tau1P4.z(),tau1P4.t());
	  TauMinusSVFit.SetPxPyPzE(tau2P4.x(),tau2P4.y(),tau2P4.z(),tau2P4.t());
	  TauPlusMTT.SetPxPyPzE(tau1P4mtt.x(),tau1P4mtt.y(),tau1P4mtt.z(),tau1P4mtt.t());
	  TauMinusMTT.SetPxPyPzE(tau2P4mtt.x(),tau2P4mtt.y(),tau2P4mtt.z(),tau2P4mtt.t());

	}
      else
	{
	  TauPlusSVFit.SetPxPyPzE(tau2P4.x(),tau2P4.y(),tau2P4.z(),tau2P4.t());
	  TauMinusSVFit.SetPxPyPzE(tau1P4.x(),tau1P4.y(),tau1P4.z(),tau1P4.t());
	  TauPlusMTT.SetPxPyPzE(tau2P4mtt.x(),tau2P4mtt.y(),tau2P4mtt.z(),tau2P4mtt.t());
	  TauMinusMTT.SetPxPyPzE(tau1P4mtt.x(),tau1P4mtt.y(),tau1P4mtt.z(),tau1P4mtt.t());
	}
      //MIXED TAUS
      if(TauMinusSVFit!=zeroLV && TauPlusSVFit!=zeroLV && TauMinusMTT!=zeroLV && TauPlusMTT!=zeroLV && TauMinusSVFit!=TauPlusSVFit && TauMinusMTT!=TauPlusMTT){
        //if(TauMinusSVFit.DeltaR(TauMinusMTT) < 0.03 && TauPlusSVFit.DeltaR(TauPlusMTT) < 0.03){
          TauPlusMixed.SetVect(TauPlusSVFit.P()*TauPlusMTT.Vect().Unit());
          TauPlusMixed.SetE(TauPlusSVFit.E());
          TauMinusMixed.SetVect(TauMinusSVFit.P()*TauMinusMTT.Vect().Unit());
          TauMinusMixed.SetE(TauMinusSVFit.E());
        //}	
	/*else{
	  TauPlusMixed = TauPlusMTT;
	  TauMinusMixed = TauMinusMTT;
	}*/
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
    m_vis_=(Tau1P4+Tau2P4).M();
    met_=PUPPIMET;

    pt_1_=Tau1P4.Pt();
    //pt_2_=Tau2P4.Pt();
    pt_tt_=(Tau1P4+Tau2P4+TLorentzVector(PUPPImetCorr_px,PUPPImetCorr_py,0,0)).Pt();
    pt_vis_=(Tau1P4+Tau2P4).Pt();
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

    boost::hash_combine(hash, Ntp->LeptonHash(Tauminus));
    boost::hash_combine(hash, Ntp->LeptonHash(Tauplus));
    hashes.push_back(hash);
    hash = 0;
    boost::hash_combine(hash, Ntp->LeptonHash(Tauplus));
    boost::hash_combine(hash, Ntp->LeptonHash(Tauminus));
    hashes.push_back(hash);

    if(a1a1){
      //GEF
      TauminusSecondaryVertex = Ntp->PFTau_secondaryVertex_pos(Tauminus);
      TauplusSecondaryVertex = Ntp->PFTau_secondaryVertex_pos(Tauplus);
      a1LVRefit_Tauminus = Ntp->PFTauRefit_PionsP4(Tauminus,0) + Ntp->PFTauRefit_PionsP4(Tauminus,1) + Ntp->PFTauRefit_PionsP4(Tauminus,2);
      a1LVRefit_Tauplus = Ntp->PFTauRefit_PionsP4(Tauplus,0) + Ntp->PFTauRefit_PionsP4(Tauplus,1) + Ntp->PFTauRefit_PionsP4(Tauplus,2);
      tauBSPrimaryVertex = GetRefittedPV(hashes, tauPrimaryVertex, PVRefitBS_X_temp ,PVRefitBS_Y_temp,PVRefitBS_Z_temp ,VertexHashBS1_temp, VertexHashBS2_temp,isRefitBS);
      solutionsBS=tauPairMomentumSolutions(TauminusSecondaryVertex-tauBSPrimaryVertex, a1LVRefit_Tauminus, a1LV_Tauminus, isMinusReal, TauplusSecondaryVertex-tauBSPrimaryVertex, a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal,isRefitBS);
      TauminusPairConstraintBS=solutionsBS.at(3);
      TauplusPairConstraintBS=solutionsBS.at(7);
      HPairConstraint= TauplusPairConstraint+TauminusPairConstraint;

      //
      HadRefitPions_minus.push_back(Ntp->PFTauRefit_PionsP4(Tauminus,0));
      HadRefitPions_minus.push_back(Ntp->PFTauRefit_PionsP4(Tauminus,1));
      HadRefitPions_minus.push_back(Ntp->PFTauRefit_PionsP4(Tauminus,2));
      HadRefitPionsCharge_minus.push_back(Ntp->PFTauRefit_PionsCharge(Tauminus, 0));
      HadRefitPionsCharge_minus.push_back(Ntp->PFTauRefit_PionsCharge(Tauminus, 1));
      HadRefitPionsCharge_minus.push_back(Ntp->PFTauRefit_PionsCharge(Tauminus, 2));

      HadRefitPions_plus.push_back(Ntp->PFTauRefit_PionsP4(Tauplus,0));
      HadRefitPions_plus.push_back(Ntp->PFTauRefit_PionsP4(Tauplus,1));
      HadRefitPions_plus.push_back(Ntp->PFTauRefit_PionsP4(Tauplus,2));
      HadRefitPionsCharge_plus.push_back(Ntp->PFTauRefit_PionsCharge(Tauplus, 0));
      HadRefitPionsCharge_plus.push_back(Ntp->PFTauRefit_PionsCharge(Tauplus, 1));
      HadRefitPionsCharge_plus.push_back(Ntp->PFTauRefit_PionsCharge(Tauplus, 2));
    }
    if(rhorho){
      HadRefitPions_minus.push_back(Ntp->ChargedDaughters_P4(Tauminus));
      HadRefitPions_minus.push_back(Tauminusvis - Ntp->ChargedDaughters_P4(Tauminus));
      HadRefitPionsCharge_minus.push_back(0);
      HadRefitPionsCharge_minus.push_back(-1);

      HadRefitPions_plus.push_back(Ntp->ChargedDaughters_P4(Tauplus));
      HadRefitPions_plus.push_back(Tauplusvis - Ntp->ChargedDaughters_P4(Tauplus));
      HadRefitPionsCharge_plus.push_back(0);
      HadRefitPionsCharge_plus.push_back(1);
    }
    if(pionpion){
      HadRefitPions_minus.push_back(Ntp->ChargedDaughters_P4(Tauminus));
      HadRefitPionsCharge_minus.push_back(-1);
      PionMinus_ref = Ntp->Daughters_vertex(Tauminus) - tauBSPrimaryVertex;

      HadRefitPions_plus.push_back(Ntp->ChargedDaughters_P4(Tauplus));
      HadRefitPionsCharge_plus.push_back(1);
      PionPlus_ref = Ntp->Daughters_vertex(Tauplus) - tauBSPrimaryVertex;
    }
    if(a1a1 || rhorho || pionpion){
      Tau1MTT = TauMinusMTT; Tau2MTT = TauPlusMTT;
      Tau1SVFit = TauMinusSVFit; Tau2SVFit = TauPlusSVFit;
      Tau1Mixed = TauMinusMixed; Tau2Mixed = TauPlusMixed;
    }
    if(plusminus){
      if(a1rho){
	HadRefitPions_plus.push_back(Ntp->PFTauRefit_PionsP4(Tauplus,0));
	HadRefitPions_plus.push_back(Ntp->PFTauRefit_PionsP4(Tauplus,1));
	HadRefitPions_plus.push_back(Ntp->PFTauRefit_PionsP4(Tauplus,2));
	HadRefitPionsCharge_plus.push_back(Ntp->PFTauRefit_PionsCharge(Tauplus, 0));
	HadRefitPionsCharge_plus.push_back(Ntp->PFTauRefit_PionsCharge(Tauplus, 1));
	HadRefitPionsCharge_plus.push_back(Ntp->PFTauRefit_PionsCharge(Tauplus, 2));

	HadRefitPions_minus.push_back(Ntp->ChargedDaughters_P4(Tauminus));
	HadRefitPions_minus.push_back(Tauminusvis - Ntp->ChargedDaughters_P4(Tauminus));
	HadRefitPionsCharge_minus.push_back(0);
	HadRefitPionsCharge_minus.push_back(-1);
      }
      if(a1pion){
        //GEF
        if(Ntp->PFTau_TrackHasMomentum(Tauminus) == true && Ntp->isPVCovAvailable() == true){
          TMatrixTSym<double> PV_cov = Ntp->PFTau_TIP_primaryVertex_cov();
          LorentzVectorParticle A1 = Ntp->PFTau_a1_lvp(Tauplus);
          TrackParticle PION = Ntp->PFTau_Track(Tauminus);
          GlobalEventFit GEF(PION,A1,METInput,tauBSPrimaryVertex,PV_cov);
          GEF.setMassConstraint(125.3);
          GEF.SetCorrectPt(false);
          GEFObject FitTaus = GEF.Fit();
          if(FitTaus.Fitconverged()){
            TauHGEF = FitTaus.getTauH().LV();
            TauPiGEF = FitTaus.getTauMu().LV();
          }
        }
   
        //
	HadRefitPions_plus.push_back(Ntp->PFTauRefit_PionsP4(Tauplus,0));
	HadRefitPions_plus.push_back(Ntp->PFTauRefit_PionsP4(Tauplus,1));
	HadRefitPions_plus.push_back(Ntp->PFTauRefit_PionsP4(Tauplus,2));
	HadRefitPionsCharge_plus.push_back(Ntp->PFTauRefit_PionsCharge(Tauplus, 0));
	HadRefitPionsCharge_plus.push_back(Ntp->PFTauRefit_PionsCharge(Tauplus, 1));
	HadRefitPionsCharge_plus.push_back(Ntp->PFTauRefit_PionsCharge(Tauplus, 2));

	HadRefitPions_minus.push_back(Ntp->ChargedDaughters_P4(Tauminus));
	HadRefitPionsCharge_minus.push_back(-1);
        PionMinus_ref = Ntp->Daughters_vertex(Tauminus) - tauBSPrimaryVertex;
      }
      if(rhopion){
	HadRefitPions_plus.push_back(Ntp->ChargedDaughters_P4(Tauplus));
	HadRefitPions_plus.push_back(Tauplusvis - Ntp->ChargedDaughters_P4(Tauplus));
	HadRefitPionsCharge_plus.push_back(0);
	HadRefitPionsCharge_plus.push_back(1);

	HadRefitPions_minus.push_back(Ntp->ChargedDaughters_P4(Tauminus));
	HadRefitPionsCharge_minus.push_back(-1);
 	PionMinus_ref = Ntp->Daughters_vertex(Tauminus) - tauBSPrimaryVertex;
      }
      Tau1MTT = TauPlusMTT; Tau2MTT = TauMinusMTT;
      Tau1SVFit = TauPlusSVFit; Tau2SVFit = TauMinusSVFit;
      Tau1Mixed = TauPlusMixed; Tau2Mixed = TauMinusMixed;
    }
    if(minusplus){
      if(a1rho){
	HadRefitPions_minus.push_back(Ntp->PFTauRefit_PionsP4(Tauminus,0));
	HadRefitPions_minus.push_back(Ntp->PFTauRefit_PionsP4(Tauminus,1));
	HadRefitPions_minus.push_back(Ntp->PFTauRefit_PionsP4(Tauminus,2));
	HadRefitPionsCharge_minus.push_back(Ntp->PFTauRefit_PionsCharge(Tauminus, 0));
	HadRefitPionsCharge_minus.push_back(Ntp->PFTauRefit_PionsCharge(Tauminus, 1));
	HadRefitPionsCharge_minus.push_back(Ntp->PFTauRefit_PionsCharge(Tauminus, 2));

	HadRefitPions_plus.push_back(Ntp->ChargedDaughters_P4(Tauplus));
	HadRefitPions_plus.push_back(Tauplusvis - Ntp->ChargedDaughters_P4(Tauplus));
	HadRefitPionsCharge_plus.push_back(0);
	HadRefitPionsCharge_plus.push_back(1);
      }
      if(a1pion){
        //GEF
        if(Ntp->PFTau_TrackHasMomentum(Tauplus) == true && Ntp->isPVCovAvailable() == true){
          TMatrixTSym<double> PV_cov = Ntp->PFTau_TIP_primaryVertex_cov();
          LorentzVectorParticle A1 = Ntp->PFTau_a1_lvp(Tauminus);
          TrackParticle PION = Ntp->PFTau_Track(Tauplus);
          GlobalEventFit GEF(PION,A1,METInput,tauBSPrimaryVertex,PV_cov);
          GEF.setMassConstraint(125.3);
          GEF.SetCorrectPt(false);
          GEFObject FitTaus = GEF.Fit();
          if(FitTaus.Fitconverged()){
            TauHGEF = FitTaus.getTauH().LV();
            TauPiGEF = FitTaus.getTauMu().LV();
          }
        }
	
	//
	HadRefitPions_minus.push_back(Ntp->PFTauRefit_PionsP4(Tauminus,0));
	HadRefitPions_minus.push_back(Ntp->PFTauRefit_PionsP4(Tauminus,1));
	HadRefitPions_minus.push_back(Ntp->PFTauRefit_PionsP4(Tauminus,2));
	HadRefitPionsCharge_minus.push_back(Ntp->PFTauRefit_PionsCharge(Tauminus, 0));
	HadRefitPionsCharge_minus.push_back(Ntp->PFTauRefit_PionsCharge(Tauminus, 1));
	HadRefitPionsCharge_minus.push_back(Ntp->PFTauRefit_PionsCharge(Tauminus, 2));

	HadRefitPions_plus.push_back(Ntp->ChargedDaughters_P4(Tauplus));
	HadRefitPionsCharge_plus.push_back(1);
 	PionPlus_ref = Ntp->Daughters_vertex(Tauplus) - tauBSPrimaryVertex;
      }
      if(rhopion){
	HadRefitPions_minus.push_back(Ntp->ChargedDaughters_P4(Tauminus));
	HadRefitPions_minus.push_back(Tauminusvis - Ntp->ChargedDaughters_P4(Tauminus));
	HadRefitPionsCharge_minus.push_back(0);
	HadRefitPionsCharge_minus.push_back(-1);

	HadRefitPions_plus.push_back(Ntp->ChargedDaughters_P4(Tauplus));
	HadRefitPionsCharge_plus.push_back(1);
        PionPlus_ref = Ntp->Daughters_vertex(Tauplus) - tauBSPrimaryVertex;
      }
      Tau1MTT = TauMinusMTT; Tau2MTT = TauPlusMTT;
      Tau1SVFit = TauMinusSVFit; Tau2SVFit = TauPlusSVFit;
      Tau1Mixed = TauMinusMixed; Tau2Mixed = TauPlusMixed;
    }
  }

  if((HadPionsTruth_minus!=HadPionsTruth_plus) && (HadPionsTruth_minus!=VectZeroLV) && (HadPionsTruth_plus!=VectZeroLV)) genPionsAreOk = true;
  if((TauMinusTruth != TauPlusTruth) && (TauMinusTruth != zeroLV) && (TauPlusTruth != zeroLV)) genTausAreOk = true;
   
  if (std::isnan(Wspin)!=true && GenMatchSelection)
    {
      if(channeltruth)
	{
	  //GEN LEVEL PLOTS
	  if(genPionsAreOk)
	    {
	      if(a1a1truth){
		if(genTausAreOk){
		  Acop_PVTruth = ScalcPVTruth.AcopAngle("a1", "a1", TauMinusTruth, HadPionsTruth_minus, HadPionsChargeTruth_minus, TauPlusTruth, HadPionsTruth_plus, HadPionsChargeTruth_plus);
		}
		Acop_DPTruth = ScalcDPTruth.AcopAngle_DP("a1", "a1", HadPionsTruth_minus, HadPionsTruth_plus);
	      }
	      if(rhorhotruth){
		if(genTausAreOk){
		  Acop_PVTruth = ScalcPVTruth.AcopAngle("rho", "rho", TauMinusTruth, HadPionsTruth_minus, HadPionsChargeTruth_minus, TauPlusTruth, HadPionsTruth_plus, HadPionsChargeTruth_plus);
		}
		Acop_DPTruth = ScalcDPTruth.AcopAngle_DP("rho", "rho", HadPionsTruth_minus, HadPionsTruth_plus);
	      }
	      if(pionpiontruth){
		if(genTausAreOk){
		  Acop_PVTruth = ScalcPVTruth.AcopAngle("pion", "pion", TauMinusTruth, HadPionsTruth_minus, HadPionsChargeTruth_minus, TauPlusTruth, HadPionsTruth_plus, HadPionsChargeTruth_plus);
		}
		Acop_IPTruth = ScalcIPTruth.AcopAngle_IP(HadPionsTruth_minus.at(0), PionMinus_refTruth, HadPionsTruth_plus.at(0), PionPlus_refTruth);
	      }
	      if(a1rhotruth){
		if(plusminustruth){
		  if(genTausAreOk){
		    Acop_PVTruth = ScalcPVTruth.AcopAngle("rho", "a1", TauMinusTruth, HadPionsTruth_minus, HadPionsChargeTruth_minus, TauPlusTruth, HadPionsTruth_plus, HadPionsChargeTruth_plus);
		  }
		  Acop_DPTruth = ScalcDPTruth.AcopAngle_DP("rho", "a1", HadPionsTruth_minus, HadPionsTruth_plus);
		}
		if(minusplustruth){
		  if(genTausAreOk){
		    Acop_PVTruth = ScalcPVTruth.AcopAngle("a1", "rho", TauMinusTruth, HadPionsTruth_minus, HadPionsChargeTruth_minus, TauPlusTruth, HadPionsTruth_plus, HadPionsChargeTruth_plus);
		  }
		  Acop_DPTruth = ScalcDPTruth.AcopAngle_DP("a1", "rho", HadPionsTruth_minus, HadPionsTruth_plus);
		}
	      }
	      if(a1piontruth){
		if(plusminustruth){
		  if(genTausAreOk){
		    Acop_PVTruth = ScalcPVTruth.AcopAngle("pion", "a1", TauMinusTruth, HadPionsTruth_minus, HadPionsChargeTruth_minus, TauPlusTruth, HadPionsTruth_plus, HadPionsChargeTruth_plus);
		    Acop_PVIPTruth = ScalcPVIPTruth.AcopAngle_PVIP("a1","pion", TauPlusTruth, 1., HadPionsTruth_plus, HadPionsChargeTruth_plus, TauMinusTruth, HadPionsTruth_minus.at(0), PionMinus_refTruth);
		  }
		  Acop_DPIPTruth = ScalcDPIPTruth.AcopAngle_DPIP("a1","pion", HadPionsTruth_plus, HadPionsTruth_minus.at(0), PionMinus_refTruth);
		}
		if(minusplustruth){
		  if(genTausAreOk){
		    Acop_PVTruth = ScalcPVTruth.AcopAngle("a1", "pion", TauMinusTruth, HadPionsTruth_minus, HadPionsChargeTruth_minus, TauPlusTruth, HadPionsTruth_plus, HadPionsChargeTruth_plus);
		    Acop_PVIPTruth = ScalcPVIPTruth.AcopAngle_PVIP("a1","pion", TauMinusTruth, -1., HadPionsTruth_minus, HadPionsChargeTruth_minus, TauPlusTruth, HadPionsTruth_plus.at(0), PionPlus_refTruth);
		  }
		  Acop_DPIPTruth = ScalcDPIPTruth.AcopAngle_DPIP("a1","pion", HadPionsTruth_minus, HadPionsTruth_plus.at(0), PionPlus_refTruth);
		}
	      }
	      if(rhopiontruth){
		if(plusminustruth){
		  if(genTausAreOk){
		    Acop_PVTruth = ScalcPVTruth.AcopAngle("pion", "rho", TauMinusTruth, HadPionsTruth_minus, HadPionsChargeTruth_minus, TauPlusTruth, HadPionsTruth_plus, HadPionsChargeTruth_plus);
		    Acop_PVIPTruth = ScalcPVIPTruth.AcopAngle_PVIP("rho","pion", TauPlusTruth, 1., HadPionsTruth_plus, HadPionsChargeTruth_plus, TauMinusTruth, HadPionsTruth_minus.at(0), PionMinus_refTruth);
		  }
		  Acop_DPIPTruth = ScalcDPIPTruth.AcopAngle_DPIP("rho","pion", HadPionsTruth_plus, HadPionsTruth_minus.at(0), PionMinus_refTruth);
		}
		if(minusplustruth){
		  if(genTausAreOk){
		    Acop_PVTruth = ScalcPVTruth.AcopAngle("rho", "pion", TauMinusTruth, HadPionsTruth_minus, HadPionsChargeTruth_minus, TauPlusTruth, HadPionsTruth_plus, HadPionsChargeTruth_plus);
		    Acop_PVIPTruth = ScalcPVIPTruth.AcopAngle_PVIP("rho","pion", TauMinusTruth, -1.,  HadPionsTruth_minus, HadPionsChargeTruth_minus, TauPlusTruth, HadPionsTruth_plus.at(0), PionPlus_refTruth);
		  }
		  Acop_DPIPTruth = ScalcDPIPTruth.AcopAngle_DPIP("rho","pion", HadPionsTruth_minus, HadPionsTruth_plus.at(0), PionPlus_refTruth);
		}
	      }
	      decayplaneAcopAngleTruth.at(t).Fill(Acop_DPTruth,Wspin);
	      polarimetricAcopAngleTruth.at(t).Fill(Acop_PVTruth,Wspin);
	      impactparameterAcopAngleTruth.at(t).Fill(Acop_IPTruth,Wspin);
	      DPIPAcopAngleTruth.at(t).Fill(Acop_DPIPTruth,Wspin);
	      PVIPAcopAngleTruth.at(t).Fill(Acop_PVIPTruth,Wspin);
	    } //gen pions
	} //gen channel

  if((HadRefitPions_minus!=HadRefitPions_plus) && (HadRefitPions_minus!=VectZeroLV) && (HadRefitPions_plus!=VectZeroLV)) recoPionsAreOk = true;
  if((TauPlusMixed!=TauMinusMixed) && (TauPlusMixed!=zeroLV) && (TauMinusMixed!=zeroLV)) mixedTausAreOk = true;
  if((TauminusPairConstraintBS!=TauplusPairConstraintBS) && (TauminusPairConstraintBS!=zeroLV) && (TauplusPairConstraintBS!=zeroLV)) recoTausAreOk = true;
  if((TauPiGEF!=TauHGEF) && (TauPiGEF!=zeroLV) && (TauHGEF!=zeroLV)) GEFTausAreOk = true;

      //RECO LEVEL PLOTS
      if(channel)
	{
	  if(recoPionsAreOk)
	    {
	      if(a1a1){
		    if(recoTausAreOk){
		      Acop_PVGEF = ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBS, HadRefitPions_plus, HadRefitPionsCharge_plus);
		    }
		if(mixedTausAreOk){
		  Acop_PV = ScalcPV.AcopAngle("a1", "a1", TauMinusMixed, HadRefitPions_minus, HadRefitPionsCharge_minus, TauPlusMixed, HadRefitPions_plus, HadRefitPionsCharge_plus);
		}
		  Acop_DP = ScalcDP.AcopAngle_DP("a1", "a1", HadRefitPions_minus, HadRefitPions_plus);
	}
		if(rhorho){
		  if(mixedTausAreOk){
		  Acop_PV = ScalcPV.AcopAngle("rho", "rho", TauMinusMixed, HadRefitPions_minus, HadRefitPionsCharge_minus, TauPlusMixed, HadRefitPions_plus, HadRefitPionsCharge_plus);
		}
		  Acop_DP = ScalcDP.AcopAngle_DP("rho", "rho", HadRefitPions_minus, HadRefitPions_plus);
		}
		if(pionpion){
		if(mixedTausAreOk){
		  Acop_PV = ScalcPV.AcopAngle("pion", "pion", TauMinusMixed, HadRefitPions_minus, HadRefitPionsCharge_minus, TauPlusMixed, HadRefitPions_plus, HadRefitPionsCharge_plus);
}
		  Acop_IP = ScalcIP.AcopAngle_IP(HadRefitPions_minus.at(0), PionMinus_ref, HadRefitPions_plus.at(0), PionPlus_ref);
		}
		if(a1rho){
		  if(plusminus){
		if(mixedTausAreOk){
		    Acop_PV = ScalcPV.AcopAngle("rho", "a1", TauMinusMixed, HadRefitPions_minus, HadRefitPionsCharge_minus, TauPlusMixed, HadRefitPions_plus, HadRefitPionsCharge_plus);
}
		    Acop_DP = ScalcDP.AcopAngle_DP("rho", "a1", HadRefitPions_minus, HadRefitPions_plus);
		  }
		  if(minusplus){
		if(mixedTausAreOk){
		    Acop_PV = ScalcPV.AcopAngle("a1", "rho", TauMinusMixed, HadRefitPions_minus, HadRefitPionsCharge_minus, TauPlusMixed, HadRefitPions_plus, HadRefitPionsCharge_plus);
}
		    Acop_DP = ScalcDP.AcopAngle_DP("a1", "rho", HadRefitPions_minus, HadRefitPions_plus);
		  }
		}
		if(a1pion){
		  if(plusminus){
		if(mixedTausAreOk){
		    Acop_PV = ScalcPV.AcopAngle("pion", "a1", TauMinusMixed, HadRefitPions_minus, HadRefitPionsCharge_minus, TauPlusMixed, HadRefitPions_plus, HadRefitPionsCharge_plus);
                    Acop_PVIP = ScalcPVIP.AcopAngle_PVIP("a1","pion", TauPlusMixed, 1., HadRefitPions_plus, HadRefitPionsCharge_plus, TauMinusMixed, HadRefitPions_minus.at(0), PionMinus_ref);
}
		if(GEFTausAreOk){
                    Acop_PVGEF = ScalcPV.AcopAngle("pion", "a1", TauPiGEF, HadRefitPions_minus, HadRefitPionsCharge_minus, TauHGEF, HadRefitPions_plus, HadRefitPionsCharge_plus);
}
		    Acop_DPIP = ScalcDPIP.AcopAngle_DPIP("a1","pion", HadRefitPions_plus, HadRefitPions_minus.at(0), PionMinus_ref);
		  }
		  if(minusplus){
		if(mixedTausAreOk){
		    Acop_PV = ScalcPV.AcopAngle("a1", "pion", TauMinusMixed, HadRefitPions_minus, HadRefitPionsCharge_minus, TauPlusMixed, HadRefitPions_plus, HadRefitPionsCharge_plus);
                    Acop_PVIP = ScalcPVIP.AcopAngle_PVIP("a1","pion", TauMinusMixed, -1.,  HadRefitPions_minus, HadRefitPionsCharge_minus, TauPlusMixed, HadRefitPions_plus.at(0), PionPlus_ref);
}
		if(GEFTausAreOk){
                    Acop_PVGEF = ScalcPV.AcopAngle("a1", "pion", TauHGEF, HadRefitPions_minus, HadRefitPionsCharge_minus, TauPiGEF, HadRefitPions_plus, HadRefitPionsCharge_plus);
}
		    Acop_DPIP = ScalcDPIP.AcopAngle_DPIP("a1","pion", HadRefitPions_minus, HadRefitPions_plus.at(0), PionPlus_ref);
		  }
		}
		if(rhopion){
		  if(plusminus){
		if(mixedTausAreOk){
		    Acop_PV = ScalcPV.AcopAngle("pion", "rho", TauMinusMixed, HadRefitPions_minus, HadRefitPionsCharge_minus, TauPlusMixed, HadRefitPions_plus, HadRefitPionsCharge_plus);
                    Acop_PVIP = ScalcPVIP.AcopAngle_PVIP("rho","pion", TauPlusMixed, 1.,  HadRefitPions_plus, HadRefitPionsCharge_plus, TauMinusMixed, HadRefitPions_minus.at(0), PionMinus_ref);
}
		    Acop_DPIP = ScalcDPIP.AcopAngle_DPIP("rho","pion", HadRefitPions_plus, HadRefitPions_minus.at(0), PionMinus_ref);
		  }
		  if(minusplus){
		if(mixedTausAreOk){
		    Acop_PV = ScalcPV.AcopAngle("rho", "pion", TauMinusMixed, HadRefitPions_minus, HadRefitPionsCharge_minus, TauPlusMixed, HadRefitPions_plus, HadRefitPionsCharge_plus);
                    Acop_PVIP = ScalcPVIP.AcopAngle_PVIP("rho","pion", TauMinusMixed, -1., HadRefitPions_minus, HadRefitPionsCharge_minus, TauPlusMixed, HadRefitPions_plus.at(0), PionPlus_ref);

}
		    Acop_DPIP = ScalcDPIP.AcopAngle_DPIP("rho","pion", HadRefitPions_minus, HadRefitPions_plus.at(0), PionPlus_ref);
		  }
		}
		decayplaneAcopAngle.at(t).Fill(Acop_DP,Wspin);
		polarimetricAcopAngle.at(t).Fill(Acop_PV,Wspin);
		polarimetricGEFAcopAngle.at(t).Fill(Acop_PVGEF,Wspin);
 		impactparameterAcopAngle.at(t).Fill(Acop_IP,Wspin);
		DPIPAcopAngle.at(t).Fill(Acop_DPIP,Wspin);
		PVIPAcopAngle.at(t).Fill(Acop_PVIP,Wspin);
	    } //pions
	} //channel

      //PULLS
      Ref.at(t).Fill(Wspin);
      if(mixedTausAreOk){
	SVfitMTTdR1.at(t).Fill(Tau1MTT.DeltaR(Tau1SVFit),w);
	SVfitMTTdR2.at(t).Fill(Tau2MTT.DeltaR(Tau2SVFit),w);

        double P1svfit = (Tau1SVFit.P()-Tau1Truth.P())/Tau1Truth.P();
        double P2svfit = (Tau2SVFit.P()-Tau2Truth.P())/Tau2Truth.P();
        double Eta1MTT = (Tau1MTT.Eta()-Tau1Truth.Eta())/Tau1Truth.Eta();
        double Phi1MTT = (Tau1MTT.Phi()-Tau1Truth.Phi())/Tau1Truth.Phi();
        double Eta2MTT = (Tau2MTT.Eta()-Tau2Truth.Eta())/Tau2Truth.Eta();
        double Phi2MTT = (Tau2MTT.Phi()-Tau2Truth.Phi())/Tau2Truth.Phi();
        double DAcopPV = (Acop_PV-Acop_PVTruth)/Acop_PVTruth;
	
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

      }

      //FastMTT
      if(Tau1MTT!=zeroLV && Tau2MTT!=zeroLV){
	DPhitau2MTT = (Tau2MTT.Phi()-Tau2Truth.Phi())/Tau2Truth.Phi();
	DEtatau2MTT = (Tau2MTT.Eta()-Tau2Truth.Eta())/Tau2Truth.Eta();
	DPtau2MTT = (Tau2MTT.P()-Tau2Truth.P())/Tau2Truth.P();
	DEtau2MTT = (Tau2MTT.E()-Tau2Truth.E())/Tau2Truth.E();
	DPhitau1MTT = (Tau1MTT.Phi()-Tau1Truth.Phi())/Tau1Truth.Phi();
	DEtatau1MTT = (Tau1MTT.Eta()-Tau1Truth.Eta())/Tau1Truth.Eta();
	DPtau1MTT = (Tau1MTT.P()-Tau1Truth.P())/Tau1Truth.P();
	DEtau1MTT = (Tau1MTT.E()-Tau1Truth.E())/Tau1Truth.E();
	DeltaPhitau2MTT.at(t).Fill(DPhitau2MTT,w);
	DeltaEtatau2MTT.at(t).Fill(DEtatau2MTT,w);
	DeltaPtau2MTT.at(t).Fill(DPtau2MTT,w);
	DeltaEtau2MTT.at(t).Fill(DEtau2MTT,w);
	DeltaPhitau1MTT.at(t).Fill(DPhitau1MTT,w);
	DeltaEtatau1MTT.at(t).Fill(DEtatau1MTT,w);
	DeltaPtau1MTT.at(t).Fill(DPtau1MTT,w);
	DeltaEtau1MTT.at(t).Fill(DEtau1MTT,w);
      }
      //SVFit
      if(Tau1SVFit!=zeroLV && Tau2SVFit!=zeroLV){
	DPhitau2SVFit = (Tau2SVFit.Phi()-Tau2Truth.Phi())/Tau2Truth.Phi();
	DEtatau2SVFit = (Tau2SVFit.Eta()-Tau2Truth.Eta())/Tau2Truth.Eta();
	DPtau2SVFit = (Tau2SVFit.P()-Tau2Truth.P())/Tau2Truth.P();
	DEtau2SVFit = (Tau2SVFit.E()-Tau2Truth.E())/Tau2Truth.E();
	DPhitau1SVFit = (Tau1SVFit.Phi()-Tau1Truth.Phi())/Tau1Truth.Phi();
	DEtatau1SVFit = (Tau1SVFit.Eta()-Tau1Truth.Eta())/Tau1Truth.Eta();
	DPtau1SVFit = (Tau1SVFit.P()-Tau1Truth.P())/Tau1Truth.P();
	DEtau1SVFit = (Tau1SVFit.E()-Tau1Truth.E())/Tau1Truth.E();
	DeltaPhitau2SVFit.at(t).Fill(DPhitau2SVFit,w);
	DeltaEtatau2SVFit.at(t).Fill(DEtatau2SVFit,w);
	DeltaPtau2SVFit.at(t).Fill(DPtau2SVFit,w);
	DeltaEtau2SVFit.at(t).Fill(DEtau2SVFit,w);
	DeltaPhitau1SVFit.at(t).Fill(DPhitau1SVFit,w);
	DeltaEtatau1SVFit.at(t).Fill(DEtatau1SVFit,w);
	DeltaPtau1SVFit.at(t).Fill(DPtau1SVFit,w);
	DeltaEtau1SVFit.at(t).Fill(DEtau1SVFit,w);
      }
      //Mixed
      if(Tau1Mixed!=zeroLV && Tau2Mixed!=zeroLV){
	DPhitau2Mixed = (Tau2Mixed.Phi()-Tau2Truth.Phi())/Tau2Truth.Phi();
	DEtatau2Mixed = (Tau2Mixed.Eta()-Tau2Truth.Eta())/Tau2Truth.Eta();
	DPtau2Mixed = (Tau2Mixed.P()-Tau2Truth.P())/Tau2Truth.P();
	DEtau2Mixed = (Tau2Mixed.E()-Tau2Truth.E())/Tau2Truth.E();
	DPhitau1Mixed = (Tau1Mixed.Phi()-Tau1Truth.Phi())/Tau1Truth.Phi();
	DEtatau1Mixed = (Tau1Mixed.Eta()-Tau1Truth.Eta())/Tau1Truth.Eta();
	DPtau1Mixed = (Tau1Mixed.P()-Tau1Truth.P())/Tau1Truth.P();
	DEtau1Mixed = (Tau1Mixed.E()-Tau1Truth.E())/Tau1Truth.E();
	DeltaPhitau2Mixed.at(t).Fill(DPhitau2Mixed,w);
	DeltaEtatau2Mixed.at(t).Fill(DEtatau2Mixed,w);
	DeltaPtau2Mixed.at(t).Fill(DPtau2Mixed,w);
	DeltaEtau2Mixed.at(t).Fill(DEtau2Mixed,w);
	DeltaPhitau1Mixed.at(t).Fill(DPhitau1Mixed,w);
	DeltaEtatau1Mixed.at(t).Fill(DEtatau1Mixed,w);
	DeltaPtau1Mixed.at(t).Fill(DPtau1Mixed,w);
	DeltaEtau1Mixed.at(t).Fill(DEtau1Mixed,w);
      }
      //GEF
      if(TauplusPairConstraintBS!=zeroLV && TauminusPairConstraintBS!=zeroLV){
	DPhitau2GEF = (TauplusPairConstraintBS.Phi()-TauPlusTruth.Phi())/TauPlusTruth.Phi();
	DEtatau2GEF = (TauplusPairConstraintBS.Eta()-TauPlusTruth.Eta())/TauPlusTruth.Eta();
	DPtau2GEF = (TauplusPairConstraintBS.P()-TauPlusTruth.P())/TauPlusTruth.P();
	DEtau2GEF = (TauplusPairConstraintBS.E()-TauPlusTruth.E())/TauPlusTruth.E();
	DPhitau1GEF = (TauminusPairConstraintBS.Phi()-TauMinusTruth.Phi())/TauMinusTruth.Phi();
	DEtatau1GEF = (TauminusPairConstraintBS.Eta()-TauMinusTruth.Eta())/TauMinusTruth.Eta();
	DPtau1GEF = (TauminusPairConstraintBS.P()-TauMinusTruth.P())/TauMinusTruth.P();
	DEtau1GEF = (TauminusPairConstraintBS.E()-TauMinusTruth.E())/TauMinusTruth.E();
	DeltaPhitau2GEF.at(t).Fill(DPhitau2GEF,w);
	DeltaEtatau2GEF.at(t).Fill(DEtatau2GEF,w);
	DeltaPtau2GEF.at(t).Fill(DPtau2GEF,w);
	DeltaEtau2GEF.at(t).Fill(DEtau2GEF,w);
	DeltaPhitau1GEF.at(t).Fill(DPhitau1GEF,w);
	DeltaEtatau1GEF.at(t).Fill(DEtatau1GEF,w);
	DeltaPtau1GEF.at(t).Fill(DPtau1GEF,w);
	DeltaEtau1GEF.at(t).Fill(DEtau1GEF,w);
      }
      if(TauPiGEF !=zeroLV && TauHGEF!=zeroLV){
        DPhitauPiGEF = (TauPiGEF.Phi()-TauPiTruth.Phi())/TauPiTruth.Phi();
        DEtatauPiGEF = (TauPiGEF.Eta()-TauPiTruth.Eta())/TauPiTruth.Eta();
        DPtauPiGEF = (TauPiGEF.P()-TauPiTruth.P())/TauPiTruth.P();
        DEtauPiGEF = (TauPiGEF.E()-TauPiTruth.E())/TauPiTruth.E();
        DPhitauHGEF = (TauHGEF.Phi()-TauHTruth.Phi())/TauHTruth.Phi();
        DEtatauHGEF = (TauHGEF.Eta()-TauHTruth.Eta())/TauHTruth.Eta();
        DPtauHGEF = (TauHGEF.P()-TauHTruth.P())/TauHTruth.P();
        DEtauHGEF = (TauHGEF.E()-TauHTruth.E())/TauHTruth.E();
        DeltaPhitauPiGEF.at(t).Fill(DPhitauPiGEF,w);
        DeltaEtatauPiGEF.at(t).Fill(DEtatauPiGEF,w);
        DeltaPtauPiGEF.at(t).Fill(DPtauPiGEF,w);
        DeltaEtauPiGEF.at(t).Fill(DEtauPiGEF,w);
        DeltaPhitauHGEF.at(t).Fill(DPhitauHGEF,w);
        DeltaEtatauHGEF.at(t).Fill(DEtatauHGEF,w);
        DeltaPtauHGEF.at(t).Fill(DPtauHGEF,w);
        DeltaEtauHGEF.at(t).Fill(DEtauHGEF,w);
      }
    } //isnan
} //do event

//  This is a function if you want to do something after the event loop
void HCPTauTau::Finish() {

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
  Selection::Finish();
}
