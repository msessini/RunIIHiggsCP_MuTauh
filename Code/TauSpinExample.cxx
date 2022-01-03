#include "TauSpinExample.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "TF1.h"
#include "TauSolver.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TauSpinerInterface.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

TauSpinExample::TauSpinExample(TString Name_, TString id_):
  Selection(Name_,id_),
  zsbins(20),
  zsmin(-0.5),
  zsmax(0.5)
{
  //  verbose=true;
}

TauSpinExample::~TauSpinExample(){
  for(unsigned int j=0; j<Npassed.size(); j++){
    Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "done." << std::endl;
}

void  TauSpinExample::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==isZtautauto3pimu)          cut.at(isZtautauto3pimu)=1;
  }

  TString hlabel;
  TString htitle; for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
  
    if(i==isZtautauto3pimu){
      title.at(i)="Is $Z\\rightarrow\\tau\\tau\\rightarrow\\mu\\pi\\pi\\pi$ MC (bool)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Is Z#rightarrow#tau#tau#rightarrow#mu#pi#pi#pi MC (bool)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_isZtautauto3pimu_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_isZtautauto3pimu_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  //mu
  mu_ExoverEtau=HConfig.GetTH1D(Name+"_mu_ExoverEtau","ExoverEtau",20,0.0,1.0,"E_{#mu}/E_{#tau}","Events");
  mu_ExoverEtau_hplus=HConfig.GetTH1D(Name+"_mu_ExoverEtau_hplus","ExoverEtau_hplus",20,0.0,1.0,"E_{#mu}/E_{#tau}|_{h^{+}}","Events");
  mu_ExoverEtau_hminus=HConfig.GetTH1D(Name+"_mu_ExoverEtau_hminus","ExoverEtau_hminus",20,0.0,1.0,"E_{#mu}/E_{#tau}|_{h^{-}}","Events");
  mu_ExoverEtau_Spin=HConfig.GetTH1D(Name+"_mu_ExoverEtau_Spin","ExoverEtau_Spin",20,0.0,1.0,"E_{#mu}/E_{#tau}|_{Spin}","Events");
  mu_ExoverEtau_UnSpin=HConfig.GetTH1D(Name+"_mu_ExoverEtau_UnSpin","ExoverEtau_UnSpin",20,0.0,1.0,"E_{#mu}/E_{#tau}|_{UnSpin}","Events");
  mu_ExoverEtau_FlipSpin=HConfig.GetTH1D(Name+"_mu_ExoverEtau_FlipSpin","ExoverEtau_FlipSpin",20,0.0,1.0,"E_{#mu}/E_{#tau}|_{FlipSpin}","Events");
  mu_WT_Spin=HConfig.GetTH1D(Name+"_mu_WT_Spin","WT_Spin",40,0.0,4.0,"WT|_{#mu}","Events");
  mu_WT_UnSpin=HConfig.GetTH1D(Name+"_mu_WT_UnSpin","WT_UnSpin",100,0.0,10.0,"1/WT|_{#mu}","Events");
  mu_WT_FlipSpin=HConfig.GetTH1D(Name+"_mu_WT_FlipSpin","WT_FlipSpin",100,0.0,10,"(2-WT)/(WT)|_{#mu}","Events");
  mu_PtRatio_hplus=HConfig.GetTH1D(Name+"_mu_PtRatio_hplus","PtRatio_hplus",20,0.0,1.0,"P_{t,#mu}/P_{t,#tau}|_{h^{+}}","Events");
  mu_PtRatio_hminus=HConfig.GetTH1D(Name+"_mu_PtRatio_hminus","PtRatio_hminus",20,0.0,1.0,"P_{t,#mu}/P_{t,#tau}|_{h^{-}}","Events");

  //pinu
  pi_ExoverEtau=HConfig.GetTH1D(Name+"_pi_ExoverEtau","ExoverEtau",20,0.0,1.0,"E_{#pi}/E_{#tau}","Events");
  pi_ExoverEtau_hplus=HConfig.GetTH1D(Name+"_pi_ExoverEtau_hplus","ExoverEtau_hplus",20,0.0,1.0,"E_{#pi}/E_{#tau}|_{h^{+}}","Events");
  pi_ExoverEtau_hminus=HConfig.GetTH1D(Name+"_pi_ExoverEtau_hminus","ExoverEtau_hminus",20,0.0,1.0,"E_{#pi}/E_{#tau}|_{h^{-}}","Events");
  pi_ExoverEtau_Spin=HConfig.GetTH1D(Name+"_pi_ExoverEtau_Spin","ExoverEtau_Spin",20,0.0,1.0,"E_{#pi}/E_{#tau}|_{Spin}","Events");
  pi_ExoverEtau_UnSpin=HConfig.GetTH1D(Name+"_pi_ExoverEtau_UnSpin","ExoverEtau_UnSpin",20,0.0,1.0,"E_{#pi}/E_{#tau}|_{UnSpin}","Events");
  pi_ExoverEtau_FlipSpin=HConfig.GetTH1D(Name+"_pi_ExoverEtau_FlipSpin","ExoverEtau_FlipSpin",20,0.0,1.0,"E_{#pi}/E_{#tau}|_{FlipSpin}","Events");

  pi_zs=HConfig.GetTH1D(Name+"_pi_zs","zs",zsbins ,zsmin,zsmax,"zs","Events");
  pi_zs_hplus=HConfig.GetTH1D(Name+"_pi_zs_hplus","zs_hplus",zsbins ,zsmin,zsmax,"zs|_{h^{+}}","Events");
  pi_zs_hminus=HConfig.GetTH1D(Name+"_pi_zs_hminus","zs_hminus",zsbins ,zsmin,zsmax,"zs|_{h^{-}}","Events");
  pi_zs_Spin=HConfig.GetTH1D(Name+"_pi_zs_Spin","zs_Spin",zsbins ,zsmin,zsmax,"zs|_{Spin}","Events");
  pi_zs_UnSpin=HConfig.GetTH1D(Name+"_pi_zs_UnSpin","zs_UnSpin",zsbins ,zsmin,zsmax,"zs|_{UnSpin}","Events");
  pi_zs_FlipSpin=HConfig.GetTH1D(Name+"_pi_zs_FlipSpin","zs_FlipSpin",zsbins ,zsmin,zsmax,"zs|_{FlipSpin}","Events");

  pi_WT_Spin=HConfig.GetTH1D(Name+"_pi_WT_Spin","WT_Spin",40,0.0,4.0,"WT|_{#pi}","Events");
  pi_WT_UnSpin=HConfig.GetTH1D(Name+"_pi_WT_UnSpin","WT_UnSpin",100,0.0,10.0,"1/WT|_{#pi}","Events");
  pi_WT_FlipSpin=HConfig.GetTH1D(Name+"_pi_WT_FlipSpin","WT_FlipSpin",100,0.0,10,"(2-WT)/(WT)|_{#pi}","Events");
  pi_PtRatio_hplus=HConfig.GetTH1D(Name+"_pi_PtRatio_hplus","PtRatio_hplus",20,0.0,1.0,"P_{t,#pi}/P_{t,#tau}|_{h^{+}}","Events");
  pi_PtRatio_hminus=HConfig.GetTH1D(Name+"_pi_PtRatio_hminus","PtRatio_hminus",20,0.0,1.0,"P_{t,#pi}/P_t,{#tau}|_{h^{-}}","Events");

  //a1nu
  a1_ExoverEtau=HConfig.GetTH1D(Name+"_a1_ExoverEtau","ExoverEtau",20,0.0,1.0,"E_{a_{1}(1260)}/E_{#tau}","Events");
  a1_ExoverEtau_hplus=HConfig.GetTH1D(Name+"_a1_ExoverEtau_hplus","ExoverEtau_hplus",20,0.0,1.0,"E_{a_{1}(1260)}/E_{#tau}|_{h^{+}}","Events");
  a1_ExoverEtau_hminus=HConfig.GetTH1D(Name+"_a1_ExoverEtau_hminus","ExoverEtau_hminus",20,0.0,1.0,"E_{a_{1}(1260)}/E_{#tau}|_{h^{-}}","Events");
  a1_ExoverEtau_Spin=HConfig.GetTH1D(Name+"_a1_ExoverEtau_Spin","ExoverEtau_Spin",20,0.0,1.0,"E_{a_{1}(1260)}/E_{#tau}|_{Spin}","Events");
  a1_ExoverEtau_UnSpin=HConfig.GetTH1D(Name+"_a1_ExoverEtau_UnSpin","ExoverEtau_UnSpin",20,0.0,1.0,"E_{a_{1}(1260)}/E_{#tau}|_{UnSpin}","Events");
  a1_ExoverEtau_FlipSpin=HConfig.GetTH1D(Name+"_a1_ExoverEtau_FlipSpin","ExoverEtau_FlipSpin",20,0.0,1.0,"E_{a_{1}(1260)}/E_{#tau}|_{FlipSpin}","Events");

  a1_Gamma=HConfig.GetTH1D(Name+"_a1_Gamma","Gamma",20,-1.0,1.0,"#gamma_{a1(1260)}","Events");
  a1_Gamma_hplus=HConfig.GetTH1D(Name+"_a1_Gamma_hplus","Gamma_hplus",20,-1.0,1.0,"#gamma_{a1(1260)}|_{h^{+}}","Events");
  a1_Gamma_hminus=HConfig.GetTH1D(Name+"_a1_Gamma_hminus","Gamma_hminus",20,-1.0,1.0,"#gamma_{a1(1260)}|_{h^{-}}","Events");
  a1_Gamma_Spin=HConfig.GetTH1D(Name+"_a1_Gamma_Spin","Gamma_Spin",20,-1.0,1.0,"#gamma_{a1(1260)}|_{Spin}","Events");
  a1_Gamma_UnSpin=HConfig.GetTH1D(Name+"_a1_Gamma_UnSpin","Gamma_UnSpin",20,-1.0,1.0,"#gamma_{a1(1260)}|_{UnSpin}","Events");
  a1_Gamma_FlipSpin=HConfig.GetTH1D(Name+"_a1_Gamma_FlipSpin","Gamma_FlipSpin",20,-1.0,1.0,"#gamma_{a1(1260)}|_{FlipSpin}","Events");

  pi_Mvis=HConfig.GetTH1D(Name+"_pi_Mvis","Mvis",20,0.0,1.0,"M_{Vis[#pi#pi]}}","Events");
  pi_Mvis_hplus=HConfig.GetTH1D(Name+"_pi_Mvis_hplus","Mvis_hplus",20,0.0,1.0,"M_{Vis[#pi#pi]}}|_{h^{+}}","Events");
  pi_Mvis_hminus=HConfig.GetTH1D(Name+"_pi_Mvis_hminus","Mvis_hminus",20,0.0,1.0,"M_{Vis[#pi#pi]}}|_{h^{-}}","Events");
  pi_Mvis_Spin=HConfig.GetTH1D(Name+"_pi_Mvis_Spin","Mvis_Spin",20,0.0,1.0,"M_{Vis[#pi#pi]}}|_{Spin}","Events");
  pi_Mvis_UnSpin=HConfig.GetTH1D(Name+"_pi_Mvis_UnSpin","Mvis_UnSpin",20,0.0,1.0,"M_{Vis[#pi#pi]}}|_{UnSpin}","Events");
  pi_Mvis_FlipSpin=HConfig.GetTH1D(Name+"_pi_Mvis_FlipSpin","Mvis_FlipSpin",20,0.0,1.0,"M_{Vis[#pi#pi]}}|_{FlipSpin}","Events");

  a1_WT_Spin=HConfig.GetTH1D(Name+"_a1_WT_Spin","WT_Spin",40,0.0,4.0,"WT|_{a_{1}(1260)}","Events");
  a1_WT_UnSpin=HConfig.GetTH1D(Name+"_a1_WT_UnSpin","WT_UnSpin",100,0.0,10.0,"1/WT|_{a_{1}(1260)}","Events");
  a1_WT_FlipSpin=HConfig.GetTH1D(Name+"_a1_WT_FlipSpin","WT_FlipSpin",100,0.0,10,"(2-WT)/(WT)|_{a_{1}(1260)}","Events");
  a1_PtRatio_hplus=HConfig.GetTH1D(Name+"_a1_PtRatio_hplus","PtRatio_hplus",20,0.0,1.0,"P_{t,a_{1}(1260)}/P_{t,#tau}|_{h^{+}}","Events");
  a1_PtRatio_hminus=HConfig.GetTH1D(Name+"_a1_PtRatio_hminus","PtRatio_hminus",20,0.0,1.0,"P_{t,a_{1}(1260)}/P_{t,#tau}|_{h^{-}}","Events");

  a1_cosbeta_hminus=HConfig.GetTH1D(Name+"_a1_cosbeta_hminus","cosbeta_hminus",20,0.0,1.0,"cos#beta|_{h^{-}}","Events");
  a1_cosbeta_hplus=HConfig.GetTH1D(Name+"_a1_cosbeta_hplus","cosbeta_hplus",20,0.0,1.0,"cos#beta|_{h^{+}}","Events");

  //rhonu
  rho_ExoverEtau=HConfig.GetTH1D(Name+"_rho_ExoverEtau","ExoverEtau",20,0.0,1.0,"E_{#rho}/E_{#tau}","Events");
  rho_ExoverEtau_hplus=HConfig.GetTH1D(Name+"_rho_ExoverEtau_hplus","ExoverEtau_hplus",20,0.0,1.0,"E_{#rho}/E_{#tau}|_{h^{+}}","Events");
  rho_ExoverEtau_hminus=HConfig.GetTH1D(Name+"_rho_ExoverEtau_hminus","ExoverEtau_hminus",20,0.0,1.0,"E_{#rho}/E_{#tau}|_{h^{-}}","Events");
  rho_ExoverEtau_Spin=HConfig.GetTH1D(Name+"_rho_ExoverEtau_Spin","ExoverEtau_Spin",20,0.0,1.0,"E_{#rho}/E_{#tau}|_{Spin}","Events");
  rho_ExoverEtau_UnSpin=HConfig.GetTH1D(Name+"_rho_ExoverEtau_UnSpin","ExoverEtau_UnSpin",20,0.0,1.0,"E_{#rho}/E_{#tau}|_{UnSpin}","Events");
  rho_ExoverEtau_FlipSpin=HConfig.GetTH1D(Name+"_rho_ExoverEtau_FlipSpin","ExoverEtau_FlipSpin",20,0.0,1.0,"E_{#rho}/E_{#tau}|_{FlipSpin}","Events");

  rho_Gamma=HConfig.GetTH1D(Name+"_rho_Gamma","Gamma",20,-1.0,1.0,"#gamma_{#rho}","Events");
  rho_Gamma_hplus=HConfig.GetTH1D(Name+"_rho_Gamma_hplus","Gamma_hplus",20,-1.0,1.0,"#gamma_{#rho}|_{h^{+}}","Events");
  rho_Gamma_hminus=HConfig.GetTH1D(Name+"_rho_Gamma_hminus","Gamma_hminus",20,-1.0,1.0,"#gamma_{#rho}|_{h^{-}}","Events");
  rho_Gamma_Spin=HConfig.GetTH1D(Name+"_rho_Gamma_Spin","Gamma_Spin",20,-1.0,1.0,"#gamma_{#rho}|_{Spin}","Events");
  rho_Gamma_UnSpin=HConfig.GetTH1D(Name+"_rho_Gamma_UnSpin","Gamma_UnSpin",20,-1.0,1.0,"#gamma_{#rho}|_{UnSpin}","Events");
  rho_Gamma_FlipSpin=HConfig.GetTH1D(Name+"_rho_Gamma_FlipSpin","Gamma_FlipSpin",20,-1.0,1.0,"#gamma_{#rho}|_{FlipSpin}","Events");

  rho_WT_Spin=HConfig.GetTH1D(Name+"_rho_WT_Spin","WT_Spin",40,0.0,4.0,"WT|_{#rho}","Events");
  rho_WT_UnSpin=HConfig.GetTH1D(Name+"_rho_WT_UnSpin","WT_UnSpin",100,0.0,10.0,"1/WT|_{#rho}","Events");
  rho_WT_FlipSpin=HConfig.GetTH1D(Name+"_rho_WT_FlipSpin","WT_FlipSpin",100,0.0,10,"(2-WT)/(WT)|_{#rho}","Events");
  rho_PtRatio_hplus=HConfig.GetTH1D(Name+"_rho_PtRatio_hplus","PtRatio_hplus",20,0.0,1.0,"P_{t,#rho}/P_{t,#tau}|_{h^{+}}","Events");
  rho_PtRatio_hminus=HConfig.GetTH1D(Name+"_rho_PtRatio_hminus","PtRatio_hminus",20,0.0,1.0,"P_{t,#rho}/P_{t,#tau}|_{h^{-}}","Events");

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
}




void  TauSpinExample::Store_ExtraDist(){

 Extradist1d.push_back(&mu_ExoverEtau);
 Extradist1d.push_back(&mu_ExoverEtau_hplus);
 Extradist1d.push_back(&mu_ExoverEtau_hminus);
 //Extradist1d.push_back(&mu_ExoverEtau_Spin);
 //Extradist1d.push_back(&mu_ExoverEtau_UnSpin);
 //Extradist1d.push_back(&mu_ExoverEtau_FlipSpin);
 //Extradist1d.push_back(&mu_WT_Spin);
 //Extradist1d.push_back(&mu_WT_UnSpin);
 //Extradist1d.push_back(&mu_WT_FlipSpin);
 Extradist1d.push_back(&mu_PtRatio_hplus);
 Extradist1d.push_back(&mu_PtRatio_hminus);

 Extradist1d.push_back(&pi_ExoverEtau);
 Extradist1d.push_back(&pi_ExoverEtau_hplus);
 Extradist1d.push_back(&pi_ExoverEtau_hminus);
 //Extradist1d.push_back(&pi_ExoverEtau_Spin);
 //Extradist1d.push_back(&pi_ExoverEtau_UnSpin);
 //Extradist1d.push_back(&pi_ExoverEtau_FlipSpin);

 Extradist1d.push_back(&pi_zs);
 Extradist1d.push_back(&pi_zs_hplus);
 Extradist1d.push_back(&pi_zs_hminus);
 //Extradist1d.push_back(&pi_zs_Spin);
 //Extradist1d.push_back(&pi_zs_UnSpin);
 //Extradist1d.push_back(&pi_zs_FlipSpin);

 //Extradist1d.push_back(&pi_WT_Spin);
 //Extradist1d.push_back(&pi_WT_UnSpin);
 //Extradist1d.push_back(&pi_WT_FlipSpin);
 Extradist1d.push_back(&pi_PtRatio_hplus);
 Extradist1d.push_back(&pi_PtRatio_hminus);
 
 Extradist1d.push_back(&a1_ExoverEtau);
 Extradist1d.push_back(&a1_ExoverEtau_hplus);
 Extradist1d.push_back(&a1_ExoverEtau_hminus);
 //Extradist1d.push_back(&a1_ExoverEtau_Spin);
 //Extradist1d.push_back(&a1_ExoverEtau_UnSpin);
 //Extradist1d.push_back(&a1_ExoverEtau_FlipSpin);

 Extradist1d.push_back(&a1_Gamma);
 Extradist1d.push_back(&a1_Gamma_hplus);
 Extradist1d.push_back(&a1_Gamma_hminus);
 //Extradist1d.push_back(&a1_Gamma_Spin);
 //Extradist1d.push_back(&a1_Gamma_UnSpin);
 //Extradist1d.push_back(&a1_Gamma_FlipSpin);

 Extradist1d.push_back(&pi_Mvis);
 Extradist1d.push_back(&pi_Mvis_hplus);
 Extradist1d.push_back(&pi_Mvis_hminus);
 //Extradist1d.push_back(&pi_Mvis_Spin);
 //Extradist1d.push_back(&pi_Mvis_UnSpin);
 //Extradist1d.push_back(&pi_Mvis_FlipSpin);

 //Extradist1d.push_back(&a1_WT_Spin);
 //Extradist1d.push_back(&a1_WT_UnSpin);
 //Extradist1d.push_back(&a1_WT_FlipSpin);
 Extradist1d.push_back(&a1_PtRatio_hplus);
 Extradist1d.push_back(&a1_PtRatio_hminus);

 //Extradist1d.push_back(&a1_cosbeta_hplus);
 //Extradist1d.push_back(&a1_cosbeta_hminus);
 
 Extradist1d.push_back(&rho_ExoverEtau);
 Extradist1d.push_back(&rho_ExoverEtau_hplus);
 Extradist1d.push_back(&rho_ExoverEtau_hminus);
 //Extradist1d.push_back(&rho_ExoverEtau_Spin);
 //Extradist1d.push_back(&rho_ExoverEtau_UnSpin);
 //Extradist1d.push_back(&rho_ExoverEtau_FlipSpin);

 Extradist1d.push_back(&rho_Gamma);
 Extradist1d.push_back(&rho_Gamma_hplus);
 Extradist1d.push_back(&rho_Gamma_hminus);
 //Extradist1d.push_back(&rho_Gamma_Spin);
 //Extradist1d.push_back(&rho_Gamma_UnSpin);
 //Extradist1d.push_back(&rho_Gamma_FlipSpin);

 Extradist1d.push_back(&rho_WT_Spin);
 //Extradist1d.push_back(&rho_WT_UnSpin);
 //Extradist1d.push_back(&rho_WT_FlipSpin);
 Extradist1d.push_back(&rho_PtRatio_hplus);
 Extradist1d.push_back(&rho_PtRatio_hminus);

}

void  TauSpinExample::doEvent(){
  unsigned int t(0);
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){t=2;}// std::cout << "failed to find id" <<std::endl; return;}
  unsigned int Boson_idx,tau_idx,tau1_idx(0),tau2_idx(0);
  value.at(isZtautauto3pimu)=0;
  pass.at(isZtautauto3pimu) = true;
  if(pass.at(isZtautauto3pimu))value.at(isZtautauto3pimu)=1;

  double wobs=1;
  double w=1;
  Logger(Logger::Verbose) << Ntp->GetMCID() << " " << Npassed.size() << " " << t << " " << Boson_idx << " " << " " << tau_idx
	    << " " << Ntp->NMCTaus() << std::endl;
  bool status=AnalysisCuts(t,w,wobs); 
  ///////////////////////////////////////////////////////////
  // Add plots
  if(status){
	  Logger(Logger::Verbose) <<"MC type: " << Ntp->GetMCID() <<std::endl;
    // reject unsupported modes
    if(Ntp->hasSignalTauDecay(PDGInfo::Z0,Boson_idx,tau1_idx,tau2_idx)){
      bool ImpTau=true;
      double jakid1=Ntp->MCTau_JAK(tau1_idx);
      double jakid2=Ntp->MCTau_JAK(tau2_idx);
      if(!(jakid1==TauDecay::JAK_PION || jakid1==TauDecay::JAK_MUON || jakid1==TauDecay::JAK_RHO_PIPI0 || jakid1==TauDecay::JAK_A1_3PI)) ImpTau=false;
      if(!(jakid2==TauDecay::JAK_PION || jakid2==TauDecay::JAK_MUON || jakid2==TauDecay::JAK_RHO_PIPI0 || jakid2==TauDecay::JAK_A1_3PI)) ImpTau=false;
      //if(jakid1!=TauDecay::JAK_MUON && jakid2!=TauDecay::JAK_MUON) return;
      if(jakid1!=jakid2) return;
      if(!ImpTau){ Logger(Logger::Warning) << "Decay modes not implemented in TauSpinner " << std::endl; return;}
    }
    else{
      Logger(Logger::Warning) << "Not a valid Z0->tau+tau- decay" <<std::endl;
      return;
    }
    double Spin_WT=Ntp->TauSpinerGet(TauSpinerInterface::Spin);
    double UnSpin_WT=1;//Ntp->TauSpinerGet(TauSpinerInterface::UnSpin);
    double FlipSpin_WT=1;//Ntp->TauSpinerGet(TauSpinerInterface::FlipSpin);
    double hplus=Ntp->TauSpinerGet(TauSpinerInterface::hplus);
    double hminus=1-hplus;//Ntp->TauSpinerGet(TauSpinerInterface::hminus);//1-hplus;
    Logger(Logger::Info) << "hplus " << hplus << " hminus " << hminus << std::endl;
    //if(Spin_WT<=0 || 2<=Spin_WT){Spin_WT=0.000001;UnSpin_WT=0.000001;FlipSpin_WT=0.000001;}
      ////////////////////////////////////////////////
      //
      // Spin Validation
      //
    if(Ntp->hasSignalTauDecay(PDGInfo::Z0,Boson_idx,TauDecay::JAK_MUON,tau_idx)){
    	Logger(Logger::Verbose) << "muon" << std::endl;
      TLorentzVector Boson_LV=Ntp->MCSignalParticle_p4(Boson_idx);
      TLorentzVector Tau_LV(0,0,0,0);
      TLorentzVector X_LV(0,0,0,0);
      for(int i=0; i<Ntp->NMCTauDecayProducts(tau_idx);i++){
        if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::tau_minus)){
          Tau_LV=Ntp->MCTauandProd_p4(tau_idx,i);
        }
        else if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::mu_plus)){
          X_LV+=Ntp->MCTauandProd_p4(tau_idx,i);
        }
      }
      if(Tau_LV.E()>0){
        mu_PtRatio_hplus.at(t).Fill(X_LV.Pt()/Tau_LV.Pt(),w*hplus*UnSpin_WT);
        mu_PtRatio_hminus.at(t).Fill(X_LV.Pt()/Tau_LV.Pt(),w*hminus*UnSpin_WT);

        Tau_LV.Boost(-Boson_LV.BoostVector());
        X_LV.Boost(-Boson_LV.BoostVector());
        Boson_LV.Boost(-Boson_LV.BoostVector());
        // Now fill results                                                                                                                                                                                                                  
        Ntp->TauSpinerSetSignal(Ntp->MCTau_charge(tau_idx));
	mu_ExoverEtau.at(t).Fill(X_LV.E()/Tau_LV.E(),w);
	mu_ExoverEtau_Spin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Spin_WT);
	mu_ExoverEtau_UnSpin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*UnSpin_WT);
	mu_ExoverEtau_FlipSpin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*FlipSpin_WT);
	mu_ExoverEtau_hplus.at(t).Fill(X_LV.E()/Tau_LV.E(),w*hplus*UnSpin_WT);
	mu_ExoverEtau_hminus.at(t).Fill(X_LV.E()/Tau_LV.E(),w*hminus*UnSpin_WT);
	mu_WT_Spin.at(t).Fill(Spin_WT,w);
	mu_WT_UnSpin.at(t).Fill(UnSpin_WT,w);
	mu_WT_FlipSpin.at(t).Fill(FlipSpin_WT,w);

	  for(int i=0;i<Ntp->NMCTauDecayProducts(tau_idx);i++){
	    TLorentzVector LV=Ntp->MCTauandProd_p4(tau_idx,i);
	    Logger(Logger::Verbose) << Ntp->MCTauandProd_pdgid(tau_idx,i) << " " << LV.Px() << " " << LV.Py() << " " << LV.Pz() << " " << LV.E() << std::endl;
	  }
	  Logger(Logger::Verbose) << "mu : x " << X_LV.E()/Tau_LV.E() << " w " << w << "spin " <<  Spin_WT << " unspin " << UnSpin_WT << " spin flip " << FlipSpin_WT << " hplus " << hplus << " hminus " << hminus  << std::endl;
      }
      Logger(Logger::Verbose)<< "muon done " <<std::endl;
    }
    if(Ntp->hasSignalTauDecay(PDGInfo::Z0,Boson_idx,tau1_idx,tau2_idx) ){
      if(Ntp->MCTau_JAK(tau1_idx)==TauDecay::JAK_PION && Ntp->MCTau_JAK(tau2_idx)==TauDecay::JAK_PION){
    Logger(Logger::Verbose)<< "pi-pi " <<std::endl;
	double x1(0),x2(0); 
	TLorentzVector Zboson=Ntp->MCSignalParticle_p4(Boson_idx);
	TLorentzVector LVTau=Ntp->MCTau_p4(tau2_idx);
	LVTau.Boost(-1*Zboson.BoostVector());
	TLorentzVector tautau=LVTau;
	TLorentzVector pipi(0,0,0,0);
	int charge=Ntp->MCTau_charge(tau2_idx);
	for(int i=0;i<Ntp->NMCTauDecayProducts(tau2_idx);i++){
	  if(abs(Ntp->MCTauandProd_pdgid(tau2_idx,i))==abs(PDGInfo::pi_plus)){
	    TLorentzVector LVpi=Ntp->MCTauandProd_p4(tau2_idx,i);
	    pipi+=LVpi;
	    LVpi.Boost(-1*Zboson.BoostVector());
	    if(charge<0){x1=LVpi.P()/LVTau.E(); }
	    else{ x2=LVpi.P()/LVTau.E();}
	  }
	}
	//
	LVTau=Ntp->MCTau_p4(tau1_idx);
	LVTau.Boost(-1*Zboson.BoostVector());
	tautau+=LVTau;
        charge=Ntp->MCTau_charge(tau1_idx);
        for(int i=0;i<Ntp->NMCTauDecayProducts(tau1_idx);i++){
          if(abs(Ntp->MCTauandProd_pdgid(tau1_idx,i))==abs(PDGInfo::pi_plus)){
            TLorentzVector LVpi=Ntp->MCTauandProd_p4(tau1_idx,i);
            pipi+=LVpi;
            LVpi.Boost(-1*Zboson.BoostVector());
            if(charge<0){x1=LVpi.P()/LVTau.E(); }
            else{ x2=LVpi.P()/LVTau.E();}
          }
        }
	//
	double myzs=-999;
	if(tautau.M()!= 0) {
	  for(int i=0;i<zsbins;i++){
	    double zslow=((double)i)*(zsmax-zsmin)/((double)zsbins)+zsmin; 
	    double zsup=((double)i+1)*(zsmax-zsmin)/((double)zsbins)+zsmin;
	    double aup=Zstoa(zsup), alow=Zstoa(zslow);
	    if(x2-x1>alow && x2-x1<aup){
	      double zs=(zsup+zslow)/2;
	      pi_zs.at(t).Fill(zs,w);
	      pi_zs_Spin.at(t).Fill(zs,w*Spin_WT);
	      pi_zs_UnSpin.at(t).Fill(zs,w*UnSpin_WT);
	      pi_zs_FlipSpin.at(t).Fill(zs,w*FlipSpin_WT);
	      pi_zs_hplus.at(t).Fill(zs,w*hplus*UnSpin_WT);
	      pi_zs_hminus.at(t).Fill(zs,w*hminus*UnSpin_WT);
	      myzs=zs;
	      break;
	    }
	  }
	  pi_Mvis.at(t).Fill(pipi.M()/tautau.M(),w);
	  pi_Mvis_Spin.at(t).Fill(pipi.M()/tautau.M(),w*Spin_WT);
	  pi_Mvis_UnSpin.at(t).Fill(pipi.M()/tautau.M(),w*UnSpin_WT);
	  pi_Mvis_FlipSpin.at(t).Fill(pipi.M()/tautau.M(),w*FlipSpin_WT);
	  pi_Mvis_hplus.at(t).Fill(pipi.M()/tautau.M(),w*hplus*UnSpin_WT);
	  pi_Mvis_hminus.at(t).Fill(pipi.M()/tautau.M(),w*hminus*UnSpin_WT);
	  Logger(Logger::Verbose) << "pipi : mtau " << pipi.M()/tautau.M() << " zs " << myzs
			       << " w " << w << "spin " <<  Spin_WT << " unspin " << UnSpin_WT << " spin flip " << FlipSpin_WT << " hplus " << hplus << " hminus " << hminus  << std::endl;
	}
      }
      Logger(Logger::Verbose)<< "pi-pi done" <<std::endl;
    }
   if(Ntp->hasSignalTauDecay(PDGInfo::Z0,Boson_idx,TauDecay::JAK_PION,tau_idx)){
	   Logger(Logger::Verbose) << "pion" << std::endl;
      TLorentzVector Boson_LV=Ntp->MCSignalParticle_p4(Boson_idx);
      TLorentzVector Tau_LV(0,0,0,0);
      TLorentzVector X_LV(0,0,0,0);
      for(int i=0; i<Ntp->NMCTauDecayProducts(tau_idx);i++){
        if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::tau_minus)){
          Tau_LV=Ntp->MCTauandProd_p4(tau_idx,i);
        }
        else if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::pi_plus) ){
          X_LV+=Ntp->MCTauandProd_p4(tau_idx,i);
        }
      }
      if(Tau_LV.E()>0){
	pi_PtRatio_hplus.at(t).Fill(X_LV.Pt()/Tau_LV.Pt(),w*hplus*UnSpin_WT);
        pi_PtRatio_hminus.at(t).Fill(X_LV.Pt()/Tau_LV.Pt(),w*hminus*UnSpin_WT);

        Tau_LV.Boost(-Boson_LV.BoostVector());
        X_LV.Boost(-Boson_LV.BoostVector());
        Boson_LV.Boost(-Boson_LV.BoostVector());
        // Now fill results                                                                                                                                                                                                                  
        Ntp->TauSpinerSetSignal(Ntp->MCTau_charge(tau_idx));
	pi_ExoverEtau.at(t).Fill(X_LV.E()/Tau_LV.E(),w);
	pi_ExoverEtau_Spin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Spin_WT);
	pi_ExoverEtau_UnSpin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*UnSpin_WT);
	pi_ExoverEtau_FlipSpin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*FlipSpin_WT);
	pi_ExoverEtau_hplus.at(t).Fill(X_LV.E()/Tau_LV.E(),w*hplus*UnSpin_WT);
	pi_ExoverEtau_hminus.at(t).Fill(X_LV.E()/Tau_LV.E(),w*hminus*UnSpin_WT);
	pi_WT_Spin.at(t).Fill(Spin_WT,w);
	pi_WT_UnSpin.at(t).Fill(UnSpin_WT,w);
	pi_WT_FlipSpin.at(t).Fill(FlipSpin_WT,w);
	Logger(Logger::Verbose) << "pi : x " << X_LV.E()/Tau_LV.E() << " w " << w << "spin " <<  Spin_WT << " unspin " << UnSpin_WT << " spin flip " << FlipSpin_WT << " hplus " << hplus << " hminus " << hminus  << std::endl;
      }
      Logger(Logger::Verbose) << "pion done " <<std::endl;
    }
   if(Ntp->hasSignalTauDecay(PDGInfo::Z0,Boson_idx,TauDecay::JAK_A1_3PI,tau_idx)){
     if(Ntp->MCTau_charge(tau_idx)==-1){
    	 Logger(Logger::Verbose) << "a1to3pi" << std::endl;
      TLorentzVector Boson_LV=Ntp->MCSignalParticle_p4(Boson_idx);
      TLorentzVector Tau_LV(0,0,0,0);
      TLorentzVector X_LV(0,0,0,0);
      std::vector<TLorentzVector> pions;
      std::vector<float> pions_charge;
      for(int i=0; i<Ntp->NMCTauDecayProducts(tau_idx);i++){
        if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::tau_minus)){
          Tau_LV=Ntp->MCTauandProd_p4(tau_idx,i);
        }
        else if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::pi_plus) ||
                abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::pi0)
                ){
          X_LV+=Ntp->MCTauandProd_p4(tau_idx,i);
	  pions.push_back(Ntp->MCTauandProd_p4(tau_idx,i));
	  pions_charge.push_back(Ntp->MCTauandProd_charge(tau_idx,i));
        }
      }
      Logger(Logger::Verbose) << "a1to3pi - 1" << std::endl;
      if(Tau_LV.E()>0){
        a1_PtRatio_hplus.at(t).Fill(X_LV.Pt()/Tau_LV.Pt(),w*hplus*UnSpin_WT);
        a1_PtRatio_hminus.at(t).Fill(X_LV.Pt()/Tau_LV.Pt(),w*hminus*UnSpin_WT);

        Tau_LV.Boost(-Boson_LV.BoostVector());
        X_LV.Boost(-Boson_LV.BoostVector());
        Boson_LV.Boost(-Boson_LV.BoostVector());
        // Now fill results               
        Logger(Logger::Verbose) << "a1to3pi - 2 " << std::endl;
        Ntp->TauSpinerSetSignal(Ntp->MCTau_charge(tau_idx));
	a1_ExoverEtau.at(t).Fill(X_LV.E()/Tau_LV.E(),w);
	a1_ExoverEtau_Spin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Spin_WT);
	a1_ExoverEtau_UnSpin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*UnSpin_WT);
	a1_ExoverEtau_FlipSpin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*FlipSpin_WT);
	a1_ExoverEtau_hplus.at(t).Fill(X_LV.E()/Tau_LV.E(),w*hplus*UnSpin_WT);
	a1_ExoverEtau_hminus.at(t).Fill(X_LV.E()/Tau_LV.E(),w*hminus*UnSpin_WT);
	a1_WT_Spin.at(t).Fill(Spin_WT,w);
	a1_WT_UnSpin.at(t).Fill(UnSpin_WT,w);
	a1_WT_FlipSpin.at(t).Fill(FlipSpin_WT,w);

	Logger(Logger::Verbose) << "a1to3pi - 3 " << std::endl;

	TLorentzVector Boson_LV=Ntp->MCSignalParticle_p4(Boson_idx);
	int charge=Ntp->MCTau_charge(tau_idx);
	TLorentzVector a1(0,0,0,0),pi(0,0,0,0);
	for(int i=0;i<Ntp->NMCTauDecayProducts(tau_idx);i++){
	  TLorentzVector LV=Ntp->MCTauandProd_p4(tau_idx,i);
	  if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::pi_plus) && Ntp->MCTauandProd_pdgid(tau_idx,i)/abs(Ntp->MCTauandProd_pdgid(tau_idx,i))!=charge){pi+=LV; a1+=LV;}
          else if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::pi0) || abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::pi_plus)){a1+=LV;}
	  else if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))!=abs(PDGInfo::nu_tau)){Logger(Logger::Info) << Ntp->MCTauandProd_pdgid(tau_idx,i) << std::endl;}
	}
	Logger(Logger::Verbose) << "a1to3pi - 3 " << std::endl;
	double gamma=2*pi.P()/a1.P()-1;
     
        a1_Gamma.at(t).Fill(gamma,w);
        a1_Gamma_Spin.at(t).Fill(gamma,w*Spin_WT);
        a1_Gamma_UnSpin.at(t).Fill(gamma,w*UnSpin_WT);
        a1_Gamma_FlipSpin.at(t).Fill(gamma,w*FlipSpin_WT);
        a1_Gamma_hplus.at(t).Fill(gamma,w*hplus*UnSpin_WT);
        a1_Gamma_hminus.at(t).Fill(gamma,w*hminus*UnSpin_WT);
        Logger(Logger::Verbose) << "a1to3pi : gamma " << gamma << " w " << w << "spin " <<  Spin_WT << " unspin " << UnSpin_WT << " spin flip " << FlipSpin_WT << " hplus " << hplus << " hminus " << hminus  << std::endl;
      }
      Logger(Logger::Verbose)<< "a1 done " <<std::endl;
     }
    }
    if(Ntp->hasSignalTauDecay(PDGInfo::Z0,Boson_idx,TauDecay::JAK_RHO_PIPI0,tau_idx)){
      if(Ntp->MCTau_charge(tau_idx)==-1){
      Logger(Logger::Verbose) << "rho" << std::endl;
      TLorentzVector Boson_LV=Ntp->MCSignalParticle_p4(Boson_idx);
      TLorentzVector Tau_LV(0,0,0,0);
      TLorentzVector X_LV(0,0,0,0);
      for(int i=0; i<Ntp->NMCTauDecayProducts(tau_idx);i++){
	  if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::tau_minus)){
	    Tau_LV=Ntp->MCTauandProd_p4(tau_idx,i);
	  }
	  else if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::pi_plus) ||
		  abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::pi0)
		  ){
	    X_LV+=Ntp->MCTauandProd_p4(tau_idx,i);
	  }
      }
      if(Tau_LV.E()>0){
        rho_PtRatio_hplus.at(t).Fill(X_LV.Pt()/Tau_LV.Pt(),w*hplus*UnSpin_WT);
        rho_PtRatio_hminus.at(t).Fill(X_LV.Pt()/Tau_LV.Pt(),w*hminus*UnSpin_WT);

	Tau_LV.Boost(-Boson_LV.BoostVector());
	X_LV.Boost(-Boson_LV.BoostVector());
	Boson_LV.Boost(-Boson_LV.BoostVector());
	// Now fill results                                                                                                                                                                                                                  
	Ntp->TauSpinerSetSignal(Ntp->MCTau_charge(tau_idx));
	rho_ExoverEtau.at(t).Fill(X_LV.E()/Tau_LV.E(),w);
	rho_ExoverEtau_Spin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*Spin_WT);
	rho_ExoverEtau_UnSpin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*UnSpin_WT);
	rho_ExoverEtau_FlipSpin.at(t).Fill(X_LV.E()/Tau_LV.E(),w*FlipSpin_WT);
	rho_ExoverEtau_hplus.at(t).Fill(X_LV.E()/Tau_LV.E(),w*hplus*UnSpin_WT);
	rho_ExoverEtau_hminus.at(t).Fill(X_LV.E()/Tau_LV.E(),w*hminus*UnSpin_WT);

        TLorentzVector Boson_LV=Ntp->MCSignalParticle_p4(Boson_idx);
        //int charge=Ntp->MCTau_charge(tau_idx);
        TLorentzVector rho(0,0,0,0),pi(0,0,0,0);
        for(int i=0;i<Ntp->NMCTauDecayProducts(tau_idx);i++){
          TLorentzVector LV=Ntp->MCTauandProd_p4(tau_idx,i);
	  //std::cout << Ntp->MCTauandProd_pdgid(tau_idx,i) << " " << LV.Px() << " " << LV.Py() << " " << LV.Pz() << " " << LV.E() << std::endl;
          if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::pi_plus)){pi+=LV; rho+=LV;}
          if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::pi0)){rho+=LV;}
        }
	double gamma=2*pi.P()/rho.P()-1;
	rho_Gamma.at(t).Fill(gamma,w);
        rho_Gamma_Spin.at(t).Fill(gamma,w*Spin_WT);
        rho_Gamma_UnSpin.at(t).Fill(gamma,w*UnSpin_WT);
        rho_Gamma_FlipSpin.at(t).Fill(gamma,w*FlipSpin_WT);
        rho_Gamma_hplus.at(t).Fill(gamma,w*hplus*UnSpin_WT);
        rho_Gamma_hminus.at(t).Fill(gamma,w*hminus*UnSpin_WT);

        Logger(Logger::Verbose) << "rho : gamma " << gamma << " w " << w << "spin " <<  Spin_WT << " unspin " << UnSpin_WT << " spin flip " << FlipSpin_WT << " hplus " << hplus << " hminus " << hminus  << " rho " << rho.E() << std::endl;
	rho_WT_Spin.at(t).Fill(Spin_WT,w);
	rho_WT_UnSpin.at(t).Fill(UnSpin_WT,w);
	rho_WT_FlipSpin.at(t).Fill(FlipSpin_WT,w);
      }
    }
      Logger(Logger::Verbose) << "rho done" << std::endl;
    }
  }
  Logger(Logger::Verbose) << "done" << std::endl;
}




void  TauSpinExample::Finish(){
  unsigned int tdata=0;
  TF1 fmu_ExoverEtau_hplus( "fmu_ExoverEtau_hplus", "5.0/3.0-3.0*(x^2.0)+(4.0/3.0)*x^3.0-1*(-1/3+3*x^2-8/3*x^3)",0,1);
  TF1 fmu_ExoverEtau_hminus("fmu_ExoverEtau_hminus","5.0/3.0-3.0*(x^2.0)+(4.0/3.0)*x^3.0+1*(-1/3+3*x^2-8/3*x^3)",0,1);
  TF1 fpi_ExoverEtau_hplus( "fpi_ExoverEtau_hplus", "x",0,1);
  TF1 fpi_ExoverEtau_hminus("fpi_ExoverEtau_hminus","1-x",0,1);
  for(int i=1;i<mu_ExoverEtau_hplus.at(tdata).GetNbinsX();i++){
    double min=mu_ExoverEtau_hplus.at(tdata).GetBinLowEdge(i);
    double max=min+mu_ExoverEtau_hplus.at(tdata).GetBinWidth(i);
    mu_ExoverEtau_hplus.at(tdata).SetBinContent(i,fmu_ExoverEtau_hplus.Integral(min,max));
    mu_ExoverEtau_hplus.at(tdata).SetBinError(i,0.00001);
    mu_ExoverEtau_hminus.at(tdata).SetBinContent(i,fmu_ExoverEtau_hminus.Integral(min,max));
    mu_ExoverEtau_hminus.at(tdata).SetBinError(i,0.00001);
    pi_ExoverEtau_hplus.at(tdata).SetBinContent(i,fpi_ExoverEtau_hplus.Integral(min,max));
    pi_ExoverEtau_hplus.at(tdata).SetBinError(i,0.00001);
    pi_ExoverEtau_hminus.at(tdata).SetBinContent(i,fpi_ExoverEtau_hminus.Integral(min,max));
    pi_ExoverEtau_hminus.at(tdata).SetBinError(i,0.00001);
  }
  Selection::Finish();
}



double  TauSpinExample::Zstoa(double zs){
  double a=1-sqrt(fabs(1.0-2*fabs(zs)));
  if(zs<0){
    a*=-1.0;
  }
  return a;
}
