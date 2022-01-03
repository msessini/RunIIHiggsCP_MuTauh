#include "TH1F.h"
#include "TH1D.h"
#include "TFile.h"

void lumi_script(){

  // what you need:
  // - JSON file from collision data
  // - pileup JSON file
  // then produce a pileup distribution from data with:
  // pileupCalc.py -i Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt --inputLumiJSON pileup_JSON_DCSONLY_190389-208686_All_2012_pixelcorr.txt --calcMode true --minBiasXsec 69400 --maxPileupBin 60 --numPileupBins 60 datapileup.root
  // also watch out for possible changes in the xsec and so forth... !

  // the output file needs to be copied/moved to $CMSSW_base/src/data/

  // ------------------- MC distribution -------------------
  TH1D* MC_Summer12 = new TH1D("MC_Summer12","MC_Summer12",60,0,60);

  // values below from https://twiki.cern.ch/twiki/bin/view/CMS/Pileup_MC_Gen_Scenarios
  Double_t Summer2012_S10[60] = {
                         2.560E-06,
                         5.239E-06,
                         1.420E-05,
                         5.005E-05,
                         1.001E-04,
                         2.705E-04,
                         1.999E-03,
                         6.097E-03,
                         1.046E-02,
                         1.383E-02,
                         1.685E-02,
                         2.055E-02,
                         2.572E-02,
                         3.262E-02,
                         4.121E-02,
                         4.977E-02,
                         5.539E-02,
                         5.725E-02,
                         5.607E-02,
                         5.312E-02,
                         5.008E-02,
                         4.763E-02,
                         4.558E-02,
                         4.363E-02,
                         4.159E-02,
                         3.933E-02,
                         3.681E-02,
                         3.406E-02,
                         3.116E-02,
                         2.818E-02,
                         2.519E-02,
                         2.226E-02,
                         1.946E-02,
                         1.682E-02,
                         1.437E-02,
                         1.215E-02,
                         1.016E-02,
                         8.400E-03,
                         6.873E-03,
                         5.564E-03,
                         4.457E-03,
                         3.533E-03,
                         2.772E-03,
                         2.154E-03,
                         1.656E-03,
                         1.261E-03,
                         9.513E-04,
                         7.107E-04,
                         5.259E-04,
                         3.856E-04,
                         2.801E-04,
                         2.017E-04,
                         1.439E-04,
                         1.017E-04,
                         7.126E-05,
                         4.948E-05,
                         3.405E-05,
                         2.322E-05,
                         1.570E-05,
                         5.005E-06};


 for(int iBin = 1; iBin < 61; iBin++){
  MC_Summer12->SetBinContent(iBin, Summer2012_S10[iBin-1]);
 } 

 // ------------------- data distribution -------------------
 TFile *_file0 = TFile::Open("datapileup.root");
 TH1D* DataPU =  (TH1D *)_file0->Get("pileup");

 DataPU->SetName("h_190456_20868");

 double integral = DataPU->Integral();
 DataPU->Scale(1./integral);

 // write root file
 TFile *f_out = new TFile("Lumi_190456_208686MC_PU_S10_andData.root","RECREATE");

 MC_Summer12->Write();
 DataPU->Write();
 f_out->Write();
 f_out->Close();
}
