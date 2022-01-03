/*
 * DataMCCorrections.cxx
 *
 *  Created on: Nov  16, 2017
 *      Author: cherepanov
 */

#include "DataMCCorrections.h"
#include "Selection_Base.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

///////////////////////////
//
// Constructor
//

DataMCCorrections::DataMCCorrections(bool load_ZPtWeights, bool load_LeptonEff, bool load_LeptonIso){

	// define which scale factors should be loaded
	// individual SF can be switched off to avoid opening files which are not necessary
	loadZPtWeights = load_ZPtWeights;
	loadLeptonEff = load_LeptonEff;
	loadLeptonIso = load_LeptonIso;

	// define location of input root files
	TString basedir = "";
	basedir = (TString)std::getenv("workdir")+"/Code/CommonFiles/";

	// Z Pt Weights
	if(loadZPtWeights){
	  ZPtWeightFile= new TFile(basedir+"weights/zpt_weights_2016.root", "READ");
	  m_zPtHist = (TH2D*)(ZPtWeightFile->Get("zptmass_histo"));
	  m_zPtHistErr = (TH2D*)(ZPtWeightFile->Get("zptmass_histo_err"));
	}
	if(loadLeptonEff){
	  TString inputRootFile = basedir+"LeptonEfficiencies/Muon/Run2016BtoH/Muon_IsoMu24_OR_TkIsoMu24_2016BtoH_eff.root";

	  TFile *fileIn = new TFile(inputRootFile, "READ");
	  // if root file not found
	  if (fileIn->IsZombie() ) { std::cout << "ERROR in ScaleFactor::init_ScaleFactor(TString inputRootFile) from NTupleMaker/src/ScaleFactor.cc : ‎File " <<inputRootFile << " does not exist. Please check. " <<std::endl; exit(1); };
	  //	std::cout<<"check 2  "<<std::endl;
	  std::string HistoBaseName = "ZMass";
	  etaBinsH =(TH1D*)(fileIn->Get("etaBinsH"));
	  
	  std::string etaLabel, GraphName;
	  //	std::cout<<"check 23  "<<std::endl;
	  int nEtaBins = etaBinsH->GetNbinsX();
	  for (int iBin=0; iBin<nEtaBins; iBin++){
	    etaLabel = etaBinsH->GetXaxis()->GetBinLabel(iBin+1);
	    GraphName = HistoBaseName+etaLabel+"_Data";
	    eff_data[etaLabel] = (TGraphAsymmErrors*)fileIn->Get(TString(GraphName)); 
	    SetAxisBins(eff_data[etaLabel]);
	    GraphName = HistoBaseName+etaLabel+"_MC";
	    eff_mc[etaLabel] = (TGraphAsymmErrors*)fileIn->Get(TString(GraphName));
	    SetAxisBins(eff_mc[etaLabel]); 
	    bool sameBinning = check_SameBinning(eff_data[etaLabel], eff_mc[etaLabel]);
	    if (!sameBinning) {std::cout<< "ERROR in ScaleFactor::init_ScaleFactor(TString inputRootFile) from LepEffInterface/src/ScaleFactor.cc . Can not proceed because ScaleFactor::check_SameBinning returned different pT binning for data and MC for eta label " << etaLabel << std::endl; exit(1); }; 
	  }
	}
	if(loadLeptonIso){
	  TString inputRootFile = basedir+"LeptonEfficiencies/Muon/Run2016BtoH/Muon_IdIso_IsoLt0p15_2016BtoH_eff.root";

	  TFile *fileIn = new TFile(inputRootFile, "READ");
	  // if root file not found
	  if (fileIn->IsZombie() ) { std::cout << "ERROR in ScaleFactor::init_ScaleFactor(TString inputRootFile) from NTupleMaker/src/ScaleFactor.cc : ‎File " <<inputRootFile << " does not exist. Please check. " <<std::endl; exit(1); };
	  //	std::cout<<"check 2  "<<std::endl;
	  std::string HistoBaseName = "ZMass";
	  etaBinsH =(TH1D*)(fileIn->Get("etaBinsH"));
	  
	  std::string etaLabel, GraphName;
	  //	std::cout<<"check 23  "<<std::endl;
	  int nEtaBins = etaBinsH->GetNbinsX();
	  for (int iBin=0; iBin<nEtaBins; iBin++){
	    etaLabel = etaBinsH->GetXaxis()->GetBinLabel(iBin+1);
	    GraphName = HistoBaseName+etaLabel+"_Data";
	    eff_data[etaLabel] = (TGraphAsymmErrors*)fileIn->Get(TString(GraphName)); 
	    SetAxisBins(eff_data[etaLabel]);
	    GraphName = HistoBaseName+etaLabel+"_MC";
	    eff_mc[etaLabel] = (TGraphAsymmErrors*)fileIn->Get(TString(GraphName));
	    SetAxisBins(eff_mc[etaLabel]); 
	    bool sameBinning = check_SameBinning(eff_data[etaLabel], eff_mc[etaLabel]);
	    if (!sameBinning) {std::cout<< "ERROR in ScaleFactor::init_ScaleFactor(TString inputRootFile) from LepEffInterface/src/ScaleFactor.cc . Can not proceed because ScaleFactor::check_SameBinning returned different pT binning for data and MC for eta label " << etaLabel << std::endl; exit(1); }; 
	  }
	}	

	
} 

DataMCCorrections::~DataMCCorrections(){
 }

///////////////////////////
//
// Muon scale factors
//

float DataMCCorrections::ZPTWeight(float genMass, float genPt){
  float zPtReweight = m_zPtHist->GetBinContent(m_zPtHist->GetXaxis()->FindBin(genMass),m_zPtHist->GetYaxis()->FindBin(genPt));
  return zPtReweight;
}

float DataMCCorrections::ZPTWeightErr(float genMass, float genPt){
  float zPtReweightErr = m_zPtHistErr->GetBinContent(m_zPtHistErr->GetXaxis()->FindBin(genMass),m_zPtHistErr->GetYaxis()->FindBin(genPt));
  return zPtReweightErr;
}

float DataMCCorrections::AgainstElectronDataMCCorrection(TLorentzVector p4, TString type){
  float scale(1);
  if(type=="AgainstElectronMVAVLoose"){
    if(fabs(p4.Eta()) < 1.460)scale=1.213;
    else if(fabs(p4.Eta()) > 1.558) scale = 1.375;
  }
  if(type=="AgainstElectronMVALoose"){
    if(fabs(p4.Eta()) < 1.460)scale=1.320;
    else if(fabs(p4.Eta()) > 1.558) scale = 1.380;
  }
  if(type=="AgainstElectronMVAMedium"){
    if(fabs(p4.Eta()) < 1.460)scale=1.323;
    else if(fabs(p4.Eta()) > 1.558) scale = 1.527;
  }
  if(type=="AgainstElectronMVATight"){
    if(fabs(p4.Eta()) < 1.460)scale=1.402;
    else if(fabs(p4.Eta()) > 1.558) scale = 1.900;
  }
  if(type=="AgainstElectronMVAVTight"){
    if(fabs(p4.Eta()) < 1.460)scale=1.207;
    else if(fabs(p4.Eta()) > 1.558) scale = 1.968;
  }
  return scale;
}
// float DataMCCorrections::LeptonTriggerEfficiencyScaleFactor(TLorentzVector leptonP4){
//   float weight = LeptonTriggerEff-> get_ScaleFactor(leptonP4.Pt(),leptonP4.Eta()); 
//   return weight;
// }


float DataMCCorrections::AgainstMuonDataMCCorrection(TLorentzVector p4, TString type){
  float scale(1);
  if(type=="AgainstMuonMVALoose3"){
    if(fabs(p4.Eta()) <0.4)scale=1.010;
    else if(fabs(p4.Eta()) > 0.4 && fabs(p4.Eta()) < 0.8) scale = 1.007;
    else if(fabs(p4.Eta()) > 0.8 && fabs(p4.Eta()) < 1.2) scale = 0.870;
    else if(fabs(p4.Eta()) > 1.2 && fabs(p4.Eta()) < 1.7) scale = 1.154;
    else if(fabs(p4.Eta()) > 1.7 && fabs(p4.Eta()) < 2.3) scale = 2.281;
  }
  if(type=="AgainstMuonMVATight3"){
    if(fabs(p4.Eta()) <0.4)scale=1.263;
    else if(fabs(p4.Eta()) > 0.4 && fabs(p4.Eta()) < 0.8) scale = 1.364;
    else if(fabs(p4.Eta()) > 0.8 && fabs(p4.Eta()) < 1.2) scale = 0.854;
    else if(fabs(p4.Eta()) > 1.2 && fabs(p4.Eta()) < 1.7) scale = 1.712;
    else if(fabs(p4.Eta()) > 1.7 && fabs(p4.Eta()) < 2.3) scale = 2.324;
  }
  return scale;
}




void DataMCCorrections::SetAxisBins(TGraphAsymmErrors* graph) {

	int NPOINTS = graph->GetN(); 
	double AXISBINS[NPOINTS+1];
	for (int i=0; i<NPOINTS; i++) { AXISBINS[i] = (graph->GetX()[i] - graph->GetErrorXlow(i)); }
	AXISBINS[NPOINTS] = (graph->GetX()[NPOINTS-1] + graph->GetErrorXhigh(NPOINTS-1));
	graph->GetXaxis()->Set(NPOINTS, AXISBINS);
	return;
}

bool DataMCCorrections::check_SameBinning(TGraphAsymmErrors* graph1, TGraphAsymmErrors* graph2){
	bool haveSameBins = false;
	int n1 = graph1->GetXaxis()->GetNbins();
	int n2 = graph2->GetXaxis()->GetNbins();
	if (n1 != n2 ) {return false;}
	else {
		haveSameBins = true;
		const int nbins = n1;
		double x1, x2;
		for (int i=0; i<nbins; i++){ 
			x1 = (graph1->GetXaxis()->GetXbins())->GetArray()[i];
			x2 = (graph2->GetXaxis()->GetXbins())->GetArray()[i]; 
			haveSameBins = haveSameBins and (x1== x2) ;
		}
	}

	return haveSameBins;
}



std::string DataMCCorrections::FindEtaLabel(double Eta){
	Eta = fabs(Eta);
	int binNumber = etaBinsH->GetXaxis()->FindFixBin(Eta);
	std::string EtaLabel = etaBinsH->GetXaxis()->GetBinLabel(binNumber);
	std::map<std::string, TGraphAsymmErrors*>::iterator it;
	it =  eff_data.find(EtaLabel);
	if ( it == eff_data.end()) { 
	std::cout << "ERROR in DataMCCorrections::get_EfficiencyData(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc : no object corresponding to eta label "<< EtaLabel << " for data " << std::endl; exit(1);
	}
	else return EtaLabel;
}


int DataMCCorrections::FindPtBin( std::map<std::string, TGraphAsymmErrors *> eff_map, std::string EtaLabel, double Pt){

        int Npoints = eff_map[EtaLabel]->GetN();
	double ptMAX = (eff_map[EtaLabel]->GetX()[Npoints-1])+(eff_map[EtaLabel]->GetErrorXhigh(Npoints-1));
	double ptMIN = (eff_map[EtaLabel]->GetX()[0])-(eff_map[EtaLabel]->GetErrorXlow(0));
	// if pt is overflow, return last pt bin
 	if (Pt >= ptMAX ) return Npoints; 
	// if pt is underflow, return nonsense number and warning
	else if (Pt < ptMIN){ 
 	std::cout<< "WARNING in DataMCCorrections::get_EfficiencyData(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc: pT too low (pt = " << Pt << "), min value is " << ptMIN << ". Returned efficiency =1. Weight will be 1. " << std::endl;
	return -99;}
	// if pt is in range
	else {return eff_map[EtaLabel]->GetXaxis()->FindFixBin(Pt);} 
	}


double DataMCCorrections::get_EfficiencyData(double pt, double eta){

        double eff;
	std::string label = FindEtaLabel(eta);

	int ptbin = FindPtBin(eff_data, label, pt); 
	if (ptbin == -99){eff =1;} // if pt is underflow 
	else eff = eff_data[label]->GetY()[ptbin-1];

	if (eff > 1.) {std::cout<< "WARNING in DataMCCorrections::get_EfficiencyData(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc: Returned efficiency in data > 1. " << std::endl;} 
	if (eff < 0 ) {std::cout<<"WARNING in DataMCCorrections::get_EfficiencyData(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc: Returned negative efficiency in data" <<std::endl;}

	return eff;
	
}


double DataMCCorrections::get_EfficiencyMC(double pt, double eta) {

	double eff;		
	std::string label = FindEtaLabel(eta);

	int ptbin = FindPtBin(eff_mc, label, pt); 
	if (ptbin == -99){eff =1;} // if pt is underflow 
	else eff= eff_mc[label]->GetY()[ptbin-1];

	if (eff > 1. ) {std::cout << "WARNING in DataMCCorrections::get_EfficiencyMC(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc : Returned efficiency in MC > 1. " << std::endl;} 		
	if (eff < 0 ) {std::cout<<"WARNING in DataMCCorrections::get_EfficiencyMC(double pt, double eta) from LepEffIntrface/src/ScaleFactor.cc : Returned negative efficiency in MC. " <<std::endl;}
	

	return eff;

}



double DataMCCorrections::get_ScaleFactor(double pt, double eta){
	
	double efficiency_data = get_EfficiencyData(pt, eta);
	double efficiency_mc = get_EfficiencyMC(pt, eta);
	double SF;

	if ( efficiency_mc != 0) {SF = efficiency_data/efficiency_mc;}
	else {
	SF=0.; std::cout << "WARNING in DataMCCorrections::get_ScaleFactor(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc : MC efficiency = 0. Scale Factor set to 0. ";
	}

	return SF;	
	
}


double DataMCCorrections::get_EfficiencyDataError(double pt, double eta){

	double eff_error;
	std::string label = FindEtaLabel(eta);
	int ptbin = FindPtBin(eff_data, label, pt); 
	if (ptbin == -99){eff_error =0.;} // if pt is underflow 
	else eff_error= eff_data[label]->GetErrorYhigh(ptbin-1); 
        // errors are supposed to be symmetric, can use GetErrorYhigh or GetErrorYlow

	double effData = get_EfficiencyData(pt,eta);
	if (eff_error > effData) eff_error = 0.5*effData;
	return eff_error;
}
	
	

double DataMCCorrections::get_EfficiencyMCError(double pt, double eta){

	double eff_error;
	std::string label = FindEtaLabel(eta);
	int ptbin = FindPtBin(eff_mc, label, pt); 
	if (ptbin == -99){eff_error =0.;} // if pt is underflow 
	else eff_error= eff_mc[label]->GetErrorYhigh(ptbin-1); 
	// errors are supposed to be symmetric, can use GetErrorYhigh or GetErrorYlow

	double effMC = get_EfficiencyMC(pt,eta);
	if (eff_error > effMC ) eff_error = 0.5*effMC;
	return eff_error;
}

double DataMCCorrections::get_ScaleFactorError(double pt, double eta){

	double SF_error = 0.;
	
	double effData = get_EfficiencyData(pt, eta);
	double effMC = get_EfficiencyMC(pt, eta);
	double errData = get_EfficiencyDataError(pt, eta);
	double errMC =  get_EfficiencyMCError(pt, eta);

	if (errData == 0) {std::cout<<"WARNING in DataMCCorrections::get_ScaleFactorError(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc: uncertainty on data point = 0, can not calculate uncertainty on scale factor. Uncertainty set to 0." << std::endl;}
	if (errMC ==0) {std::cout<<"WARNING in DataMCCorrections::get_ScaleFactorError(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc: uncertainty on MC = 0, can not calculate uncerttainty on scale factor. Uncertainty set to 0." << std::endl;}
	if (effData ==0) {std::cout<<"WARNING in DataMCCorrections::get_ScaleFactorError(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc: efficiency in data = 0, can not calculate uncertainty on scale factor. Uncertainty set to 0." << std::endl;}
	if (effMC ==0) {std::cout<<"WARNING in DataMCCorrections::get_ScaleFactorError(double pt, double eta) from LepEffInterface/src/ScaleFactor.cc: efficiency in MC = 0, can not calculate uncertainty on scale factor. Uncertainty set to 0." << std::endl;}
	else {	
	SF_error = pow((errData/effData),2) + pow((errMC/effMC),2);
	SF_error = pow(SF_error, 0.5)*(effData/effMC);
	}
	return SF_error;
}
	

