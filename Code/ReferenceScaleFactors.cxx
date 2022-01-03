/*
 * ReferenceScaleFactors.cxx
 *
 *  Created on: Aug 13, 2014
 *      Author: nehrkorn
 */

#include "ReferenceScaleFactors.h"
#include "Selection_Base.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

///////////////////////////
//
// Constructor
//

ReferenceScaleFactors::ReferenceScaleFactors(int runType, bool load_ElectronID, bool load_EMuTriggerEff, bool load_HiggsPtWeights){

	// define which scale factors should be loaded
	// individual SF can be switched off to avoid opening files which are not necessary
	loadElectronID = load_ElectronID;
	loadEMuTriggerEff = load_EMuTriggerEff;
	loadHiggsPtWeights = load_HiggsPtWeights;

	// define location of input root files
	TString basedir = "";
	if(runType==Selection_Base::GRID){
		basedir = (TString)std::getenv("PWD")+"/Code/CommonFiles/";
	}
	else if(runType==Selection_Base::Local){
		basedir = (TString)std::getenv("workdir")+"/Code/CommonFiles/";
	}

	// Electron Id's
	if(loadElectronID){
		// open root files
		ETrigIdEffFile = new TFile(basedir+"ElectronEfficiencies_Run2012ReReco_53X_Trig.root", "READ");
		ENonTrigIdEffFile = new TFile(basedir+"ElectronEfficiencies_Run2012ReReco_53X_NonTrig.root", "READ");
		ERecoEffFile = new TFile(basedir+"Electrons_ScaleFactors_Reco_8TeV.root", "READ");
		// load histograms
		ElectronTrigEff = (TH2D*)(ETrigIdEffFile->Get("electronsDATAMCratio_FO_ID_ISO"));
		ElectronNonTrigEff = (TH2D*)(ENonTrigIdEffFile->Get("h_electronScaleFactor_IdIsoSip"));
		ElectronRecoEff = (TH2D*)(ERecoEffFile->Get("h_electronScaleFactor_RECO"));
	}

	// Trigger efficiencies
	if(loadEMuTriggerEff){
		// open root files
		HWW_TriggerEfficiencies = new TFile(basedir+"TriggerEfficienciesWW_TH1D.root", "READ");
		// load histograms
		HiggsWW_EMu_SingleEle15 = (TH1D*)(HWW_TriggerEfficiencies->Get("SingleEle15"));
		HiggsWW_EMu_SingleEle25 = (TH1D*)(HWW_TriggerEfficiencies->Get("SingleEle25"));
		HiggsWW_EMu_DoubleEleLead15 = (TH1D*)(HWW_TriggerEfficiencies->Get("DoubleEleLead15"));
		HiggsWW_EMu_DoubleEleLead25 = (TH1D*)(HWW_TriggerEfficiencies->Get("DoubleEleLead25"));
		HiggsWW_EMu_DoubleEleTrail15 = (TH1D*)(HWW_TriggerEfficiencies->Get("DoubleEleTrail15"));
		HiggsWW_EMu_DoubleEleTrail25 = (TH1D*)(HWW_TriggerEfficiencies->Get("DoubleEleTrail25"));
		HiggsWW_EMu_SingleMu08 = (TH1D*)(HWW_TriggerEfficiencies->Get("SingleMu08"));
		HiggsWW_EMu_SingleMu12 = (TH1D*)(HWW_TriggerEfficiencies->Get("SingleMu12"));
		HiggsWW_EMu_SingleMu21 = (TH1D*)(HWW_TriggerEfficiencies->Get("SingleMu21"));
		HiggsWW_EMu_SingleMu25 = (TH1D*)(HWW_TriggerEfficiencies->Get("SingleMu25"));
		HiggsWW_EMu_DoubleMuLead12 = (TH1D*)(HWW_TriggerEfficiencies->Get("DoubleMuLead12"));
		HiggsWW_EMu_DoubleMuLead21 = (TH1D*)(HWW_TriggerEfficiencies->Get("DoubleMuLead21"));
		HiggsWW_EMu_DoubleMuLead25 = (TH1D*)(HWW_TriggerEfficiencies->Get("DoubleMuLead25"));
		HiggsWW_EMu_DoubleMuTrail12 = (TH1D*)(HWW_TriggerEfficiencies->Get("DoubleMuTrail12"));
		HiggsWW_EMu_DoubleMuTrail21 = (TH1D*)(HWW_TriggerEfficiencies->Get("DoubleMuTrail21"));
		HiggsWW_EMu_DoubleMuTrail25 = (TH1D*)(HWW_TriggerEfficiencies->Get("DoubleMuTrail25"));
	}

	// Higgs Pt Weights
	if(loadHiggsPtWeights){
		// loop over masses
		for (unsigned mass = 100; mass <= 160; mass += 5){
			// open root files
			// files are taken from
			HiggsPtWeightFiles.insert(
					std::pair<unsigned, TFile*>( mass,
							new TFile(basedir + "HiggsPtWeights/weight_ptH_"+ TString::UItoa(mass,10) + "_8TeV.root", "READ")));

			// load histograms
			HiggsPtWeightsNominal.insert(
					std::pair<unsigned, TH1D*>( mass,
							(TH1D*)(HiggsPtWeightFiles[mass]->Get("Nominal"))));
			HiggsPtWeightsDown.insert(
					std::pair<unsigned, TH1D*>( mass,
							(TH1D*)(HiggsPtWeightFiles[mass]->Get("Down"))));
			HiggsPtWeightsUp.insert(
					std::pair<unsigned, TH1D*>( mass,
							(TH1D*)(HiggsPtWeightFiles[mass]->Get("Up"))));
		}
	}
}

ReferenceScaleFactors::~ReferenceScaleFactors(){
	if(loadElectronID){
		delete ETrigIdEffFile;
		delete ENonTrigIdEffFile;
		delete ERecoEffFile;
	}
	if(loadEMuTriggerEff){
		delete HWW_TriggerEfficiencies;
	}
	if(loadHiggsPtWeights){
	    for(std::map<unsigned, TFile*>::iterator itr = HiggsPtWeightFiles.begin(); itr != HiggsPtWeightFiles.end(); itr++){
	        delete itr->second;
	    }
	}
}

///////////////////////////
//
// Muon scale factors
//

double ReferenceScaleFactors::MuonIdTight2012(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	double eff = 1.;
	if(pt>=10. && pt<20.){
		if(eta<0.9) eff = 0.970;
		if(eta>=0.9 && eta<1.2) eff = 1.002;
		if(eta>=1.2 && eta<2.1) eff = 1.018;
		if(eta>=2.1 && eta<2.4) eff = 1.005;
	}
	if(pt>=20. && pt<25.){
		if(eta<0.9) eff = 0.989;
		if(eta>=0.9 && eta<1.2) eff = 0.994;
		if(eta>=1.2 && eta<2.1) eff = 1.000;
		if(eta>=2.1 && eta<2.4) eff = 0.998;
	}
	if(pt>=25. && pt<30.){
		if(eta<0.9) eff = 0.992;
		if(eta>=0.9 && eta<1.2) eff = 0.995;
		if(eta>=1.2 && eta<2.1) eff = 0.998;
		if(eta>=2.1 && eta<2.4) eff = 0.996;
	}
	if(pt>=30. && pt<35.){
		if(eta<0.9) eff = 0.993;
		if(eta>=0.9 && eta<1.2) eff = 0.993;
		if(eta>=1.2 && eta<2.1) eff = 0.997;
		if(eta>=2.1 && eta<2.4) eff = 1.001;
	}
	if(pt>=35. && pt<40.){
		if(eta<0.9) eff = 0.994;
		if(eta>=0.9 && eta<1.2) eff = 0.992;
		if(eta>=1.2 && eta<2.1) eff = 0.996;
		if(eta>=2.1 && eta<2.4) eff = 0.993;
	}
	if(pt>=40. && pt<50.){
		if(eta<0.9) eff = 0.992;
		if(eta>=0.9 && eta<1.2) eff = 0.992;
		if(eta>=1.2 && eta<2.1) eff = 0.996;
		if(eta>=2.1 && eta<2.4) eff = 0.995;
	}
	if(pt>=50. && pt<60.){
		if(eta<0.9) eff = 0.991;
		if(eta>=0.9 && eta<1.2) eff = 0.995;
		if(eta>=1.2 && eta<2.1) eff = 0.995;
		if(eta>=2.1 && eta<2.4) eff = 0.994;
	}
	if(pt>=60. && pt<90.){
		if(eta<0.9) eff = 0.989;
		if(eta>=0.9 && eta<1.2) eff = 0.990;
		if(eta>=1.2 && eta<2.1) eff = 0.992;
		if(eta>=2.1 && eta<2.4) eff = 0.989;
	}
	if(pt>=90. && pt<140.){
		if(eta<0.9) eff = 1.004;
		if(eta>=0.9 && eta<1.2) eff = 1.009;
		if(eta>=1.2 && eta<2.1) eff = 1.023;
		if(eta>=2.1 && eta<2.4) eff = 1.060;
	}
	if(pt>=140. && pt<300.){
		if(eta<0.9) eff = 1.019;
		if(eta>=0.9 && eta<1.2) eff = 1.011;
		if(eta>=1.2 && eta<2.1) eff = 0.975;
		if(eta>=2.1 && eta<2.4) eff = 0.891;
	}
	return eff;
}

double ReferenceScaleFactors::MuonIdUncTight2012(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	double err = 0.;
	if(pt>=10. && pt<20.){
		if(eta<0.9) err = 0.005;
		if(eta>=0.9 && eta<1.2) err = 0.007;
		if(eta>=1.2 && eta<2.1) err = 0.004;
		if(eta>=2.1 && eta<2.4) err = 0.007;
	}
	if(pt>=20. && pt<25.){
		if(eta<0.9) err = 0.0016;
		if(eta>=0.9 && eta<1.2) err = 0.0025;
		if(eta>=1.2 && eta<2.1) err = 0.0014;
		if(eta>=2.1 && eta<2.4) err = 0.003;
	}
	if(pt>=25. && pt<30.){
		if(eta<0.9) err = 0.0008;
		if(eta>=0.9 && eta<1.2) err = 0.0014;
		if(eta>=1.2 && eta<2.1) err = 0.0009;
		if(eta>=2.1 && eta<2.4) err = 0.0018;
	}
	if(pt>=30. && pt<35.){
		if(eta<0.9) err = 0.0005;
		if(eta>=0.9 && eta<1.2) err = 0.001;
		if(eta>=1.2 && eta<2.1) err = 0.0007;
		if(eta>=2.1 && eta<2.4) err = 0.0015;
	}
	if(pt>=35. && pt<40.){
		if(eta<0.9) err = 0.0004;
		if(eta>=0.9 && eta<1.2) err = 0.0008;
		if(eta>=1.2 && eta<2.1) err = 0.0006;
		if(eta>=2.1 && eta<2.4) err = 0.0013;
	}
	if(pt>=40. && pt<50.){
		if(eta<0.9) err = 0.0003;
		if(eta>=0.9 && eta<1.2) err = 0.0005;
		if(eta>=1.2 && eta<2.1) err = 0.00024;
		if(eta>=2.1 && eta<2.4) err = 0.001;
	}
	if(pt>=50. && pt<60.){
		if(eta<0.9) err = 0.0007;
		if(eta>=0.9 && eta<1.2) err = 0.0012;
		if(eta>=1.2 && eta<2.1) err = 0.0009;
		if(eta>=2.1 && eta<2.4) err = 0.0026;
	}
	if(pt>=60. && pt<90.){
		if(eta<0.9) err = 0.001;
		if(eta>=0.9 && eta<1.2) err = 0.0019;
		if(eta>=1.2 && eta<2.1) err = 0.0015;
		if(eta>=2.1 && eta<2.4) err = 0.005;
	}
	if(pt>=90. && pt<140.){
		if(eta<0.9) err = 0.003;
		if(eta>=0.9 && eta<1.2) err = 0.0062;
		if(eta>=1.2 && eta<2.1) err = 0.0054;
		if(eta>=2.1 && eta<2.4) err = 0.013;
	}
	if(pt>=140. && pt<300.){
		if(eta<0.9) err = 0.017;
		if(eta>=0.9 && eta<1.2) err = 0.03;
		if(eta>=1.2 && eta<2.1) err = 0.03;
		if(eta>=2.1 && eta<2.4) err = 0.14;
	}
	return err;
}

double ReferenceScaleFactors::MuonIsoTight2012(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	double eff = 1.;
	if(pt>=10. && pt<20.){
		if(eta<0.9) eff = 0.940;
		if(eta>=0.9 && eta<1.2) eff = 0.948;
		if(eta>=1.2 && eta<2.1) eff = 0.972;
		if(eta>=2.1 && eta<2.4) eff = 1.117;
	}
	if(pt>=20. && pt<25.){
		if(eta<0.9) eff = 0.977;
		if(eta>=0.9 && eta<1.2) eff = 0.986;
		if(eta>=1.2 && eta<2.1) eff = 0.990;
		if(eta>=2.1 && eta<2.4) eff = 1.116;
	}
	if(pt>=25. && pt<30.){
		if(eta<0.9) eff = 0.996;
		if(eta>=0.9 && eta<1.2) eff = 1.000;
		if(eta>=1.2 && eta<2.1) eff = 1.003;
		if(eta>=2.1 && eta<2.4) eff = 1.097;
	}
	if(pt>=30. && pt<35.){
		if(eta<0.9) eff = 0.993;
		if(eta>=0.9 && eta<1.2) eff = 1.000;
		if(eta>=1.2 && eta<2.1) eff = 1.004;
		if(eta>=2.1 && eta<2.4) eff = 1.075;
	}
	if(pt>=35. && pt<40.){
		if(eta<0.9) eff = 0.994;
		if(eta>=0.9 && eta<1.2) eff = 0.999;
		if(eta>=1.2 && eta<2.1) eff = 1.002;
		if(eta>=2.1 && eta<2.4) eff = 1.061;
	}
	if(pt>=40. && pt<50.){
		if(eta<0.9) eff = 0.994;
		if(eta>=0.9 && eta<1.2) eff = 0.999;
		if(eta>=1.2 && eta<2.1) eff = 1.001;
		if(eta>=2.1 && eta<2.4) eff = 1.034;
	}
	if(pt>=50. && pt<60.){
		if(eta<0.9) eff = 0.996;
		if(eta>=0.9 && eta<1.2) eff = 0.998;
		if(eta>=1.2 && eta<2.1) eff = 1.000;
		if(eta>=2.1 && eta<2.4) eff = 1.025;
	}
	if(pt>=60. && pt<90.){
		if(eta<0.9) eff = 0.999;
		if(eta>=0.9 && eta<1.2) eff = 0.999;
		if(eta>=1.2 && eta<2.1) eff = 1.001;
		if(eta>=2.1 && eta<2.4) eff = 1.015;
	}
	if(pt>=90. && pt<140.){
		if(eta<0.9) eff = 1.000;
		if(eta>=0.9 && eta<1.2) eff = 1.001;
		if(eta>=1.2 && eta<2.1) eff = 0.999;
		if(eta>=2.1 && eta<2.4) eff = 1.008;
	}
	if(pt>=140. && pt<300.){
		if(eta<0.9) eff = 0.999;
		if(eta>=0.9 && eta<1.2) eff = 1.002;
		if(eta>=1.2 && eta<2.1) eff = 0.996;
		if(eta>=2.1 && eta<2.4) eff = 1.011;
	}
	return eff;
}

double ReferenceScaleFactors::MuonIsoUncTight2012(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	double err = 0.;
	if(pt>=10. && pt<20.){
		if(eta<0.9) err = 0.0043;
		if(eta>=0.9 && eta<1.2) err = 0.0048;
		if(eta>=1.2 && eta<2.1) err = 0.0023;
		if(eta>=2.1 && eta<2.4) err = 0.0055;
	}
	if(pt>=20. && pt<25.){
		if(eta<0.9) err = 0.0022;
		if(eta>=0.9 && eta<1.2) err = 0.0031;
		if(eta>=1.2 && eta<2.1) err = 0.0015;
		if(eta>=2.1 && eta<2.4) err = 0.0041;
	}
	if(pt>=25. && pt<30.){
		if(eta<0.9) err = 0.0012;
		if(eta>=0.9 && eta<1.2) err = 0.002;
		if(eta>=1.2 && eta<2.1) err = 0.001;
		if(eta>=2.1 && eta<2.4) err = 0.0028;
	}
	if(pt>=30. && pt<35.){
		if(eta<0.9) err = 0.0008;
		if(eta>=0.9 && eta<1.2) err = 0.0014;
		if(eta>=1.2 && eta<2.1) err = 0.0008;
		if(eta>=2.1 && eta<2.4) err = 0.0022;
	}
	if(pt>=35. && pt<40.){
		if(eta<0.9) err = 0.00054;
		if(eta>=0.9 && eta<1.2) err = 0.00092;
		if(eta>=1.2 && eta<2.1) err = 0.00056;
		if(eta>=2.1 && eta<2.4) err = 0.00017;
	}
	if(pt>=40. && pt<50.){
		if(eta<0.9) err = 0.00027;
		if(eta>=0.9 && eta<1.2) err = 0.00044;
		if(eta>=1.2 && eta<2.1) err = 0.00026;
		if(eta>=2.1 && eta<2.4) err = 0.00092;
	}
	if(pt>=50. && pt<60.){
		if(eta<0.9) err = 0.0005;
		if(eta>=0.9 && eta<1.2) err = 0.0008;
		if(eta>=1.2 && eta<2.1) err = 0.00048;
		if(eta>=2.1 && eta<2.4) err = 0.00017;
	}
	if(pt>=60. && pt<90.){
		if(eta<0.9) err = 0.00058;
		if(eta>=0.9 && eta<1.2) err = 0.001;
		if(eta>=1.2 && eta<2.1) err = 0.0006;
		if(eta>=2.1 && eta<2.4) err = 0.0021;
	}
	if(pt>=90. && pt<140.){
		if(eta<0.9) err = 0.0011;
		if(eta>=0.9 && eta<1.2) err = 0.0018;
		if(eta>=1.2 && eta<2.1) err = 0.0012;
		if(eta>=2.1 && eta<2.4) err = 0.0045;
	}
	if(pt>=140. && pt<300.){
		if(eta<0.9) err = 0.0026;
		if(eta>=0.9 && eta<1.2) err = 0.005;
		if(eta>=1.2 && eta<2.1) err = 0.0028;
		if(eta>=2.1 && eta<2.4) err = 0.018;
	}
	return err;
}

double ReferenceScaleFactors::TrackingEfficiency2012(TLorentzVector vect){
	double eta = vect.Eta();
	double eff = 1.;
	if(eta>=-2.4 && eta<-2.2) eff = 0.9861;
	if(eta>=-2.2 && eta<-2.0) eff = 0.9901;
	if(eta>=-2.0 && eta<-1.8) eff = 0.9946;
	if(eta>=-1.8 && eta<-1.6) eff = 0.9966;
	if(eta>=-1.6 && eta<-1.4) eff = 0.9964;
	if(eta>=-1.4 && eta<-1.2) eff = 0.9970;
	if(eta>=-1.2 && eta<-1.0) eff = 0.9973;
	if(eta>=-1.0 && eta<-0.8) eff = 0.9977;
	if(eta>=-0.8 && eta<-0.6) eff = 0.9980;
	if(eta>=-0.6 && eta<-0.4) eff = 0.9981;
	if(eta>=-0.4 && eta<-0.2) eff = 0.9974;
	if(eta>=-0.2 && eta<0.0) eff = 0.9967;
	if(eta>=0.0 && eta<0.2) eff = 0.9960;
	if(eta>=0.2 && eta<0.4) eff = 0.9978;
	if(eta>=0.4 && eta<0.6) eff = 0.9977;
	if(eta>=0.6 && eta<0.8) eff = 0.9975;
	if(eta>=0.8 && eta<1.0) eff = 0.9976;
	if(eta>=1.0 && eta<1.2) eff = 0.9964;
	if(eta>=1.2 && eta<1.4) eff = 0.9970;
	if(eta>=1.4 && eta<1.6) eff = 0.9948;
	if(eta>=1.6 && eta<1.8) eff = 0.9977;
	if(eta>=1.8 && eta<2.0) eff = 0.9974;
	if(eta>=2.0 && eta<2.2) eff = 0.9917;
	if(eta>=2.2 && eta<2.4) eff = 0.9806;
	return eff;
}

double ReferenceScaleFactors::TrackingEfficiencyUnc2012(TLorentzVector vect){
	double eta = vect.Eta();
	double eff = 0.;
	if(eta>=-2.4 && eta<-2.2) eff = 0.0012;
	if(eta>=-2.2 && eta<-2.0) eff = 0.0005;
	if(eta>=-2.0 && eta<-1.8) eff = 0.0003;
	if(eta>=-1.8 && eta<-1.6) eff = 0.0003;
	if(eta>=-1.6 && eta<-1.4) eff = 0.0003;
	if(eta>=-1.4 && eta<-1.2) eff = 0.0003;
	if(eta>=-1.2 && eta<-1.0) eff = 0.0002;
	if(eta>=-1.0 && eta<-0.8) eff = 0.0001;
	if(eta>=-0.8 && eta<-0.6) eff = 0.0001;
	if(eta>=-0.6 && eta<-0.4) eff = 0.0001;
	if(eta>=-0.4 && eta<-0.2) eff = 0.0001;
	if(eta>=-0.2 && eta<0.0) eff = 0.0001;
	if(eta>=0.0 && eta<0.2) eff = 0.0002;
	if(eta>=0.2 && eta<0.4) eff = 0.0001;
	if(eta>=0.4 && eta<0.6) eff = 0.0001;
	if(eta>=0.6 && eta<0.8) eff = 0.0001;
	if(eta>=0.8 && eta<1.0) eff = 0.0002;
	if(eta>=1.0 && eta<1.2) eff = 0.0002;
	if(eta>=1.2 && eta<1.4) eff = 0.0004;
	if(eta>=1.4 && eta<1.6) eff = 0.0003;
	if(eta>=1.6 && eta<1.8) eff = 0.0002;
	if(eta>=1.8 && eta<2.0) eff = 0.0002;
	if(eta>=2.0 && eta<2.2) eff = 0.0004;
	if(eta>=2.2 && eta<2.4) eff = 0.0013;
	return eff;
}

//source: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorkingSummer2013#Muon_ID_Isolation_EMu_Channel
double ReferenceScaleFactors::HiggsTauTau_EMu_Id_Mu(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	double eff = 1.;
	if(pt>10 && pt<=15){
		if(eta>=0 && eta<0.8){
			eff = 0.9771;
		}else if(eta>=0.8 && eta<1.2){
			eff = 0.9746;
		}else if(eta>=1.2 && eta<1.6){
			eff = 0.9644;
		}else if(eta>=1.6 && eta<2.1){
			eff = 0.9891;
		}
	}else if(pt>15 && pt<=20){
		if(eta>=0 && eta<0.8){
			eff = 0.9548;
		}else if(eta>=0.8 && eta<1.2){
			eff = 0.9701;
		}else if(eta>=1.2 && eta<1.6){
			eff = 0.9766;
		}else if(eta>=1.6 && eta<2.1){
			eff = 0.9892;
		}
	}else if(pt>20 && pt<=25){
		if(eta>=0 && eta<0.8){
			eff = 0.9648;
		}else if(eta>=0.8 && eta<1.2){
			eff = 0.9836;
		}else if(eta>=1.2 && eta<1.6){
			eff = 0.9820;
		}else if(eta>=1.6 && eta<2.1){
			eff = 0.9909;
		}
	}else if(pt>25 && pt<=30){
		if(eta>=0 && eta<0.8){
			eff = 0.9676;
		}else if(eta>=0.8 && eta<1.2){
			eff = 0.9817;
		}else if(eta>=1.2 && eta<1.6){
			eff = 0.9886;
		}else if(eta>=1.6 && eta<2.1){
			eff = 0.9883;
		}
	}else if(pt>30 && pt<=35){
		if(eta>=0 && eta<0.8){
			eff = 0.9883;
		}else if(eta>=0.8 && eta<1.2){
			eff = 0.9833;
		}else if(eta>=1.2 && eta<1.6){
			eff = 0.9910;
		}else if(eta>=1.6 && eta<2.1){
			eff = 0.9900;
		}
	}else if(pt>35){
		if(eta>=0 && eta<0.8){
			eff = 0.9826;
		}else if(eta>=0.8 && eta<1.2){
			eff = 0.9841;
		}else if(eta>=1.2 && eta<1.6){
			eff = 0.9900;
		}else if(eta>=1.6 && eta<2.1){
			eff = 0.9886;
		}
	}
	return eff;
}

double ReferenceScaleFactors::HiggsTauTau_EMu_IdUnc_Mu(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	double err = 0.;
	if(pt>10 && pt<=15){
		if(eta>=0 && eta<0.8){
			err = 0.0107;
		}else if(eta>=0.8 && eta<1.2){
			err = 0.0091;
		}else if(eta>=1.2 && eta<1.6){
			err = 0.0078;
		}else if(eta>=1.6 && eta<2.1){
			err = 0.0080;
		}
	}else if(pt>15 && pt<=20){
		if(eta>=0 && eta<0.8){
			err = 0.0046;
		}else if(eta>=0.8 && eta<1.2){
			err = 0.0049;
		}else if(eta>=1.2 && eta<1.6){
			err = 0.0049;
		}else if(eta>=1.6 && eta<2.1){
			err = 0.0049;
		}
	}else if(pt>20 && pt<=25){
		if(eta>=0 && eta<0.8){
			err = 0.0023;
		}else if(eta>=0.8 && eta<1.2){
			err = 0.0030;
		}else if(eta>=1.2 && eta<1.6){
			err = 0.0029;
		}else if(eta>=1.6 && eta<2.1){
			err = 0.0030;
		}
	}else if(pt>25 && pt<=30){
		if(eta>=0 && eta<0.8){
			err = 0.0012;
		}else if(eta>=0.8 && eta<1.2){
			err = 0.0018;
		}else if(eta>=1.2 && eta<1.6){
			err = 0.0018;
		}else if(eta>=1.6 && eta<2.1){
			err = 0.0019;
		}
	}else if(pt>30 && pt<=35){
		if(eta>=0 && eta<0.8){
			err = 0.0008;
		}else if(eta>=0.8 && eta<1.2){
			err = 0.0012;
		}else if(eta>=1.2 && eta<1.6){
			err = 0.0013;
		}else if(eta>=1.6 && eta<2.1){
			err = 0.0016;
		}
	}else if(pt>35){
		if(eta>=0 && eta<0.8){
			err = 0.0005;
		}else if(eta>=0.8 && eta<1.2){
			err = 0.0004;
		}else if(eta>=1.2 && eta<1.6){
			err = 0.0003;
		}else if(eta>=1.6 && eta<2.1){
			err = 0.0004;
		}
	}
	return err;
}

//source: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorkingSummer2013#Muon_ID_Isolation_MuTau_Channel
double ReferenceScaleFactors::HiggsTauTau_MuTau_Id_Mu(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	double weight = 1.;
	if(20 <= pt && pt < 30){
		if(eta>=0 && eta<0.8) weight = 0.9818;
		else if(eta>=0.8 && eta<1.2) weight = 0.9829;
		else if(eta>=1.2 && eta<2.1) weight = 0.9869;
	}
	else if(pt >= 30){
		if(eta>=0 && eta<0.8) weight = 0.9852;
		else if(eta>=0.8 && eta<1.2) weight = 0.9852;
		else if(eta>=1.2 && eta<2.1) weight = 0.9884;
	}
	return weight;
}

double ReferenceScaleFactors::HiggsTauTau_MuTau_IdUnc_Mu(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	double err = 0;
	if(20 <= pt && pt < 30){
		if(eta>=0 && eta<0.8) err = 0.0005;
		else if(eta>=0.8 && eta<1.2) err = 0.0009;
		else if(eta>=1.2 && eta<2.1) err = 0.0007;
	}
	else if(pt >= 30){
		if(eta>=0 && eta<0.8) err = 0.0001;
		else if(eta>=0.8 && eta<1.2) err = 0.0002;
		else if(eta>=1.2 && eta<2.1) err = 0.0001;
	}
	return err;
}
//source: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorkingSummer2013#Muon_ID_Isolation_MuTau_Channel
double ReferenceScaleFactors::HiggsTauTau_MuTau_Iso_Mu(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	double weight = 1.;
	if(20 <= pt && pt < 30){
		if(eta>=0 && eta<0.8) weight = 0.9494;
		else if(eta>=0.8 && eta<1.2) weight = 0.9835;
		else if(eta>=1.2 && eta<2.1) weight = 0.9923;
	}
	else if(pt >= 30){
		if(eta>=0 && eta<0.8) weight = 0.9883;
		else if(eta>=0.8 && eta<1.2) weight = 0.9937;
		else if(eta>=1.2 && eta<2.1) weight = 0.9996;
	}
	return weight;
}

double ReferenceScaleFactors::HiggsTauTau_MuTau_IsoUnc_Mu(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	double err = 0;
	if(20 <= pt && pt < 30){
		if(eta>=0 && eta<0.8) err = 0.0015;
		else if(eta>=0.8 && eta<1.2) err = 0.0020;
		else if(eta>=1.2 && eta<2.1) err = 0.0013;
	}
	else if(pt >= 30){
		if(eta>=0 && eta<0.8) err = 0.0003;
		else if(eta>=0.8 && eta<1.2) err = 0.0004;
		else if(eta>=1.2 && eta<2.1) err = 0.0005;
	}
	return err;
}

///////////////////////////
//
// Electron scale factors
//

// Electrons need eta from supercluster -> not using four vector as argument
double ReferenceScaleFactors::ElectronReconstruction2012(double Et, double Eta){
	if(!loadElectronID){Logger(Logger::Error) << "Electron ID not loaded in ReferenceScaleFactors." << std::endl; return -999;}
	double et = Et;
	double eff = 1.;
	if(fabs(Eta)<2.5){
		if(et>=200) et = 199.;
		eff = ElectronRecoEff->GetBinContent(ElectronRecoEff->FindFixBin(et,Eta));
	}
	return eff;
}

double ReferenceScaleFactors::ElectronReconstructionUnc2012(double Et, double Eta){
	if(!loadElectronID){Logger(Logger::Error) << "Electron ID not loaded in ReferenceScaleFactors." << std::endl; return -999;}
	double et = Et;
	double err = 0.;
	if(fabs(Eta)<2.5){
		if(et>=200) et = 199.;
		err = ElectronRecoEff->GetBinError(ElectronRecoEff->FindFixBin(et,Eta));
	}
	return err;
}

double ReferenceScaleFactors::ElectronIdTrig2012(double Et, double Eta){
	if(!loadElectronID){Logger(Logger::Error) << "Electron ID not loaded in ReferenceScaleFactors." << std::endl; return -999;}
	double et = Et;
	double eta = fabs(Eta);
	double eff = 1.;
	if(eta<2.5){
		if(et>=200) et = 199;
		eff = ElectronTrigEff->GetBinContent(ElectronTrigEff->FindFixBin(eta,et));
	}
	return eff;
}

double ReferenceScaleFactors::ElectronIdTrigUnc2012(double Et, double Eta){
	if(!loadElectronID){Logger(Logger::Error) << "Electron ID not loaded in ReferenceScaleFactors." << std::endl; return -999;}
	double et = Et;
	double eta = fabs(Eta);
	double err = 0.;
	if(eta<2.5){
		if(et>=200) et = 199;
		err = ElectronTrigEff->GetBinError(ElectronTrigEff->FindFixBin(eta,et));
	}
	return err;
}

double ReferenceScaleFactors::ElectronIdNonTrig2012(double Et, double Eta){
	if(!loadElectronID){Logger(Logger::Error) << "Electron ID not loaded in ReferenceScaleFactors." << std::endl; return -999;}
	double et = Et;
	double eff = 1.;
	if(fabs(Eta)<2.5){
		if(et>=200) et = 199;
		eff = ElectronNonTrigEff->GetBinContent(ElectronNonTrigEff->FindFixBin(et,Eta));
	}
	return eff;
}

double ReferenceScaleFactors::ElectronIdNonTrigUnc2012(double Et, double Eta){
	if(!loadElectronID){Logger(Logger::Error) << "Electron ID not loaded in ReferenceScaleFactors." << std::endl; return -999;}
	double et = Et;
	double err = 0.;
	if(fabs(Eta)<2.5){
		if(et>=200) et = 199;
		err = ElectronNonTrigEff->GetBinError(ElectronNonTrigEff->FindFixBin(et,Eta));
	}
	return err;
}

double ReferenceScaleFactors::ElectronEmbedding2012(double Et, double Eta){
	double et = Et;
	double eta = fabs(Eta);
	double eff = 1.;
	double xAxis[10] = {10,15,20,25,30,40,55,70,100,200};
	double yAxis[4] = {0,0.8,1.479,2.5};
	TH2D hPtEtaSFL("hPtEtaSFL","",9,xAxis,3,yAxis);
	hPtEtaSFL.SetBinContent(12,0.81);
	hPtEtaSFL.SetBinContent(13,0.91);
	hPtEtaSFL.SetBinContent(14,0.95);
	hPtEtaSFL.SetBinContent(15,0.96);
	hPtEtaSFL.SetBinContent(16,0.97);
	hPtEtaSFL.SetBinContent(17,0.98);
	hPtEtaSFL.SetBinContent(18,0.99);
	hPtEtaSFL.SetBinContent(19,0.98);
	hPtEtaSFL.SetBinContent(20,0.99);
	hPtEtaSFL.SetBinContent(21,0.98);
	hPtEtaSFL.SetBinContent(23,0.78);
	hPtEtaSFL.SetBinContent(24,0.89);
	hPtEtaSFL.SetBinContent(25,0.92);
	hPtEtaSFL.SetBinContent(26,0.94);
	hPtEtaSFL.SetBinContent(27,0.94);
	hPtEtaSFL.SetBinContent(28,0.97);
	hPtEtaSFL.SetBinContent(29,0.97);
	hPtEtaSFL.SetBinContent(30,0.99);
	hPtEtaSFL.SetBinContent(31,1.00);
	hPtEtaSFL.SetBinContent(32,1.00);
	hPtEtaSFL.SetBinContent(34,0.46);
	hPtEtaSFL.SetBinContent(35,0.66);
	hPtEtaSFL.SetBinContent(36,0.73);
	hPtEtaSFL.SetBinContent(37,0.80);
	hPtEtaSFL.SetBinContent(38,0.83);
	hPtEtaSFL.SetBinContent(39,0.86);
	hPtEtaSFL.SetBinContent(40,0.88);
	hPtEtaSFL.SetBinContent(41,0.91);
	hPtEtaSFL.SetBinContent(42,0.93);
	hPtEtaSFL.SetBinContent(43,1.00);

	if(et>199.99) et = 199.9;
	if(eta>2.49) eta = 2.49;
	if(et<10.) eff = 0.;
	eff = hPtEtaSFL.GetBinContent(hPtEtaSFL.FindFixBin(et,eta));

	return eff;
}

double ReferenceScaleFactors::HiggsTauTau_EMu_Id_E(double Et, double Eta){
	double et = Et;
	double eta = fabs(Eta);
	double eff = 1.;
	if(et>10 && et<=15){
		if(eta>=0 && eta<0.8){
			eff = 0.7654;
		}else if(eta>=0.8 && eta<1.5){
			eff = 0.7693;
		}else if(eta>=1.5 && eta<2.3){
			eff = 0.5719;
		}
	}else if(et>15 && et<=20){
		if(eta>=0 && eta<0.8){
			eff = 0.8394;
		}else if(eta>=0.8 && eta<1.5){
			eff = 0.8457;
		}else if(eta>=1.5 && eta<2.3){
			eff = 0.7024;
		}
	}else if(et>20 && et<=25){
		if(eta>=0 && eta<0.8){
			eff = 0.8772;
		}else if(eta>=0.8 && eta<1.5){
			eff = 0.8530;
		}else if(eta>=1.5 && eta<2.3){
			eff = 0.7631;
		}
	}else if(et>25 && et<=30){
		if(eta>=0 && eta<0.8){
			eff = 0.9006;
		}else if(eta>=0.8 && eta<1.5){
			eff = 0.8874;
		}else if(eta>=1.5 && eta<2.3){
			eff = 0.8092;
		}
	}else if(et>30 && et<=35){
		if(eta>=0 && eta<0.8){
			eff = 0.9261;
		}else if(eta>=0.8 && eta<1.5){
			eff = 0.9199;
		}else if(eta>=1.5 && eta<2.3){
			eff = 0.8469;
		}
	}else if(et>35){
		if(eta>=0 && eta<0.8){
			eff = 0.9514;
		}else if(eta>=0.8 && eta<1.5){
			eff = 0.9445;
		}else if(eta>=1.5 && eta<2.3){
			eff = 0.9078;
		}
	}
	return eff;
}

double ReferenceScaleFactors::HiggsTauTau_EMu_IdUnc_E(double Et, double Eta){
	double et = Et;
	double eta = fabs(Eta);
	double err = 0.;
	if(et>10 && et<=15){
		if(eta>=0 && eta<0.8){
			err = 0.0149;
		}else if(eta>=0.8 && eta<1.5){
			err = 0.0164;
		}else if(eta>=1.5 && eta<2.3){
			err = 0.0131;
		}
	}else if(et>15 && et<=20){
		if(eta>=0 && eta<0.8){
			err = 0.0045;
		}else if(eta>=0.8 && eta<1.5){
			err = 0.0061;
		}else if(eta>=1.5 && eta<2.3){
			err = 0.0075;
		}
	}else if(et>20 && et<=25){
		if(eta>=0 && eta<0.8){
			err = 0.0023;
		}else if(eta>=0.8 && eta<1.5){
			err = 0.0039;
		}else if(eta>=1.5 && eta<2.3){
			err = 0.0061;
		}
	}else if(et>25 && et<=30){
		if(eta>=0 && eta<0.8){
			err = 0.0018;
		}else if(eta>=0.8 && eta<1.5){
			err = 0.0017;
		}else if(eta>=1.5 && eta<2.3){
			err = 0.0024;
		}
	}else if(et>30 && et<=35){
		if(eta>=0 && eta<0.8){
			err = 0.0007;
		}else if(eta>=0.8 && eta<1.5){
			err = 0.0010;
		}else if(eta>=1.5 && eta<2.3){
			err = 0.0024;
		}
	}else if(et>35){
		if(eta>=0 && eta<0.8){
			err = 0.0002;
		}else if(eta>=0.8 && eta<1.5){
			err = 0.0003;
		}else if(eta>=1.5 && eta<2.3){
			err = 0.0007;
		}
	}
	return err;
}

///////////////////////////
//
// Tau scale factors
//

///////////////////////////
//
// Trigger scale factors
//

// single lepton trigger
double ReferenceScaleFactors::IsoMu24_eta2p1(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	double eff = 1.;
	if(pt>25. && pt<=30){
		if(eta<0.9){
			eff = 0.9837;
		}else if(eta>=0.9 && eta<1.2){
			eff = 0.9684;
		}else if(eta>=1.2 && eta<2.1){
			eff = 1.0052;
		}
	}else if(pt>30. && pt<=35.){
		if(eta<0.9){
			eff = 0.9841;
		}else if(eta>=0.9 && eta<1.2){
			eff = 0.9654;
		}else if(eta>=1.2 && eta<2.1){
			eff = 1.0014;
		}
	}else if(pt>35. && pt<=40.){
		if(eta<0.9){
			eff = 0.9839;
		}else if(eta>=0.9 && eta<1.2){
			eff = 0.9670;
		}else if(eta>=1.2 && eta<2.1){
			eff = 0.9962;
		}
	}else if(pt>40. && pt<=50.){
		if(eta<0.9){
			eff = 0.9835;
		}else if(eta>=0.9 && eta<1.2){
			eff = 0.9667;
		}else if(eta>=1.2 && eta<2.1){
			eff = 0.9943;
		}
	}else if(pt>50. && pt<=60){
		if(eta<0.9){
			eff = 0.9843;
		}else if(eta>=0.9 && eta<1.2){
			eff = 0.9627;
		}else if(eta>=1.2 && eta<2.1){
			eff = 0.9905;
		}
	}else if(pt>60. && pt<=90){
		if(eta<0.9){
			eff = 0.9847;
		}else if(eta>=0.9 && eta<1.2){
			eff = 0.9595;
		}else if(eta>=1.2 && eta<2.1){
			eff = 0.9883;
		}
	}else if(pt>90. && pt<=140.){
		if(eta<0.9){
			eff = 0.9809;
		}else if(eta>=0.9 && eta<1.2){
			eff = 0.9644;
		}else if(eta>=1.2 && eta<2.1){
			eff = 0.9819;
		}
	}else if(pt>140. && pt<=500.){
		if(eta<0.9){
			eff = 0.9804;
		}else if(eta>=0.9 && eta<1.2){
			eff = 0.9713;
		}else if(eta>=1.2 && eta<2.1){
			eff = 0.9942;
		}
	}
	return eff;
}

double ReferenceScaleFactors::IsoMu24_eta2p1_unc(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	double err = 0.;
	if(pt>25. && pt<=30){
		if(eta<0.9){
			err = 0.0010;
		}else if(eta>=0.9 && eta<1.2){
			err = 0.0026;
		}else if(eta>=1.2 && eta<2.1){
			err = 0.0018;
		}
	}else if(pt>30. && pt<=35.){
		if(eta<0.9){
			err = 0.0006;
		}else if(eta>=0.9 && eta<1.2){
			err = 0.0019;
		}else if(eta>=1.2 && eta<2.1){
			err = 0.0014;
		}
	}else if(pt>35. && pt<=40.){
		if(eta<0.9){
			err = 0.0004;
		}else if(eta>=0.9 && eta<1.2){
			err = 0.0014;
		}else if(eta>=1.2 && eta<2.1){
			err = 0.0011;
		}
	}else if(pt>40. && pt<=50.){
		if(eta<0.9){
			err = 0.0003;
		}else if(eta>=0.9 && eta<1.2){
			err = 0.0005;
		}else if(eta>=1.2 && eta<2.1){
			err = 0.0007;
		}
	}else if(pt>50. && pt<=60){
		if(eta<0.9){
			err = 0.0006;
		}else if(eta>=0.9 && eta<1.2){
			err = 0.0019;
		}else if(eta>=1.2 && eta<2.1){
			err = 0.0016;
		}
	}else if(pt>60. && pt<=90){
		if(eta<0.9){
			err = 0.0009;
		}else if(eta>=0.9 && eta<1.2){
			err = 0.0030;
		}else if(eta>=1.2 && eta<2.1){
			err = 0.0025;
		}
	}else if(pt>90. && pt<=140.){
		if(eta<0.9){
			err = 0.0031;
		}else if(eta>=0.9 && eta<1.2){
			err = 0.0101;
		}else if(eta>=1.2 && eta<2.1){
			err = 0.0081;
		}
	}else if(pt>140. && pt<=500.){
		if(eta<0.9){
			err = 0.0081;
		}else if(eta>=0.9 && eta<1.2){
			err = 0.0254;
		}else if(eta>=1.2 && eta<2.1){
			err = 0.0276;
		}
	}
	return err;
}

//trigger turn-on parameterization from https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2012#ETau_MuTau_trigger_turn_on_Joshu
//m is pt
double ReferenceScaleFactors::Efficiency(double m, double m0, double sigma, double alpha, double n, double norm){
	const double sqrtPiOver2 = 1.2533141373;
	const double sqrt2 = 1.4142135624;
	double sig = fabs((double) sigma);
	double t = (m - m0)/sig;
	if(alpha < 0) t = -t;
	double absAlpha = fabs(alpha/sig);
	double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
	double b = absAlpha - n/absAlpha;
	double ApproxErf;
	double arg = absAlpha / sqrt2;
	if (arg > 5.) ApproxErf = 1;
	else if (arg < -5.) ApproxErf = -1;
	else ApproxErf = TMath::Erf(arg);
	double leftArea = (1 + ApproxErf) * sqrtPiOver2;
	double rightArea = ( a * 1/TMath::Power(absAlpha - b,n-1)) / (n - 1);
	double area = leftArea + rightArea;
	if( t <= absAlpha ){
		arg = t / sqrt2;
		if(arg > 5.) ApproxErf = 1;
		else if (arg < -5.) ApproxErf = -1;
		else ApproxErf = TMath::Erf(arg);
		return norm * (1 + ApproxErf) * sqrtPiOver2 / area;
	}
	else{
		return norm * (leftArea + a * (1/TMath::Power(t-b,n-1) - 1/TMath::Power(absAlpha - b,n-1)) / (1 - n)) / area;
	}
}

///////////////////////////
//
// Final state mu+tau
//

//source: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorkingSummer2013#Electron_Muon_Tau_Trigger
double ReferenceScaleFactors::HiggsTauTau_MuTau_Trigger_Mu_Eff_Data(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = vect.Eta();
	double eff = 1.;
	if(pt >= 20){
		if(eta < -1.2){
			eff = Efficiency(pt, 15.9977, 7.64004e-05, 6.4951e-08, 1.57403, 0.865325);
		}
		else if (-1.2 <= eta && eta < -0.8){
			eff = Efficiency(pt, 17.3974, 0.804001, 1.47145, 1.24295, 0.928198);
		}
		else if (-0.8 <= eta && eta < 0.){
			eff = Efficiency(pt, 16.4307, 0.226312, 0.265553, 1.55756, 0.974462);
		}
		else if (0. <= eta && eta < 0.8){
			eff = Efficiency(pt, 17.313, 0.662731, 1.3412, 1.05778, 1.26624);
		}
		else if (0.8 <= eta && eta < 1.2){
			eff = Efficiency(pt, 16.9966, 0.550532, 0.807863, 1.55402, 0.885134);
		}
		else if (eta >= 1.2){
			eff = Efficiency(pt, 15.9962, 0.000106195, 4.95058e-08, 1.9991, 0.851294);
		}
	}
	return eff;
}

double ReferenceScaleFactors::HiggsTauTau_MuTau_Trigger_Mu_Eff_MC(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = vect.Eta();
	double eff = 1.;
	if(pt >= 20){
		if(eta < -1.2){
			eff = Efficiency(pt, 16.0051, 2.45144e-05, 4.3335e-09, 1.66134, 0.87045);
		}
		else if (-1.2 <= eta && eta < -0.8){
			eff = Efficiency(pt, 17.3135, 0.747636, 1.21803, 1.40611, 0.934983);
		}
		else if (-0.8 <= eta && eta < 0.){
			eff = Efficiency(pt, 15.9556, 0.0236127, 0.00589832, 1.75409, 0.981338);
		}
		else if (0. <= eta && eta < 0.8){
			eff = Efficiency(pt, 15.9289, 0.0271317, 0.00448573, 1.92101, 0.978625);
		}
		else if (0.8 <= eta && eta < 1.2){
			eff = Efficiency(pt, 16.5678, 0.328333, 0.354533, 1.67085, 0.916992);
		}
		else if (eta >= 1.2){
			eff = Efficiency(pt, 15.997, 7.90069e-05, 4.40036e-08, 1.66272, 0.884502);
		}
	}
	return eff;
}

double ReferenceScaleFactors::HiggsTauTau_MuTau_Trigger_Mu_ScaleMCtoData(TLorentzVector vect){
	double Data_eff =  HiggsTauTau_MuTau_Trigger_Mu_Eff_Data(vect);
	double MC_eff   =  HiggsTauTau_MuTau_Trigger_Mu_Eff_MC(vect);
	return Data_eff/MC_eff;
}

//NOTE: the set of trigger efficiencies below labelled http://benitezj.web.cern.ch/benitezj/Summer13Studies/TauTrigger/muTauABCD_June30/results.txt was the file we finally used for the legacy result: Without Tau ES corrections
double ReferenceScaleFactors::HiggsTauTau_MuTau_Trigger_Tau_Eff_Data(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = vect.Eta();
	double fabs_eta = fabs(eta);
	double eff = 1.;
	if(pt >= 20){
		if(fabs_eta < 1.5){
			eff = Efficiency(pt, 18.604910, 0.276042, 0.137039, 2.698437, 0.940721);
		}
		if(fabs_eta >= 1.5){
			eff = Efficiency(pt, 18.701715, 0.216523, 0.148111, 2.245081, 0.895320);
		}
	}
	return eff;
}

double ReferenceScaleFactors::HiggsTauTau_MuTau_Trigger_Tau_Eff_MC(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = vect.Eta();
	double fabs_eta = fabs(eta);
	double eff = 1.;
	if(pt >= 20){
		if(fabs_eta < 1.5){
			eff = Efficiency(pt, 18.532997, 1.027880, 2.262950, 1.003322, 5.297292);
		}
		if(fabs_eta >= 1.5){
			eff = Efficiency(pt, 18.212782, 0.338119, 0.122828, 12.577926, 0.893975);
		}
	}
	return eff;
}

double ReferenceScaleFactors::HiggsTauTau_MuTau_Trigger_Tau_ScaleMCtoData(TLorentzVector vect){
	double Data_eff =  HiggsTauTau_MuTau_Trigger_Tau_Eff_Data(vect);
	double MC_eff   =  HiggsTauTau_MuTau_Trigger_Tau_Eff_MC(vect);
	return Data_eff/MC_eff;
}

///////////////////////////
//
// Final state e+mu
//

double ReferenceScaleFactors::HiggsWW_EMu_Trigger(TLorentzVector mu_vect, double e_et, double e_eta, TString path){
	double eff = 1.;
	if(path.Contains("Mu17_Ele8")){
		eff = HiggsWW_EMu_SingleMu(mu_vect) + (1-HiggsWW_EMu_SingleMu(mu_vect))*HiggsWW_EMu_SingleEle(e_et, e_eta)
				+ (1-HiggsWW_EMu_SingleMu(mu_vect))*(1-HiggsWW_EMu_SingleEle(e_et, e_eta))*
				(HiggsWW_EMu_DoubleMuLeading(mu_vect)*HiggsWW_EMu_DoubleEleTrailing(e_et, e_eta) + (1-HiggsWW_EMu_DoubleMuLeading(mu_vect)*HiggsWW_EMu_DoubleEleTrailing(e_et, e_eta))*HiggsWW_EMu_DoubleEleLeading(e_et, e_eta)*HiggsWW_EMu_DoubleMuTrailing(mu_vect));
	}
	if(path.Contains("Mu8_Ele17")){
		eff = HiggsWW_EMu_SingleEle(e_et, e_eta) + (1-HiggsWW_EMu_SingleEle(e_et, e_eta))*HiggsWW_EMu_SingleMu(mu_vect)
				+ (1-HiggsWW_EMu_SingleEle(e_et, e_eta))*(1-HiggsWW_EMu_SingleMu(mu_vect))*
				(HiggsWW_EMu_DoubleEleLeading(e_et, e_eta)*HiggsWW_EMu_DoubleMuTrailing(mu_vect) + (1-HiggsWW_EMu_DoubleEleLeading(e_et, e_eta)*HiggsWW_EMu_DoubleMuTrailing(mu_vect))*HiggsWW_EMu_DoubleMuLeading(mu_vect)*HiggsWW_EMu_DoubleEleTrailing(e_et, e_eta));
	}
	return eff;
}

double ReferenceScaleFactors::HiggsWW_EMu_SingleEle(double Et, double Eta){
	if(!loadEMuTriggerEff){Logger(Logger::Error) << "EMu trigger efficiency not loaded in ReferenceScaleFactors." << std::endl; return -999;}
	double et = Et;
	double eta = fabs(Eta);
	double eff = 1.;
	if(et>=60.)et=55.;
	if(eta<1.5) eff = HiggsWW_EMu_SingleEle15->GetBinContent(HiggsWW_EMu_SingleEle15->FindFixBin(et));
	if(eta>=1.5 && eta<2.5) eff = HiggsWW_EMu_SingleEle25->GetBinContent(HiggsWW_EMu_SingleEle25->FindFixBin(et));
	return eff;
}

double ReferenceScaleFactors::HiggsWW_EMu_DoubleEleLeading(double Et, double Eta){
	if(!loadEMuTriggerEff){Logger(Logger::Error) << "EMu trigger efficiency not loaded in ReferenceScaleFactors." << std::endl; return -999;}
	double et = Et;
	double eta = fabs(Eta);
	double eff = 1.;
	if(et>=60.)et=55.;
	if(eta<1.5) eff = HiggsWW_EMu_DoubleEleLead15->GetBinContent(HiggsWW_EMu_DoubleEleLead15->FindFixBin(et));
	if(eta>=1.5 && eta<2.5) eff = HiggsWW_EMu_DoubleEleLead25->GetBinContent(HiggsWW_EMu_DoubleEleLead25->FindFixBin(et));
	return eff;
}

double ReferenceScaleFactors::HiggsWW_EMu_DoubleEleTrailing(double Et, double Eta){
	if(!loadEMuTriggerEff){Logger(Logger::Error) << "EMu trigger efficiency not loaded in ReferenceScaleFactors." << std::endl; return -999;}
	double et = Et;
	double eta = fabs(Eta);
	double eff = 1.;
	if(et>=60.)et=55.;
	if(eta<1.5) eff = HiggsWW_EMu_DoubleEleTrail15->GetBinContent(HiggsWW_EMu_DoubleEleTrail15->FindFixBin(et));
	if(eta>=1.5 && eta<2.5) eff = HiggsWW_EMu_DoubleEleTrail25->GetBinContent(HiggsWW_EMu_DoubleEleTrail25->FindFixBin(et));
	return eff;
}

double ReferenceScaleFactors::HiggsWW_EMu_SingleMu(TLorentzVector vect){
	if(!loadEMuTriggerEff){Logger(Logger::Error) << "EMu trigger efficiency not loaded in ReferenceScaleFactors." << std::endl; return -999;}
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	double eff = 1.;
	if(pt>=60.)pt=55.;
	if(eta<0.8) eff = HiggsWW_EMu_SingleMu08->GetBinContent(HiggsWW_EMu_SingleMu08->FindFixBin(pt));
	if(eta>=0.8 && eta<1.2) eff = HiggsWW_EMu_SingleMu12->GetBinContent(HiggsWW_EMu_SingleMu12->FindFixBin(pt));
	if(eta>=1.2 && eta<2.1) eff = HiggsWW_EMu_SingleMu21->GetBinContent(HiggsWW_EMu_SingleMu21->FindFixBin(pt));
	if(eta>=2.1 && eta<2.5) eff = HiggsWW_EMu_SingleMu25->GetBinContent(HiggsWW_EMu_SingleMu25->FindFixBin(pt));
	return eff;
}

double ReferenceScaleFactors::HiggsWW_EMu_DoubleMuLeading(TLorentzVector vect){
	if(!loadEMuTriggerEff){Logger(Logger::Error) << "EMu trigger efficiency not loaded in ReferenceScaleFactors." << std::endl; return -999;}
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	double eff = 1.;
	if(pt>=60.)pt=55.;
	if(eta<1.2) eff = HiggsWW_EMu_DoubleMuLead12->GetBinContent(HiggsWW_EMu_DoubleMuLead12->FindFixBin(pt));
	if(eta>=1.2 && eta<2.1) eff = HiggsWW_EMu_DoubleMuLead21->GetBinContent(HiggsWW_EMu_DoubleMuLead21->FindFixBin(pt));
	if(eta>=2.1 && eta<2.5) eff = HiggsWW_EMu_DoubleMuLead25->GetBinContent(HiggsWW_EMu_DoubleMuLead25->FindFixBin(pt));
	return eff;
}

double ReferenceScaleFactors::HiggsWW_EMu_DoubleMuTrailing(TLorentzVector vect){
	if(!loadEMuTriggerEff){Logger(Logger::Error) << "EMu trigger efficiency not loaded in ReferenceScaleFactors." << std::endl; return -999;}
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	double eff = 1.;
	if(pt>=60.)pt=55.;
	if(eta<1.2) eff = HiggsWW_EMu_DoubleMuTrail12->GetBinContent(HiggsWW_EMu_DoubleMuTrail12->FindFixBin(pt));
	if(eta>=1.2 && eta<2.1) eff = HiggsWW_EMu_DoubleMuTrail21->GetBinContent(HiggsWW_EMu_DoubleMuTrail21->FindFixBin(pt));
	if(eta>=2.1 && eta<2.5) eff = HiggsWW_EMu_DoubleMuTrail25->GetBinContent(HiggsWW_EMu_DoubleMuTrail25->FindFixBin(pt));
	return eff;
}

double ReferenceScaleFactors::HiggsTauTau_EMu_Trigger_Mu(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	if(pt>10 && pt<=15){
		if(eta>=0 && eta<0.8){
			return 0.9829;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9745;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9943;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9158;
		}
	}else if(pt>15 && pt<=20){
		if(eta>=0 && eta<0.8){
			return 0.9850;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9852;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9743;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9333;
		}
	}else if(pt>20 && pt<=25){
		if(eta>=0 && eta<0.8){
			return 0.9951;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9610;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9716;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9459;
		}
	}else if(pt>25 && pt<=30){
		if(eta>=0 && eta<0.8){
			return 0.9869;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9779;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9665;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9501;
		}
	}else if(pt>30 && pt<=35){
		if(eta>=0 && eta<0.8){
			return 0.9959;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9881;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9932;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9391;
		}
	}else if(pt>35){
		if(eta>=0 && eta<0.8){
			return 0.9986;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9540;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9549;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9386;
		}
	}
	return 1.;
}

double ReferenceScaleFactors::HiggsTauTau_EMu_TriggerUnc_Mu(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	if(pt>10 && pt<=15){
		if(eta>=0 && eta<0.8){
			return 0.0058;
		}else if(eta>=0.8 && eta<1.2){
			return 0.0124;
		}else if(eta>=1.2 && eta<1.6){
			return 0.0164;
		}else if(eta>=1.6 && eta<2.1){
			return 0.0176;
		}
	}else if(pt>15 && pt<=20){
		if(eta>=0 && eta<0.8){
			return 0.0056;
		}else if(eta>=0.8 && eta<1.2){
			return 0.0171;
		}else if(eta>=1.2 && eta<1.6){
			return 0.0179;
		}else if(eta>=1.6 && eta<2.1){
			return 0.0162;
		}
	}else if(pt>20 && pt<=25){
		if(eta>=0 && eta<0.8){
			return 0.0060;
		}else if(eta>=0.8 && eta<1.2){
			return 0.0116;
		}else if(eta>=1.2 && eta<1.6){
			return 0.0141;
		}else if(eta>=1.6 && eta<2.1){
			return 0.0159;
		}
	}else if(pt>25 && pt<=30){
		if(eta>=0 && eta<0.8){
			return 0.0074;
		}else if(eta>=0.8 && eta<1.2){
			return 0.0187;
		}else if(eta>=1.2 && eta<1.6){
			return 0.0184;
		}else if(eta>=1.6 && eta<2.1){
			return 0.0251;
		}
	}else if(pt>30 && pt<=35){
		if(eta>=0 && eta<0.8){
			return 0.0085;
		}else if(eta>=0.8 && eta<1.2){
			return 0.0227;
		}else if(eta>=1.2 && eta<1.6){
			return 0.0271;
		}else if(eta>=1.6 && eta<2.1){
			return 0.0307;
		}
	}else if(pt>35){
		if(eta>=0 && eta<0.8){
			return 0.0087;
		}else if(eta>=0.8 && eta<1.2){
			return 0.0165;
		}else if(eta>=1.2 && eta<1.6){
			return 0.0211;
		}else if(eta>=1.6 && eta<2.1){
			return 0.0209;
		}
	}
	return 0.;
}

double ReferenceScaleFactors::HiggsTauTau_EMu_Trigger_E(double Et, double Eta){
	double et = Et;
	double eta = fabs(Eta);
	if(et>10 && et<=15){
		if(eta>=0 && eta<0.8){
			return 0.9548;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9015;
		}else if(eta>=1.5 && eta<2.3){
			return 0.9017;
		}
	}else if(et>15 && et<=20){
		if(eta>=0 && eta<0.8){
			return 0.9830;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9672;
		}else if(eta>=1.5 && eta<2.3){
			return 0.9463;
		}
	}else if(et>20 && et<=25){
		if(eta>=0 && eta<0.8){
			return 0.9707;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9731;
		}else if(eta>=1.5 && eta<2.3){
			return 0.9691;
		}
	}else if(et>25 && et<=30){
		if(eta>=0 && eta<0.8){
			return 0.9768;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9870;
		}else if(eta>=1.5 && eta<2.3){
			return 0.9727;
		}
	}else if(et>30 && et<=35){
		if(eta>=0 && eta<0.8){
			return 1.0047;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9891;
		}else if(eta>=1.5 && eta<2.3){
			return 0.9858;
		}
	}else if(et>35){
		if(eta>=0 && eta<0.8){
			return 1.0063;
		}else if(eta>=0.8 && eta<1.5){
			return 1.0047;
		}else if(eta>=1.5 && eta<2.3){
			return 1.0015;
		}
	}
	return 1.;
}

double ReferenceScaleFactors::HiggsTauTau_EMu_TriggerUnc_E(double Et, double Eta){
	double et = Et;
	double eta = fabs(Eta);
	if(et>10 && et<=15){
		if(eta>=0 && eta<0.8){
			return 0.0197;
		}else if(eta>=0.8 && eta<1.5){
			return 0.0205;
		}else if(eta>=1.5 && eta<2.3){
			return 0.0470;
		}
	}else if(et>15 && et<=20){
		if(eta>=0 && eta<0.8){
			return 0.0115;
		}else if(eta>=0.8 && eta<1.5){
			return 0.0113;
		}else if(eta>=1.5 && eta<2.3){
			return 0.0212;
		}
	}else if(et>20 && et<=25){
		if(eta>=0 && eta<0.8){
			return 0.0087;
		}else if(eta>=0.8 && eta<1.5){
			return 0.0083;
		}else if(eta>=1.5 && eta<2.3){
			return 0.0149;
		}
	}else if(et>25 && et<=30){
		if(eta>=0 && eta<0.8){
			return 0.0084;
		}else if(eta>=0.8 && eta<1.5){
			return 0.0083;
		}else if(eta>=1.5 && eta<2.3){
			return 0.0162;
		}
	}else if(et>30 && et<=35){
		if(eta>=0 && eta<0.8){
			return 0.0100;
		}else if(eta>=0.8 && eta<1.5){
			return 0.0111;
		}else if(eta>=1.5 && eta<2.3){
			return 0.0112;
		}
	}else if(et>35){
		if(eta>=0 && eta<0.8){
			return 0.0078;
		}else if(eta>=0.8 && eta<1.5){
			return 0.0073;
		}else if(eta>=1.5 && eta<2.3){
			return 0.0135;
		}
	}
	return 0.;
}

// Higgs pT reweighting
double ReferenceScaleFactors::HiggsPtWeight(TLorentzVector vect, int mass, TString shift){
	// define which scale to use. Default is "nominal".
	std::map<unsigned, TH1D*> weightMap;
	if (shift == "nominal")		weightMap = HiggsPtWeightsNominal;
	else if (shift == "down")	weightMap = HiggsPtWeightsDown;
	else if (shift == "up")		weightMap = HiggsPtWeightsUp;
	else {
		Logger(Logger::Error) << "ERROR: shift of type " << shift.Data() << " not known for Higgs pT weights." << std::endl;
		return -999;
	}

	// do some checks
	if( weightMap.find(mass) == weightMap.end() ){
		Logger(Logger::Error) << "Requested mass " << mass << " not available for Higgs pT reweighting." << std::endl;
		return -999;
	}

	// read weight from histogram
	return weightMap[mass]->GetBinContent( weightMap[mass]->FindFixBin( vect.Pt() ) );
}

double ReferenceScaleFactors::TauTriggerTurnOnPars(TLorentzVector vect, TString type, int dm ){
  std::vector<double> out; // "m_{0}" "sigma"  "alpha"   "n"  "norm"


  int dmindex;
  if(dm==0) dmindex=0;
  if(dm==1) dmindex=1;
  if(dm==10)dmindex=2;



 
  float  data_m0[3]       = {38.73447611695629, 36.690947837524625,   42.83389974719572};
  float  data_sigma[3]    = {7.407431850780567,  4.615715187723146,   6.618279442098442};
  float  data_alpha[3]    = {13.773094478855745,6.837452816365016,    5.586768802617782};
  float  data_n[3]        = {1.2242791989806183,  1.3923658691165495, 138.11526635111798};
  float  data_norm[3]    = {0.999999999045921, 0.9999998778570275,    0.9846027163460338} ;



  float mc_m0[3]     = { 38.40062512717636, 36.19161339284503,40.034392442679504};
  float mc_sigma[3]  = {  6.880123201583417, 4.354633288490193,5.013762507971054};
  float mc_alpha[3]  = { 11.139690096268433, 5.5715514431807485,2.0072577361956396};
  float mc_n[3]      = {  1.5069036953587036,  1.6738765374606683,123.76902846216875};
  float mc_norm[3]   = {  0.999999986561927,  0.9999999897484987,0.9912887339725038};

  double pt = vect.Pt();

  double eff(1);
  if(type.Contains("data"))     eff = Efficiency(pt, data_m0[dmindex], data_sigma[dmindex], data_alpha[dmindex], data_n[dmindex], data_norm[dmindex]);	
  if(type.Contains("mc"))     eff = Efficiency(pt, mc_m0[dmindex], mc_sigma[dmindex], mc_alpha[dmindex], mc_n[dmindex], mc_norm[dmindex]);	
  return eff;
}


 
double ReferenceScaleFactors::DiTauTrigger2016_ScaleMCtoData(TLorentzVector vect, int dm){
  double Data_eff =  TauTriggerTurnOnPars(vect,"data",dm);
  double MC_eff   =  TauTriggerTurnOnPars(vect,"mc",dm);
  return Data_eff/MC_eff;
}
