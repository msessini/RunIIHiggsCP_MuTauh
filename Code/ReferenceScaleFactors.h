/*
 * ReferenceScaleFactors.h
 *
 *  Created on: Aug 13, 2014
 *      Author: nehrkorn
 */

#ifndef REFERENCESCALEFACTORS_H_
#define REFERENCESCALEFACTORS_H_

#include <map>
#include "TFile.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TLorentzVector.h"

class ReferenceScaleFactors {

public:

	ReferenceScaleFactors(int runType, bool load_ElectronID = true, bool load_EMuTriggerEff = true, bool load_HiggsPtWeights = true);
	virtual ~ReferenceScaleFactors();

	// Muon Id scale factors
	double MuonIdTight2012(TLorentzVector vect);
	double MuonIdUncTight2012(TLorentzVector vect);
	double MuonIsoTight2012(TLorentzVector vect);
	double MuonIsoUncTight2012(TLorentzVector vect);
	double TrackingEfficiency2012(TLorentzVector vect);
	double TrackingEfficiencyUnc2012(TLorentzVector vect);
	double HiggsTauTau_EMu_Id_Mu(TLorentzVector vect);
	double HiggsTauTau_EMu_IdUnc_Mu(TLorentzVector vect);
	double HiggsTauTau_MuTau_Id_Mu(TLorentzVector vect);
	double HiggsTauTau_MuTau_IdUnc_Mu(TLorentzVector vect);
	double HiggsTauTau_MuTau_Iso_Mu(TLorentzVector vect);
	double HiggsTauTau_MuTau_IsoUnc_Mu(TLorentzVector vect);

	// Electron Id scale factors
	double ElectronReconstruction2012(double Et, double Eta);
	double ElectronReconstructionUnc2012(double Et, double Eta);
	double ElectronIdTrig2012(double Et, double Eta);
	double ElectronIdTrigUnc2012(double Et, double Eta);
	double ElectronIdNonTrig2012(double Et, double Eta);
	double ElectronIdNonTrigUnc2012(double Et, double Eta);
	double ElectronEmbedding2012(double Et, double Eta);
	double HiggsTauTau_EMu_Id_E(double Et, double Eta);
	double HiggsTauTau_EMu_IdUnc_E(double Et, double Eta);

	// Tau scale factors

	//
	// Trigger efficiency scale factors
	//

	// ++++++++++ Single lepton ++++++++++
	double IsoMu24_eta2p1(TLorentzVector vect);
	double IsoMu24_eta2p1_unc(TLorentzVector vect);
	// ++++++++++ Final state: e+mu ++++++++++
	// use this for full trigger efficiency
	double HiggsWW_EMu_Trigger(TLorentzVector mu_vect, double e_et, double e_eta, TString path);
	// efficiencies of individual trigger legs
	double HiggsWW_EMu_SingleEle(double Et, double Eta);
	double HiggsWW_EMu_DoubleEleLeading(double Et, double Eta);
	double HiggsWW_EMu_DoubleEleTrailing(double Et, double Eta);
	double HiggsWW_EMu_SingleMu(TLorentzVector vect);
	double HiggsWW_EMu_DoubleMuLeading(TLorentzVector vect);
	double HiggsWW_EMu_DoubleMuTrailing(TLorentzVector vect);
	double HiggsTauTau_EMu_Trigger_Mu(TLorentzVector vect);
	double HiggsTauTau_EMu_TriggerUnc_Mu(TLorentzVector vect);
	double HiggsTauTau_EMu_Trigger_E(double Et, double Eta);
	double HiggsTauTau_EMu_TriggerUnc_E(double Et, double Eta);

	// Final state: mu+tau
	double Efficiency(double m, double m0, double sigma, double alpha, double n, double norm);
	double HiggsTauTau_MuTau_Trigger_Mu_Eff_Data(TLorentzVector vect);
	double HiggsTauTau_MuTau_Trigger_Mu_Eff_MC(TLorentzVector vect);
	double HiggsTauTau_MuTau_Trigger_Mu_ScaleMCtoData(TLorentzVector vect);
	double HiggsTauTau_MuTau_Trigger_Tau_Eff_Data(TLorentzVector vect);
	double HiggsTauTau_MuTau_Trigger_Tau_Eff_MC(TLorentzVector vect);
	double HiggsTauTau_MuTau_Trigger_Tau_ScaleMCtoData(TLorentzVector vect);

	// Higgs pT reweighting
	double HiggsPtWeight(TLorentzVector vect, int mass, TString shift = "nominal");




	//Final state tau+tau
	double TauTriggerTurnOnPars(TLorentzVector vect, TString type, int dm );
 	double DiTauTrigger2016_ScaleMCtoData(TLorentzVector vect, int dm);
 

private:
	// booleans to switch on/off individual scale factors
	bool loadElectronID;
	bool loadEMuTriggerEff;
	bool loadHiggsPtWeights;

	//
	// Root files for scale factors
	//

	// Electron Id's
	TFile* ETrigIdEffFile;
	TFile* ENonTrigIdEffFile;
	TFile* ERecoEffFile;
	// Trigger efficiencies
	TFile* HWW_TriggerEfficiencies;
	// Higgs pT reweighting
	std::map<unsigned, TFile*> HiggsPtWeightFiles;

	//
	// Histograms for scale factors
	//

	// Electron Id's
	TH2D* ElectronTrigEff;
	TH2D* ElectronNonTrigEff;
	TH2D* ElectronRecoEff;
	// Trigger efficiencies
	TH1D* HiggsWW_EMu_SingleEle15;
	TH1D* HiggsWW_EMu_SingleEle25;
	TH1D* HiggsWW_EMu_DoubleEleLead15;
	TH1D* HiggsWW_EMu_DoubleEleLead25;
	TH1D* HiggsWW_EMu_DoubleEleTrail15;
	TH1D* HiggsWW_EMu_DoubleEleTrail25;
	TH1D* HiggsWW_EMu_SingleMu08;
	TH1D* HiggsWW_EMu_SingleMu12;
	TH1D* HiggsWW_EMu_SingleMu21;
	TH1D* HiggsWW_EMu_SingleMu25;
	TH1D* HiggsWW_EMu_DoubleMuLead12;
	TH1D* HiggsWW_EMu_DoubleMuLead21;
	TH1D* HiggsWW_EMu_DoubleMuLead25;
	TH1D* HiggsWW_EMu_DoubleMuTrail12;
	TH1D* HiggsWW_EMu_DoubleMuTrail21;
	TH1D* HiggsWW_EMu_DoubleMuTrail25;
	// Higgs pT reweighting
	std::map<unsigned, TH1D*> HiggsPtWeightsNominal;
	std::map<unsigned, TH1D*> HiggsPtWeightsDown;
	std::map<unsigned, TH1D*> HiggsPtWeightsUp;
};



#endif /* REFERENCESCALEFACTORS_H_ */
