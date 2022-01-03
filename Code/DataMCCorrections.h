/*
 * DataMCCorrections.cxx
 *
 *  Created on: Nov  16, 2017
 *      Author: cherepanov
 */


#ifndef DataMCCorrections_H_
#define DataMCCorrections_H_

#include <map>
#include "TFile.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "ScaleFactor.h"
class DataMCCorrections {

public:

  DataMCCorrections( bool load_ZPtWeights = true, bool load_LeptonEff=true, bool load_LeptonIso=false);
	virtual ~DataMCCorrections();


	// Higgs/Z pT reweighting
	double HiggsPtWeight(TLorentzVector vect, int mass, TString shift = "nominal");
	float ZPTWeight(float genMass, float genPt);
	float ZPTWeightErr(float genMass, float genPt);

	float AgainstElectronDataMCCorrection(TLorentzVector p4, TString type);
	float AgainstMuonDataMCCorrection(TLorentzVector p4, TString type);

	float LeptonTriggerEfficiencyScaleFactor(TLorentzVector leptonP4);
	double get_EfficiencyData(double, double); //pt, eta
	double get_EfficiencyMC(double, double);
	double get_ScaleFactor(double, double); 
	double get_EfficiencyDataError(double, double);
	double get_EfficiencyMCError(double, double);
	double get_ScaleFactorError(double, double);
	
private:
	// booleans to switch on/off individual scale factors
	bool loadZPtWeights;
	bool loadLeptonEff;
	bool loadLeptonIso;

	//
	// Root files for scale factors
	//

	// Z pT reweighting
	TFile* ZPtWeightFile;

	//
	// Histograms for scale factors
	//

	// Z pT reweighting
	TH2D* m_zPtHist;
	TH2D* m_zPtHistErr;

	//Lepton Corrections
	ScaleFactor *LeptonTriggerEff;
	TH1D* etaBinsH;
	std::map<std::string, TGraphAsymmErrors *> eff_data;
	std::map<std::string, TGraphAsymmErrors *> eff_mc;


	void  SetAxisBins(TGraphAsymmErrors*);
	bool  check_SameBinning(TGraphAsymmErrors*, TGraphAsymmErrors*);
	std::string FindEtaLabel(double);
        int FindPtBin( std::map<std::string, TGraphAsymmErrors *>, std::string, double);


};



#endif /* REFERENCESCALEFACTORS_H_ */
