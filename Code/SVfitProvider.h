/*
 * SVfitProvider.h
 *
 *  Created on: Jan 19, 2015
 *      Author: cherepanov
 *
 *  Interface to SVfit mass reconstruction framework
 *	Holds input information for SVfit, executes SVfit and provides results
 */

#ifndef SVFITPROVIDER_H_
#define SVFITPROVIDER_H_

#include "Objects.h"
#include "TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h"
#include "DataFormats/SVFitObject.h"
#include "TLorentzVector.h"

class Ntuple_Controller;

class SVfitProvider {
public:
	// Constructor to be used in analysis:
	// Use default CMS leptons as input (electron, muon, hadronic tau)
	SVfitProvider(Ntuple_Controller* const Ntp, objects::MET& met, TString typeLep1, int idxLep1, TString typeLep2, int idxLep2,
			int verbosity = 1, double scaleLep1 = 1 , double scaleLep2 = 1 );

	// Constructor to be used in analysis:
	// Use default one default CMS lepton (electron, muon, hadronic tau) and a fully reconstructed 3prong tau
	SVfitProvider(Ntuple_Controller* const Ntp, objects::MET& met, TString typeLep1, int idxLep1, TLorentzVector lvec3ProngTau,
			int verbosity/* =1 */, double scaleLep1 /* =1 */, double scaleLep2 /* =1 */);

	virtual ~SVfitProvider();

	// run SVfit algorithm and create  SVfitObject
	SVFitObject runAndMakeObject();

	// execute SVfit
	void run();

	// create SVfitObject
	SVFitObject makeObject(const SVfitStandaloneAlgorithm* svfitAlgo, TString fitMethod);

	// conversions
	static svFitStandalone::LorentzVector convert_p4Vect(const TLorentzVector& in);

	// getters and setters
	const Ntuple_Controller* get_ntp() const {return ntp_;}
	void set_ntp(Ntuple_Controller* ntp) {this->ntp_ = ntp; isSetup_ = false;}

	const objects::MET& get_inputMet() const {return inputMet_;}
	void set_inputMet(const objects::MET& inputMet) {inputMet_ = inputMet; isSetup_ = false;}

	const std::vector<svFitStandalone::MeasuredTauLepton>& get_inputTauLeptons() const {return inputTauLeptons_;}
	void set_inputTauLeptons(const std::vector<svFitStandalone::MeasuredTauLepton>& inputTauLeptons) {inputTauLeptons_ = inputTauLeptons; isSetup_ = false;}

	int get_verbosity() const {return verbosity_;}
	void set_verbosity(int verbosity) {this->verbosity_ = verbosity;}

	TString get_fitMethod() const {return fitMethod_;}
	void set_fitMethod(TString method) {this->fitMethod_ = method;}

	bool get_addLogM() const {return addLogM_;}
	void set_addLogM(bool addLogM) {addLogM_ = addLogM; isSetup_ = false;}

	unsigned get_maxObjFunctionCalls() const {return maxObjFunctionCalls_;}
	void set_maxObjFunctionCalls(unsigned maxObjFunctionCalls) {maxObjFunctionCalls_ = maxObjFunctionCalls; isSetup_ = false;}

	float get_metPower() const {return metPower_;}
	void set_metPower(float metPower) {metPower_ = metPower; isSetup_ = false;}

	// access to SVfit object which holds the result
	const SVfitStandaloneAlgorithm* result() const {return svFitAlgo_;}

private:
	// input information
	Ntuple_Controller* ntp_;
	objects::MET inputMet_;
	std::vector<svFitStandalone::MeasuredTauLepton> inputTauLeptons_;
	int verbosity_;
	TString fitMethod_;
	// SVfit configuration
	bool addLogM_;
	float metPower_;
	int maxObjFunctionCalls_;

	bool isSetup_;

	// output information
	SVfitStandaloneAlgorithm* svFitAlgo_;

	void addMeasuredLepton(TString type, int index, double energyScale = 1);
	void addFullReco3ProngTau(TLorentzVector lv);
	void createSvFitAlgo();
};

#endif /* SVFITPROVIDER_H_ */
