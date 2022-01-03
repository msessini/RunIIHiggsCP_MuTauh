/*
 * SVfitProvider.cxx
 *
 *  Created on: Jan 19, 2015
 *      Author: cherepanov
 */

#include "Ntuple_Controller.h"
#include "SVfitProvider.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

void SVfitProvider::createSvFitAlgo() {
	svFitAlgo_ = new SVfitStandaloneAlgorithm(inputTauLeptons_, inputMet_.ex(), inputMet_.ey(), inputMet_.significanceMatrix(), verbosity_);
	svFitAlgo_->addLogM(addLogM_);
	if (maxObjFunctionCalls_ >= 0) svFitAlgo_->maxObjFunctionCalls(maxObjFunctionCalls_);
	if (metPower_ >= 0)	svFitAlgo_->metPower(metPower_);
	isSetup_ = true;
}

// Constructor to be used in analysis:
// Use default CMS leptons as input (electron, muon, hadronic tau)
SVfitProvider::SVfitProvider(Ntuple_Controller* const Ntp, objects::MET& met, TString typeLep1, int idxLep1, TString typeLep2, int idxLep2,
		int verbosity/* =1 */, double scaleLep1 /* =1 */, double scaleLep2 /* =1 */){
	ntp_ = Ntp;
	inputMet_ = met;

	addMeasuredLepton(typeLep1, idxLep1, scaleLep1);
	addMeasuredLepton(typeLep2, idxLep2, scaleLep2);

	verbosity_ = verbosity;

	// use default settings for SVfit
	// can be switched by individual functions
	fitMethod_ = "MarkovChain";
	addLogM_ = false;
	maxObjFunctionCalls_ = -1;
	metPower_ = -1;

	createSvFitAlgo();
}

// Constructor to be used in analysis:
// Use default one default CMS lepton (electron, muon, hadronic tau) and a fully reconstructed 3prong tau
SVfitProvider::SVfitProvider(Ntuple_Controller* const Ntp, objects::MET& met, TString typeLep1, int idxLep1, TLorentzVector lvec3ProngTau,
		int verbosity/* =1 */, double scaleLep1 /* =1 */, double scaleLep2 /* =1 */){
	ntp_ = Ntp;
	inputMet_ = met;

	addMeasuredLepton(typeLep1, idxLep1, scaleLep1);
	addFullReco3ProngTau(lvec3ProngTau);

	verbosity_ = verbosity;

	// use default settings for SVfit
	// can be switched by individual functions
	fitMethod_ = "MarkovChain";
	addLogM_ = false;
	maxObjFunctionCalls_ = -1;
	metPower_ = -1;

	createSvFitAlgo();

}

SVfitProvider::~SVfitProvider() {
	delete svFitAlgo_;
}

// create SVfitObject (includes running of SVfit algorithm)
SVFitObject SVfitProvider::runAndMakeObject(){
	run();
	SVFitObject res = makeObject(result(), get_fitMethod());
	return res;
}

// execute SVfit
void SVfitProvider::run(){
	if ( !isSetup_ )
		createSvFitAlgo();

	if(fitMethod_ == "MarkovChain")
		svFitAlgo_->integrateMarkovChain();
	else if(fitMethod_ == "Vegas")
		svFitAlgo_->integrateVEGAS();
	else if(fitMethod_ == "Minuit")
		svFitAlgo_->fit();
	else{
		Logger(Logger::Warning) << "Method " << fitMethod_ << " not available for SVfit. Using MarkovChain integration..." << std::endl;
		svFitAlgo_->integrateMarkovChain();
	}

}

// constructor from SVfitStandaloneAlgorithm
SVFitObject SVfitProvider::makeObject(const SVfitStandaloneAlgorithm* svfitAlgo, TString fitMethod){
	SVFitObject obj;

	obj.fitMethod_ = fitMethod;

	obj.mass_ = svfitAlgo->mass();
	obj.massUncert_ = svfitAlgo->massUncert();

	obj.pt_ = svfitAlgo->pt();
	obj.ptUncert_ = svfitAlgo->ptUncert();
	obj.eta_ = svfitAlgo->eta();
	obj.etaUncert_ = svfitAlgo->etaUncert();
	obj.phi_ = svfitAlgo->phi();
	obj.phiUncert_ = svfitAlgo->phiUncert();

	obj.massLmax_ = svfitAlgo->massLmax();
	obj.ptLmax_ = svfitAlgo->ptLmax();
	obj.etaLmax_ = svfitAlgo->etaLmax();
	obj.phiLmax_ = svfitAlgo->phiLmax();

	obj.fittedTauLeptons_ = svfitAlgo->fittedTauLeptons();
	obj.measuredTauLeptons_ = svfitAlgo->measuredTauLeptons();
	obj.fittedMET_ = svfitAlgo->fittedMET();
	obj.measuredMET_ = svfitAlgo->measuredMET();

	obj.tauCorr_ = ntp_->GetTauCorrections();
	obj.muonCorr_= ntp_->GetMuonCorrections();
	obj.elecCorr_= ntp_->GetElecCorrections();
	obj.metType_ = inputMet_.metType();
	obj.addLogM_ = addLogM_;
	obj.maxObjFunctionCalls_ = maxObjFunctionCalls_;
	obj.metPower_ = metPower_;


	// check if SVfit gives a valid solution
	// Have to use const_cast here, as isValidSolution() is a non-const getter function. Please cry quietly.
	if ( const_cast<SVfitStandaloneAlgorithm*>(svfitAlgo) ->isValidSolution() )
		obj.valid_ = true;
	else
		obj.valid_ = false;

	return obj;
}

// convert TLorentzVector into ROOT::Math::LorentzVector
svFitStandalone::LorentzVector SVfitProvider::convert_p4Vect(const TLorentzVector& in){
	return svFitStandalone::LorentzVector(in.X(), in.Y(), in.Z(), in.T());
}

void SVfitProvider::addMeasuredLepton(TString type, int index, double energyScale /* = 1 */){
	svFitStandalone::MeasuredTauLepton lep;
	type.ToLower();
	if (type == "mu"){
		lep = svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToMuDecay, energyScale * ntp_->Daughters_P4(index).Pt(), ntp_->Daughters_P4(index).Eta(), ntp_->Daughters_P4(index).Phi(), energyScale * ntp_->Daughters_P4(index).M());
	}
	else if (type == "ele"){
		lep = svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToElecDecay, energyScale * ntp_->Daughters_P4(index).Pt(), ntp_->Daughters_P4(index).Eta(), ntp_->Daughters_P4(index).Phi(), energyScale * ntp_->Daughters_P4(index).M());
	}
	else if (type == "tau"){
		lep = svFitStandalone::MeasuredTauLepton(svFitStandalone::kTauToHadDecay, energyScale * ntp_->TauP4_Corrected(index).Pt(), ntp_->TauP4_Corrected(index).Eta(), ntp_->TauP4_Corrected(index).Phi(), energyScale * ntp_->TauP4_Corrected(index).M());
	}
	else
		Logger(Logger::Error) << "Object type " << type << " not implemented in SVfitProvider." << std::endl;

	inputTauLeptons_.push_back(lep);
	isSetup_ = false;
}

void SVfitProvider::addFullReco3ProngTau(TLorentzVector lv){
	svFitStandalone::MeasuredTauLepton lep(svFitStandalone::kPrompt, lv.Pt(), lv.Eta(), lv.Phi(), lv.M());
	inputTauLeptons_.push_back(lep);
	isSetup_ = false;
}
