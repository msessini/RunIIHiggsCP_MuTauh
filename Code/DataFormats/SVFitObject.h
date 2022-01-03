#ifndef SVFitObject_H
#define SVFitObject_H

// Object class to hold information about SVfit result
// to be used for storage of SVfit results on dCache.
// This class is basically a slimmed down version of SVfitStandaloneAlgorithm.

#include "Math/LorentzVector.h"
#include "Math/Vector3D.h"
#include "Math/PxPyPzE4D.h"
#include "TString.h"

class SVFitObject{

	// give SVfitProvider full access to SVFitObject class
	// and use SVfitProvider::makeObject(const SVfitStandaloneAlgorithm*, TString fitMethod)
	// to create an instance of SVFitObject
	friend class SVfitProvider;

 public:
	// default constructor, creates an invalid object
	SVFitObject(): valid_(false) {};
	virtual ~SVFitObject() {};

	// There is no proper constructor for this class here on purpose,
	// in order to keep this data container class independent from other
	// code. Otherwise CINT has problems with that.
	// To creat an object from SVfitStandaloneAlgorithm object, have a look at
	// SVfitProvider::makeObject(const SVfitStandaloneAlgorithm*, TString fitMethod)

	// validity MUST be checked before using
	bool isValid() const {return valid_;}

	// compare with other SVFitObject (allow for numerical differences)
	bool operator==(const SVFitObject& other){
		if ( fabs(mass_ - other.get_mass()) > 0.001 ) 				return false;
		if ( fabs(massUncert_ - other.get_massUncert()) > 0.001 )	return false;
		return true;
	}

	//// getter functions
	// kinematics of di-tau system
	double get_mass() const {return mass_;}
	double get_massUncert() const {return massUncert_;}
	double get_pt() const {return pt_;}
	double get_ptUncert() const {return ptUncert_;}
	double get_eta() const {return eta_;}
	double get_etaUncert() const {return etaUncert_;}
	double get_phi() const {return phi_;}
	double get_phiUncert() const {return phiUncert_;}

	// SVFit method
	TString get_fitMethod() const {return fitMethod_;}

	// maximum of likelihood scan
	double get_massLmax() const {return massLmax_;}
	double get_phiLmax() const {return phiLmax_;}
	double get_etaLmax() const {return etaLmax_;}
	double get_ptLmax() const {return ptLmax_;}

	// typedefs (similar to SVFitStandalone definitions)
	typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> > Vector;
	typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

	// original and re-fitted MET and tau leptons
	const SVFitObject::Vector& get_fittedMET() const {return fittedMET_;}
	const std::vector<SVFitObject::LorentzVector>& get_fittedTauLeptons() const {return fittedTauLeptons_;}
	const SVFitObject::Vector& get_measuredMET() const {return measuredMET_;}
	const std::vector<SVFitObject::LorentzVector>& get_measuredTauLeptons() const {return measuredTauLeptons_;}

	// definition of setup used to create SVFit mass
	const TString& get_elecCorr() const {return elecCorr_;}
	const TString& get_muonCorr() const {return muonCorr_;}
	const TString& get_tauCorr() const {return tauCorr_;}
	const TString& get_metType() const {return metType_;}
	bool get_addLogM() const {return addLogM_;}
	int get_maxObjFunctionCalls() const {return maxObjFunctionCalls_;}
	int get_metPower() const {return metPower_;}


private:
	// Validity of SVfit object
	bool valid_;

	// method used to obtain diTau mass
	// can be "unknown", "MarkovChain", "Vegas", "Minuit"
	TString fitMethod_;

	// mass of the di-tau system
	double mass_;
	double massUncert_;

	// kinematics of the di-tau system (available only if Markov-Chain was used)
	double pt_;
	double ptUncert_;
	double eta_;
	double etaUncert_;
	double phi_;
	double phiUncert_;

	// maxima of likelihood scan
	double massLmax_;
	double ptLmax_;
	double etaLmax_;
	double phiLmax_;

	// 4-vectors of fitted tau leptons
	std::vector<LorentzVector> fittedTauLeptons_;
	// 4-vectors of measured tau leptons
	std::vector<LorentzVector> measuredTauLeptons_;
	// return spacial vector of the fitted MET
	Vector fittedMET_;
	// return spacial vector of the measured MET
	Vector measuredMET_;

	// information on the setup used to produce SVFit mass
	TString tauCorr_;
	TString muonCorr_;
	TString elecCorr_;
	TString metType_;
	bool addLogM_;
	int maxObjFunctionCalls_;
	int metPower_;

};

#endif 
