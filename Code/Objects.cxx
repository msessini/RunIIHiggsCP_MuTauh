/*
 * Objects.cxx
 *
 *  Created on: May  16, 2017
 *      Author: Cherepanov
 */

#include "Objects.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "Ntuple_Controller.h"
#include "TVector2.h"

namespace objects {

// default constructor
MET::MET() {}

// constructor of MET object from Ntuple_Controller and type of MET
MET::MET(Ntuple_Controller* const Ntp, TString met_type) {
	metType_ = met_type;

	Init(Ntp);

}

void MET::Init(Ntuple_Controller* const Ntp){
	if(metType_ == "Uncorr"){
	   et_             = Ntp->MET();
	   phi_            = Ntp->METphi();
	   ex_             = Ntp->MET()*cos(phi_);
	   ey_             = Ntp->MET()*sin(phi_);
	   significance_   = Ntp->PFMETsignif();
	   significanceXX_ = Ntp->PFMETCov00();
	   significanceXY_ = Ntp->PFMETCov01();
	   significanceYY_ = Ntp->PFMETCov11();
	   hasSignificance_= true;
	}
	// if(metType_ == "CorrT0rt"){   //  give here other type of met
	//    et_             = Ntp->MET_CorrT0rt_et();
	//    phi_            = Ntp->MET_CorrT0rt_phi();
	//    ex_             = Ntp->MET_CorrT0rt_ex();
	//    ey_             = Ntp->MET_CorrT0rt_ey();
	//    hasSignificance_= false;
	// }

	else
		Logger(Logger::Error) << "MET type " << metType_ << " unknown!" << std::endl;
}

MET::~MET() {
}

void MET::subtractNeutrino(LorentzVectorParticle neutrino){
	metType_ += "NeutrinoSubtracted";

	TVectorD metvec(2), nuvec(2);
	metvec[0] = ex_;
	metvec[1] = ey_;
	nuvec[0] = neutrino.LV().Px();
	nuvec[1] = neutrino.LV().Py();
	TVectorD diff = metvec - nuvec;
	ex_ = diff[0];
	ey_ = diff[1];
	et_ = sqrt(ex_*ex_ + ey_*ey_);
	TVector2 temp(ex_, ey_);
	phi_ = temp.Phi();

	significanceXX_ = significanceXX_ + neutrino.Covariance(LorentzVectorParticle::px,LorentzVectorParticle::px);
	significanceXY_ = significanceXY_ + neutrino.Covariance(LorentzVectorParticle::px,LorentzVectorParticle::py);
	significanceYY_ = significanceYY_ + neutrino.Covariance(LorentzVectorParticle::py,LorentzVectorParticle::py);
	TMatrixD metmat = significanceMatrix();
	if (fabs(metmat.Determinant())>0.00001){
		metmat.Invert();
		significance_ = diff * (metmat * diff);
	}
	else{
		Logger(Logger::Warning) << "Cannot update significance of MET object." << std::endl;
		significance_ = -1;
		hasSignificance_ = false;
	}
}

TMatrixD MET::significanceMatrix() const {
	TMatrixD covMET(2,2);
	covMET[0][0] = significanceXX_;
	covMET[1][0] = significanceXY_;
	covMET[0][1] = significanceXY_;
	covMET[1][1] = significanceYY_;
	return covMET;
}

Vector3D MET::met3D() const {
	return Vector3D(ex_, ey_, 0.);
}

} /* namespace objects */
