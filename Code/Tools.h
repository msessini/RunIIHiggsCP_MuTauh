#ifndef Tools_h
#define Tools_h

#include "TMath.h"
#include "TLorentzVector.h"

namespace Tools {

  inline double DeltaPhi(double phi1, double phi2){
    double dphi=fabs(phi1-phi2);
    if (dphi>TMath::Pi())   dphi=2*TMath::Pi()-dphi;
    double sign=1;
    if(phi1-phi2<0) sign=-1;
    return dphi*sign;
  }
  inline double DeltaPhi(TLorentzVector LV1, TLorentzVector LV2){
    return  DeltaPhi(LV1.Phi(),LV2.Phi());
  }

  inline double DeltaEta(double eta1, double eta2){
    double deta=fabs(eta1-eta2);
    return deta;
  }
  inline double DeltaEta(TLorentzVector LV1, TLorentzVector LV2){
    return DeltaEta(LV1.Eta(),LV2.Eta());
  }


  inline double dr(double phi1, double eta1,double phi2, double eta2){
    return sqrt(pow(DeltaEta(eta1,eta2),2.0)+pow(DeltaPhi(phi1,phi2),2.0));
  }
  
  inline double dr(TLorentzVector LV1, TLorentzVector LV2){
    return LV1.DeltaR(LV2);
  }

 
}
#endif
