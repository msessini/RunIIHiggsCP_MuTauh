#include "TauSolver.h"
#include "PDG_Var.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

TauSolver::TauSolver(TVector3 Tau_, TLorentzVector a1_):
  Tau(Tau_),
  a1(a1_)
{
  if(Tau.Mag()!=0)Tau*=1/Tau.Mag();
}

TauSolver::~TauSolver(){

}

void TauSolver::SolvebyRotation(TLorentzVector &Tau1,TLorentzVector &Tau2,TLorentzVector &nu1,TLorentzVector &nu2, int mode){
  TLorentzVector a1rot=a1;
  double phi(Tau.Phi()),theta(Tau.Theta());
  a1rot.RotateZ(-phi);
  a1rot.RotateY(-theta);
  double Enu1(0), Enu2(0), Ea1(a1.E()),Pz(a1rot.Pz()),Pt(a1rot.Pt()),ma1(a1.M());
  if(mode==PZ)PzSolver(Enu1,Enu2,Ea1,ma1,Pz,Pt);
  if(mode==E)ESolver(Enu1,Enu2,Ea1,ma1,Pz,Pt);
  TLorentzVector Neutrino1(-a1rot.Px(),-a1rot.Py(),sqrt(Enu1*Enu1-Pt*Pt),Enu1);
  Neutrino1.RotateY(theta);
  Neutrino1.RotateZ(phi);
  Tau1=a1+Neutrino1;
  nu1=Neutrino1;
  TLorentzVector Neutrino2(-a1rot.Px(),-a1rot.Py(),sqrt(Enu2*Enu2-Pt*Pt),Enu2);
  Neutrino2.RotateY(theta);
  Neutrino2.RotateZ(phi);
  Tau2=a1+Neutrino2;
  nu2=Neutrino2;
}

void TauSolver::SolvebyProjection(TLorentzVector &Tau1,TLorentzVector &Tau2,TLorentzVector &nu1,TLorentzVector &nu2, int mode){
  double Enu1(0), Enu2(0), Ea1(a1.E()),Pz(0),Pt(0),ma1(a1.M());
  Pz=Tau.Dot(a1.Vect());
  TVector3 Pz_vec=Tau;
  Pz_vec*=Pz;
  TVector3 Pt_vec=a1.Vect();
  Pt_vec-=Pz_vec;
  Pt=Pt_vec.Mag();
  if(mode==PZ) PzSolver(Enu1,Enu2,Ea1,ma1,Pz,Pt);
  if(mode==E) ESolver(Enu1,Enu2,Ea1,ma1,Pz,Pt);
  double Pnuz1=0;if(Enu1-Pt>0)Pnuz1=sqrt(Enu1*Enu1-Pt*Pt);
  double Pnuz2=0;if(Enu2-Pt>0)Pnuz1=sqrt(Enu2*Enu2-Pt*Pt);
  TLorentzVector Neutrino1(Pnuz1*Tau.Px()-Pt_vec.Px(),Pnuz1*Tau.Py()-Pt_vec.Py(),Pnuz1*Tau.Pz()-Pt_vec.Pz(),Enu1);
  Tau1=a1+Neutrino1;
  nu1=Neutrino1;
  TLorentzVector Neutrino2(Pnuz2*Tau.Px()-Pt_vec.Px(),Pnuz2*Tau.Py()-Pt_vec.Py(),Pnuz2*Tau.Pz()-Pt_vec.Pz(),Enu2);
  Tau1=a1+Neutrino2;
  nu2=Neutrino1;
}

void TauSolver::quadratic(double &x_plus,double &x_minus,double a, double b, double c){
  double R=b*b-4*a*c;
  if(R<0){R=0;}
  x_minus=(-b-sqrt(R))/(2.0*a);
  x_plus=(-b+sqrt(R))/(2.0*a);
}

void TauSolver::ESolver(double &Enu1,double &Enu2,double Ea1,double ma1, double Pz, double Pt){
  double mtau(PDG_Var::Tau_mass());
  double a=1-Ea1*Ea1/(Pz*Pz);
  double K=(mtau*mtau-ma1*ma1-2*Pt*Pt);
  double b=K*Ea1/(Pz*Pz);
  double c=-(Pt*Pt+K*K/(4*Pz*Pz));
  quadratic(Enu1,Enu2,a,b,c);
}

void TauSolver::PzSolver(double &Enu1,double &Enu2,double Ea1,double ma1, double Pz, double Pt){
  double mtau(PDG_Var::Tau_mass());
  double a=1-Pz*Pz/(Ea1*Ea1);
  double K=(mtau*mtau-ma1*ma1-2.0*Pt*Pt);
  double b=-K*Pz/(Ea1*Ea1);
  double c=Pt*Pt-K*K/(4.0*Ea1*Ea1);
  double Pnuz1(0),Pnuz2(0);
  quadratic(Pnuz1,Pnuz2,a,b,c);
  Enu1=sqrt(Pnuz1*Pnuz1+Pt*Pt);
  Enu2=sqrt(Pnuz2*Pnuz2+Pt*Pt);
  if(Pnuz1<0)Enu1=Pt;
  if(Pnuz2<0)Enu2=Pt;
}





///////////////////////////////////////////////////////////////////////////////////////
//
// Computed Euler Angles for Tau: See Kuhn and Mirkes Phys. Lett. B286 (1992) pg 381
//  
//
bool TauSolver::EulerAnglesfor3prong(std::vector<TLorentzVector> Particle, std::vector<float> Charge,float &cosbeta, float &gamma, bool sortbymass, bool sortbyq1q2crossq3){
  TLorentzVector Q_LV(0,0,0,0),q1_LV(0,0,0,0),q2_LV(0,0,0,0),q3_LV(0,0,0,0);
  // Begin sorting input parameters q1 and q2 are sorted in s' frame
  if(Particle.size()!=Charge.size() || Charge.size()!=3){
    return false;
  }
  double chargesum=0;
  for(unsigned int i=0; i<Charge.size();i++){
    chargesum+=Charge.at(i);
  }
  bool foundq1=false;
  for(unsigned int i=0; i<Charge.size();i++){
    Q_LV+=Particle.at(i);
    if(chargesum!=Charge.at(i)){ 
      q3_LV=Particle.at(i); 
      Logger(Logger::Debug) << "Found q3: " << i << " " << Charge.at(i) << " "
	   << Particle.at(i).Px() << " " << Particle.at(i).Py() << " " << Particle.at(i).Pz() << endl;
    }
    else if(!foundq1){
      q1_LV=Particle.at(i);foundq1=true;
      Logger(Logger::Debug) << "Found q1: " << i << " " << Charge.at(i) << " "
           << Particle.at(i).Px() << " " << Particle.at(i).Py() << " " << Particle.at(i).Pz() << endl;
    }
    else{ 
      q2_LV=Particle.at(i);
      Logger(Logger::Debug) << "Found q1: " << i << " " << Charge.at(i) << " "
           << Particle.at(i).Px() << " " << Particle.at(i).Py() << " " << Particle.at(i).Pz() << endl;
    }
  }
  Logger(Logger::Debug) << "Q" << Q_LV.M() << endl;
  Logger(Logger::Debug) << "Q direction (noboost or rot) " << Q_LV.Px() << " " <<  Q_LV.Py() << " " << Q_LV.Pz() << endl;
  if(Q_LV.M()>=PDG_Var::Tau_mass())return false;
  // Rotate and boost system  to enter s' frame
  double phi(Q_LV.Phi()),theta(Q_LV.Theta());
  Q_LV.RotateZ(-phi);
  Q_LV.RotateY(-theta);
  
  Logger(Logger::Debug) << "Q direction (noboost - check rotation) " << Q_LV.Px() << " " <<  Q_LV.Py() << " " << Q_LV.Pz() << endl;

  q1_LV.RotateZ(-phi);
  q1_LV.RotateY(-theta);
  
  q2_LV.RotateZ(-phi);
  q2_LV.RotateY(-theta);
  
  q3_LV.RotateZ(-phi);
  q3_LV.RotateY(-theta);
  
  q1_LV.Boost(-Q_LV.BoostVector());
  q2_LV.Boost(-Q_LV.BoostVector());
  q3_LV.Boost(-Q_LV.BoostVector());
    
  // sort q1 and q2 either by mass or track by momentum in hhh system
  if(sortbymass){
    TLorentzVector q13_LV=q3_LV;q13_LV+=q1_LV;
    TLorentzVector q23_LV=q3_LV;q23_LV+=q2_LV;
    TLorentzVector tmp;
    if(q13_LV.M()<q23_LV.M()){
      tmp=q1_LV;
      q1_LV=q2_LV;
      q2_LV=tmp;
    }
  }
  else{
    if(q1_LV.P()>q2_LV.P()){
      TLorentzVector tmp;
      tmp=q1_LV;
      q1_LV=q2_LV;
      q2_LV=tmp;
    }
  }
  
  // now compute nl,np
  TVector3 nl(-Q_LV.Px(),-Q_LV.Py(),-Q_LV.Pz()); 
  TVector3 q1(q1_LV.Px(),q1_LV.Py(),q1_LV.Pz()); 
  TVector3 q2(q2_LV.Px(),q2_LV.Py(),q2_LV.Pz()); 
  TVector3 q3(q3_LV.Px(),q3_LV.Py(),q3_LV.Pz());
  TLorentzVector Q_a1frame=q1_LV; Q_a1frame+=q2_LV+q3_LV;
  Logger(Logger::Debug) <<  "Q in a1 rest frame (checks boost is correct): " << Q_a1frame.Px() << " " <<  Q_a1frame.Py() << " " <<  Q_a1frame.Pz() << " " <<  Q_a1frame.M() << endl;
  Logger(Logger::Debug) << "nl direction (checks nl direction - should be along z) " << -Q_LV.Px() << " " << -Q_LV.Py() << " " << -Q_LV.Pz() << endl;
  Logger(Logger::Debug) << "mag (before normilization) " << q1.Mag() << " " << q2.Mag() << " " << q3.Mag() << " " <<endl;
  if(nl.Mag()>0 && q1.Mag()>0 &&q2.Mag()>0 && q3.Mag()>0 ){
    nl*=1/nl.Mag();
    q1*=1/q1.Mag();
    q2*=1/q2.Mag();
    q3*=1/q3.Mag();
    TVector3 nperp=q1.Cross(q2);
    nperp*=1/nperp.Mag();
    if(sortbyq1q2crossq3){
      TVector3 q12=q1;
      q12+=q2;
      TVector3 nperp12=q12.Cross(q2);
      if(nperp12.Dot(nperp)<0)nperp*=-1;
    }
    Logger(Logger::Debug) << "norm. mag (check normilization) " << q1.Mag() << " " << q2.Mag() << " " << q3.Mag() << " " << nperp.Mag() <<endl;
    TVector3 nlcrossnperp=nl.Cross(nperp);
    
    cosbeta=nl.Dot(nperp);
    Logger(Logger::Debug) << "Dot products (show nperp is really normal to plane): " << q1.Dot(nperp) << " " << q2.Dot(nperp) << " " << q3.Dot(nperp) << " " << nl.Dot(nperp)  << " " << cosbeta << endl;
    double sine_gamma=nlcrossnperp.Dot(q3)/nlcrossnperp.Mag();
    double cosine_gamma=-nl.Dot(q3)/nlcrossnperp.Mag();
    gamma=atan2(sine_gamma,cosine_gamma);
    return true;
  }
 return false;
}
