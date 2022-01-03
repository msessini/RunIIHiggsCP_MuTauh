#ifndef TauSolver_h
#define TauSolver_h

#include "Selection.h"
#include <vector>
#include "TVector3.h"
#include "TLorentzVector.h"

class TauSolver {

 public:
  TauSolver(TVector3 Tau_, TLorentzVector a1_);
  virtual ~TauSolver();

  void SolvebyRotation(TLorentzVector &Tau1,TLorentzVector &Tau2,TLorentzVector &nu1,TLorentzVector &nu2, int mode);
  void SolvebyProjection(TLorentzVector &Tau1,TLorentzVector &Tau2,TLorentzVector &nu1,TLorentzVector &nu2, int mode);
  bool EulerAnglesfor3prong(std::vector<TLorentzVector> Particle, std::vector<float> Charge,float &cosbeta, float &gamma, bool sortbymass=false, bool sortbyq1q2crossq3=true);

  enum SolutionType{E,PZ};

 private:
  void quadratic(double &x_plus,double &x_minus,double a, double b, double c);
  void ESolver(double &Enu1,double &Enu2,double Ea1,double ma1, double Pz, double Pt);
  void PzSolver(double &Enu1,double &Enu2,double Ea1,double ma1, double Pz, double Pt);

  TVector3 Tau;
  TLorentzVector a1;
};
#endif
