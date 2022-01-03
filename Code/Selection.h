#ifndef Selection_h
#define Selection_h

#include "Selection_Base.h"
#include "Ntuple_Controller.h"
#include "TH1D.h"
#include "TH2D.h"
#include <vector>
#include "HistoConfig.h"

class Selection : public Selection_Base {

 public:
  Selection(TString Name_, TString id_);
  virtual ~Selection();

  virtual void  Finish();
  virtual void  LoadResults(std::vector<TString> files);

  virtual bool Passed();
  virtual bool NMinusL(int a, int b=-1, int c=-1, int d=-1, int e=-1);
  virtual bool NMinus1(int a);
  virtual bool NMinus2(int a, int b);

  std::vector<TH1D>&                 Get_Npassed(){return Npassed;}
  std::vector<TH1D>&                 Get_Npassed_noweight(){return Npassed_noweight;}
  std::vector<std::vector<TH1D> >&   Get_Nminus1(){return Nminus1;}
  std::vector<std::vector<TH1D> >&   Get_Nminus0(){return Nminus0;}
  std::vector<std::vector<TH1D> >&   Get_Nminus1dist(){return Nminus1dist;}
  std::vector<std::vector<TH1D> >&   Get_Accumdist(){return Accumdist;}
  std::vector<std::vector<TH1D>* >&  Get_Extradist1d(){return Extradist1d;}
  std::vector<std::vector<TH2D>* >&  Get_Extradist2d(){return Extradist2d;}
  std::vector<std::vector<TH3F>* >&  Get_Extradist3d(){return Extradist3d;}

  virtual double Compute(double thisdata,double thissignal, double thissignalTotal, double thisbkg,double data,double signal,
			 double signalTotal, double bkg);
  virtual void EvaluateSystematics(Selection_Base* &selectionsys, double w);
  static TString splitString(const std::string &s, char delim, std::string splitpoint);

  void Save(TString fName);

 protected:
  virtual bool AnalysisCuts(int t,double w,double wobjs=1.0);
  virtual void Store_ExtraDist()=0;
  virtual void ConfigureHistograms();
  void SetDataIdx(unsigned int d){data=d;}
  unsigned int DataIdx(){return data;}
  virtual void ResetEvent();
  virtual void ScaleAllHistOfType(unsigned int t,float);
  void suppressDrawingHistOfType(unsigned int t);
  bool passAllBut(std::vector<unsigned int> index);
  bool passAllBut(unsigned int i_cut);
  bool passAllUntil(unsigned int i_cut);
  double scaleFactorToLumi(unsigned int id);

  HistoConfig HConfig;


  std::vector<Long64_t> types;
  std::vector<TString> legend;
  std::vector<TH1D> Npassed; //[type]
  std::vector<TH1D> Npassed_noweight; //[type]
  std::vector<float> CrossSectionandAcceptance; //[type]
  std::vector<float> kFactor; //[type]
  std::vector<int> colour; //[type]
  std::vector<float> nevents_noweight_default; // to keep track of number of event

  std::vector<TString> title;
  std::vector<std::vector<TH1D> >  Nminus1;      //[cut n-1][type]
  std::vector<std::vector<TH1D> >  Nminus0;      //[cut n-0][type]
  std::vector<std::vector<TH1D> >  Nminus1dist;  //[distindx][type]
  std::vector<std::vector<TH1D> >  Accumdist;    //[distindx][type]
  std::vector<std::vector<TH1D>* >  Extradist1d;
  std::vector<std::vector<TH2D>* >  Extradist2d;
  std::vector<std::vector<TH3F>* > Extradist3d;
  std::vector<float> value; 
  std::vector<bool>  pass; 
  std::vector<float> cut; 
  std::vector<std::vector<float> > dist; 
  std::vector<bool>   distindx;

  int NGoodFiles;
  int NBadFiles;
  std::vector<TString>  ListofBadFiles;

  bool histsAreScaled; // info if histograms have been scaled already or not

 private:
  bool isStored;
  unsigned int data;

};
#endif
