#ifndef Selection_Base_h
#define Selection_Base_h

#include "Ntuple_Controller.h"
#include "TString.h"
#include <vector>

class Selection_Base {

 public:
  Selection_Base(TString Name_, TString id_);
  virtual ~Selection_Base();
  virtual void  Event(){
    if(isNtp){
      Ntp->SetSysID(sysid);
      ResetEvent();
      doEvent();
    }
  };


  virtual void  Configure()=0;
  virtual void  Finish()=0;
  virtual void  LoadResults(std::vector<TString> files)=0;
  virtual bool Passed()=0;

  virtual TString Get_Analysis(){return Analysis;}
  virtual TString Get_Name(){return Name;};
  virtual TString Get_SysType(){return systype;};
  virtual void Set_Ntuple(Ntuple_Controller *Ntp_){Ntp=Ntp_;isNtp=true;sysid=Ntp->SetupSystematics(systype);}
  virtual void EvaluateSystematics(Selection_Base* &selectionsys, double w)=0;

  enum SelectionMode {ANALYSIS,RECONSTRUCT};
  enum RunType {GRID,Local};
  enum Verbosity {FATAL,WARNING,VERBOSE,DEBUG};
  void SetMode(int mode_);
  void SetRunType(int runtype_);
  void SetDetail(bool i){doDetails=i;}
  void SetLumi(double l){Lumi=l;}
  double GetLumi(){return Lumi;}

 protected:
  virtual void doEvent()=0;
  virtual void ResetEvent()=0;

  Ntuple_Controller *Ntp;
  double GeV;
  TString Analysis;
  TString Name;
  TString systype;
  int sysid;
  bool isNtp;
  int mode;
  int runtype;
  bool doDetails;
  double Lumi;
};
#endif
