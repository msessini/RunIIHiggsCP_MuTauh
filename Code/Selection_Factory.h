#ifndef Selection_Factory_h
#define Selection_Factory_h

#include "Selection_Base.h"

class Selection_Factory {

 public:
  Selection_Factory();
  virtual ~Selection_Factory();

  Selection_Base* Factory(TString Analysis, TString UncertType,int mode, int runtype,double lumi);
    
};
#endif



