#ifndef MyTest_h
#define MyTest_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "SVFitStorage.h"

class MyTest : public Selection {

 public:
  MyTest(TString Name_, TString id_);
  virtual ~MyTest();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0,PrimeVtx,NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables and Histos
  std::vector<TH1D> DaughtersPt;
  std::vector<TH1D> MissingTEnergy;
  std::vector<TH1D> NumVertices;


};
#endif
