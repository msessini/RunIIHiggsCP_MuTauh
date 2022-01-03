#ifndef Example_h
#define Example_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "SVFitStorage.h"

class Example : public Selection {

 public:
  Example(TString Name_, TString id_);
  virtual ~Example();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0,PrimeVtx,NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables

  std::vector<TH1D> NVtx;
  std::vector<TH1D> NGoodVtx;
  std::vector<TH1D> NTrackperVtx;

};
#endif
