#ifndef SkimConfig_h
#define SkimConfig_h

#include <vector>
#include "TString.h"
#include "TH1D.h"

class SkimConfig {

 public:
  SkimConfig();
  virtual ~SkimConfig();

  bool Load();
  bool Load(TString Name_);
  void CheckNEvents(std::vector<Long64_t> ids, std::vector<float> nevts);
  void SaveEfficiency(TString Name,std::vector<Long64_t> ids,std::vector<TH1D> NPassed, std::vector<TH1D> NPassed_noweight);
  void ApplySkimEfficiency(std::vector<Long64_t> ids,std::vector<TH1D> &NPassed, std::vector<TH1D> &NPassed_noweight);
  double GetNEvents(Long64_t id);

 private:
  bool CovertToHistoFormat();

  static std::vector<Long64_t>	   SkimIDs;
  static std::vector<float>        NEvents;
  static std::vector<float>        NEventsErr;
  static std::vector<float>        NEvents_sel;
  static std::vector<float>        NEventsErr_sel;
  static std::vector<float>        NEvents_noweight;
  static std::vector<float>        NEvents_noweight_sel;
  static bool loaded;
  static bool converted;
};
#endif
