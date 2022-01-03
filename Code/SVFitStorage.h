#ifndef SVFitStorage_h
#define SVFitStorage_h

#include <vector>
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "TString.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TChainIndex.h"
#include "TFile.h"
#include "SVFitObject.h"
#include "DataStorage.h"
#include "TBranch.h"

class SVFitStorage : public DataStorage {
 public:
  SVFitStorage();
  ~SVFitStorage();

  void Configure(TString datasetName, TString suffix = "");

  void SaveTree();
  void SaveEvent(Int_t RunNumber, Int_t LumiNumber, Int_t EventNumber, SVFitObject* svfit);
  // obtain SVFitObject from Tree. Make sure to test validity of object
  SVFitObject* GetEvent(UInt_t RunNumber, UInt_t LumiNumber, UInt_t EventNumber);
  
  bool isConfigured(){return isConfigured_;}

 private:
  void LoadTree();
  bool isTreeInFile(TString fileName);

  TFile *outfile_;
  TTree *outtree_;
  TChain *intree_;
  TTreeIndex *index_;
  
  TString treeName_;
  TString suffix_; // optional identifier for modifications (e.g. systematics)
  TString storageFileName_;

  UInt_t RunNumber_;
  UInt_t LumiNumber_;
  UInt_t EventNumber_;
  SVFitObject *svfit_;

  TBranch *b_RunNumber_;
  TBranch *b_LumiNumber_;
  TBranch *b_EventNumber_;
  TBranch *b_svfit_;

  bool isConfigured_;
  bool intreeLoaded_;
};
#endif
