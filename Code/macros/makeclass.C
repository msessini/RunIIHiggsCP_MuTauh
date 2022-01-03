void makeclass(TString file){

  TFile *_file0 = TFile::Open(file);
  _file0->Cd("HTauTauTree");
  TTree* treePtr = (TTree*)_file0 ->Get("HTauTauTree/HTauTauTree");
  //   TChain c("HTauTauTree");
  //   c.Add(file);
   treePtr->MakeClass("NtupleReader");
}
