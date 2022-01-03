#ifndef Tables_h
#define Tables_h


#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "Riostream.h"
#include <vector>
#include <fstream>
#include <istream>
#include <strstream>

class Tables {

 public:
  Tables(TString name);
  ~Tables();
  void MakeNEventsTable(std::vector<TH1D> histo,std::vector<TString> names);
  void MakeEffTable(std::vector<TH1D> histo, std::vector<TString> names,float Lumi,std::vector<float> CrossSectionandAcceptance);
  void AddPlots(std::vector<TString> names);
  void GeneratePDF();

 private:
  TString Name;
  bool plotCutOptimization;

};
#endif
