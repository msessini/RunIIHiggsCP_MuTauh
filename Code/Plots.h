#ifndef Plots_h
#define Plots_h

#include <string.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include<vector>


class Plots {
 public:
 Plots();
 ~Plots();

 enum PLOTLABEL {none=0,cmsInternal,cmsPrivate,cmsPreliminary,cmsPublic};

 enum PLOTTYPE {cmsStyle1,cmsStyle2};

 void Set_Plot_Type(TString style, TString label);
 void SaveHistograms(TString File, std::vector<TString> HistogramNames);

 void CMSStyle1();
 void CMSStyle2();
 void CMSPrivateLabel(Double_t x,Double_t y,bool Preliminary,Color_t color);
 void CMSLabel(Double_t x,Double_t y,Color_t color);

 void Plot1D(std::vector<TH1D> histo,std::vector<int> colour,std::vector<TString> legend);
 void Plot1D(std::vector<std::vector<TH1D> > histo,std::vector<int> colour,std::vector<TString> legend);
 void Plot2D(std::vector<TH2D>  histo,std::vector<int> colour,std::vector<TString> legend);
 void Plot3D(std::vector<TH3F>  histo,std::vector<int> colour,std::vector<TString> legend);
 void Make_Figure(TString name,TString cap);
 void Plot1DSignificance(std::vector<TH1D> histo, bool gt,bool lt,std::vector<int> colour,std::vector<TString> legend);
 void Plot1Dsigtobkg(std::vector<TH1D> histo, bool gt,bool lt,std::vector<int> colour,std::vector<TString> legend);
 void Plot1D_DataMC_Compare(std::vector<TH1D> histo,std::vector<int> colour,std::vector<TString> legend);

 private:
 static TString File_;
 static std::vector<TString> HistogramNames_;
 static int plotLabel;
 bool doscale;
 bool doprofiles;
 bool dooneprofile;
};
#endif
