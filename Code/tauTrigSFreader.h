#ifndef tauTrigSFreader_h
#define tauTrigSFreader_h


#include <string>
#include "TFile.h"
#include "TF1.h"

class tauTrigSFreader {
    public:
        tauTrigSFreader(std::string Isotype);
        ~tauTrigSFreader();
        double getSF(double pt, int decayMode);
    private:
        TFile* fIn_;
        
        // 3 TF1 per decay mode 0, 1, and 10
        TF1* fMC_[3];
        TF1* fDataBG_[3];
        TF1* fDataH_[3];
        
        double fracBG_;
        double fracH_;

        double xmin_;
        double xmax_;
};


#endif
