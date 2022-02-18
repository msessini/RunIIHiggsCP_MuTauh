import argparse
from ROOT import TFile, TH1D
from TH1maker import make_TH1
from dicts import hist_dict
 
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--inputFile', help='Input File Name', required=True)
    args = parser.parse_args()
    
    fname = args.inputFile
    tfile = TFile.Open(fname)

    test = make_TH1(tfile,'GEF','tauH','P')
