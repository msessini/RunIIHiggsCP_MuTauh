import argparse
from ROOT import TFile
from plotter import plot

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--oddFile', help='Odd Input File Name', required=True)
    parser.add_argument('-e', '--evenFile', help='Odd Input File Name', required=True)
    parser.add_argument('-c', '--channel', help='Available channels are : pipi, rhorho, a1a1, a1pi, a1rho, pirho, a1mu, rhomu, pimu', required=True)
    parser.add_argument('-p', '--process', help='Production process', default='ggfH')
    parser.add_argument('-y', '--year', help='Year (2016,2017,2018)', required=True)
    parser.add_argument('--genLevel', dest='genLevel', action='store_true')
    parser.set_defaults(genLevel=False)
    args = parser.parse_args()

    evenfile = TFile(args.evenFile)
    oddfile = TFile(args.oddFile)
    channel = args.channel
    process = args.process
    genlevel = args.genLevel
    year = args.year

    plot(evenfile,'even',channel,process,year,output=channel+'_even_reco',genlevel=False)
    plot(evenfile,'even',channel,process,year,output=channel+'_even_gen',genlevel=True)
    plot(oddfile,'odd',channel,process,year,output=channel+'_odd_reco',genlevel=False)
    plot(oddfile,'odd',channel,process,year,output=channel+'_odd_gen',genlevel=True)
