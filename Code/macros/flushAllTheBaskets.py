#!/usr/bin/env python

import argparse
import ROOT

parser = argparse.ArgumentParser(description='Prepare ROOT files for hadd: Loop through all Sets and call FlushBasket on the tree.')
parser.add_argument('nSets', help="Amount of Sets you have in your directory.", type=int)
args = parser.parse_args()

for i in range(args.nSets):
    set = i+1
    print "Flushing Set", set
    
    file = ROOT.TFile("Set_"+str(set)+"/SKIMMED_NTUP.root","update")
    tree = file.Get("t")
    tree.FlushBaskets()
    tree.Write()
    file.Close()
