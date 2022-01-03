#!/usr/bin/env python

import argparse
import ROOT
import sys

# allow accessing additional data types
ROOT.gROOT.ProcessLine('.L Loader.C+')
#ROOT.gInterpreter.GenerateDictionary("std::vector<std::vector<float>>", "vector;vector;float")

parser = argparse.ArgumentParser(description='Check for each object type if all vectors in Ntuple have the same size.')
parser.add_argument('headerFile', help="Location of the header code file. This is used to determine the vector names.")
parser.add_argument('sourceFile', help="Location of the source code file.")
parser.add_argument('NtupleFile', help="Location of the input root file containing the Ntuple.")
parser.add_argument('--treeName', dest="treeName", default="t",help="Name of the tree in the Ntuple file.")
parser.add_argument('--verbose', action='store_true',help="Print more stuff.")
args = parser.parse_args()

if (args.verbose):
    print "header file: ", args.headerFile
    print "source file: ", args.sourceFile

header = open(args.headerFile,'r')
source = open(args.sourceFile,'r')

# create a dictionary to hold the lists of variable names
variables = {"rest":[],'Vtx':[], 'Muon':[], 'PFTau':[], 'Electron':[], 'PFJet':[], 'MET':[], 'Track':[], 'MC_':[], 'MCTau':[],'MCSignalParticle':[]}

# analyze source file, to make sure branches are created for each variable
lBranches = []
for line in source:
    # looks for branches
    if "Branch" in line:
        if line.strip().startswith('//'):
            continue
        b = line.split("&")[-1]
        b = b.strip() # remove linebreak
        b = b[0:-2] # remove the );
        lBranches.append(b)

# analyze header file and get all defined vectors
varsNotInTree = set()
for line in header:
    if line.strip().startswith("std::vector"):
        var = line.strip().split(">")[-1]
        if (not ")" in var) and (not "myManual" in var): # don't consider function definitions
            var = var.split(";")[0].strip()
            # loop over object types
            foundCategory = False
            for key in variables.keys():
                if (not foundCategory) and (var.startswith(key)):
                    if not var in lBranches:
                        varsNotInTree.add(var)
                        continue
                    variables[key].append(var)
                    foundCategory = True
            if not foundCategory:
                variables['rest'].append(var)

print "\nThe following categories of variables are considered:"
for key in variables.keys():
    if (args.verbose):
        print "\n++++++++++", key
        print variables[key]
    else:
        print "   * ", key

print "\nThe variables in the rest category are not analysed. Please take care of them yourself:"
print variables["rest"]


## now load the ROOT file
ntupleFile = ROOT.TFile(args.NtupleFile, "READ")
tree = ntupleFile.Get(args.treeName)
nevts = tree.GetEntries()

# container to hold all bad variables
badVars = set()
emptyVars = set()

# check all categories individually
for key in variables.keys():
    if key == "rest":
        continue # do not analyze the rest category
    print "\nStart check of category", key
    # loop over events
    for event in tree:
        # loop over variables
        foundSizes = set()
        for var in variables[key]:
            if var in badVars:
                continue
            functionCall = "event."+var+".size()"
            size = eval(functionCall)
            if size == 0:
                emptyVars.add(var)
            foundSizes.add(size)
            #print var, "has size", size
            if len(foundSizes) > 1:
                badVars.add(var)
                foundSizes.discard(size)
                #print "WARNING:", var, "is of different length than other variables!"

print "\n\n\nTest has been performed on", nevts, "events"
if len(varsNotInTree) != 0:
    print "There are variables defined, which are not stored in the Ntuple (no branch defined):"
    for v in sorted(varsNotInTree):
        print "   * ", v
if len(emptyVars) != 0:
    print "The following variables seem to be empty (at least sometimes):"
    for v in sorted(emptyVars):
        print "   * ", v
if len(badVars) == 0:
    print "No inconsistency has been found."
else:
    print " !!! WARNING !!! \nThe following variables seem not to be filled correctly:"
    for v in sorted(badVars):
        print "   * ", v

