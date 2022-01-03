#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='Check if all vectors saved in ROOT branches are cleared somewhere.')
parser.add_argument('headerFile', help="Location of the header code file. Can be identical to source code file.")
parser.add_argument('sourceFile', help="Location of the source code file.")
args = parser.parse_args()

print "header file: ", args.headerFile
print "source file: ", args.sourceFile

header = open(args.headerFile,'r')
source = open(args.sourceFile,'r')

lVar = []
# analyze header file and get all defined vectors
for line in header:
    if "vector" in line:
        var = line.split()[-1]
        var = var[0:-1] # remove the ;
        if not ")" in var: # don't consider function definitions
            lVar.append(var)



# analyze the source file
lBranches = []
lClean = []

# loop through file
for line in source:
    # looks for branches
    if "Branch" in line:
        b = line.split("&")[-1]
        b = b[0:-4] # remove the ); and linebreak
        lBranches.append(b)
    # look for vector clearing
    if "clear()" in line:
        c = line.split(".clear()")[0]
        c = c.strip()
        if not "//" in c:
            lClean.append(c)

toClean = set(lVar).intersection(lBranches)

notCleaned = [x for x in toClean if x not in lClean]

if notCleaned != []:
    print "The following vectors are defined and stored in vectors, but NOT CLEANED:"
    for i in notCleaned:
        print i
else:
    print "Don't worry: All vectors are cleared."

