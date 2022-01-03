#!/usr/bin/env python

import pickle
from optparse import OptionParser

ptval = 0

parser = OptionParser()
parser.add_option("--pt",help="pt range to use",type="int",dest="pt")
(options,args) = parser.parse_args()

ptval = options.pt

# muon id efficiencies
f = open('MuonEfficiencies_ISO_Run_2012ReReco_53X.pkl','r')
dict = pickle.load(f)

eta_id = ['ptabseta<0.9','ptabseta0.9-1.2','ptabseta1.2-2.1','ptabseta2.1-2.4']
pt_id = ['10_20','20_25','25_30','30_35','35_40','40_50','50_60','60_90','90_140','140_300']

for i in range(len(eta_id)):
	if ptval > len(pt_id)-1:
		break
	print pt_id[ptval]+","+eta_id[i]#+"\n"
	print dict['combRelIsoPF04dBeta<012_Tight'][eta_id[i]][pt_id[ptval]]['data/mc']['efficiency_ratio']
	print dict['combRelIsoPF04dBeta<012_Tight'][eta_id[i]][pt_id[ptval]]['data/mc']['err_hi']
	print dict['combRelIsoPF04dBeta<012_Tight'][eta_id[i]][pt_id[ptval]]['data/mc']['err_low']
	print "\n"

# isomu24_eta2p1 efficiencies
g = open('SingleMuonTriggerEfficiencies_eta2p1_Run2012ABCD_v5trees.pkl','r')
dict = pickle.load(g)

eta_trig = ['PT_ABSETA_Barrel_0to0p9','PT_ABSETA_Transition_0p9to1p2','PT_ABSETA_Endcaps_1p2to2p1']
pt_trig = ['25_30','30_35','35_40','40_50','50_60','60_90','90_140','140_500']

for j in range(len(eta_trig)):
	if ptval > len(pt_trig)-1:
		break
	print pt_trig[ptval]+","+eta_trig[j]#+"\n"
	print dict['IsoMu24']['TightID_IsodB'][eta_trig[j]][pt_trig[ptval]]['data/mc']['efficiency_ratio']
	print dict['IsoMu24']['TightID_IsodB'][eta_trig[j]][pt_trig[ptval]]['data/mc']['err_hi']
	print dict['IsoMu24']['TightID_IsodB'][eta_trig[j]][pt_trig[ptval]]['data/mc']['err_low']
	print "\n"
