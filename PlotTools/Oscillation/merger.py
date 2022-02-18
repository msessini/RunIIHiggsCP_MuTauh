from ROOT import TFile,TMath,TH1D

processes = [
	'ggf',
	'vbf'
]

def merge(file,channel,CPstate):

  hist = []

  hcptautau_default_polarimetricAcopAngleTruthMC_all_tautau = TH1D('hcptautau_default_polarimetricAcopAngleTruthMC_all_tautau',' hcptautau_default_polarimetricAcopAngleTruthMC_all_tautau',5,0,2*TMath.Pi())
  hcptautau_default_decayplaneAcopAngleTruthMC_all_tautau = TH1D('hcptautau_default_decayplaneAcopAngleTruthMC_all_tautau',' hcptautau_default_decayplaneAcopAngleTruthMC_all_tautau',5,0,2*TMath.Pi())
  hcptautau_default_impactparameterAcopAngleTruthMC_all_tautau = TH1D('hcptautau_default_impactparameterAcopAngleTruthMC_all_tautau',' hcptautau_default_impactparameterAcopAngleTruthMC_all_tautau',5,0,2*TMath.Pi())
  hcptautau_default_DPIPAcopAngleTruthMC_all_tautau = TH1D('hcptautau_default_DPIPAcopAngleTruthMC_all_tautau',' hcptautau_default_DPIPAcopAngleTruthMC_all_tautau',5,0,2*TMath.Pi())
  hcptautau_default_PVIPAcopAngleTruthMC_all_tautau = TH1D('hcptautau_default_PVIPAcopAngleTruthMC_all_tautau',' hcptautau_default_PVIPAcopAngleTruthMC_all_tautau',5,0,2*TMath.Pi())

  hcptautau_default_polarimetricAcopAngleMC_all_tautau = TH1D('hcptautau_default_polarimetricAcopAngleMC_all_tautau',' hcptautau_default_polarimetricAcopAngleMC_all_tautau',5,0,2*TMath.Pi())
  hcptautau_default_polarimetricGEFAcopAngleMC_all_tautau = TH1D('hcptautau_default_polarimetricGEFAcopAngleMC_all_tautau',' hcptautau_default_polarimetricGEFAcopAngleMC_all_tautau',5,0,2*TMath.Pi())
  hcptautau_default_decayplaneAcopAngleMC_all_tautau = TH1D('hcptautau_default_decayplaneAcopAngleMC_all_tautau',' hcptautau_default_decayplaneAcopAngleMC_all_tautau',5,0,2*TMath.Pi())
  hcptautau_default_impactparameterAcopAngleMC_all_tautau = TH1D('hcptautau_default_impactparameterAcopAngleMC_all_tautau',' hcptautau_default_impactparameterAcopAngleMC_all_tautau',5,0,2*TMath.Pi())
  hcptautau_default_DPIPAcopAngleMC_all_tautau = TH1D('hcptautau_default_DPIPAcopAngleMC_all_tautau',' hcptautau_default_DPIPAcopAngleMC_all_tautau',5,0,2*TMath.Pi())
  hcptautau_default_PVIPAcopAngleMC_all_tautau = TH1D('hcptautau_default_PVIPAcopAngleMC_all_tautau',' hcptautau_default_PVIPAcopAngleMC_all_tautau',5,0,2*TMath.Pi())

  tfile = TFile.Open(file)
  
  for key in tfile.GetListOfKeys():
    h = key.ReadObj()
    if h.ClassName() == 'TH1D' and 'AcopAngleTruth' in h.GetName() and any(process in h.GetName() for process in processes):
      if 'polarimetric' in h.GetName():
        hcptautau_default_polarimetricAcopAngleTruthMC_all_tautau.Add(tfile.Get(h.GetName()))
      if 'decayplane' in h.GetName():
	hcptautau_default_decayplaneAcopAngleTruthMC_all_tautau.Add(tfile.Get(h.GetName()))
      if 'impactparameter' in h.GetName():
	hcptautau_default_impactparameterAcopAngleTruthMC_all_tautau.Add(tfile.Get(h.GetName()))
      if 'DPIP' in h.GetName():
	hcptautau_default_DPIPAcopAngleTruthMC_all_tautau.Add(tfile.Get(h.GetName()))
      if 'PVIP' in h.GetName():
        hcptautau_default_PVIPAcopAngleTruthMC_all_tautau.Add(tfile.Get(h.GetName()))

    if h.ClassName() == 'TH1D' and 'AcopAngleMC' in h.GetName() and any(process in h.GetName() for process in processes):
      if 'polarimetricAcop' in h.GetName():
        hcptautau_default_polarimetricAcopAngleMC_all_tautau.Add(tfile.Get(h.GetName()))
      if 'polarimetricGEF' in h.GetName():
        hcptautau_default_polarimetricGEFAcopAngleMC_all_tautau.Add(tfile.Get(h.GetName()))
      if 'decayplane' in h.GetName():
        hcptautau_default_decayplaneAcopAngleMC_all_tautau.Add(tfile.Get(h.GetName()))
      if 'impactparameter' in h.GetName():
        hcptautau_default_impactparameterAcopAngleMC_all_tautau.Add(tfile.Get(h.GetName()))
      if 'DPIP' in h.GetName():
        hcptautau_default_DPIPAcopAngleMC_all_tautau.Add(tfile.Get(h.GetName()))
      if 'PVIP' in h.GetName():
        hcptautau_default_PVIPAcopAngleMC_all_tautau.Add(tfile.Get(h.GetName()))

  tfile.Close()

  hist.append(hcptautau_default_polarimetricAcopAngleTruthMC_all_tautau)
  hist.append(hcptautau_default_decayplaneAcopAngleTruthMC_all_tautau)
  hist.append(hcptautau_default_impactparameterAcopAngleTruthMC_all_tautau)
  hist.append(hcptautau_default_DPIPAcopAngleTruthMC_all_tautau)
  hist.append(hcptautau_default_PVIPAcopAngleTruthMC_all_tautau)
  hist.append(hcptautau_default_polarimetricAcopAngleMC_all_tautau)
  hist.append(hcptautau_default_polarimetricGEFAcopAngleMC_all_tautau)
  hist.append(hcptautau_default_decayplaneAcopAngleMC_all_tautau)
  hist.append(hcptautau_default_impactparameterAcopAngleMC_all_tautau)
  hist.append(hcptautau_default_DPIPAcopAngleMC_all_tautau)
  hist.append(hcptautau_default_PVIPAcopAngleMC_all_tautau)

  mergedfile = TFile('merged_'+channel+'_'+CPstate+'.root','RECREATE')

  for h in hist:
      h.Write()

  mergedfile.Close()

  return mergedfile.GetName()




