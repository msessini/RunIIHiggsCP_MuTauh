processes = [
  'ggf',
  'vbf'
]

hist={
  'MTT': {
	'tau1': {'Eta': [], 'Phi': [], 'P': [], 'E': []},
	'tau2': {'Eta': [], 'Phi': [], 'P': [], 'E': []}
  },
  'SVFit': {
	'tau1': {'Eta': [], 'Phi': [], 'P': [], 'E': []},
	'tau2': {'Eta': [], 'Phi': [], 'P': [], 'E': []}
  },
  'Mixed': {
	'tau1': {'Eta': [], 'Phi': [], 'P': [], 'E': []},
	'tau2': {'Eta': [], 'Phi': [], 'P': [], 'E': []}
  },
  'GEF': {
	'tauH': {'Eta': [], 'Phi': [], 'P': [], 'E': []}
  }
}

def hist_dict(tfile):
    for key in tfile.GetListOfKeys():
      h = key.ReadObj()
      if h.ClassName() == 'TH1D':
        if 'Delta' in h.GetName() and any(process in h.GetName() for process in processes):
          if 'MTT' in h.GetName():
            if 'tau1' in h.GetName():
              if 'Etatau1' in h.GetName():
                hist['MTT']['tau1']['Eta'].append(h.GetName())
              if 'Phitau1' in h.GetName():
                hist['MTT']['tau1']['Phi'].append(h.GetName())
              if 'Ptau1' in h.GetName():
                hist['MTT']['tau1']['P'].append(h.GetName())
              if 'Etau1' in h.GetName():
                hist['MTT']['tau1']['E'].append(h.GetName())
            if 'tau2' in h.GetName():
              if 'Etatau2' in h.GetName():
                hist['MTT']['tau2']['Eta'].append(h.GetName())
              if 'Phitau2' in h.GetName():
                hist['MTT']['tau2']['Phi'].append(h.GetName())
              if 'Ptau2' in h.GetName():
                hist['MTT']['tau2']['P'].append(h.GetName())
              if 'Etau2' in h.GetName():
                hist['MTT']['tau2']['E'].append(h.GetName())
          if 'SVFit' in h.GetName():
            if 'tau1' in h.GetName():
              if 'Etatau1' in h.GetName():
                hist['SVFit']['tau1']['Eta'].append(h.GetName())
              if 'Phitau1' in h.GetName():
                hist['SVFit']['tau1']['Phi'].append(h.GetName())
              if 'Ptau1' in h.GetName():
                hist['SVFit']['tau1']['P'].append(h.GetName())
              if 'Etau1' in h.GetName():
                hist['SVFit']['tau1']['E'].append(h.GetName())
            if 'tau2' in h.GetName():
              if 'Etatau2' in h.GetName():
                hist['SVFit']['tau2']['Eta'].append(h.GetName())
              if 'Phitau2' in h.GetName():
                hist['SVFit']['tau2']['Phi'].append(h.GetName())
              if 'Ptau2' in h.GetName():
                hist['SVFit']['tau2']['P'].append(h.GetName())
              if 'Etau2' in h.GetName():
                hist['SVFit']['tau2']['E'].append(h.GetName())
          if 'Mixed' in h.GetName():
            if 'tau1' in h.GetName():
              if 'Etatau1' in h.GetName():
                hist['Mixed']['tau1']['Eta'].append(h.GetName())
              if 'Phitau1' in h.GetName():
                hist['Mixed']['tau1']['Phi'].append(h.GetName())
              if 'Ptau1' in h.GetName():
                hist['Mixed']['tau1']['P'].append(h.GetName())
              if 'Etau1' in h.GetName():
                hist['Mixed']['tau1']['E'].append(h.GetName())
            if 'tau2' in h.GetName():
              if 'Etatau2' in h.GetName():
                hist['Mixed']['tau2']['Eta'].append(h.GetName())
              if 'Phitau2' in h.GetName():
                hist['Mixed']['tau2']['Phi'].append(h.GetName())
              if 'Ptau2' in h.GetName():
                hist['Mixed']['tau2']['P'].append(h.GetName())
              if 'Etau2' in h.GetName():
                hist['Mixed']['tau2']['E'].append(h.GetName())
          if 'GEF' in h.GetName():
            if 'tauH' in h.GetName():
              if 'EtatauH' in h.GetName():
                hist['GEF']['tauH']['Eta'].append(h.GetName())
              if 'PhitauH' in h.GetName():
                hist['GEF']['tauH']['Phi'].append(h.GetName())
              if 'PtauH' in h.GetName():
                hist['GEF']['tauH']['P'].append(h.GetName())
              if 'EtauH' in h.GetName():
                hist['GEF']['tauH']['E'].append(h.GetName())
    return hist
