from ROOT import TH1D, TFile
from dicts import hist_dict

hist_style = [
    {'color': 1, 'marker': 21, 'width': 1, 'size': 0.8},
    {'color': 8, 'marker': 26, 'width': 1, 'size': 0.8},
    {'color': 9, 'marker': 25, 'width': 1, 'size': 0.8},
]

def make_TH1(tfile,algo,tau,obs):
  if obs == 'P' or obs == 'E':
    x = 1
  if obs == 'Phi' or obs == 'Eta':
    x = 0.2
  th1 = TH1D(' ',' ',50,-x,x)
  hist_list = hist_dict(tfile)
  for i, hist in enumerate(hist_list[algo][tau][obs]):
    th1.Add(tfile.Get(hist))

  return th1
  
