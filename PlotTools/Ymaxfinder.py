from histograms import hist_dict

def findYmax(file,channel,process,genlevel):
  Ymax = 0.
  for hist in hist_dict(file,channel,process,genlevel):
    if(hist['TH1'].Integral()>0):
      hist['TH1'].Scale(1/hist['TH1'].Integral())
    Y = hist['TH1'].GetMaximum(1.)
    if(Y>Ymax):
      Ymax = Y

  return 1.8*Ymax
