from histograms import hist_dict, create_title
from legend import config_legend
from saver import create_dir, save
from Ymaxfinder import findYmax
from ROOT import THStack, TCanvas, TLegend

def plot(file,CPstate,channel,process,year,output,genlevel):
    legend = TLegend(0.2, 0.7, 0.6, 0.9)
    config_legend(legend,CPstate,genlevel,year)
    canvas = TCanvas()
    canvas.SetGrid()
    Ymax = findYmax(file,channel,process,genlevel)
    stack = THStack()
    for i, hist in enumerate(hist_dict(file,channel,process,genlevel)):
        hist['TH1'].SetLineColor(hist['style']['color'])
        hist['TH1'].SetLineWidth(hist['style']['width'])
        hist['TH1'].SetMarkerColor(hist['style']['color'])
        hist['TH1'].SetMarkerStyle(hist['style']['marker'])
        hist['TH1'].SetMarkerSize(hist['style']['width'])
        if(hist['TH1'].Integral()>0):
          hist['TH1'].Scale(1/hist['TH1'].Integral())
        hist['TH1'].GetYaxis().SetRangeUser(0,Ymax)
        hist['TF1'].SetLineStyle(hist['fit']['style'])
        hist['TF1'].SetLineColor(hist['fit']['color'])
        hist['TF1'].SetLineWidth(hist['fit']['width'])
        hist['TH1'].Fit(hist['TF1'])
        stack.Add(hist['TH1'])
        a = hist['TF1'].GetParameter(0)
        b = hist['TF1'].GetParameter(1)
        if(a and b != 0):
          ratio = "{:.3f}".format(abs(a/b))
        else: ratio = str(0)
        legend.AddEntry(hist['TH1'],hist['legend'] + ', |a/b| = ' + ratio,'lep')

    stack.SetTitle(create_title(channel,process))
    stack.Draw('NOSTACK')
    stack.GetXaxis().SetTitle(r'\phi_{CP} (rad)')
    stack.GetYaxis().SetTitle('a.u')
    legend.Draw()

    dir_name = channel
    subdir_name = channel+'_'+process
    save(canvas,dir_name+'/'+subdir_name+'/'+output)
