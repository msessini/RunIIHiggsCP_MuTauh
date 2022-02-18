#from histograms import hist_dict, create_title
#from legend import config_legend
from saver import create_dir, save
from TH1maker import make_TH1
from ROOT import THStack, TCanvas, TLegend

def plot(tfile,channel,output):
    legend = TLegend(0.2, 0.7, 0.6, 0.9)
    #config_legend(legend,CPstate,genlevel,year)
    canvas = TCanvas()
    #canvas.SetGrid()
    #Ymax = findYmax(file,channel,process,genlevel)
    stack = THStack()
    th1P = make(
    for i, hist in enumerate(hist_dict(file,channel,process,genlevel)):
        hist['TH1'].SetLineColor(hist['style']['color'])
        hist['TH1'].SetLineWidth(hist['style']['width'])
        hist['TH1'].SetMarkerColor(hist['style']['color'])
        hist['TH1'].SetMarkerStyle(hist['style']['marker'])
        hist['TH1'].SetMarkerSize(hist['style']['width'])
        if(hist['TH1'].Integral()>0):
          hist['TH1'].Scale(1/hist['TH1'].Integral())
        hist['TH1'].GetYaxis().SetRangeUser(0,Ymax)
        stack.Add(hist['TH1'])
        #legend.AddEntry(hist['TH1'],hist['legend'] + ', |a/b| = ' + ratio,'lep')

    stack.SetTitle(create_title(channel,process))
    stack.Draw('NOSTACK')
    #stack.GetXaxis().SetTitle(r'\phi_{CP} (rad)')
    #stack.GetYaxis().SetTitle('a.u')
    #legend.Draw()

    dir_name = test
    save(canvas,dir_name+'/'+output)
