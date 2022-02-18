def config_legend(legend,CPstate,genlevel,year):
    if(genlevel):
        legend.SetHeader('CP ' + CPstate + ', gen. level, ' + year,'C')
    else:
        legend.SetHeader('CP ' + CPstate + ', reco. level, '+ year,'C')
    legend.SetBorderSize(1)
    legend.SetLineColor(1)
    legend.SetTextSize(0.03)
    legend.SetTextFont(42)
