def config_legend(legend,CPstate,genlevel):
    if(genlevel):
        legend.SetHeader('CP ' + CPstate + ', gen. level','C')
    else:
        legend.SetHeader('CP ' + CPstate + ', reco. level','C')
    legend.SetBorderSize(1)
    legend.SetLineColor(1)
    legend.SetTextSize(0.03)
    legend.SetTextFont(42)
