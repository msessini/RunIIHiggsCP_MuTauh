from CPmethods import build_methods, methods_translator
from fit import fit_style, fit_function

hist_style = [
    {'color': 1, 'marker': 21, 'width': 1, 'size': 0.8},
    {'color': 8, 'marker': 26, 'width': 1, 'size': 0.8},
    {'color': 9, 'marker': 25, 'width': 1, 'size': 0.8},
    {'color': 13, 'marker': 27, 'width': 1, 'size': 0.8}
]

ch_translator = {
    "A1A1": r'a_{1} + a_{1}',
    "RHORHO": r'\rho + \rho',
    "PIONPION": r'\pi + \pi',
    "A1RHO": r'a_{1} + \rho',
    "A1PION": r'a_{1} + \pi',
    "RHOPION": r'\pi + \rho',
    "A1MU": r'\a_{1} + \mu',
    "RHOMU": r'\rho + \mu',
    "PIONMU": r'\pi + \mu'
}

def create_names(channel,process,genlevel=False):
    methods = build_methods(channel,genlevel)
    if genlevel==True:
        return [('hcptautau_default_' + method + 'AcopAngleTruthMC_' + process + '_tautau') for method in methods]
    else:
        return [('hcptautau_default_' + method + 'AcopAngleMC_' + process + '_tautau') for method in methods]

def create_title(channel,process):
    if(process == 'all'):
      return r'H \rightarrow \tau\tau \rightarrow ' + ch_translator[channel]
    else:
      return process + r' \rightarrow \tau\tau \rightarrow ' + ch_translator[channel]

def hist_dict(file,channel,process,genlevel):
    return [
        {"name": name,
         "TH1": file.Get(name),
         "style": hist_style[i],
         "legend": methods_translator[method],
         "fit": fit_style[i],
         "TF1": fit_function
        }
    for i,(name, method) in enumerate(zip(create_names(channel,process,genlevel),build_methods(channel,genlevel)))
    ]

